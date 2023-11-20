#include <argp.h>
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "binRange.h"
#include "basicBed.h"
#include "bigWig.h"
#include "sqlNum.h"
#include "obscure.h"
#include "localmem.h"
#include "dystring.h"
#include "cirTree.h"
#include "sig.h"
#include "zlibFace.h"
#include "bPlusTree.h"
#include "bbiFile.h"
#include "bwgInternal.h"
#include "sam.h"

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };

typedef unsigned int unitSize;

#define MAXCOUNT (unitSize)~0
#define MAXMESSAGE "Overflow of overlap counts. Max is %lu.  Recompile with bigger unitSize or use -max option"
#define INCWOVERFLOW(countArray,x) if(countArray[x] == MAXCOUNT) {if(!doMax) errAbort(MAXMESSAGE,(unsigned long)MAXCOUNT);} else countArray[x]++

const char *argp_program_version = "ChipSeqHandler 0.1";

const char *argp_program_bug_address = "<http://wang.wustl.edu>";

boolean doMax = FALSE;   /* if overlap count will overflow, just keep max */
boolean doZero = FALSE;  /* add blocks with 0 counts */
boolean doBed12 = FALSE;  /* expect bed12 and process block by block */
boolean doOutBounds = FALSE;  /* output min/max to stderr */
unitSize overMin = ~1;
unitSize overMax = 0;

static int blockSize = 256;
static int itemsPerSlot = 1024;
static boolean doCompress = TRUE;

struct sectionItem
/* An item in a section of a bedGraph. */
    {
    bits32 start, end;			/* Position in chromosome, half open. */
    float val;				/* Single precision value. */
    };

void writeSections(struct bbiChromUsage *usageList, struct lineFile *lf, 
	int itemsPerSlot, struct bbiBoundsArray *bounds, int sectionCount, FILE *f,
	int resTryCount, int resScales[], int resSizes[], 
	boolean doCompress, bits32 *retMaxSectionSize)
/* Read through lf, chunking it into sections that get written to f.  Save info
 * about sections in bounds. */
{
int maxSectionSize = 0;
struct bbiChromUsage *usage = usageList;
int itemIx = 0, sectionIx = 0;
bits32 reserved32 = 0;
UBYTE reserved8 = 0;
struct sectionItem items[itemsPerSlot];
struct sectionItem *lastB = NULL;
bits32 resEnds[resTryCount];
int resTry;
for (resTry = 0; resTry < resTryCount; ++resTry)
    resEnds[resTry] = 0;
struct dyString *stream = dyStringNew(0);

/* remove initial browser and track lines */
lineFileRemoveInitialCustomTrackLines(lf);

for (;;)
    {
    /* Get next line of input if any. */
    char *row[5];
    int rowSize = lineFileChopNext(lf, row, ArraySize(row));

    /* Figure out whether need to output section. */
    boolean sameChrom = FALSE;
    if (rowSize > 0)
	sameChrom = sameString(row[0], usage->name);
    if (itemIx >= itemsPerSlot || rowSize == 0 || !sameChrom)
        {
	/* Figure out section position. */
	bits32 chromId = usage->id;
	bits32 sectionStart = items[0].start;
	bits32 sectionEnd = items[itemIx-1].end;

	/* Save section info for indexing. */
	assert(sectionIx < sectionCount);
	struct bbiBoundsArray *section = &bounds[sectionIx++];
	section->offset = ftell(f);
	section->range.chromIx = chromId;
	section->range.start = sectionStart;
	section->range.end = sectionEnd;

	/* Output section header to stream. */
	dyStringClear(stream);
	UBYTE type = bwgTypeBedGraph;
	bits16 itemCount = itemIx;
	dyStringWriteOne(stream, chromId);			// chromId
	dyStringWriteOne(stream, sectionStart);		// start
	dyStringWriteOne(stream, sectionEnd);	// end
	dyStringWriteOne(stream, reserved32);		// itemStep
	dyStringWriteOne(stream, reserved32);		// itemSpan
	dyStringWriteOne(stream, type);			// type
	dyStringWriteOne(stream, reserved8);			// reserved
	dyStringWriteOne(stream, itemCount);			// itemCount

	/* Output each item in section to stream. */
	int i;
	for (i=0; i<itemIx; ++i)
	    {
	    struct sectionItem *item = &items[i];
	    dyStringWriteOne(stream, item->start);
	    dyStringWriteOne(stream, item->end);
	    dyStringWriteOne(stream, item->val);
	    }

	/* Save stream to file, compressing if need be. */
	if (stream->stringSize > maxSectionSize)
	    maxSectionSize = stream->stringSize;
	if (doCompress)
	    {
	    size_t maxCompSize = zCompBufSize(stream->stringSize);
	    char compBuf[maxCompSize];
	    int compSize = zCompress(stream->string, stream->stringSize, compBuf, maxCompSize);
	    mustWrite(f, compBuf, compSize);
	    }
	else
	    mustWrite(f, stream->string, stream->stringSize);


	/* If at end of input we are done. */
	if (rowSize == 0)
	    break;

	/* Set up for next section. */
	itemIx = 0;

	if (!sameChrom)
	    {
	    usage = usage->next;
	    assert(usage != NULL);
            if (!sameString(row[0], usage->name))
                errAbort("read %s, expecting %s on line %d in file %s\n", 
                    row[0], usage->name, lf->lineIx, lf->fileName);
	    assert(sameString(row[0], usage->name));
	    lastB = NULL;
	    for (resTry = 0; resTry < resTryCount; ++resTry)
		resEnds[resTry] = 0;
	    }
	}

    /* Parse out input. */
    lineFileExpectWords(lf, 4, rowSize);
    bits32 start = lineFileNeedNum(lf, row, 1);
    bits32 end = lineFileNeedNum(lf, row, 2);
    float val = lineFileNeedDouble(lf, row, 3);

    /* Verify that inputs meets our assumption - that it is a sorted bedGraph file. */
    if (start > end)
        errAbort("Start (%u) after end (%u) line %d of %s", start, end, lf->lineIx, lf->fileName);
    if (lastB != NULL)
        {
	if (lastB->start > start)
	    errAbort("BedGraph not sorted on start line %d of %s", lf->lineIx, lf->fileName);
	if (lastB->end > start)
	    errAbort("Overlapping regions in bedGraph line %d of %s", lf->lineIx, lf->fileName);
	}


    /* Do zoom counting. */
    for (resTry = 0; resTry < resTryCount; ++resTry)
        {
	bits32 resEnd = resEnds[resTry];
	if (start >= resEnd)
	    {
	    resSizes[resTry] += 1;
	    resEnds[resTry] = resEnd = start + resScales[resTry];
	    }
	while (end > resEnd)
	    {
	    resSizes[resTry] += 1;
	    resEnds[resTry] = resEnd = resEnd + resScales[resTry];
	    }
	}

    /* Save values in output array. */
    struct sectionItem *b = &items[itemIx];
    b->start = start;
    b->end = end;
    b->val = val;
    lastB = b;
    itemIx += 1;
    }
assert(sectionIx == sectionCount);

*retMaxSectionSize = maxSectionSize;
}

static struct bbiSummary *writeReducedOnceReturnReducedTwice(struct bbiChromUsage *usageList, 
	struct lineFile *lf, bits32 initialReduction, bits32 initialReductionCount, 
	int zoomIncrement, int blockSize, int itemsPerSlot, boolean doCompress,
	struct lm *lm, FILE *f, bits64 *retDataStart, bits64 *retIndexStart,
	struct bbiSummaryElement *totalSum)
/* Write out data reduced by factor of initialReduction.  Also calculate and keep in memory
 * next reduction level.  This is more work than some ways, but it keeps us from having to
 * keep the first reduction entirely in memory. */
{
struct bbiSummary *twiceReducedList = NULL;
bits32 doubleReductionSize = initialReduction * zoomIncrement;
struct bbiChromUsage *usage = usageList;
struct bbiSummary oneSummary, *sum = NULL;
struct bbiBoundsArray *boundsArray, *boundsPt, *boundsEnd;
boundsPt = AllocArray(boundsArray, initialReductionCount);
boundsEnd = boundsPt + initialReductionCount;

*retDataStart = ftell(f);
writeOne(f, initialReductionCount);
boolean firstRow = TRUE;

struct bbiSumOutStream *stream = bbiSumOutStreamOpen(itemsPerSlot, f, doCompress);

/* remove initial browser and track lines */
lineFileRemoveInitialCustomTrackLines(lf);

for (;;)
    {
    /* Get next line of input if any. */
    char *row[5];
    int rowSize = lineFileChopNext(lf, row, ArraySize(row));

    /* Output last section and break if at end of file. */
    if (rowSize == 0 && sum != NULL)
	{
	bbiOutputOneSummaryFurtherReduce(sum, &twiceReducedList, doubleReductionSize, 
		&boundsPt, boundsEnd, usage->size, lm, stream);
	break;
	}

    /* Parse out row. */
    char *chrom = row[0];
    bits32 start = sqlUnsigned(row[1]);
    bits32 end = sqlUnsigned(row[2]);
    float val = sqlFloat(row[3]);

    /* Update total summary stuff. */
    bits32 size = end-start;
    if (firstRow)
	{
        totalSum->validCount = size;
	totalSum->minVal = totalSum->maxVal = val;
	totalSum->sumData = val*size;
	totalSum->sumSquares = val*val*size;
	firstRow = FALSE;
	}
    else
        {
	totalSum->validCount += size;
	if (val < totalSum->minVal) totalSum->minVal = val;
	if (val > totalSum->maxVal) totalSum->maxVal = val;
	totalSum->sumData += val*size;
	totalSum->sumSquares += val*val*size;
	}

    /* If new chromosome output existing block. */
    if (differentString(chrom, usage->name))
        {
	usage = usage->next;
	bbiOutputOneSummaryFurtherReduce(sum, &twiceReducedList, doubleReductionSize,
		&boundsPt, boundsEnd, usage->size, lm, stream);
	sum = NULL;
	}

    /* If start past existing block then output it. */
    else if (sum != NULL && sum->end <= start)
	{
	bbiOutputOneSummaryFurtherReduce(sum, &twiceReducedList, doubleReductionSize, 
		&boundsPt, boundsEnd, usage->size, lm, stream);
	sum = NULL;
	}

    /* If don't have a summary we're working on now, make one. */
    if (sum == NULL)
        {
	oneSummary.chromId = usage->id;
	oneSummary.start = start;
	oneSummary.end = start + initialReduction;
	if (oneSummary.end > usage->size) oneSummary.end = usage->size;
	oneSummary.minVal = oneSummary.maxVal = val;
	oneSummary.sumData = oneSummary.sumSquares = 0.0;
	oneSummary.validCount = 0;
	sum = &oneSummary;
	}
    
    /* Deal with case where might have to split an item between multiple summaries.  This
     * loop handles all but the final affected summary in that case. */
    while (end > sum->end)
        {
	verbose(3, "Splitting start %d, end %d, sum->start %d, sum->end %d\n", start, end, sum->start, sum->end);
	/* Fold in bits that overlap with existing summary and output. */
	bits32 overlap = rangeIntersection(start, end, sum->start, sum->end);
	sum->validCount += overlap;
	if (sum->minVal > val) sum->minVal = val;
	if (sum->maxVal < val) sum->maxVal = val;
	sum->sumData += val * overlap;
	sum->sumSquares += val*val * overlap;
	bbiOutputOneSummaryFurtherReduce(sum, &twiceReducedList, doubleReductionSize, 
		&boundsPt, boundsEnd, usage->size, lm, stream);
	size -= overlap;

	/* Move summary to next part. */
	sum->start = start = sum->end;
	sum->end = start + initialReduction;
	if (sum->end > usage->size) sum->end = usage->size;
	sum->minVal = sum->maxVal = val;
	sum->sumData = sum->sumSquares = 0.0;
	sum->validCount = 0;
	}

    /* Add to summary. */
    sum->validCount += size;
    if (sum->minVal > val) sum->minVal = val;
    if (sum->maxVal < val) sum->maxVal = val;
    sum->sumData += val * size;
    sum->sumSquares += val*val * size;
    }
bbiSumOutStreamClose(&stream);

/* Write out 1st zoom index. */
int indexOffset = *retIndexStart = ftell(f);
assert(boundsPt == boundsEnd);
cirTreeFileBulkIndexToOpenFile(boundsArray, sizeof(boundsArray[0]), initialReductionCount,
    blockSize, itemsPerSlot, NULL, bbiBoundsArrayFetchKey, bbiBoundsArrayFetchOffset, 
    indexOffset, f);

freez(&boundsArray);
slReverse(&twiceReducedList);
return twiceReducedList;
}

void bedGraphToBigWig(char *inName, char *chromSizes, char *outName)
/* bedGraphToBigWig - Convert a bedGraph program to bigWig.. */
{
verboseTimeInit();
struct lineFile *lf = lineFileOpen(inName, TRUE);
struct hash *chromSizesHash = bbiChromSizesFromFile(chromSizes);
verbose(2, "%d chroms in %s\n", chromSizesHash->elCount, chromSizes);
int minDiff = 0, i;
double aveSize = 0;
bits64 bedCount = 0;
bits32 uncompressBufSize = 0;
struct bbiChromUsage *usageList = bbiChromUsageFromBedFile(lf, chromSizesHash, &minDiff, &aveSize, &bedCount);
verboseTime(2, "pass1");
verbose(2, "%d chroms in %s\n", slCount(usageList), inName);

/* Write out dummy header, zoom offsets. */
FILE *f = mustOpen(outName, "wb");
bbiWriteDummyHeader(f);
bbiWriteDummyZooms(f);

/* Write out dummy total summary. */
struct bbiSummaryElement totalSum;
ZeroVar(&totalSum);
bits64 totalSummaryOffset = ftell(f);
bbiSummaryElementWrite(f, &totalSum);

/* Write out chromosome/size database. */
bits64 chromTreeOffset = ftell(f);
bbiWriteChromInfo(usageList, blockSize, f);

/* Set up to keep track of possible initial reduction levels. */
int resTryCount = 10, resTry;
int resIncrement = 4;
int resScales[resTryCount], resSizes[resTryCount];
int res = minDiff * 2;
if (res > 0)
    {
    for (resTry = 0; resTry < resTryCount; ++resTry)
	{
	resSizes[resTry] = 0;
	resScales[resTry] = res;
	res *= resIncrement;
	}
    }
else
    resTryCount = 0;

/* Write out primary full resolution data in sections, collect stats to use for reductions. */
bits64 dataOffset = ftell(f);
bits64 sectionCount = bbiCountSectionsNeeded(usageList, itemsPerSlot);
writeOne(f, sectionCount);
struct bbiBoundsArray *boundsArray;
AllocArray(boundsArray, sectionCount);
lineFileRewind(lf);
bits32 maxSectionSize = 0;
writeSections(usageList, lf, itemsPerSlot, boundsArray, sectionCount, f,
	resTryCount, resScales, resSizes, doCompress, &maxSectionSize);
verboseTime(2, "pass2");

/* Write out primary data index. */
bits64 indexOffset = ftell(f);
cirTreeFileBulkIndexToOpenFile(boundsArray, sizeof(boundsArray[0]), sectionCount,
    blockSize, 1, NULL, bbiBoundsArrayFetchKey, bbiBoundsArrayFetchOffset, 
    indexOffset, f);
verboseTime(2, "index write");

/* Declare arrays and vars that track the zoom levels we actually output. */
bits32 zoomAmounts[bbiMaxZoomLevels];
bits64 zoomDataOffsets[bbiMaxZoomLevels];
bits64 zoomIndexOffsets[bbiMaxZoomLevels];
int zoomLevels = 0;

/* Write out first zoomed section while storing in memory next zoom level. */
if (minDiff > 0)
    {
    bits64 dataSize = indexOffset - dataOffset;
    int maxReducedSize = dataSize/2;
    int initialReduction = 0, initialReducedCount = 0;

    /* Figure out initialReduction for zoom. */
    for (resTry = 0; resTry < resTryCount; ++resTry)
	{
	bits64 reducedSize = resSizes[resTry] * sizeof(struct bbiSummaryOnDisk);
	if (doCompress)
	    reducedSize /= 2;	// Estimate!
	if (reducedSize <= maxReducedSize)
	    {
	    initialReduction = resScales[resTry];
	    initialReducedCount = resSizes[resTry];
	    break;
	    }
	}
    verbose(2, "initialReduction %d, initialReducedCount = %d\n", 
    	initialReduction, initialReducedCount);

    if (initialReduction > 0)
        {
	struct lm *lm = lmInit(0);
	int zoomIncrement = 4;
	lineFileRewind(lf);
	struct bbiSummary *rezoomedList = writeReducedOnceReturnReducedTwice(usageList, 
		lf, initialReduction, initialReducedCount,
		resIncrement, blockSize, itemsPerSlot, doCompress, lm, 
		f, &zoomDataOffsets[0], &zoomIndexOffsets[0], &totalSum);
	verboseTime(2, "writeReducedOnceReturnReducedTwice");
	zoomAmounts[0] = initialReduction;
	zoomLevels = 1;

	int zoomCount = initialReducedCount;
	int reduction = initialReduction * zoomIncrement;
	while (zoomLevels < bbiMaxZoomLevels)
	    {
	    int rezoomCount = slCount(rezoomedList);
	    if (rezoomCount >= zoomCount)
	        break;
	    zoomCount = rezoomCount;
	    zoomDataOffsets[zoomLevels] = ftell(f);
	    zoomIndexOffsets[zoomLevels] = bbiWriteSummaryAndIndex(rezoomedList, 
	    	blockSize, itemsPerSlot, doCompress, f);
	    zoomAmounts[zoomLevels] = reduction;
	    ++zoomLevels;
	    reduction *= zoomIncrement;
	    rezoomedList = bbiSummarySimpleReduce(rezoomedList, reduction, lm);
	    }
	lmCleanup(&lm);
	verboseTime(2, "further reductions");
	}

    }

/* Figure out buffer size needed for uncompression if need be. */
if (doCompress)
    {
    int maxZoomUncompSize = itemsPerSlot * sizeof(struct bbiSummaryOnDisk);
    uncompressBufSize = max(maxSectionSize, maxZoomUncompSize);
    }

/* Go back and rewrite header. */
rewind(f);
bits32 sig = bigWigSig;
bits16 version = bbiCurrentVersion;
bits16 summaryCount = zoomLevels;
bits16 reserved16 = 0;
bits32 reserved32 = 0;
bits64 reserved64 = 0;

/* Write fixed header */
writeOne(f, sig);
writeOne(f, version);
writeOne(f, summaryCount);
writeOne(f, chromTreeOffset);
writeOne(f, dataOffset);
writeOne(f, indexOffset);
writeOne(f, reserved16);	// fieldCount
writeOne(f, reserved16);	// definedFieldCount
writeOne(f, reserved64);	// autoSqlOffset
writeOne(f, totalSummaryOffset);
writeOne(f, uncompressBufSize);
for (i=0; i<2; ++i)
    writeOne(f, reserved32);
assert(ftell(f) == 64);

/* Write summary headers with data. */
verbose(2, "Writing %d levels of zoom\n", zoomLevels);
for (i=0; i<zoomLevels; ++i)
    {
    verbose(3, "zoomAmounts[%d] = %d\n", i, (int)zoomAmounts[i]);
    writeOne(f, zoomAmounts[i]);
    writeOne(f, reserved32);
    writeOne(f, zoomDataOffsets[i]);
    writeOne(f, zoomIndexOffsets[i]);
    }
/* Write rest of summary headers with no data. */
for (i=zoomLevels; i<bbiMaxZoomLevels; ++i)
    {
    writeOne(f, reserved32);
    writeOne(f, reserved32);
    writeOne(f, reserved64);
    writeOne(f, reserved64);
    }

/* Write total summary. */
fseek(f, totalSummaryOffset, SEEK_SET);
bbiSummaryElementWrite(f, &totalSum);

/* Write end signature. */
fseek(f, 0L, SEEK_END);
writeOne(f, sig);

lineFileClose(&lf);
carefulClose(&f);
}


/* definitions of structures*/

//struct hold contens from sam line
struct sam {
    unsigned int start, end, length;
    char *chr;
    char *name;
    char *seq;
    uint32_t qual:8;
    char strand;
};

/* definitions of functions */

static void outputCounts(unitSize *counts, char *chrom, unsigned size, FILE *f){
if (size == 0)
    errAbort("got 0 for size of chrom %s\n", chrom);

if (doOutBounds)
    {
    if (counts[0] < overMin)
	overMin = counts[0];
    if (counts[0] > overMax)
	overMax = counts[0];
    }

int ii;
int prevValue = counts[0];
int startPoint = 0;
for(ii=1; ii < size; ii++)
    {
    if (doOutBounds)
	{
	if (counts[ii] < overMin)
	    overMin = counts[ii];
	if (counts[ii] > overMax)
	    overMax = counts[ii];
	}
    if (counts[ii] != prevValue)
	{
	if (doZero || (prevValue != 0))
	    fprintf(f, "%s\t%u\t%u\t%u\n", chrom, startPoint, ii, prevValue);
	startPoint = ii;
	prevValue = counts[ii];
	}
    }

if (doZero || (prevValue != 0))
    fprintf(f, "%s\t%u\t%u\t%u\n", chrom, startPoint, ii, prevValue);
}

static void bedItemOverlapCount(struct hash *chromHash, char *infile, char *outfile){
unsigned maxChromSize = 0;
unitSize *counts = (unitSize *)NULL;
FILE *f = mustOpen(outfile, "w");
struct hashCookie hc = hashFirst(chromHash);
struct hashEl *hel;
while( (hel = hashNext(&hc)) != NULL) {
    unsigned num = (unsigned) ptToInt(hel->val);
    maxChromSize = max(num, maxChromSize);
}
verbose(2,"#\tmaxChromSize: %u\n", maxChromSize);
if (maxChromSize < 1)
    errAbort("maxChromSize is zero ?");

/*	Allocate just once for the largest chrom and reuse this array */
counts = needHugeMem(sizeof(unitSize) * maxChromSize);

/*	Reset the array to be zero to be reused */
memset((void *)counts, 0, sizeof(unitSize)*(size_t)maxChromSize);

unsigned chromSize = 0;
char *prevChrom = (char *)NULL;
boolean outputToDo = FALSE;
struct hash *seenHash = newHash(5);

    struct lineFile *bf = lineFileOpen(infile , TRUE);
    struct bed *bed = (struct bed *)NULL;
    char *row[12];
    int numFields = doBed12 ? 12 : 3;

    while (lineFileNextRow(bf,row, numFields))
	{
	int i;
	bed = bedLoadN(row, numFields);

	verbose(3,"#\t%s\t%d\t%d\n",bed->chrom,bed->chromStart, bed->chromEnd);

	if (prevChrom && differentWord(bed->chrom,prevChrom)) // End a chr
	    {
	    verbose(2,"#\tchrom %s done, size %d\n", prevChrom, chromSize);
	    if (outputToDo)
		outputCounts(counts, prevChrom, chromSize, f);
	    outputToDo = FALSE;
	    memset((void *)counts, 0,
		sizeof(unitSize)*(size_t)maxChromSize); /* zero counts */
	    freez(&prevChrom); 
	    // prevChrom is now NULL so it will be caught by next if!
	    }
	if ((char *)NULL == prevChrom)  // begin a chr
	    {
	    if (hashLookup(seenHash, bed->chrom))
		errAbort("ERROR:input file not sorted. %s seen before on line %d\n",
		    bed->chrom, bf->lineIx);

	    hashAdd(seenHash, bed->chrom, NULL);
	    prevChrom = cloneString(bed->chrom);
	    chromSize = hashIntVal(chromHash, prevChrom);
	    verbose(2,"#\tchrom %s starting, size %d\n", prevChrom,chromSize);
	    }
	if (bed->chromEnd > chromSize)
	    {
	    // check for circular chrM
	    if (doBed12 || bed->chromStart>=chromSize 
		|| differentWord(bed->chrom,"chrM")) 
		{
		warn("ERROR: %s\t%d\t%d", bed->chrom, bed->chromStart,
		bed->chromEnd);
		errAbort("chromEnd > chromSize ?  %d > %d", 
		    bed->chromEnd,chromSize);
		}

	    for (i = bed->chromStart; i < chromSize; ++i)
		INCWOVERFLOW(counts,i);
	    for (i = 0; i < (bed->chromEnd - chromSize); ++i)
		INCWOVERFLOW(counts,i);
	    }
	else if (doBed12)
	    {
	    int *starts = bed->chromStarts;
	    int *sizes = bed->blockSizes;
	    int *endStarts = &bed->chromStarts[bed->blockCount];

	    for(; starts < endStarts; starts++, sizes++)
		{
		unsigned int end = *starts + *sizes + bed->chromStart;
		for (i = *starts + bed->chromStart; i < end; ++i)
		    INCWOVERFLOW(counts,i);
		}
	    }
	else
	    {
	    for (i = bed->chromStart; i < bed->chromEnd; ++i)
		INCWOVERFLOW(counts, i);
	    }
	outputToDo = TRUE;
	bedFree(&bed); // plug the memory leak
	}

    lineFileClose(&bf);
    // Note, next file could be on same chr!

if (outputToDo)
    outputCounts(counts, prevChrom, chromSize, f);

if (doOutBounds)
    fprintf(stderr, "min %lu max %lu\n", (unsigned long)overMin, (unsigned long)overMax);

verbose(2,"#\tchrom %s done, size %d\n", prevChrom, chromSize);
carefulClose(&f);
freeMem(counts);
freez(&prevChrom);
// hashFreeWithVals(&chromHash, freez);
freeHash(&seenHash);
}

static int itemRgbColumn(char *column9){
int itemRgb = 0;
/*  Allow comma separated list of rgb values here   */
char *comma = strchr(column9, ',');
if (comma)
    {
    if (-1 == (itemRgb = bedParseRgb(column9)))
	errAbort("ERROR: expecting r,g,b specification, "
		    "found: '%s'", column9);
    }
else
    itemRgb = sqlUnsigned(column9);
return itemRgb;
}

char *get_filename_without_ext(char *filename) {
    char *s;
    s = malloc(strlen(filename) + 1);
    strcpy(s, filename);
    char *dot = strrchr(s, '.');
    if(!dot || dot == s) return s;
    *dot = '\0';
    return s;
}

char *get_filename_ext(char *filename) {
    char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

double cal_rpkm (unsigned int reads_count, unsigned int total_length, unsigned int mapped_reads_num) {
    return reads_count / (mapped_reads_num * 1e-9 * total_length);
}

struct arguments {
    char *args[1];
    int Sam;
    unsigned int extend;
    char *output, *sizef, *db;
};

static struct argp_option options[] = 
{
    {"db", 'd', "dataBase", 0, "database assembly, such as hg19, hg18, mm9 and rn4 (default: hg19)"},
    {"sizef", 's', "sizeFile", 0, "chromosome size file"},
    {"extend", 'e', "extendLength", 0, "Extend length (default: 150"},
    {"Sam", 'S', 0, 0, "Input is a SAM file (default: off)"},
    {"output", 'o', "outputBase", 0, "Base name for output files (default: the basename of input file)"},
    {0}
};
static char args_doc[] = "<bam/sam alignment file>";
static char doc[] = "\nChipSeqHandler program.\n\nPlease note a chromosome size file was required, size file path for assembly hg19, hg18, mm9 and rn4 was included in the program.\nYou can always use the -s option to sepcify your own size file whether your assembly was included or not.\nIf your assembly was not in the list above, using -s specify a size file, -d could be anything then.\nAlso noticed that: if reads mapped to the chromosomes which didn't existed in size file, this type of reads will be discarded.\n";

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = state->input;

    switch (key)
    {
        case 'd':
            arguments->db = arg;
            break;
        case 's':
            arguments->sizef = arg;
            break;
        case 'S':
            arguments->Sam = 1;
            break;
        case 'e':
            arguments->extend = (unsigned int)strtol(arg, NULL, 0);
            break;
        case 'o':
            arguments->output = arg;
            break;
        case ARGP_KEY_ARG:
            if( state->arg_num >= 1)
            {
                argp_usage(state);
            }
            arguments->args[state->arg_num] = arg;
            break;
        case ARGP_KEY_END:
            if (state->arg_num < 1)
            {
                argp_usage(state);
            }
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

struct sam *fetch_sa (const bam1_t *b, void *data){
    struct sam *s = malloc(sizeof(struct sam));
    samfile_t *fp = (samfile_t *) data;
    uint32_t *cigar = bam1_cigar(b);
    const bam1_core_t *c = &b->core;
    int i, l, j;
    uint8_t *ss = bam1_seq(b);
    s->name = cloneString(bam1_qname(b));
    s->chr = cloneString("*");
    s->start = 0;
    s->end = 0;
    s->length = 0;
    s->qual = 0;
    s->strand = '*';
    s->seq = needLargeMem(c->l_qseq + 1);
    for (j = 0; j < c->l_qseq; ++j) s->seq[j] = bam_nt16_rev_table[bam1_seqi(ss, j)];
    s->seq[c->l_qseq] = '\0';
    if (b->core.tid < 0) return s;
    for (i = l = 0; i < c->n_cigar; ++i) {
        int op = cigar[i]&0xf;
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
            l += cigar[i]>>4;
    }
    s->chr = cloneString(fp->header->target_name[c->tid]);
    s->start = c->pos;
    s->end = c->pos + l;
    s->length = l;
    s->qual = c->qual;
    s->strand = (c->flag&BAM_FREVERSE) ? '-' : '+';
    return s;
}

struct bed *sam2bed(struct sam *sam) {
    struct bed *ret;
    char *chr;
    AllocVar(ret);
    if (startsWith("GL", sam->chr)) {
        return 0;
    } else if (sameWord(sam->chr, "MT")) {
        chr = cloneString("chrM");
    } else if (!startsWith("chr", sam->chr)) {
        chr = catTwoStrings("chr", sam->chr);
    } else {
        chr = cloneString(sam->chr);
    }
    ret->chrom = cloneString(chr);
    ret->chromStart = sam->start;
    ret->chromEnd = sam->end;
    ret->name = cloneString(sam->seq);
    ret->score = 0;
    ret->strand[0] = sam->strand;
    ret->thickStart = 0;
    ret->thickEnd = 0;
    ret->itemRgb = itemRgbColumn((sam->strand == '-') ? "0,0,255" : "255,0,0");
    freeMem(chr);
    return ret;
}

int writesam2bed(struct sam *sam, struct hash *hash, FILE *f){
    char *chr;
    unsigned int end;
    if (startsWith("GL", sam->chr)) {
        return 0;
    } else if (sameWord(sam->chr, "MT")) {
        chr = cloneString("chrM");
    } else if (!startsWith("chr", sam->chr)) {
        chr = catTwoStrings("chr", sam->chr);
    } else {
        chr = cloneString(sam->chr);
    }
    end = (unsigned int) (hashIntValDefault(hash, chr, 2) - 1);
    if (end == 1)
        return 0;
    end = min(end, sam->end);
    fprintf(f, "%s\t%u\t%u\t%s\t%d\t%c\t%d\t%d\t%s\n", chr, sam->start, end, sam->seq, 0, sam->strand, 0, 0, (sam->strand == '-') ? "0,0,255" : "255,0,0");
    freeMem(chr);
    return 0;
}

unsigned long long int removeBedDup(char *infile, char *outfile) {
    char *preChrom = "empty", *line;
    int preStart = 0, lineSize;
    unsigned long long int i = 0;
    struct lineFile *lf = NULL;
    struct bedLine *bl;
    FILE *f = mustOpen(outfile, "w");
    lf = lineFileOpen(infile, TRUE);
    while (lineFileNext(lf, &line, &lineSize)) {
        if (line[0] == '#')
            continue;
        int num = chopByWhite(line, NULL, 0);
        if (num != 9)
            continue;
        bl = bedLineNew(line);
        if (sameWord(bl->chrom, preChrom) && (bl->chromStart == preStart))
            continue;
        fprintf(f, "%s\t%s\n", bl->chrom, bl->line);
        preChrom = cloneString(bl->chrom);
        preStart = bl->chromStart;
        i++;
    }
    carefulClose(&f);
    return i;
}

int extendBed(struct hash *hash, int extend, char *infile, char *outfile){
    struct bed *bedList = NULL, *bed;
    FILE *f  = mustOpen(outfile, "w");
    bedList = bedLoadAll(infile);
    for (bed = bedList; bed != NULL; bed = bed->next){
        int len = bed->chromEnd - bed->chromStart;
        if (len >= extend){
            warn("* Warning: read length %d longer than extend length %d, no need for extending", len, extend);
            carefulClose(&f);
            bedFreeList(&bedList);
            return 1;
        }
        int tlen = hashIntVal(hash, bed->chrom);
        if (bed->strand[0] == '+'){
            bed->chromEnd += (extend - len);
            if (bed->chromEnd >= (tlen - 1))
                bed->chromEnd = tlen - 1;
        } else {
            bed->chromStart -= (extend - len);
            if (bed->chromStart <= 0)
                bed->chromStart = 0;
        }
        bedOutputN(bed, 9, f, '\t', '\n');
    }
    carefulClose(&f);
    bedFreeList(&bedList);
    return 0;
}

void freeSam(struct sam *s){
    freeMem(s->chr);
    freeMem(s->name);
    freeMem(s->seq);
    free(s);
}

void sortBedfile(char *bedfile) {
    struct lineFile *lf = NULL;
    FILE *f = NULL;
    struct bedLine *blList = NULL, *bl;
    char *line;
    int lineSize;

    lf = lineFileOpen(bedfile, TRUE);
    while (lineFileNext(lf, &line, &lineSize)){
        if (line[0] == '#')
            continue;
        bl = bedLineNew(line);
        slAddHead(&blList, bl);
    }
    lineFileClose(&lf);

    slSort(&blList, bedLineCmp);

    f = mustOpen(bedfile, "w");
    for (bl = blList; bl != NULL; bl = bl->next){
        fprintf(f, "%s\t%s\n", bl->chrom, bl->line);
        if (ferror(f)){
    	    perror("Writing error\n");
	    errAbort("%s is truncated, sorry.", bedfile);
	}
    }
    carefulClose(&f);
    bedLineFreeList(&blList);
}

void writeReport(char *outfile, unsigned long long int *cnt){
    FILE *f = mustOpen(outfile, "w");
    fprintf(f, " Total reads: %llu\n", cnt[0]);
    fprintf(f, "Mapped reads: %llu\n", cnt[1]);
    fprintf(f, "  Used reads: %llu\n", cnt[2]);
    fprintf(f, "Unique reads: %llu\n", cnt[3]);
    carefulClose(&f);
}

unsigned long long int *samFile2nodupExtbedFile(char *samfile, char *bedfile, struct hash *hash, int isSam, unsigned int extend) {
    samfile_t *samfp;
    char *chr = NULL, *prn, *key;
    FILE *outBed;
    unsigned int start, end, cend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 4);
    unsigned long long int mapped_reads_num = 0, reads_num = 0, reads_used = 0, unique_reads = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
    struct hashEl *hel, *he;
    boolean doExtend = TRUE;
    int outputWarn = 0;
    //uint32_t *cigar;
    bam1_core_t *c;
    uint8_t *s;
    //int i, l, j;
    int j;
    if (isSam) {
        if ( (samfp = samopen(samfile, "r", 0)) == 0) {
            fprintf(stderr, "Fail to open SAM file %s\n", samfile);
            errAbort("Error\n");
        }
    } else {
        if ( (samfp = samopen(samfile, "rb", 0)) == 0) {
            fprintf(stderr, "Fail to open BAM file %s\n", samfile);
            errAbort("Error\n");
        }
    }
    outBed = mustOpen(bedfile, "w");
    prn = cloneString("empty");
    bam1_t *b = bam_init1();
    while ( samread(samfp, b) >= 0) {
        //cigar = bam1_cigar(b);
        c = &b->core;
        s = bam1_seq(b);
        if ( sameString (bam1_qname(b), prn)) 
            continue;
        reads_num++;
        if ((reads_num % 10000) == 0)
            fprintf(stderr, "\r* Processed reads: %llu", reads_num);
        prn = cloneString(bam1_qname(b));
        if (b->core.tid < 0)
            continue;
        mapped_reads_num++;
        //for (i = l = 0; i < c->n_cigar; ++i) {
        //    int op = cigar[i]&0xf;
        //    if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
        //        l += cigar[i]>>4;
        //}
        //change chr name to chr1, chr2 ...
        if (startsWith("GL", samfp->header->target_name[c->tid])) {
            continue;
        } else if (sameWord(samfp->header->target_name[c->tid], "MT")) {
            chr = cloneString("chrM");
        } else if (!startsWith("chr", samfp->header->target_name[c->tid])) {
            chr = catTwoStrings("chr", samfp->header->target_name[c->tid]);
        } else {
            chr = cloneString(samfp->header->target_name[c->tid]);
        }
        he = hashLookup(nochr, chr);
        if (he != NULL)
            continue;
        cend = (unsigned int) (hashIntValDefault(hash, chr, 2) - 1);
        if (cend == 1){
            hashAddInt(nochr, chr, 1);
            warn("* Warning: reads mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
            continue;
        }
        reads_used++;
        start = (unsigned int) c->pos;
        //remove dup first
        if (asprintf(&key, "%s:%u", chr, start) < 0)
            errAbort("Mem ERROR");
        hel = hashLookup(dup, key);
        if (hel == NULL) {
            hashAddInt(dup, key, 1);
        } else {
            continue;
        }
        unique_reads++;
        //extend
        end = min(cend, (unsigned int)c->pos + c->l_qseq);
        if ( (unsigned int)c->l_qseq >= extend){
            doExtend = FALSE;
            outputWarn++;
            if (outputWarn == 1)
                warn("* Warning: read length %d longer than extend length %u, do not extend", c->l_qseq, extend);
        }
        if (doExtend){
            if (c->flag&BAM_FREVERSE){
                start = end - extend;
                if (start < 0)
                    start = 0;
            } else {
                end = start + extend;
                end = min(cend, end);
            }
        }
        fprintf(outBed, "%s\t%u\t%u\t", chr, start, end);
        //fprintf(outBed, "%s\t", bam1_qname(b));
        for (j = 0; j < c->l_qseq; ++j) fprintf(outBed, "%c", bam_nt16_rev_table[bam1_seqi(s, j)]);
        //fprintf(outBed, "\t%d\t%c\t%d\t%d\t%s\n", 0, (c->flag&BAM_FREVERSE) ? '-' : '+', 0, 0, (c->flag&BAM_FREVERSE) ? "0,0,255" : "255,0,0");
        fprintf(outBed, "\t%d\t%c\n", 0, (c->flag&BAM_FREVERSE) ? '-' : '+');
    }
    fprintf(stderr, "\r* Processed reads: %llu", reads_num);
    bam_destroy1(b);
    samclose(samfp);
    freeMem(prn);
    freeMem(chr);
    free(key);
    freeHash(&nochr);
    freeHash(&dup);
    carefulClose(&outBed);
    cnt[0] = reads_num;
    cnt[1] = mapped_reads_num;
    cnt[2] = reads_used;
    cnt[3] = unique_reads;
    return cnt;
}

unsigned long long int *samFile2nodupExtbedFile1(char *samfile, char *bedfile, struct hash *hash, int isSam, unsigned int extend) {
    samfile_t *samfp;
    char chr[100], prn[500], key[100];
    FILE *outBed;
    unsigned int start, end, cend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 4);
    unsigned long long int mapped_reads_num = 0, reads_num = 0, reads_used = 0, unique_reads = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
    boolean doExtend = TRUE;
    int outputWarn = 0;
    if (isSam) {
        if ( (samfp = samopen(samfile, "r", 0)) == 0) {
            fprintf(stderr, "Fail to open SAM file %s\n", samfile);
            errAbort("Error\n");
        }
    } else {
        if ( (samfp = samopen(samfile, "rb", 0)) == 0) {
            fprintf(stderr, "Fail to open BAM file %s\n", samfile);
            errAbort("Error\n");
        }
    }
    outBed = mustOpen(bedfile, "w");
    strcpy(prn, "empty");
    bam1_t *b;
    bam_header_t *h;
    int8_t *buf;
    int max_buf;
    h = samfp->header;
    b = bam_init1();
    buf = 0;
    max_buf = 0;
    while ( samread(samfp, b) >= 0) {
        if ( sameString (bam1_qname(b), prn)) 
            continue;
        reads_num++;
        if ((reads_num % 10000) == 0)
            fprintf(stderr, "\r* Processed reads: %llu", reads_num);
        strcpy(prn, bam1_qname(b));
        if (b->core.tid < 0)
            continue;
        mapped_reads_num++;
        //change chr name to chr1, chr2 ...
        if (startsWith("GL", h->target_name[b->core.tid])) {
            continue;
        } else if (sameWord(h->target_name[b->core.tid], "MT")) {
            strcpy(chr,"chrM");
        } else if (!startsWith("chr", h->target_name[b->core.tid])) {
            strcpy(chr, "chr");
            strcat(chr, h->target_name[b->core.tid]);
        } else {
            strcpy(chr, h->target_name[b->core.tid]);
        }
        struct hashEl *he = hashLookup(nochr, chr);
        if (he != NULL)
            continue;
        cend = (unsigned int) (hashIntValDefault(hash, chr, 2) - 1);
        if (cend == 1){
            hashAddInt(nochr, chr, 1);
            warn("* Warning: reads mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
            continue;
        }
        reads_used++;
        int i, qlen = b->core.l_qseq;
        uint8_t *seq;
        start = (unsigned int) b->core.pos;
        //remove dup first
        //strcpy(key, chr);
        //strcat(key, ":");
        //strcat(key, (char*)start);
        if (sprintf(key, "%s:%u", chr, start) < 0)
            errAbort("Mem ERROR");
        struct hashEl *hel = hashLookup(dup, key);
        if (hel == NULL) {
            hashAddInt(dup, key, 1);
        } else {
            continue;
        }
        unique_reads++;
        //extend
        end = min(cend, (unsigned int)b->core.pos + qlen);
        if ( (unsigned int)qlen >= extend){
            doExtend = FALSE;
            outputWarn++;
            if (outputWarn == 1)
                warn("* Warning: read length %d longer than extend length %u, do not extend", qlen, extend);
        }
        if (doExtend){
            if (b->core.flag&BAM_FREVERSE){
                start = end - extend;
                if (start < 0)
                    start = 0;
            } else {
                end = start + extend;
                end = min(cend, end);
            }
        }
        fprintf(outBed, "%s\t%u\t%u\t", chr, start, end);

        if (max_buf < qlen + 1 ) {
            max_buf = qlen + 1;
            kroundup32(max_buf);
            buf = realloc(buf, max_buf);
        }
        buf[qlen] = 0;
        seq = bam1_seq(b);
        for (i = 0; i < qlen; ++i)
            buf[i] = bam1_seqi(seq, i);
        if (b->core.flag & 16) {
            for (i = 0; i < qlen>>1; ++i){
                int8_t t = seq_comp_table[buf[qlen - 1 - i]];
                buf[qlen - 1 - i] = seq_comp_table[buf[i]];
                buf[i] = t;
            }
            if (qlen&1) buf[i] = seq_comp_table[buf[i]];
        }
        for (i = 0; i < qlen; ++i)
            buf[i] = bam_nt16_rev_table[buf[i]];
        fprintf(outBed, "%s", (char*)buf);
        fprintf(outBed, "\t%d\t%c\n", 0, (b->core.flag&BAM_FREVERSE) ? '-' : '+');
    }
    fprintf(stderr, "\r* Processed reads: %llu", reads_num);
    samclose(samfp);
    free(buf);
    bam_destroy1(b);
    //bam_header_destroy(h);
    //freeMem(prn);
    //freeMem(chr);
    freeHash(&nochr);
    freeHash(&dup);
    carefulClose(&outBed);
    cnt[0] = reads_num;
    cnt[1] = mapped_reads_num;
    cnt[2] = reads_used;
    cnt[3] = unique_reads;
    return cnt;
}

/* main funtion */
int main (int argc, char *argv[]) {
    
    char *output, *chr_size_file = NULL, *outReportfile, *outExtfile, *outbedGraphfile, *outbigWigfile;
    unsigned long long int *cnt;
    struct arguments arguments;
    time_t start_time, end_time;
    start_time = time(NULL);

    /* set default parameters*/
    arguments.output = NULL;
    arguments.db = "hg19";
    arguments.extend = 150;
    arguments.sizef = NULL;
    arguments.Sam = 0;
    /* parameter magic*/
    argp_parse(&argp, argc, argv, 0, 0, &arguments);
    
    if (arguments.sizef == NULL){
        if (sameWord(arguments.db, "hg19")){
            chr_size_file = "/home/comp/twlab/twang/twlab-shared/genomes/hg19/hg19_chrom_sizes_tab";
        } else if (sameWord(arguments.db, "hg18")){
            chr_size_file = "/home/comp/twlab/twang/twlab-shared/genomes/hg18/hg18_chrom_sizes";
        } else if (sameWord(arguments.db, "mm9")){
            chr_size_file = "/home/comp/twlab/twang/twlab-shared/genomes/mm9/mm9_chrom_sizes";
        } else if (sameWord(arguments.db, "rn4")){
            chr_size_file = "/home/comp/twlab/mxie/Rat/seq/Rat_chrom_sizes";
        }
    } else{
        chr_size_file = cloneString(arguments.sizef);
    }

    if (chr_size_file == NULL)
        errAbort("A chromosome size file was required, specify it by -s option.\n");

    char *sam_file = arguments.args[0];
    if(arguments.output) {
        output = arguments.output;
    } else {
        //output = get_filename_without_ext(basename(sam_file));
        //if (asprintf(&output, get_filename_without_ext(basename(sam_file))) < 0)
        //    errAbort("Preparing output wrong");
        output = cloneString(get_filename_without_ext(basename(sam_file)));
    }
    //if(asprintf(&outBedfile, "%s.original.bed", output) < 0)
    //    errAbort("Mem Error.\n");
    //if(asprintf(&outFilterfile, "%s.filter.bed", output) < 0)
    //    errAbort("Mem Error.\n");
    if(asprintf(&outExtfile, "%s.extended.bed", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outbedGraphfile, "%s.extended.bedGraph", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outbigWigfile, "%s.bigWig", output) < 0)
        errAbort("Mem Error.\n");
    if (asprintf(&outReportfile, "%s.report", output) < 0)
        errAbort("Preparing output wrong");
    
    struct hash *hash = hashNameIntFile(chr_size_file);
    
    //sam file to bed file
    fprintf(stderr, "* Start to parse the SAM/BAM file ...\n");
    cnt = samFile2nodupExtbedFile1(sam_file, outExtfile, hash, arguments.Sam, arguments.extend);
    //sort
    //fprintf(stderr, "\n* Sorting\n");
    //bedSortFile(outBedfile, outBedfile);

    //remove dup
    //fprintf(stderr, "* Removing duplication\n");
    //uniqueBed = removeBedDup(outBedfile, outFilterfile);

    //extend and write extend bed
    //fprintf(stderr, "* Extending to %d and writing extended bed\n", arguments.extend);
    //int extendWarn = extendBed(hash, arguments.extend, outFilterfile, outExtfile);
    //if (extendWarn == 1)
    //    outExtfile = cloneString(outFilterfile);
    //
    //if (extendWarn != 1){
        //sort extend bed
    //    fprintf(stderr, "* Sorting extended bed\n");
    //    bedSortFile(outExtfile, outExtfile);
    //}

    //sort extend bed
    fprintf(stderr, "\n* Sorting extended bed\n");
    sortBedfile(outExtfile);
    
    //bedItemOverlap step
    fprintf(stderr, "* Generating bedGraph\n");
    bedItemOverlapCount(hash, outExtfile, outbedGraphfile);

    //generate bigWig
    fprintf(stderr, "* Generating bigWig\n");
    //bigWigFileCreate(outbedGraphfile, chr_size_file, 256, 1024, 0, 1, outbigWigfile);
    bedGraphToBigWig(outbedGraphfile, chr_size_file, outbigWigfile);

    //write report file
    fprintf(stderr, "* Preparing report file\n");
    writeReport(outReportfile, cnt);
    
    //cleaning
    hashFree(&hash);
    free(outExtfile);
    free(outbedGraphfile);
    free(outbigWigfile);
    free(outReportfile);
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}
