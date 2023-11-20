//#include <stdbool.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <unistd.h>
#include "generic.h"

int stat_usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "Obtain alignment statistics for each repeat subfamily, family and class.\n\n");
    fprintf(stderr, "Usage:   iteres stat [options] <chromosome size file> <repeat size file> <rmsk.txt> <bam/sam alignment file1,file2,file3...>\n\n");
    fprintf(stderr, "Options: -S       input is SAM [off]\n");
    fprintf(stderr, "         -Q       unique reads mapping Quality threshold [10]\n");
    fprintf(stderr, "         -c       coverage threshold for overlapping [0.0001]\n");
    fprintf(stderr, "         -x       discard multi-reads if mapped to different subfamily [on]\n");
    fprintf(stderr, "         -N       normalized by number of (0: reads in repeats, 1: non-redundant reads, 2: mapped reads, 3: total reads) [0])\n");
    fprintf(stderr, "         -U       unique reads normalized by number of (0: unique mapped reads in repeats, 1: unique mapped reads, 2: total reads) [0])\n");
    fprintf(stderr, "         -R       remove redundant reads [off]\n");
    fprintf(stderr, "         -T       treat 1 paired-end read as 2 single-end reads [off]\n");
    fprintf(stderr, "         -D       discard if only one end mapped in a paired end reads [off]\n");
    fprintf(stderr, "         -w       keep the wiggle file [off]\n");
    fprintf(stderr, "         -B       output bed file of mapped reads [off]\n");
    fprintf(stderr, "         -V       output bed file of unique mapped reads [off]\n");
    fprintf(stderr, "         -C       add 'chr' string as prefix of reference sequence [off]\n");
    fprintf(stderr, "         -E       extend reads to represent fragment [150], specify 0 if want no extension\n");
    fprintf(stderr, "         -I       Insert length threshold [500]\n");
    fprintf(stderr, "         -o       output prefix [basename of input without extension]\n");
    fprintf(stderr, "         -W       input the length of the cage-seq window\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    return 1;
}

/* main stat function */
int main_stat (int argc, char *argv[]) {
    char *output, *outReport, *outWig, *outWigPlus, *outWigMinus, *outbigWig, *outWigUniq, *outWigUniqPlus, *outWigUniqMinus, *outbigWigUniq, *outStat, *outFam, *outCla, *row[100], *samfilecopy;
    unsigned long long int *cnt;
    int optSam = 0, optkeepWig = 0, c, optDup = 0, optaddChr = 0, optDis = 0, optTreat = 0, optBed = 0,
        optBedUniq = 0, optdiffSubfam = 1, optcagewindow = 0;
    unsigned int optQual = 10, optNorm = 0, optisize = 500, optNorm2 = 0, optExt = 150;
//    bool cageFlag = false;
    int cageFlag = 0;
    float optCov = 0.0001;
    char *optoutput = NULL;
    char *outBed = NULL;
    char *outBedUniq = NULL;
    time_t start_time, end_time;
    struct hash *hashRmsk = newHash(0);
    struct hash *hashRep = newHash(0);
    struct hash *hashFam = newHash(0);
    struct hash *hashCla = newHash(0);
    start_time = time(NULL);
    while ((c = getopt(argc, argv, "SQ:c:xN:U:RTDwBVCo:E:I:h?W:")) >= 0) {
        switch (c) {
            case 'S': optSam = 1; break;
            case 'Q': optQual = (unsigned int)strtol(optarg, 0, 0); break;
            case 'c': optCov = atof(optarg); break;
            case 'x': optdiffSubfam = 0; break;
            case 'N': optNorm = (unsigned int)strtol(optarg, 0, 0); break;
            case 'U': optNorm2 = (unsigned int)strtol(optarg, 0, 0); break;
            case 'R': optDup = 1; break;
            case 'T': optTreat = 1; break;
            case 'D': optDis = 1; break;
            case 'w': optkeepWig = 1; break;
            case 'B': optBed = 1; break;
            case 'V': optBedUniq = 1; break;
            case 'C': optaddChr = 1; break;
            case 'E': optExt = (unsigned int)strtol(optarg, 0, 0); break;
            case 'I': optisize = (unsigned int)strtol(optarg, 0, 0); break;
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return stat_usage(); break;
            case 'W': if (optarg == NULL) {
                    // Handle -W without parameters
                    optcagewindow = 50;
//                    cageFlag = true;
                    cageFlag = 1;
                } else {
                    optcagewindow = (unsigned int)strtol(optarg, 0, 0);
//                    cageFlag = true;
                    cageFlag = 1;
                } break;
            default: return 1;
        }
    }

    printf("optwindow length: %d, bool: %d, optind: %d\n", optcagewindow, cageFlag, optind);

    if (optind + 4 > argc)
        return stat_usage();

    char *chr_size_file = argv[optind];
    char *rep_size_file = argv[optind+1];
    char *rmsk_file = argv[optind+2];
    char *sam_file = argv[optind+3];
    
    samfilecopy = cloneString(sam_file);
    int numFields = chopByChar(samfilecopy, ',', row, ArraySize(row));
    fprintf(stderr, "* Provided %i BAM/SAM file(s)\n", numFields);
    
    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(row[0])));
    }

    printf("%s\n", output);
    
    if(asprintf(&outWig, "%s.iteres.wig", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outWigPlus, "%s_+.iteres.wig", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outWigMinus, "%s_-.iteres.wig", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outbigWig, "%s.iteres.bigWig", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outWigUniq, "%s.iteres.unique.wig", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outWigUniqPlus, "%s_+.iteres.unique.wig", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outWigUniqMinus, "%s_-.iteres.unique.wig", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outbigWigUniq, "%s.iteres.unique.bigWig", output) < 0)
        errAbort("Mem Error.\n");
    if (asprintf(&outReport, "%s.iteres.report", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outStat, "%s.iteres.subfamily.stat", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outFam, "%s.iteres.family.stat", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outCla, "%s.iteres.class.stat", output) < 0)
        errAbort("Preparing output wrong");

    if(optBed){
        if(asprintf(&outBed, "%s.iteres.bed", output) < 0)
            errAbort("Mem Error.\n");
    }
    if(optBedUniq){
        if(asprintf(&outBedUniq, "%s.iteres.unique.bed", output) < 0)
            errAbort("Mem Error.\n");
    }
    
    int nindex = 0;
    if (optNorm == 0){
        nindex = 9;
    } else if (optNorm == 1){
        nindex = 8;
    } else if (optNorm == 2) {
        nindex = 6;    
    } else if (optNorm == 3) {
        nindex = 0;    
    } else{
        errAbort("Wrong normalization method specified");
    }
    
    int nindex2 = 0;
    if (optNorm2 == 0){
        nindex2 = 10;
    } else if (optNorm2 == 1){
        nindex2 = 7;
    } else if (optNorm2 == 2){
        nindex2 = 0;
    } else{
        errAbort("Wrong normalization method specified");
    }

    struct hash *chrHash = hashNameIntFile(chr_size_file);
    struct hash *repHash = hashNameIntFile(rep_size_file);
    fprintf(stderr, "* Parsing the rmsk file\n");
    rmsk2binKeeperHash(rmsk_file, chrHash, repHash, &hashRmsk, &hashRep, &hashFam, &hashCla, 0, "ALL");
	FILE *ftest1 = mustOpen("asddsk.txt", "w");    
    //sam file
    fprintf(stderr, "* Parsing the SAM/BAM file\n");
    //if (optPair){
    //    cnt = PEsamFile2nodupRepbedFile(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, optQual, 0, optDup, optaddChr, optisize);
    //} else {
    //    cnt = samFile2nodupRepbedFile(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, optQual, 0, optDup, optaddChr);
    //}

    if(cageFlag == 1)
        cnt = samFiles2nodupRepbedFileNewCage(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, 255, 0,
                                          optDup, optaddChr, optDis, optisize, optExt, optCov, optTreat, outBed, outBedUniq,
                                          optdiffSubfam, optcagewindow);
    else
        cnt = samFiles2nodupRepbedFileNew(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, optQual, 0,
                                          optDup, optaddChr, optDis, optisize, optExt, optCov, optTreat, outBed, outBedUniq,
                                          optdiffSubfam);
    fprintf(stderr, "* Writing stats and Wig file\n");
//    writeWigandStat(hashRep, hashFam, hashCla, outStat, outWig, outFam, outCla, outWigUniq, cnt[nindex], cnt[nindex2]);
    if(cageFlag == 1){
	writeWigandStatCage(hashRep, hashFam, hashCla, outStat, outWigPlus, outWigMinus, outFam, outCla, outWigUniqPlus, outWigUniqMinus, cnt[nindex], cnt[nindex2]);
    }
    writeWigandStat(hashRep, hashFam, hashCla, outStat, outWig, outFam, outCla, outWigUniq, cnt[nindex], cnt[nindex2]);

    fprintf(stderr, "* Generating bigWig files\n");
    bigWigFileCreate(outWig, rep_size_file, 256, 1024, 0, 1, outbigWig);
    bigWigFileCreate(outWigUniq, rep_size_file, 256, 1024, 0, 1, outbigWigUniq);

    //write report file
    fprintf(stderr, "* Preparing report file\n");
    writeReport(outReport, cnt, optQual, "ALL");
    
    if (!optkeepWig){
        unlink(outWig);
        unlink(outWigUniq);
    }
    //cleaning
    hashFree(&chrHash);
    hashFree(&repHash);
    hashFree(&hashRmsk);
    hashFree(&hashRep);
    hashFree(&hashFam);
    hashFree(&hashCla);
    free(outWig);
    free(outbigWig);
    free(outWigUniq);
    free(outbigWigUniq);
    free(outReport);
    free(outStat);
    free(outFam);
    free(outCla);
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}

