#include "generic.h"
#include <stdio.h>

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };

//Debug Tools
typedef struct Node {
    int number;
    int count;
    struct Node* next;
} Node;

Node* head = NULL;

// Function to search for a number in the list and update its count
Node* searchAndUpdate(Node* head, int number) {
    Node* current = head;
    while (current != NULL) {
        if (current->number == number) {
            current->count++;
            return head;
        }
        current = current->next;
    }

    // Number not found, create a new node and add it to the list
    Node* newNode = (Node*)malloc(sizeof(Node));
    newNode->number = number;
    newNode->count = 1;
    newNode->next = head;
    return newNode;
}

// Function to display the occurrence count of numbers in the list
void displayList(Node* head) {
    Node* current = head;
    printf("Number\tCount\n");
    while (current != NULL) {
        printf("%d\t%d\n", current->number, current->count);
        current = current->next;
    }
}

// Function to free the memory allocated for the list
void freeList(Node* head) {
    Node* current = head;
    while (current != NULL) {
        Node* next = current->next;
        free(current);
        current = next;
    }
}

/* definitions of functions */
//The iteres functions.
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

bool is_file(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISREG(buf.st_mode);
}

bool is_dir(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISDIR(buf.st_mode);
}

double cal_rpkm (unsigned long long int reads_count, unsigned long long int total_length, unsigned long long int mapped_reads_num) {
    return reads_count / (mapped_reads_num * 1e-9 * total_length);
}

double cal_rpm (unsigned long long int reads_count, unsigned long long int mapped_reads_num) {
    return reads_count / (mapped_reads_num * 1e-6 );
}

struct lineFile *lineFileOpen2(char *fileName, bool zTerm){
/* Open up a lineFile or die trying. */
if (is_dir(fileName))
    errAbort("Error: %s is a directory not a file", fileName);
struct lineFile *lf = lineFileMayOpen(fileName, zTerm);
if (lf == NULL)
    errAbort("Couldn't open %s , %s", fileName, strerror(errno));
return lf;
}

void writeReport(char *outfile, unsigned long long int *cnt, unsigned int mapQ, char *subfam){
    FILE *f = mustOpen(outfile, "w");
    fprintf(f, "total reads (pair): %llu\n", cnt[0]);
    //fprintf(f, "read ends 1: %llu\n", cnt[0]);
    //fprintf(f, "read ends 2: %llu\n", cnt[1]);
    //fprintf(f, "mapped read ends 1: %llu\n", cnt[2]);
    //fprintf(f, "mapped read ends 2: %llu\n", cnt[3]);
    //fprintf(f, "used read ends 1: %llu\n", cnt[4]);
    //fprintf(f, "used read ends 2: %llu\n", cnt[5]);
    fprintf(f, "mappable reads (pair): %llu\n", cnt[6]);
    //fprintf(f, "non-redundant mappable reads (pair): %llu\n", cnt[8]);
    fprintf(f, "uniquely mapped reads (pair) (mapQ >= %u): %llu\n", mapQ, cnt[7]);
    fprintf(f, "non-redundant uniquely mapped reads (pair): %llu\n", cnt[11]);
    fprintf(f, "mapped reads (pair) overlap with repeats but discarded due to mapped to different subfamilies: %llu\n", cnt[12]);
    fprintf(f, "mapped reads (pair) overlap with [%s] repeats: %llu\n", subfam, cnt[9]);
    fprintf(f, "uniquely mapped reads (pair) overlap with [%s] repeats: %llu\n", subfam, cnt[10]);
    carefulClose(&f);
}

void writeWigandStat(struct hash *hash, struct hash *hash1, struct hash *hash2, char *of1, char *of2, char *of3, char *of4, char *of5, unsigned long long int reads_num, unsigned long int reads_num_unique){
    FILE *f1 = mustOpen(of1, "w");
    FILE *f2 = mustOpen(of2, "w");
    FILE *f5 = mustOpen(of5, "w");

    unsigned int m;
    struct hashEl *helr;
    struct hashCookie cookier = hashFirst(hash);
    // Write the stat file
    fprintf(f1, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#subfamily", "family", "class", "consensus_length", "reads_count", "unique_reads_count", "total_length", "genome_count", "all_reads_RPKM", "all_reads_RPM", "unique_reads_RPKM", "unique_reads_RPM");
    while ( (helr = hashNext(&cookier)) != NULL ) {
        struct rep *or = (struct rep *) (helr->val);
        fprintf(f1, "%s\t%s\t%s\t%u\t%llu\t%llu\t%llu\t%llu\t%.3f\t%.3f\t%.3f\t%.3f\n", or->name, or->fname, or->cname, or->length, or->read_count, or->read_count_unique, or->total_length, or->genome_count, cal_rpkm(or->read_count, or->total_length, reads_num), cal_rpm(or->read_count, reads_num), cal_rpkm(or->read_count_unique, or->total_length, reads_num_unique), cal_rpm(or->read_count_unique, reads_num_unique));
        // Write the wig file
        if (or->length != 0){
            fprintf(f2, "fixedStep chrom=%s start=1 step=1 span=1\n", or->name);
            fprintf(f5, "fixedStep chrom=%s start=1 step=1 span=1\n", or->name);
            for (m = 0; m < or->length; m++){ 
                fprintf(f2, "%u\n", (or->bp_total)[m]);
                fprintf(f5, "%u\n", (or->bp_total_unique)[m]);
            }
        }
    }
    carefulClose(&f2);
    carefulClose(&f1);
    carefulClose(&f5);

    FILE *f3 = mustOpen(of3, "w");
    struct hashEl *hel3;
    struct hashCookie cookier3 = hashFirst(hash1);
    fprintf(f3, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#family", "class", "reads_count", "unique_reads_count", "total_length", "genome_count", "all_reads_RPKM", "all_reads_RPM", "unique_reads_RPKM", "unique_reads_RPM");
    while ( (hel3 = hashNext(&cookier3)) != NULL) {
        struct repfam *or = (struct repfam *) hel3->val;
        fprintf(f3, "%s\t%s\t%llu\t%llu\t%llu\t%llu\t%.3f\t%.3f\t%.3f\t%.3f\n", or->fname, or->cname, or->read_count, or->read_count_unique, or->total_length, or->genome_count, cal_rpkm(or->read_count, or->total_length, reads_num), cal_rpm(or->read_count, reads_num), cal_rpkm(or->read_count_unique, or->total_length, reads_num_unique), cal_rpm(or->read_count_unique, reads_num_unique));
    }
    carefulClose(&f3);
    FILE *f4 = mustOpen(of4, "w");
    struct hashEl *hel4;
    struct hashCookie cookier4 = hashFirst(hash2);
    fprintf(f4, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",  "#class", "reads_count", "unique_reads_count", "total_length", "genome_count", "all_reads_RPKM", "all_reads_RPM", "unique_reads_RPKM", "unique_reads_RPM");
    while ( (hel4 = hashNext(&cookier4)) != NULL) {
        struct repcla *or = (struct repcla *) hel4->val;
        fprintf(f4, "%s\t%llu\t%llu\t%llu\t%llu\t%.3f\t%.3f\t%.3f\t%.3f\n", or->cname, or->read_count, or->read_count_unique, or->total_length, or->genome_count, cal_rpkm(or->read_count, or->total_length, reads_num), cal_rpm(or->read_count, reads_num), cal_rpkm(or->read_count_unique, or->total_length, reads_num_unique), cal_rpm(or->read_count_unique, reads_num_unique));
    }
    carefulClose(&f4);
}

void writeWigandStatCage(struct hash *hash, struct hash *hash1, struct hash *hash2, char *of1, char *of2, char *of6, char *of3, char *of4, char *of5, char *of7, unsigned long long int reads_num, unsigned long int reads_num_unique){
    fprintf(stderr, "Test stats and Wig file %s. \n", of1);
    FILE *ftest = fopen("asample.txt", "w");
    fprintf(stderr, "Test example.stat \n");
    FILE *f1;
    f1  = mustOpen(of1, "w");
    FILE *f2 = mustOpen(of2, "w"); // outWig_ALL_+
    FILE *f5 = mustOpen(of5, "w"); // outWig_Uniq_+
    FILE *f6 = mustOpen(of6, "w"); // outWig_ALL_-
    FILE *f7 = mustOpen(of7, "w"); // outWig_Uniq_-
	
    unsigned int m;
    struct hashEl *helr;
    struct hashCookie cookier = hashFirst(hash);
    // Write the stat file
    fprintf(f1, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#subfamily", "family", "class", "consensus_length", "reads_count", "unique_reads_count", "total_length", "genome_count", "all_reads_RPKM", "all_reads_RPM", "unique_reads_RPKM", "unique_reads_RPM");
    while ( (helr = hashNext(&cookier)) != NULL ) {
        struct rep *or = (struct rep *) (helr->val);
        fprintf(f1, "%s\t%s\t%s\t%u\t%llu\t%llu\t%llu\t%llu\t%.3f\t%.3f\t%.3f\t%.3f\n", or->name, or->fname, or->cname, or->length, or->read_count, or->read_count_unique, or->total_length, or->genome_count, cal_rpkm(or->read_count, or->total_length, reads_num), cal_rpm(or->read_count, reads_num), cal_rpkm(or->read_count_unique, or->total_length, reads_num_unique), cal_rpm(or->read_count_unique, reads_num_unique));
        // Write the wig file
        if (or->length != 0){
            fprintf(f6, "fixedStep chrom=%s start=1 step=1 span=1\n", or->name); // outWig_ALL_-
            fprintf(f7, "fixedStep chrom=%s start=1 step=1 span=1\n", or->name); // outWig_Uniq_-
            for (m = 0; m < or->length; m++){
                fprintf(f6, "%u\n", (or->bp_total_minus)[m]);
                fprintf(f7, "%u\n", (or->bp_total_unique_minus)[m]);
            }

            fprintf(f2, "fixedStep chrom=%s start=1 step=1 span=1\n", or->name); // outWig_ALL_+
            fprintf(f5, "fixedStep chrom=%s start=1 step=1 span=1\n", or->name); // outWig_Uniq_+
            for (m = 0; m < or->length; m++){
                fprintf(f2, "%u\n", (or->bp_total_plus)[m]);
                fprintf(f5, "%u\n", (or->bp_total_unique_plus)[m]);
            }
        }
    }
    carefulClose(&f2);
    carefulClose(&f1);
    carefulClose(&f5);

    FILE *f3 = mustOpen(of3, "w");
    struct hashEl *hel3;
    struct hashCookie cookier3 = hashFirst(hash1);
    fprintf(f3, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#family", "class", "reads_count", "unique_reads_count", "total_length", "genome_count", "all_reads_RPKM", "all_reads_RPM", "unique_reads_RPKM", "unique_reads_RPM");
    while ( (hel3 = hashNext(&cookier3)) != NULL) {
        struct repfam *or = (struct repfam *) hel3->val;
        fprintf(f3, "%s\t%s\t%llu\t%llu\t%llu\t%llu\t%.3f\t%.3f\t%.3f\t%.3f\n", or->fname, or->cname, or->read_count, or->read_count_unique, or->total_length, or->genome_count, cal_rpkm(or->read_count, or->total_length, reads_num), cal_rpm(or->read_count, reads_num), cal_rpkm(or->read_count_unique, or->total_length, reads_num_unique), cal_rpm(or->read_count_unique, reads_num_unique));
    }
    carefulClose(&f3);
    FILE *f4 = mustOpen(of4, "w");
    struct hashEl *hel4;
    struct hashCookie cookier4 = hashFirst(hash2);
    fprintf(f4, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",  "#class", "reads_count", "unique_reads_count", "total_length", "genome_count", "all_reads_RPKM", "all_reads_RPM", "unique_reads_RPKM", "unique_reads_RPM");
    while ( (hel4 = hashNext(&cookier4)) != NULL) {
        struct repcla *or = (struct repcla *) hel4->val;
        fprintf(f4, "%s\t%llu\t%llu\t%llu\t%llu\t%.3f\t%.3f\t%.3f\t%.3f\n", or->cname, or->read_count, or->read_count_unique, or->total_length, or->genome_count, cal_rpkm(or->read_count, or->total_length, reads_num), cal_rpm(or->read_count, reads_num), cal_rpkm(or->read_count_unique, or->total_length, reads_num_unique), cal_rpm(or->read_count_unique, reads_num_unique));
    }
    carefulClose(&f4);
}

void MREwriteWigandStat(struct hash *hash, struct hash *hash1, struct hash *hash2, char *of1, char *of2, char *of3, char *of4){
    FILE *f1 = mustOpen(of1, "w");
    FILE *f2 = mustOpen(of2, "w");
    unsigned int m;
    struct hashEl *helr;
    struct hashCookie cookier = hashFirst(hash);
    fprintf(f1, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#subfamily", "family", "class", "consensus_length", "covered_CpG_sites", "CpG_total_score", "total_length", "genome_count");
    while ( (helr = hashNext(&cookier)) != NULL ) {
        struct rep *or = (struct rep *) (helr->val);
        fprintf(f1, "%s\t%s\t%s\t%u\t%u\t%.4f\t%llu\t%llu\n", or->name, or->fname, or->cname, or->length, or->cpgCount, or->cpgTotalScore, or->total_length, or->genome_count);
        if (or->length != 0){
            fprintf(f2, "fixedStep chrom=%s start=1 step=1 span=1\n", or->name);
            for (m = 0; m < or->length; m++){ 
                fprintf(f2, "%.4f\n", (or->cpgScore)[m]);
            }
        }
    }
    carefulClose(&f2);
    carefulClose(&f1);
    FILE *f3 = mustOpen(of3, "w");
    struct hashEl *hel3;
    struct hashCookie cookier3 = hashFirst(hash1);
    fprintf(f3, "%s\t%s\t%s\t%s\t%s\t%s\n", "#family", "class", "covered_CpG_sites", "CpG_total_score", "total_length", "genome_count");
    while ( (hel3 = hashNext(&cookier3)) != NULL) {
        struct repfam *or = (struct repfam *) hel3->val;
        fprintf(f3, "%s\t%s\t%u\t%.4f\t%llu\t%llu\n", or->fname, or->cname, or->cpgCount, or->cpgTotalScore, or->total_length, or->genome_count);
    }
    carefulClose(&f3);
    FILE *f4 = mustOpen(of4, "w");
    struct hashEl *hel4;
    struct hashCookie cookier4 = hashFirst(hash2);
    fprintf(f4, "%s\t%s\t%s\t%s\t%s\n", "#class", "covered_CpG_sites", "CpG_total_score", "total_length", "genome_count");
    while ( (hel4 = hashNext(&cookier4)) != NULL) {
        struct repcla *or = (struct repcla *) hel4->val;
        fprintf(f4, "%s\t%u\t%.4f\t%llu\t%llu\n", or->cname, or->cpgCount, or->cpgTotalScore, or->total_length, or->genome_count);
    }
    carefulClose(&f4);
}

unsigned long long int *samFile2nodupRepbedFile(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr) {
    samfile_t *samfp;
    char chr[100], prn[500], key[100], strand;
    unsigned int start, end, cend, rstart, rend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 5);
    unsigned long long int mapped_reads_num = 0, reads_num = 0, reads_used = 0, unique_reads = 0, repeat_reads = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
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
    strcpy(prn, "empty");
    bam1_t *b;
    bam_header_t *h;
    h = samfp->header;
    b = bam_init1();
    while ( samread(samfp, b) >= 0) {
        //if ( sameString (bam1_qname(b), prn)) 
        //    continue;
        reads_num++;
        if ((reads_num % 10000) == 0)
            fprintf(stderr, "\r* Processed reads: %llu", reads_num);
        //strcpy(prn, bam1_qname(b));
        //if (b->core.tid < 0)
        if (b->core.flag & BAM_FUNMAP)
            continue;
        if (b->core.qual < mapQ)
            continue;
        mapped_reads_num++;
        //change chr name to chr1, chr2 ...
        strcpy(chr, h->target_name[b->core.tid]);
        if (addChr){
            if (startsWith("GL", h->target_name[b->core.tid])) {
                continue;
            } else if (sameWord(h->target_name[b->core.tid], "MT")) {
                strcpy(chr,"chrM");
            } else if (!startsWith("chr", h->target_name[b->core.tid])) {
                strcpy(chr, "chr");
                strcat(chr, h->target_name[b->core.tid]);
            }
        }
        struct hashEl *he = hashLookup(nochr, chr);
        if (he != NULL)
            continue;
        cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
        if (cend == 1){
            hashAddInt(nochr, chr, 1);
            warn("* Warning: reads mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
            continue;
        }
        reads_used++;
        start = (unsigned int) b->core.pos;
        int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
        end = min(cend, (unsigned int)tmpend);
        strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
        //remove dup first
        if (rmDup == 1){
            if (sprintf(key, "%s:%u:%u:%c", chr, start, end, strand) < 0)
                errAbort("Mem ERROR");
            struct hashEl *hel = hashLookup(dup, key);
            if (hel == NULL) {
                hashAddInt(dup, key, 1);
            } else {
                continue;
            }
        }
        unique_reads++;
        //transfer coordinates
        int i, j;
        unsigned int qlen = end - start;
        struct binElement *hitList = NULL, *hit;
        struct hashEl *hel2 = hashLookup(hashRmsk, chr);
        if (hel2 != NULL) {
            struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
            hitList = binKeeperFind(bs2, start, end);
            if(hitList != NULL) {
                for (hit = hitList; hit !=NULL; hit = hit->next) {
                    struct rmsk *ss = (struct rmsk *) hit->val;
                    if (filter == 0){
                        struct hashEl *hel3 = hashLookup(hashRep, ss->name);
                        if (hel3 != NULL){
                            struct rep *rs = (struct rep *) hel3->val;
                            rs->read_count++;
                            if (rs->length != 0){
                                rstart = start - ss->start;
                                rstart = (rstart < 0) ? 0 : rstart;
                                rend = rstart + qlen;
                                rend = (rend < ss->end) ? rend : ss->end;
                                for (i = rstart; i < rend; i++) {
                                    j = i + ss->consensus_start;
                                    if (j >= ss->consensus_end) {
                                        break;
                                    }
                                    if (j >= rs->length) {
                                        break;
                                    }
                                    (rs->bp_total)[j]++;
                                }
                            }
                        }
                        //fill hashFam
                        struct hashEl *hel4 = hashLookup(hashFam, ss->fname);
                        if (hel4 != NULL) {
                            struct repfam *fs = (struct repfam *) hel4->val;
                            fs->read_count++;
                        }
                        //fill hashCla
                        struct hashEl *hel5 = hashLookup(hashCla, ss->cname);
                        if (hel5 != NULL) {
                            struct repcla *cs = (struct repcla *) hel5->val;
                            cs->read_count++;
                        }
                    } else {
                        slNameAddHead(&(ss->sl), bam1_qname(b));
                    }
                    break;
                }
                repeat_reads++;
                slFreeList(hitList);
            }
        }
    }
    fprintf(stderr, "\r* Processed reads: %llu\n", reads_num);
    samclose(samfp);
    bam_destroy1(b);
    freeHash(&nochr);
    freeHash(&dup);
    cnt[0] = reads_num;
    cnt[1] = mapped_reads_num;
    cnt[2] = reads_used;
    cnt[3] = unique_reads;
    cnt[4] = repeat_reads;
    return cnt;
}

float getCov(unsigned int aStart, unsigned int aEnd, unsigned int start, unsigned int end){
    float overlap = positiveRangeIntersection((int)aStart, (int)aEnd, (int)start, (int)end);
    float denominator = (float)(aEnd - aStart);
    float cov = (denominator == 0) ? 0.0 : overlap/denominator;
    return cov;
}

int mapped2diffSubfam(struct hash *hashRmsk, char *subfam, int nm, char *ahstring, int qlen){
    char *row[100], *row2[4];
    int i, nm2, start, end, num2;
    struct binElement *hitList = NULL, *hit;
    //an example string
    //XA:Z:chr9,-69070599,36M,1;chr20,-29616939,36M,1;chr9,+68450226,36M,2;chrUn_gl000219,+94584,36M,2;chrUn_gl000241,+31782,36M,2;
    //fprintf(stderr, "row %s\n", ahstring);
    int numFields = chopByChar(ahstring, ';', row, ArraySize(row));
    for(i=0; i<numFields; i++){
        if (strlen(row[i]) > 0){
            //fprintf(stderr, "row[i] %s\n", row[i]);
            num2 = chopByChar(row[i], ',', row2, ArraySize(row2));
            //if (num2 != 4){
            //    fprintf(stderr, "num2 %i\n", num2);
            //    exit(2);
            //}
            assert(num2 == 4);
            nm2 = (int) strtol(row2[3], 0, 0);
            if (nm2 <= nm){ // hmm..., FIXME?
                start = abs((int)strtol(row2[1], 0, 0));
                end = start + qlen; //FIXME
                struct hashEl *hel2 = hashLookup(hashRmsk,row2[0]);
                if (hel2 != NULL) {
                    struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
                    hitList = binKeeperFind(bs2, start, end);
                    if(hitList != NULL) {
                        for (hit = hitList; hit !=NULL; hit = hit->next) {
                            struct rmsk *ss = (struct rmsk *) hit->val;
                            if (!sameWord(ss->name, subfam)){
                                return 1;
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}

//samFile, no s, this function is for the filter.
unsigned long long int *samFile2nodupRepbedFileNew(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, float minCoverage, int treat, char *outbed, char *outbed_unique, int diffSubfam) {
    samfile_t *samfp;
    FILE *outbed_f = NULL, *outbed_unique_f = NULL;
    char chr[100], key[100], strand, ahstring[2000];
    int nm;
    unsigned int start, end, cend, rstart, rend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 13);
    //unsigned long long int mapped_reads_num = 0, reads_num = 0, reads_used = 0, unique_reads = 0, repeat_reads = 0;
    unsigned long long int read_end1 = 0, read_end2 = 0;
    unsigned long long int read_end1_mapped = 0, read_end2_mapped = 0;
    unsigned long long int read_end1_used = 0, read_end2_used = 0;
    unsigned long long int reads_nonredundant = 0;
    unsigned long long int reads_nonredundant_unique = 0;
    unsigned long long int reads_mapped = 0;
    unsigned long long int reads_mapped_unique = 0;
    unsigned long long int reads_repeat = 0;
    unsigned long long int reads_repeat_unique = 0;
    unsigned long long int reads_diff_subfam = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
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
    //strcpy(prn, "empty");
    if (outbed != NULL)
        outbed_f = mustOpen(outbed, "w");
    if (outbed_unique != NULL)
        outbed_unique_f = mustOpen(outbed_unique, "w");
    bam1_t *b;
    bam_header_t *h;
    h = samfp->header;
    b = bam_init1();
    while ( samread(samfp, b) >= 0) {
        //if ( sameString (bam1_qname(b), prn)) 
        //    continue;
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1++;
            }else{
                if(treat)
                    read_end1++;
                else
                    read_end2++;
            }
        }else{
            read_end1++;
        }
        if (((read_end1 + read_end2) % 10000) == 0)
            fprintf(stderr, "\r* Processed read ends: %llu", (read_end1 + read_end2));
        //strcpy(prn, bam1_qname(b));
        //if (b->core.tid < 0)
        if (b->core.flag & BAM_FUNMAP)
            continue;
        //if (b->core.qual < mapQ)
        //    continue;
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_mapped++;
            }else{
                if (treat)
                    read_end1_mapped++;
                else
                    read_end2_mapped++;
            }
        }else{
            read_end1_mapped++;
        }
        //change chr name to chr1, chr2 ...
        strcpy(chr, h->target_name[b->core.tid]);
        if (addChr){
            if (startsWith("GL", h->target_name[b->core.tid])) {
                continue;
            } else if (sameWord(h->target_name[b->core.tid], "MT")) {
                strcpy(chr,"chrM");
            } else if (!startsWith("chr", h->target_name[b->core.tid])) {
                strcpy(chr, "chr");
                strcat(chr, h->target_name[b->core.tid]);
            }
        }
        //check Ref reads mapped to existed in chromosome size file or not
        struct hashEl *he = hashLookup(nochr, chr);
        if (he != NULL)
            continue;
        cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
        if (cend == 1){
            hashAddInt(nochr, chr, 1);
            warn("* Warning: read ends mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
            continue;
        }
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_used++;
            }else{
                if (treat)
                    read_end1_used++;
                else
                    read_end2_used++;
            }
        }else{
            read_end1_used++;
        }
        //get mapping location for paired-end or single-end
        if (treat){
            reads_mapped++;
            if (b->core.qual >= mapQ)
                reads_mapped_unique++;
            start = (unsigned int) b->core.pos;
            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
            end = min(cend, (unsigned int)tmpend);
            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
            if (extension) {
                if (strand == '+'){
                    end = min(start + extension, cend);
                }else{
                    if (end < extension)
                        start = 0;
                    else
                        start = end - extension;
                    //start = max(end - extension, 0);
                }
            }

        }else{
        if (b->core.flag & BAM_FPAIRED) {
            if (!(b->core.flag & BAM_FMUNMAP)){
                if (b->core.flag & BAM_FREAD1){
                    if (abs(b->core.isize) > iSize || b->core.isize == 0){
                        continue;
                    }else{
                        reads_mapped++;
                        if (b->core.qual >= mapQ)
                            reads_mapped_unique++;
                        if (b->core.isize > 0){
                            start = (unsigned int) b->core.pos;
                            strand = '+';
                            int tmpend = start + b->core.isize;
                            end = min(cend, (unsigned int)tmpend);
                        }else{
                            start = (unsigned int) b->core.mpos;
                            strand = '-';
                            int tmpend = start - b->core.isize;
                            end = min(cend, (unsigned int)tmpend);
                        }
                    }
                }else{
                    continue;
                }
            }else{
                if (discardWrongEnd){
                    continue;
                }else{
                    reads_mapped++;
                    if (b->core.qual >= mapQ)
                        reads_mapped_unique++;
                    start = (unsigned int) b->core.pos;
                    int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
                    end = min(cend, (unsigned int)tmpend);
                    strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
                    if (extension) {
                        if (strand == '+'){
                            end = min(start + extension, cend);
                        }else{
                            if (end < extension)
                                start = 0;
                            else
                                start = end - extension;
                            //start = max(end - extension, 0);
                        }
                    }
                }
            }
        }else{
            reads_mapped++;
            if (b->core.qual >= mapQ)
                reads_mapped_unique++;
            start = (unsigned int) b->core.pos;
            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
            end = min(cend, (unsigned int)tmpend);
            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
            if (extension) {
                if (strand == '+'){
                    end = min(start + extension, cend);
                }else{
                    if (end < extension)
                        start = 0;
                    else
                        start = end - extension;
                    //start = max(end - extension, 0);
                }
            }
        }
    }
        //remove dup first
        if (rmDup){
            //redundant only useful for unique reads
            if (b->core.qual >= mapQ){
                if (sprintf(key, "%s:%u:%u:%c", chr, start, end, strand) < 0)
                    errAbort("Mem ERROR");
            }
            struct hashEl *hel = hashLookup(dup, key);
            if (hel == NULL) {
                hashAddInt(dup, key, 1);
            } else {
                continue;
            }
        }
        //reads_nonredundant++;
        if (b->core.qual >= mapQ)
            reads_nonredundant_unique++;

        //output bed
        if (outbed_f){
            fprintf(outbed_f, "%s\t%u\t%u\t%s\t%i\t%c", chr, start, end, bam1_qname(b), b->core.qual, strand);
            if(bam_aux_get(b, "XA")){
                fprintf(outbed_f, "\t%i\t%s", bam_aux2i(bam_aux_get(b, "NM")), bam_aux2Z(bam_aux_get(b, "XA")) );
            }
            fprintf(outbed_f, "\n");
        }
        if (outbed_unique_f){
            if(b->core.qual >= mapQ){
                fprintf(outbed_unique_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, end, bam1_qname(b), b->core.qual, strand);
            }
        }

        //transfer coordinates
        int i, j;
        int index = 0, tindex = 0;
        float coverage = 0.0, tcoverage = 0.0;
        unsigned int qlen = end - start;
        struct rmsk *ss = NULL;
        struct binElement *hitList = NULL, *hit;
        struct hashEl *hel2 = hashLookup(hashRmsk, chr);
        if (hel2 != NULL) {
            struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
            hitList = binKeeperFind(bs2, start, end);
            if(hitList != NULL) {
                for (hit = hitList; hit !=NULL; hit = hit->next) {
                    index++;
                    struct rmsk *sss = (struct rmsk *) hit->val;
                    float cov = getCov(start, end, sss->start, sss->end);
                    //fprintf(stderr, "coverage: %.2f\n", cov);
                    if (cov > coverage){
                        tindex = index;
                        tcoverage = cov;
                    }
                    coverage = cov;
                }
                if (tcoverage < minCoverage)
                    continue;
                index = 0;
                for (hit = hitList; hit !=NULL; hit = hit->next) {
                    index++;
                    if (index == tindex){
                        ss = (struct rmsk *) hit->val;
                        break;
                    }
                }
                //filter reads mapped to differetn subfam
                if (diffSubfam){
                    //filter reads mapped to different subfamiles with same NM
                    if(bam_aux_get(b, "XA")){
                        strcpy(ahstring, bam_aux2Z(bam_aux_get(b, "XA")));
                        nm = bam_aux2i(bam_aux_get(b, "NM"));
                        if (mapped2diffSubfam(hashRmsk, ss->name, nm, ahstring, (int)qlen)){
                            reads_diff_subfam++;
                            continue;
                        }
                    }
                }
                if (filter == 0){
                    struct hashEl *hel3 = hashLookup(hashRep, ss->name);
                    if (hel3 != NULL){
                        struct rep *rs = (struct rep *) hel3->val;
                        rs->read_count++;
                        if (b->core.qual >= mapQ)
                            rs->read_count_unique++;
                        if (rs->length != 0){
                            rstart = start - ss->start;
                            rstart = (rstart < 0) ? 0 : rstart;
                            rend = rstart + qlen;
                            rend = (rend < ss->end) ? rend : ss->end;
                            for (i = rstart; i < rend; i++) {
                                j = i + ss->consensus_start;
                                if (j >= ss->consensus_end) {
                                    break;
                                }
                                if (j >= rs->length) {
                                    break;
                                }
                                (rs->bp_total)[j]++;
                                if (b->core.qual >= mapQ)
                                    (rs->bp_total_unique)[j]++;
                            }
                        }
                    }
                    //fill hashFam
                    struct hashEl *hel4 = hashLookup(hashFam, ss->fname);
                    if (hel4 != NULL) {
                        struct repfam *fs = (struct repfam *) hel4->val;
                        fs->read_count++;
                        if (b->core.qual >= mapQ)
                            fs->read_count_unique++;
                    }
                    //fill hashCla
                    struct hashEl *hel5 = hashLookup(hashCla, ss->cname);
                    if (hel5 != NULL) {
                        struct repcla *cs = (struct repcla *) hel5->val;
                        cs->read_count++;
                        if (b->core.qual >= mapQ)
                            cs->read_count_unique++;
                    }
                } else {
                    slNameAddHead(&(ss->sl), bam1_qname(b));
                    if (b->core.qual >= mapQ)
                        slNameAddHead(&(ss->sl_unique), bam1_qname(b));
                }
                reads_repeat++;
                if (b->core.qual >= mapQ)
                    reads_repeat_unique++;
                slFreeList(hitList);
            }
        }
    }
    fprintf(stderr, "\r* Processed read ends: %llu\n", (read_end1 + read_end2));
    samclose(samfp);
    bam_destroy1(b);
    freeHash(&nochr);
    freeHash(&dup);
    if (outbed_f)
        carefulClose(&outbed_f);
    if (outbed_unique_f)
        carefulClose(&outbed_unique_f);
    cnt[0] = read_end1;
    cnt[1] = read_end2;
    cnt[2] = read_end1_mapped;
    cnt[3] = read_end2_mapped;
    cnt[4] = read_end1_used;
    cnt[5] = read_end2_used;
    cnt[6] = reads_mapped;
    cnt[7] = reads_mapped_unique;
    cnt[8] = reads_nonredundant;
    cnt[9] = reads_repeat;
    cnt[10] = reads_repeat_unique;
    cnt[11] = reads_nonredundant_unique;
    cnt[12] = reads_diff_subfam;
    return cnt;
}

//This is Cage-specific for the filter.
unsigned long long int *samFile2nodupRepbedFileNewCage(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, float minCoverage, int treat, char *outbed, char *outbed_unique, int diffSubfam, int optcagewindow) {
    samfile_t *samfp;
    FILE *outbed_f = NULL, *outbed_unique_f = NULL;
    char chr[100], key[100], strand, ahstring[2000];
    int nm;
    unsigned int start, end, cend, rstart, rend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 13);
    //unsigned long long int mapped_reads_num = 0, reads_num = 0, reads_used = 0, unique_reads = 0, repeat_reads = 0;
    unsigned long long int read_end1 = 0, read_end2 = 0;
    unsigned long long int read_end1_mapped = 0, read_end2_mapped = 0;
    unsigned long long int read_end1_used = 0, read_end2_used = 0;
    unsigned long long int reads_nonredundant = 0;
    unsigned long long int reads_nonredundant_unique = 0;
    unsigned long long int reads_mapped = 0;
    unsigned long long int reads_mapped_unique = 0;
    unsigned long long int reads_repeat = 0;
    unsigned long long int reads_repeat_unique = 0;
    unsigned long long int reads_diff_subfam = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
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
    //strcpy(prn, "empty");
    if (outbed != NULL)
        outbed_f = mustOpen(outbed, "w");
    if (outbed_unique != NULL)
        outbed_unique_f = mustOpen(outbed_unique, "w");
    bam1_t *b;
    bam_header_t *h;
    h = samfp->header;
    b = bam_init1();
    while ( samread(samfp, b) >= 0) {
        //if ( sameString (bam1_qname(b), prn))
        //    continue;
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1++;
            }else{
                if(treat)
                    read_end1++;
                else
                    read_end2++;
            }
        }else{
            read_end1++;
        }
        if (((read_end1 + read_end2) % 10000) == 0)
            fprintf(stderr, "\r* Processed read ends: %llu", (read_end1 + read_end2));
        //strcpy(prn, bam1_qname(b));
        //if (b->core.tid < 0)
        if (b->core.flag & BAM_FUNMAP)
            continue;
        //if (b->core.qual < mapQ)
        //    continue;
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_mapped++;
            }else{
                if (treat)
                    read_end1_mapped++;
                else
                    read_end2_mapped++;
            }
        }else{
            read_end1_mapped++;
        }
        //change chr name to chr1, chr2 ...
        strcpy(chr, h->target_name[b->core.tid]);
        if (addChr){
            if (startsWith("GL", h->target_name[b->core.tid])) {
                continue;
            } else if (sameWord(h->target_name[b->core.tid], "MT")) {
                strcpy(chr,"chrM");
            } else if (!startsWith("chr", h->target_name[b->core.tid])) {
                strcpy(chr, "chr");
                strcat(chr, h->target_name[b->core.tid]);
            }
        }
        //check Ref reads mapped to existed in chromosome size file or not
        struct hashEl *he = hashLookup(nochr, chr);
        if (he != NULL)
            continue;
        cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
        if (cend == 1){
            hashAddInt(nochr, chr, 1);
            warn("* Warning: read ends mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
            continue;
        }
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_used++;
            }else{
                if (treat)
                    read_end1_used++;
                else
                    read_end2_used++;
            }
        }else{
            read_end1_used++;
        }
        //get mapping location for paired-end or single-end
//        if (treat){
//            reads_mapped++;
//            if (b->core.qual >= mapQ)
//                reads_mapped_unique++;
//            start = (unsigned int) b->core.pos;
//            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
//            end = min(cend, (unsigned int)tmpend);
//            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
//            if (extension) {
//                if (strand == '+'){
//                    end = min(start + extension, cend);
//                }else{
//                    if (end < extension)
//                        start = 0;
//                    else
//                        start = end - extension;
//                    //start = max(end - extension, 0);
//                }
//            }
//
//        }else{
//            if (b->core.flag & BAM_FPAIRED) {
//                if (!(b->core.flag & BAM_FMUNMAP)){
//                    if (b->core.flag & BAM_FREAD1){
//                        if (abs(b->core.isize) > iSize || b->core.isize == 0){
//                            continue;
//                        }else{
//                            reads_mapped++;
//                            if (b->core.qual >= mapQ)
//                                reads_mapped_unique++;
//                            if (b->core.isize > 0){
//                                start = (unsigned int) b->core.pos;
//                                strand = '+';
//                                int tmpend = start + b->core.isize;
//                                end = min(cend, (unsigned int)tmpend);
//                            }else{
//                                start = (unsigned int) b->core.mpos;
//                                strand = '-';
//                                int tmpend = start - b->core.isize;
//                                end = min(cend, (unsigned int)tmpend);
//                            }
//                        }
//                    }else{
//                        continue;
//                    }
//                }else{
//                    if (discardWrongEnd){
//                        continue;
//                    }else{
//                        reads_mapped++;
//                        if (b->core.qual >= mapQ)
//                            reads_mapped_unique++;
//                        start = (unsigned int) b->core.pos;
//                        int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
//                        end = min(cend, (unsigned int)tmpend);
//                        strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
//                        if (extension) {
//                            if (strand == '+'){
//                                end = min(start + extension, cend);
//                            }else{
//                                if (end < extension)
//                                    start = 0;
//                                else
//                                    start = end - extension;
//                                //start = max(end - extension, 0);
//                            }
//                        }
//                    }
//                }
//            }else{
//                reads_mapped++;
//                if (b->core.qual >= mapQ)
//                    reads_mapped_unique++;
//                start = (unsigned int) b->core.pos;
//                int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
//                end = min(cend, (unsigned int)tmpend);
//                strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
//                if (extension) {
//                    if (strand == '+'){
//                        end = min(start + extension, cend);
//                    }else{
//                        if (end < extension)
//                            start = 0;
//                        else
//                            start = end - extension;
//                        //start = max(end - extension, 0);
//                    }
//                }
//            }
//        }

        //Cage-seq pipeline, extract the 5' end.
        if(optcagewindow >= 0){
            reads_mapped++;
            if (b->core.qual >= mapQ)
                reads_mapped_unique++;
            start = (unsigned int) b->core.pos;
            end = start + 1 + optcagewindow > cend ? cend : start + 1 + optcagewindow;
            if(optcagewindow == 0){
                end = start + 1;
            } else{
                start = optcagewindow < start ? start - optcagewindow : 0;
                end = start + 1 + optcagewindow > cend ? cend : start + 1 + optcagewindow;
            }
            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
//            if (strand == '+'){
//                //head = searchAndUpdate(head, b->core.qual);
//                printf("The strand is +\n");
//            } else{
//                printf("The strand is -\n");
//            }
        }

        //remove dup first
        if (rmDup){
            //redundant only useful for unique reads
            if (b->core.qual >= mapQ){
                if (sprintf(key, "%s:%u:%u:%c", chr, start, end, strand) < 0)
                    errAbort("Mem ERROR");
            }
            struct hashEl *hel = hashLookup(dup, key);
            if (hel == NULL) {
                hashAddInt(dup, key, 1);
            } else {
                continue;
            }
        }
        //reads_nonredundant++;
        if (b->core.qual >= mapQ)
            reads_nonredundant_unique++;

        //output bed
        if (outbed_f){
            fprintf(outbed_f, "%s\t%u\t%u\t%s\t%i\t%c", chr, start, end, bam1_qname(b), b->core.qual, strand);
            if(bam_aux_get(b, "XA")){
                fprintf(outbed_f, "\t%i\t%s", bam_aux2i(bam_aux_get(b, "NM")), bam_aux2Z(bam_aux_get(b, "XA")) );
            }
            fprintf(outbed_f, "\n");
        }
        if (outbed_unique_f){
            if(b->core.qual >= mapQ){
                fprintf(outbed_unique_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, end, bam1_qname(b), b->core.qual, strand);
            }
        }

        //transfer coordinates
        int i, j;
        int index = 0, tindex = 0;
        float coverage = 0.0, tcoverage = 0.0;
        unsigned int qlen = end - start;
        struct rmsk *ss = NULL;
        struct binElement *hitList = NULL, *hit;
        struct hashEl *hel2 = hashLookup(hashRmsk, chr);
        if (hel2 != NULL) {
            struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
            hitList = binKeeperFind(bs2, start, end);
            if(hitList != NULL) {
                for (hit = hitList; hit !=NULL; hit = hit->next) {
                    index++;
                    struct rmsk *sss = (struct rmsk *) hit->val;
                    float cov = getCov(start, end, sss->start, sss->end);
                    //fprintf(stderr, "coverage: %.2f\n", cov);
                    if (cov > coverage){
                        tindex = index;
                        tcoverage = cov;
                    }
                    coverage = cov;
                }
                if (tcoverage < minCoverage)
                    continue;
                index = 0;
                for (hit = hitList; hit !=NULL; hit = hit->next) {
                    index++;
                    if (index == tindex){
                        ss = (struct rmsk *) hit->val;
                        break;
                    }
                }
                //filter reads mapped to differetn subfam
                if (diffSubfam){
                    //filter reads mapped to different subfamiles with same NM
                    if(bam_aux_get(b, "XA")){
                        strcpy(ahstring, bam_aux2Z(bam_aux_get(b, "XA")));
                        nm = bam_aux2i(bam_aux_get(b, "NM"));
                        if (mapped2diffSubfam(hashRmsk, ss->name, nm, ahstring, (int)qlen)){
                            reads_diff_subfam++;
                            continue;
                        }
                    }
                }
                if (filter == 0){
                    struct hashEl *hel3 = hashLookup(hashRep, ss->name);
                    if (hel3 != NULL){
                        struct rep *rs = (struct rep *) hel3->val;
                        rs->read_count++;
                        if (b->core.qual >= mapQ)
                            rs->read_count_unique++;
                        if (rs->length != 0){
                            rstart = start - ss->start;
                            rstart = (rstart < 0) ? 0 : rstart;
                            rend = rstart + qlen;
                            rend = (rend < ss->end) ? rend : ss->end;
                            for (i = rstart; i < rend; i++) {
                                j = i + ss->consensus_start;
                                if (j >= ss->consensus_end) {
                                    break;
                                }
                                if (j >= rs->length) {
                                    break;
                                }
                                (rs->bp_total)[j]++;
                                if (b->core.qual >= mapQ)
                                    (rs->bp_total_unique)[j]++;
                            }
                        }
                    }
                    //fill hashFam
                    struct hashEl *hel4 = hashLookup(hashFam, ss->fname);
                    if (hel4 != NULL) {
                        struct repfam *fs = (struct repfam *) hel4->val;
                        fs->read_count++;
                        if (b->core.qual >= mapQ)
                            fs->read_count_unique++;
                    }
                    //fill hashCla
                    struct hashEl *hel5 = hashLookup(hashCla, ss->cname);
                    if (hel5 != NULL) {
                        struct repcla *cs = (struct repcla *) hel5->val;
                        cs->read_count++;
                        if (b->core.qual >= mapQ)
                            cs->read_count_unique++;
                    }
                } else {
                    slNameAddHead(&(ss->sl), bam1_qname(b));
                    if (b->core.qual >= mapQ)
                        slNameAddHead(&(ss->sl_unique), bam1_qname(b));
                }
                reads_repeat++;
                if (b->core.qual >= mapQ)
                    reads_repeat_unique++;
                slFreeList(hitList);
            }
        }
    }
    fprintf(stderr, "\r* Processed read ends: %llu\n", (read_end1 + read_end2));
    samclose(samfp);
    bam_destroy1(b);
    freeHash(&nochr);
    freeHash(&dup);
    if (outbed_f)
        carefulClose(&outbed_f);
    if (outbed_unique_f)
        carefulClose(&outbed_unique_f);
    cnt[0] = read_end1;
    cnt[1] = read_end2;
    cnt[2] = read_end1_mapped;
    cnt[3] = read_end2_mapped;
    cnt[4] = read_end1_used;
    cnt[5] = read_end2_used;
    cnt[6] = reads_mapped;
    cnt[7] = reads_mapped_unique;
    cnt[8] = reads_nonredundant;
    cnt[9] = reads_repeat;
    cnt[10] = reads_repeat_unique;
    cnt[11] = reads_nonredundant_unique;
    cnt[12] = reads_diff_subfam;
    return cnt;
}

//support many bam files for stat
unsigned long long int *samFiles2nodupRepbedFileNew(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, float minCoverage, int treat, char *outbed, char *outbed_unique, int diffSubfam) {
    FILE *outbed_f = NULL, *outbed_unique_f = NULL;
    char chr[100], key[100], strand, ahstring[2000], *row[100];
    int nm, fi;
    unsigned int start, end, cend, rstart, rend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 13);
    //unsigned long long int mapped_reads_num = 0, reads_num = 0, reads_used = 0, unique_reads = 0, repeat_reads = 0;
    unsigned long long int read_end1 = 0, read_end2 = 0;
    unsigned long long int read_end1_mapped = 0, read_end2_mapped = 0;
    unsigned long long int read_end1_used = 0, read_end2_used = 0;
    unsigned long long int reads_nonredundant = 0;
    unsigned long long int reads_nonredundant_unique = 0;
    unsigned long long int reads_mapped = 0;
    unsigned long long int reads_mapped_unique = 0;
    unsigned long long int reads_repeat = 0;
    unsigned long long int reads_repeat_unique = 0;
    unsigned long long int reads_diff_subfam = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
    if (outbed != NULL){
        outbed_f = mustOpen(outbed, "w");
    }
    if (outbed_unique != NULL){
        outbed_unique_f = mustOpen(outbed_unique, "w");
    }
    //process bam
    int numFields = chopByChar(samfile, ',', row, ArraySize(row));
    for(fi = 0; fi < numFields; fi++){
        fprintf(stderr, "\n* Processing %s\n", row[fi]);
        samfile_t *samfp;
        if (isSam) {
            if ( (samfp = samopen(row[fi], "r", 0)) == 0) {
                fprintf(stderr, "Fail to open SAM file %s\n", samfile);
                errAbort("Error\n");
            }
        } else {
            if ( (samfp = samopen(row[fi], "rb", 0)) == 0) {
                fprintf(stderr, "Fail to open BAM file %s\n", samfile);
                errAbort("Error\n");
            }
        }
        //strcpy(prn, "empty");
        bam1_t *b;
        bam_header_t *h;
        h = samfp->header;
        b = bam_init1();
        while ( samread(samfp, b) >= 0) {
            //if ( sameString (bam1_qname(b), prn))
            //    continue;
            if (b->core.flag & BAM_FPAIRED) {
                if (b->core.flag & BAM_FREAD1){
                    read_end1++;
                }else{
                    if(treat)
                        read_end1++;
                    else
                        read_end2++;
                }
            }else{
                read_end1++;
            }
            if (((read_end1 + read_end2) % 100000) == 0)
                fprintf(stderr, "\r* Processed read ends: %llu", (read_end1 + read_end2));
            //strcpy(prn, bam1_qname(b));
            //if (b->core.tid < 0)
            if (b->core.flag & BAM_FUNMAP)
                continue;
            //if (b->core.qual < mapQ)
            //    continue;
            if (b->core.flag & BAM_FPAIRED) {
                if (b->core.flag & BAM_FREAD1){
                    read_end1_mapped++;
                }else{
                    if (treat)
                        read_end1_mapped++;
                    else
                        read_end2_mapped++;
                }
            }else{
                read_end1_mapped++;
            }
            //change chr name to chr1, chr2 ...
            strcpy(chr, h->target_name[b->core.tid]);
            if (addChr){
                if (startsWith("GL", h->target_name[b->core.tid])) {
                    continue;
                } else if (sameWord(h->target_name[b->core.tid], "MT")) {
                    strcpy(chr,"chrM");
                } else if (!startsWith("chr", h->target_name[b->core.tid])) {
                    strcpy(chr, "chr");
                    strcat(chr, h->target_name[b->core.tid]);
                }
            }
            //check Ref reads mapped to existed in chromosome size file or not
            struct hashEl *he = hashLookup(nochr, chr);
            if (he != NULL)
                continue;
            cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
            if (cend == 1){
                hashAddInt(nochr, chr, 1);
                warn("* Warning: read ends mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
                continue;
            }
            if (b->core.flag & BAM_FPAIRED) {
                if (b->core.flag & BAM_FREAD1){
                    read_end1_used++;
                }else{
                    if (treat)
                        read_end1_used++;
                    else
                        read_end2_used++;
                }
            }else{
                read_end1_used++;
            }
            //get mapping location for paired-end or single-end
            if (treat){
                reads_mapped++;
                if (b->core.qual >= mapQ)
                    reads_mapped_unique++;
                start = (unsigned int) b->core.pos;
                int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
                end = min(cend, (unsigned int)tmpend);
                strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
                if (extension) {
                    if (strand == '+'){
                        end = min(start + extension, cend);
                    }else{
                        if (end < extension)
                            start = 0;
                        else
                            start = end - extension;
                        //start = max(end - extension, 0);
                    }
                }

            }else{
                if (b->core.flag & BAM_FPAIRED) {
                    if (!(b->core.flag & BAM_FMUNMAP)){
                        if (b->core.flag & BAM_FREAD1){
                            if (abs(b->core.isize) > iSize || b->core.isize == 0){
                                continue;
                            }else{
                                reads_mapped++;
                                if (b->core.qual >= mapQ)
                                    reads_mapped_unique++;
                                if (b->core.isize > 0){
                                    start = (unsigned int) b->core.pos;
                                    strand = '+';
                                    int tmpend = start + b->core.isize;
                                    end = min(cend, (unsigned int)tmpend);
                                }else{
                                    start = (unsigned int) b->core.mpos;
                                    strand = '-';
                                    int tmpend = start - b->core.isize;
                                    end = min(cend, (unsigned int)tmpend);
                                }

                            }
                        }else{
                            continue;
                        }
                    }else{
                        if (discardWrongEnd){
                            continue;
                        }else{
                            reads_mapped++;
                            if (b->core.qual >= mapQ)
                                reads_mapped_unique++;
                            start = (unsigned int) b->core.pos;
                            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
                            end = min(cend, (unsigned int)tmpend);
                            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
                            if (extension) {
                                if (strand == '+'){
                                    end = min(start + extension, cend);
                                }else{
                                    if (end < extension)
                                        start = 0;
                                    else
                                        start = end - extension;
                                    //start = max(end - extension, 0);
                                }
                            }
                        }
                    }
                }else{
                    reads_mapped++;
                    if (b->core.qual >= mapQ)
                        reads_mapped_unique++;
                    start = (unsigned int) b->core.pos;
                    int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
                    end = min(cend, (unsigned int)tmpend);
                    strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
                    if (extension) {
                        if (strand == '+'){
                            end = min(start + extension, cend);
                        }else{
                            if (end < extension)
                                start = 0;
                            else
                                start = end - extension;
                            //start = max(end - extension, 0);
                        }
                    }
                }
            }
            //remove dup first
            if (rmDup){
                //redundant only useful for unique reads
                if (b->core.qual >= mapQ){
                    if (sprintf(key, "%s:%u:%u:%c", chr, start, end, strand) < 0)
                        errAbort("Mem ERROR");
                }
                struct hashEl *hel = hashLookup(dup, key);
                if (hel == NULL) {
                    hashAddInt(dup, key, 1);
                } else {
                    continue;
                }
            }
            //reads_nonredundant++;
            if (b->core.qual >= mapQ)
                reads_nonredundant_unique++;

            //output bed
            if (outbed_f){
                fprintf(outbed_f, "%s\t%u\t%u\t%s\t%i\t%c", chr, start, end, bam1_qname(b), b->core.qual, strand);
                if(bam_aux_get(b, "XA")){
                    fprintf(outbed_f, "\t%i\t%s", bam_aux2i(bam_aux_get(b, "NM")), bam_aux2Z(bam_aux_get(b, "XA")) );
                }
                fprintf(outbed_f, "\n");
            }
            if (outbed_unique_f){
                if(b->core.qual >= mapQ){
                    fprintf(outbed_unique_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, end, bam1_qname(b), b->core.qual, strand);
                }
            }

            //transfer coordinates
            int i, j;
            int index = 0, tindex = 0;
            float coverage = 0.0, tcoverage = 0.0;
            unsigned int qlen = end - start;
            struct rmsk *ss = NULL;
            struct binElement *hitList = NULL, *hit;
            struct hashEl *hel2 = hashLookup(hashRmsk, chr);
            if (hel2 != NULL) {
                struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
                hitList = binKeeperFind(bs2, start, end);
                if(hitList != NULL) {
                    for (hit = hitList; hit !=NULL; hit = hit->next) {
                        index++;
                        struct rmsk *sss = (struct rmsk *) hit->val;
                        float cov = getCov(start, end, sss->start, sss->end);
                        //fprintf(stderr, "coverage: %.2f\n", cov);
                        if (cov > coverage){
                            tindex = index;
                            tcoverage = cov;
                        }
                        coverage = cov;
                    }
                    if (tcoverage < minCoverage)
                        continue;
                    index = 0;
                    for (hit = hitList; hit !=NULL; hit = hit->next) {
                        index++;
                        if (index == tindex){
                            ss = (struct rmsk *) hit->val;
                            break;
                        }
                    }
                    //filter reads mapped to differetn subfam
                    if (diffSubfam){
                        //filter reads mapped to different subfamiles with same NM
                        if(bam_aux_get(b, "XA")){
                            strcpy(ahstring, bam_aux2Z(bam_aux_get(b, "XA")));
                            nm = bam_aux2i(bam_aux_get(b, "NM"));
                            if (mapped2diffSubfam(hashRmsk, ss->name, nm, ahstring, (int)qlen)){
                                reads_diff_subfam++;
                                continue;
                            }
                        }
                    }
                    if (filter == 0){
                        struct hashEl *hel3 = hashLookup(hashRep, ss->name);
                        if (hel3 != NULL){
                            struct rep *rs = (struct rep *) hel3->val;
                            rs->read_count++;
                            if (b->core.qual >= mapQ)
                                rs->read_count_unique++;
                            if (rs->length != 0){
                                rstart = start - ss->start;
                                rstart = (rstart < 0) ? 0 : rstart;
                                rend = rstart + qlen;
                                rend = (rend < ss->end) ? rend : ss->end;
                                for (i = rstart; i < rend; i++) {
                                    j = i + ss->consensus_start;
                                    if (j >= ss->consensus_end) {
                                        break;
                                    }
                                    if (j >= rs->length) {
                                        break;
                                    }
                                    (rs->bp_total)[j]++;
                                    if (b->core.qual >= mapQ)
                                        (rs->bp_total_unique)[j]++;
                                }
                            }
                        }
                        //fill hashFam
                        struct hashEl *hel4 = hashLookup(hashFam, ss->fname);
                        if (hel4 != NULL) {
                            struct repfam *fs = (struct repfam *) hel4->val;
                            fs->read_count++;
                            if (b->core.qual >= mapQ)
                                fs->read_count_unique++;
                        }
                        //fill hashCla
                        struct hashEl *hel5 = hashLookup(hashCla, ss->cname);
                        if (hel5 != NULL) {
                            struct repcla *cs = (struct repcla *) hel5->val;
                            cs->read_count++;
                            if (b->core.qual >= mapQ)
                                cs->read_count_unique++;
                        }
                    } else {
                        slNameAddHead(&(ss->sl), bam1_qname(b));
                        if (b->core.qual >= mapQ)
                            slNameAddHead(&(ss->sl_unique), bam1_qname(b));
                    }
                    reads_repeat++;
                    if (b->core.qual >= mapQ)
                        reads_repeat_unique++;
                    slFreeList(hitList);
                }
            }
        }
        fprintf(stderr, "\r* Processed read ends: %llu\n", (read_end1 + read_end2));
        samclose(samfp);
        bam_destroy1(b);
    }
    //process bam ends
    freeHash(&nochr);
    freeHash(&dup);
    if (outbed_f)
        carefulClose(&outbed_f);
    if (outbed_unique_f)
        carefulClose(&outbed_unique_f);
    cnt[0] = read_end1;
    cnt[1] = read_end2;
    cnt[2] = read_end1_mapped;
    cnt[3] = read_end2_mapped;
    cnt[4] = read_end1_used;
    cnt[5] = read_end2_used;
    cnt[6] = reads_mapped;
    cnt[7] = reads_mapped_unique;
    cnt[8] = reads_nonredundant;
    cnt[9] = reads_repeat;
    cnt[10] = reads_repeat_unique;
    cnt[11] = reads_nonredundant_unique;
    cnt[12] = reads_diff_subfam;
    return cnt;
}

// Cage-seq specific function for stat
unsigned long long int *samFiles2nodupRepbedFileNewCage(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, float minCoverage, int treat, char *outbed, char *outbed_unique, int diffSubfam, int optcagewindow) {
    FILE *outbed_f = NULL, *outbed_unique_f = NULL;
    char chr[100], key[100], strand, ahstring[2000], *row[100];
    int nm, fi;
    unsigned int start, end, cend, rstart, rend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 17);
    //unsigned long long int mapped_reads_num = 0, reads_num = 0, reads_used = 0, unique_reads = 0, repeat_reads = 0;
    unsigned long long int read_end1 = 0, read_end2 = 0;
    unsigned long long int read_end1_mapped = 0, read_end2_mapped = 0;
    unsigned long long int read_end1_used = 0, read_end2_used = 0;
    unsigned long long int reads_nonredundant = 0;
    unsigned long long int reads_nonredundant_unique = 0;
    unsigned long long int reads_mapped = 0;
    unsigned long long int reads_mapped_unique = 0;
    unsigned long long int reads_repeat = 0, reads_repeat_minus = 0, reads_repeat_plus = 0;
    unsigned long long int reads_repeat_unique = 0, reads_repeat_unique_minus = 0, reads_repeat_unique_plus = 0;

    unsigned long long int reads_diff_subfam = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
    if (outbed != NULL){
        outbed_f = mustOpen(outbed, "w");
    }
    if (outbed_unique != NULL){
        outbed_unique_f = mustOpen(outbed_unique, "w");
    }
    //process bam
    int numFields = chopByChar(samfile, ',', row, ArraySize(row));
    for(fi = 0; fi < numFields; fi++){
        fprintf(stderr, "\n* Processing %s\n", row[fi]);
        samfile_t *samfp;
        if (isSam) {
            if ( (samfp = samopen(row[fi], "r", 0)) == 0) {
                fprintf(stderr, "Fail to open SAM file %s\n", samfile);
                errAbort("Error\n");
            }
        } else {
            if ( (samfp = samopen(row[fi], "rb", 0)) == 0) {
                fprintf(stderr, "Fail to open BAM file %s\n", samfile);
                errAbort("Error\n");
            }
        }
        //strcpy(prn, "empty");
        bam1_t *b;
        bam_header_t *h;
        h = samfp->header;
        b = bam_init1();
        while (samread(samfp, b) >= 0) {
            //if ( sameString (bam1_qname(b), prn)) 
            //    continue;
            if (b->core.flag & BAM_FPAIRED) {
                if (b->core.flag & BAM_FREAD1){
                    read_end1++;
                }else{
                    if(treat)
                        read_end1++;
                    else
                        read_end2++;
                }
            }else{
                read_end1++;
            }

            if (((read_end1 + read_end2) % 100000) == 0)
                fprintf(stderr, "\r* Processed read ends: %llu", (read_end1 + read_end2));
            //strcpy(prn, bam1_qname(b));
            //if (b->core.tid < 0)
            if (b->core.flag & BAM_FUNMAP)
                continue;
            //if (b->core.qual < mapQ)
            //    continue;
            if (b->core.flag & BAM_FPAIRED) {
                if (b->core.flag & BAM_FREAD1){
                    read_end1_mapped++;
                }else{
                    if (treat)
                        read_end1_mapped++;
                    else
                        read_end2_mapped++;
                }
            }else{
                read_end1_mapped++;
            }
            //change chr name to chr1, chr2 ...
            strcpy(chr, h->target_name[b->core.tid]);
            if (addChr){
                if (startsWith("GL", h->target_name[b->core.tid])) {
                    continue;
                } else if (sameWord(h->target_name[b->core.tid], "MT")) {
                    strcpy(chr,"chrM");
                } else if (!startsWith("chr", h->target_name[b->core.tid])) {
                    strcpy(chr, "chr");
                    strcat(chr, h->target_name[b->core.tid]);
                }
            }
            //check Ref reads mapped to existed in chromosome size file or not
            struct hashEl *he = hashLookup(nochr, chr);
            if (he != NULL)
                continue;
            cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
            if (cend == 1){
                hashAddInt(nochr, chr, 1);
                warn("* Warning: read ends mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
                continue;
            }
            if (b->core.flag & BAM_FPAIRED) {
                if (b->core.flag & BAM_FREAD1){
                    read_end1_used++;
                }else{
                    if (treat)
                        read_end1_used++;
                    else
                        read_end2_used++;
                }
            }else{
                read_end1_used++;
            }

            //get mapping location for paired-end or single-end
//            if (treat){
//                reads_mapped++;
//                if (b->core.qual >= mapQ)
//                    reads_mapped_unique++;
//                //STAR output reads may be different. 255.
//                start = (unsigned int) b->core.pos;
//                int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
//                end = min(cend, (unsigned int)tmpend);
//                strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
//                if (extension) {
//                    if (strand == '+'){
//                        end = min(start + extension, cend);
//                    }else{
//                        if (end < extension)
//                            start = 0;
//                        else
//                            start = end - extension;
//                        //start = max(end - extension, 0);
//                    }
//                }
//
//            }else{
//                if (b->core.flag & BAM_FPAIRED) {
//                    if (!(b->core.flag & BAM_FMUNMAP)){
//                        if (b->core.flag & BAM_FREAD1){
//                            if (abs(b->core.isize) > iSize || b->core.isize == 0){
//                                continue;
//                            }else{
//                                reads_mapped++;
//                                if (b->core.qual >= mapQ)
//                                    reads_mapped_unique++;
//                                if (b->core.isize > 0){
//                                    start = (unsigned int) b->core.pos;
//                                    strand = '+';
//                                    int tmpend = start + b->core.isize;
//                                    end = min(cend, (unsigned int)tmpend);
//                                }else{
//                                    start = (unsigned int) b->core.mpos;
//                                    strand = '-';
//                                    int tmpend = start - b->core.isize;
//                                    end = min(cend, (unsigned int)tmpend);
//                                }
//                            }
//                        }else{
//                            continue;
//                        }
//                    }else{
//                        if (discardWrongEnd){
//                            continue;
//                        }else{
//                            reads_mapped++;
//                            if (b->core.qual >= mapQ)
//                                reads_mapped_unique++;
//                            start = (unsigned int) b->core.pos;
//                            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
//                            end = min(cend, (unsigned int)tmpend);
//                            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
//                            if (extension) {
//                                if (strand == '+'){
//                                    end = min(start + extension, cend);
//                                }else{
//                                    if (end < extension)
//                                        start = 0;
//                                    else
//                                        start = end - extension;
//                                    //start = max(end - extension, 0);
//                                }
//                            }
//                        }
//                    }
//                }else{
//                    reads_mapped++;
//                    if (b->core.qual >= mapQ)
//                        reads_mapped_unique++;
//                    start = (unsigned int) b->core.pos;
//                    int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
//                    end = min(cend, (unsigned int)tmpend);
//                    strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
//                    if (extension) {
//                        if (strand == '+'){
//                            end = min(start + extension, cend);
//                        }else{
//                            if (end < extension)
//                                start = 0;
//                            else
//                                start = end - extension;
//                                printf("end < extension, end: %d, extension:%d ", end, extension);
//                            //start = max(end - extension, 0);
//                        }
//                    }
//                }
//            }

            //Cage-seq pipeline, extract the 5' end.
            if(optcagewindow >= 0){
                reads_mapped++;
                if (b->core.qual >= mapQ)
                    reads_mapped_unique++;
                start = (unsigned int) b->core.pos;
                if(optcagewindow == 0){
                    end = start + 1;
                } else{
                    end = start + 1 + optcagewindow > cend ? cend : start + 1 + optcagewindow;
                    start = optcagewindow < start ? start - optcagewindow : 0;
                }
                strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
//                if (strand == '+'){
//                    //head = searchAndUpdate(head, b->core.qual);
//                    printf("The strand is +\n");
//                } else{
//                    printf("The strand is -\n");
//                }
            }

//            printf("mapQ: %d, b->core.qual: %d\n", mapQ, b->core.qual);
//            printf("5' loci: %d, start: %d, end: %d, optcagewindow: %d, chromosome end: %d\n",(unsigned int) b->core.pos, start, end, optcagewindow, cend);

            //remove dup first
            if (rmDup){
                //redundant only useful for unique reads
                if (b->core.qual >= mapQ){
                    if (sprintf(key, "%s:%u:%u:%c", chr, start, end, strand) < 0)
                        errAbort("Mem ERROR");
                }
                struct hashEl *hel = hashLookup(dup, key);
                if (hel == NULL) {
                    hashAddInt(dup, key, 1);
                } else {
                    continue;
                }
            }
            //reads_nonredundant++;
            if (b->core.qual >= mapQ)
                reads_nonredundant_unique++;

            //output bed
            if (outbed_f){
                fprintf(outbed_f, "%s\t%u\t%u\t%s\t%i\t%c", chr, start, end, bam1_qname(b), b->core.qual, strand);
                if(bam_aux_get(b, "XA")){
                    fprintf(outbed_f, "\t%i\t%s", bam_aux2i(bam_aux_get(b, "NM")), bam_aux2Z(bam_aux_get(b, "XA")) );
                }
                fprintf(outbed_f, "\n");
            }
            if (outbed_unique_f){
                if(b->core.qual >= mapQ){
                    fprintf(outbed_unique_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, end, bam1_qname(b), b->core.qual, strand);
                }
            }

            //transfer coordinates
            int i, j;
            int index = 0, tindex = 0;
            float coverage = 0.0, tcoverage = 0.0;
            unsigned int qlen = end - start;
            struct rmsk *ss = NULL;
            struct binElement *hitList = NULL, *hit;
            struct hashEl *hel2 = hashLookup(hashRmsk, chr);
            if (hel2 != NULL) {
                struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
                hitList = binKeeperFind(bs2, start, end);
                if(hitList != NULL) {
                    for (hit = hitList; hit !=NULL; hit = hit->next) {
                        index++;
                        struct rmsk *sss = (struct rmsk *) hit->val;
                        float cov = getCov(start, end, sss->start, sss->end);
                        //fprintf(stderr, "coverage: %.2f\n", cov);
                        if (cov > coverage){
                            tindex = index;
                            tcoverage = cov;
                        }
                        coverage = cov;
                    }
                    if (tcoverage < minCoverage)
                        continue;
                    index = 0;
                    for (hit = hitList; hit !=NULL; hit = hit->next) {
                        index++;
                        if (index == tindex){
                            ss = (struct rmsk *) hit->val;
                            break;
                        }
                    }
                    //filter reads mapped to differetn subfam
                    if (diffSubfam){
                        //filter reads mapped to different subfamiles with same NM
                        if(bam_aux_get(b, "XA")){
                            strcpy(ahstring, bam_aux2Z(bam_aux_get(b, "XA")));
                            nm = bam_aux2i(bam_aux_get(b, "NM"));
                            if (mapped2diffSubfam(hashRmsk, ss->name, nm, ahstring, (int)qlen)){
                                reads_diff_subfam++;
                                continue;
                            }
                        }
                    }
                    if (filter == 0){
                        struct hashEl *hel3 = hashLookup(hashRep, ss->name);
                        if (hel3 != NULL){
                            struct rep *rs = (struct rep *) hel3->val;
                            rs->read_count++;
//                            printf("rs->read_count: %llu\n", rs->read_count);
                            if (b->core.qual >= mapQ)
                                rs->read_count_unique++;
                            if (rs->length != 0){
                                rstart = start - ss->start;
                                rstart = (rstart < 0) ? 0 : rstart;
                                rend = rstart + qlen;
                                rend = (rend < ss->end) ? rend : ss->end;
                                for (i = rstart; i < rend; i++) {
                                    j = i + ss->consensus_start;
                                    if (j >= ss->consensus_end) {
                                        break;
                                    }
                                    if (j >= rs->length) {
                                        break;
                                    }
                                    (rs->bp_total)[j]++;
                                    if(strand == '+'){
                                        (rs->bp_total_plus)[j]++;
                                    } else if (strand == '-'){
                                        (rs->bp_total_minus)[j]++;
                                    }
                                    if (b->core.qual >= mapQ){
                                        (rs->bp_total_unique)[j]++;
                                        if(strand == '+'){
                                            (rs->bp_total_unique_plus)[j]++;
                                        } else if (strand == '-'){
                                            (rs->bp_total_unique_minus)[j]++;
                                        }
                                    }
                                }
                            }
                        }
                        //fill hashFam
                        struct hashEl *hel4 = hashLookup(hashFam, ss->fname);
                        if (hel4 != NULL) {
                            struct repfam *fs = (struct repfam *) hel4->val;
                            fs->read_count++;
                            if (b->core.qual >= mapQ)
                                fs->read_count_unique++;
                        }
                        //fill hashCla
                        struct hashEl *hel5 = hashLookup(hashCla, ss->cname);
                        if (hel5 != NULL) {
                            struct repcla *cs = (struct repcla *) hel5->val;
                            cs->read_count++;
                            if (b->core.qual >= mapQ)
                                cs->read_count_unique++;
                        }
                    } else {
                        slNameAddHead(&(ss->sl), bam1_qname(b));
                        if (b->core.qual >= mapQ)
                            slNameAddHead(&(ss->sl_unique), bam1_qname(b));
                    }
                    reads_repeat++;
                    if (strand == '-'){
                        reads_repeat_minus++;
                    } else if (strand == '+'){
                        reads_repeat_plus++;
                    }
                    if (b->core.qual >= mapQ){
                        reads_repeat_unique++;
                        if (strand == '-'){
                            reads_repeat_unique_minus++;
                        } else if (strand == '+'){
                            reads_repeat_unique_plus++;
                        }
                    }
//                    printf("reads_repeat: %llu, reads_repeat_plus: %llu, reads_repeat_minus: %llu,\n"
//                           "reads_repeat_unique: %llu, reads_repeat_unique_plus: %llu, reads_repeat_unique_minus: %llu\n"
//                           "The strand is %c\n\n",
//                           reads_repeat, reads_repeat_plus, reads_repeat_minus,
//                           reads_repeat_unique, reads_repeat_unique_plus, reads_repeat_unique_minus,
//                           strand);
                    slFreeList(hitList);
                }
            }
        }
        fprintf(stderr, "\r* Processed read ends: %llu\n", (read_end1 + read_end2));
        samclose(samfp);
        bam_destroy1(b);
    }
    //process bam ends
    freeHash(&nochr);
    freeHash(&dup);
    if (outbed_f)
        carefulClose(&outbed_f);
    if (outbed_unique_f)
        carefulClose(&outbed_unique_f);
    cnt[0] = read_end1;
    cnt[1] = read_end2;
    cnt[2] = read_end1_mapped;
    cnt[3] = read_end2_mapped;
    cnt[4] = read_end1_used;
    cnt[5] = read_end2_used;
    cnt[6] = reads_mapped;
    cnt[7] = reads_mapped_unique;
    cnt[8] = reads_nonredundant;
    cnt[9] = reads_repeat;
    cnt[10] = reads_repeat_unique;
    cnt[11] = reads_nonredundant_unique;
    cnt[12] = reads_diff_subfam;

    cnt[13] = reads_repeat_minus;
    cnt[14] = reads_repeat_plus;
    cnt[15] = reads_repeat_unique_minus;
    cnt[16] = reads_repeat_unique_plus;

    return cnt;
}

void cpgBedGraphOverlapRepeat(char *cpgBedGraphFile, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int filter) {
    struct lineFile *infileStream = lineFileOpen2(cpgBedGraphFile, TRUE);
    char *row[20], *line;
    unsigned int start, end, rstart, rend, cpgInRepeat = 0, cpglines = 0;
    double score = 0;
    while ( lineFileNextReal(infileStream, &line)) {
        int numFields = chopByWhite(line, row, ArraySize(row));
        if (numFields < 4)
            errAbort("file %s doesn't appear to be in bedGraph format. At least 4 fields required, got %d", cpgBedGraphFile, numFields);
        cpglines++;
        start = (unsigned int)strtol(row[1], NULL, 0);
        end = (unsigned int)strtol(row[2], NULL, 0);
        score = strtod(row[3], NULL);
        //transfer coordinates
        int i, j;
        struct rmsk *ss = NULL;
        struct binElement *hitList = NULL, *hit;
        struct hashEl *hel2 = hashLookup(hashRmsk, row[0]);
        if (hel2 != NULL) {
            struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
            hitList = binKeeperFind(bs2, start, end);
            if(hitList != NULL) {
                for (hit = hitList; hit !=NULL; hit = hit->next) {
                    ss = (struct rmsk *) hit->val;
                    break;
                }
                if (filter){
                    ss->cpgCount++;
                    ss->cpgTotalScore += score;
                }else{
                    struct hashEl *hel3 = hashLookup(hashRep, ss->name);
                    if (hel3 != NULL){
                        struct rep *rs = (struct rep *) hel3->val;
                        rs->cpgCount++;
                        rs->cpgTotalScore += score;
                        if (rs->length != 0){
                            rstart = start - ss->start;
                            rstart = (rstart < 0) ? 0 : rstart;
                            rend = rstart + 2; //CpG site
                            rend = (rend < ss->end) ? rend : ss->end;
                            for (i = rstart; i < rend; i++) {
                                j = i + ss->consensus_start;
                                if (j >= ss->consensus_end) {
                                    break;
                                }
                                if (j >= rs->length) {
                                    break;
                                }
                                (rs->cpgScore)[j] += score;
                            }
                        }
                    }
                    //fill hashFam
                    struct hashEl *hel4 = hashLookup(hashFam, ss->fname);
                    if (hel4 != NULL) {
                        struct repfam *fs = (struct repfam *) hel4->val;
                        fs->cpgCount++;
                        fs->cpgTotalScore += score;
                    }
                    //fill hashCla
                    struct hashEl *hel5 = hashLookup(hashCla, ss->cname);
                    if (hel5 != NULL) {
                        struct repcla *cs = (struct repcla *) hel5->val;
                        cs->cpgCount++;
                        cs->cpgTotalScore += score;
                    }
                }
                slFreeList(hitList);
                cpgInRepeat++;
            }
        }
    }
    fprintf(stderr, "* Processed CpG sites: %u\n", cpglines);
    fprintf(stderr, "* CpG sites in Repeats: %u\n", cpgInRepeat);
    lineFileClose(&infileStream);
}

unsigned long long int *sam2bed(char *samfile, char *outbed, struct hash *chrHash, int isSam, unsigned int mapQ, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, int treat) {
    samfile_t *samfp;
    FILE *outbed_f = mustOpen(outbed, "w");
    char chr[100], key[100], strand;
    unsigned int start, end, cend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 10);
    unsigned long long int read_end1 = 0, read_end2 = 0;
    unsigned long long int read_end1_mapped = 0, read_end2_mapped = 0;
    unsigned long long int read_end1_used = 0, read_end2_used = 0;
    unsigned long long int reads_nonredundant = 0;
    unsigned long long int reads_nonredundant_unique = 0;
    unsigned long long int reads_mapped = 0;
    unsigned long long int reads_mapped_unique = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
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
    bam1_t *b;
    bam_header_t *h;
    h = samfp->header;
    b = bam_init1();
    int8_t *buf;
    int max_buf;
    buf = 0;
    max_buf = 0;
    uint8_t *seq;
    while ( samread(samfp, b) >= 0) {
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1++;
            }else{
                if(treat)
                    read_end1++;
                else
                    read_end2++;
            }
        }else{
            read_end1++;
        }
        if (((read_end1 + read_end2) % 10000) == 0)
            fprintf(stderr, "\r* Processed read ends: %llu", (read_end1 + read_end2));
        if (b->core.flag & BAM_FUNMAP)
            continue;
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_mapped++;
            }else{
                if (treat)
                    read_end1_mapped++;
                else
                    read_end2_mapped++;
            }
        }else{
            read_end1_mapped++;
        }
        //change chr name to chr1, chr2 ...
        strcpy(chr, h->target_name[b->core.tid]);
        if (addChr){
            if (startsWith("GL", h->target_name[b->core.tid])) {
                continue;
            } else if (sameWord(h->target_name[b->core.tid], "MT")) {
                strcpy(chr,"chrM");
            } else if (!startsWith("chr", h->target_name[b->core.tid])) {
                strcpy(chr, "chr");
                strcat(chr, h->target_name[b->core.tid]);
            }
        }
        //check Ref reads mapped to existed in chromosome size file or not
        struct hashEl *he = hashLookup(nochr, chr);
        if (he != NULL)
            continue;
        cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
        if (cend == 1){
            hashAddInt(nochr, chr, 1);
            warn("* Warning: read ends mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
            continue;
        }
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_used++;
            }else{
                if (treat)
                    read_end1_used++;
                else
                    read_end2_used++;
            }
        }else{
            read_end1_used++;
        }
        //get mapping location for paired-end or single-end
        if (treat){
            reads_mapped++;
            if (b->core.qual >= mapQ)
                reads_mapped_unique++;
            start = (unsigned int) b->core.pos;
            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
            end = min(cend, (unsigned int)tmpend);
            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
            if (extension) {
                if (strand == '+'){
                    end = min(start + extension, cend);
                }else{
                    if (end < extension)
                        start = 0;
                    else
                        start = end - extension;
                    //start = max(end - extension, 0);
                }
            }

        }else{
        if (b->core.flag & BAM_FPAIRED) {
            if (!(b->core.flag & BAM_FMUNMAP)){
                if (b->core.flag & BAM_FREAD1){
                    if (abs(b->core.isize) > iSize || b->core.isize == 0){
                        continue;
                    }else{
                        reads_mapped++;
                        if (b->core.qual >= mapQ)
                            reads_mapped_unique++;
                        if (b->core.isize > 0){
                            start = (unsigned int) b->core.pos;
                            strand = '+';
                            int tmpend = start + b->core.isize;
                            end = min(cend, (unsigned int)tmpend);
                        }else{
                            start = (unsigned int) b->core.mpos;
                            strand = '-';
                            int tmpend = start - b->core.isize;
                            end = min(cend, (unsigned int)tmpend);
                        }
                
                    }
                }else{
                    continue;
                }
            }else{
                if (discardWrongEnd){
                    continue;
                }else{
                    reads_mapped++;
                    if (b->core.qual >= mapQ)
                        reads_mapped_unique++;
                    start = (unsigned int) b->core.pos;
                    int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
                    end = min(cend, (unsigned int)tmpend);
                    strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
                    if (extension) {
                        if (strand == '+'){
                            end = min(start + extension, cend);
                        }else{
                            if (end < extension)
                                start = 0;
                            else
                                start = end - extension;
                            //start = max(end - extension, 0);
                        }
                    }
                }
            }
        }else{
            reads_mapped++;
            if (b->core.qual >= mapQ)
                reads_mapped_unique++;
            start = (unsigned int) b->core.pos;
            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
            end = min(cend, (unsigned int)tmpend);
            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
            if (extension) {
                if (strand == '+'){
                    end = min(start + extension, cend);
                }else{
                    if (end < extension)
                        start = 0;
                    else
                        start = end - extension;
                }
            }
        }
    }
        //remove dup or not
        if (rmDup){
            if (sprintf(key, "%s:%u:%u:%c", chr, start, end, strand) < 0)
                errAbort("Mem ERROR");
            struct hashEl *hel = hashLookup(dup, key);
            if (hel == NULL) {
                hashAddInt(dup, key, 1);
            } else {
                continue;
            }
        }
        reads_nonredundant++;
        if (b->core.qual >= mapQ)
            reads_nonredundant_unique++;
        //output bed
        int i, qlen = b->core.l_qseq;
        if(b->core.qual >= mapQ){
            fprintf(outbed_f, "%s\t%u\t%u\t", chr, start, end);
            //print read sequence
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
            fprintf(outbed_f, "%s", (char*)buf);

            fprintf(outbed_f, "\t%i\t%c\n", b->core.qual, strand);
        }
    }
    fprintf(stderr, "\r* Processed read ends: %llu\n", (read_end1 + read_end2));
    samclose(samfp);
    free(buf);
    bam_destroy1(b);
    freeHash(&nochr);
    freeHash(&dup);
    carefulClose(&outbed_f);
    cnt[0] = read_end1;
    cnt[1] = read_end2;
    cnt[2] = read_end1_mapped;
    cnt[3] = read_end2_mapped;
    cnt[4] = read_end1_used;
    cnt[5] = read_end2_used;
    cnt[6] = reads_mapped;
    cnt[7] = reads_mapped_unique;
    cnt[8] = reads_nonredundant;
    cnt[9] = reads_nonredundant_unique;
    return cnt;
}

unsigned long long int *PEsamFile2nodupRepbedFile(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, unsigned int iSize) {
    samfile_t *samfp;
    char chr[100], prn[500], key[100], strand;
    unsigned int start, end, cend, rstart, rend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 5);
    unsigned long long int mapped_reads_num = 0, reads_num = 0, reads_used = 0, unique_reads = 0, repeat_reads = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
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
    strcpy(prn, "empty");
    bam1_t *b[2];
    int curr, has_prev;
    bam_header_t *h;
    h = samfp->header;
    b[0] = bam_init1();
    b[1] = bam_init1();
    curr = 0;
    has_prev = 0;
    while ( samread(samfp, b[curr]) >= 0) {
        
	bam1_t *cur = b[curr], *pre = b[1-curr];
	if (has_prev) {
	    if (strcmp(bam1_qname(cur), bam1_qname(pre)) == 0) { // identical pair name
		has_prev = 0;
	    }
	} else has_prev = 1;
	curr = 1 - curr;
        
        if (!has_prev) {
            if ( sameString (bam1_qname(cur), prn)) 
                continue;
        }
        reads_num++;
        if ((reads_num % 10000) == 0)
            fprintf(stderr, "\r* Processed reads pair: %llu", reads_num);
        if (!has_prev)
            strcpy(prn, bam1_qname(cur));
        if (has_prev){
            //if (cur->core.tid < 0 || pre->core.tid < 0)
            if (cur->core.flag & BAM_FUNMAP || pre->core.flag & BAM_FUNMAP)
                continue;
        }
        if (has_prev){
            if (cur->core.qual < mapQ || pre->core.qual < mapQ)
                continue;
        }
        mapped_reads_num++;
        if (has_prev) {
            if (cur->core.tid != pre->core.tid)
                continue;
        }
        if (has_prev) {
            if (abs(pre->core.isize) > iSize || pre->core.isize == 0)
            //if (abs(pre->core.isize) > iSize)
                continue;
            strcpy(chr, h->target_name[cur->core.tid]);
            //change chr name to chr1, chr2 ...
            if (addChr){
                if (startsWith("GL", h->target_name[cur->core.tid])) {
                    continue;
                } else if (sameWord(h->target_name[cur->core.tid], "MT")) {
                    strcpy(chr,"chrM");
                } else if (!startsWith("chr", h->target_name[cur->core.tid])) {
                    strcpy(chr, "chr");
                    strcat(chr, h->target_name[cur->core.tid]);
                }
            }
            struct hashEl *he = hashLookup(nochr, chr);
            if (he != NULL)
                continue;
            cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
            if (cend == 1){
                hashAddInt(nochr, chr, 1);
                warn("* Warning: reads mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
                continue;
            }
            reads_used++;
            int tmpend;
            if (pre->core.isize > 0){
                start = (unsigned int) pre->core.pos;
                strand = '+';
                tmpend = start + pre->core.isize;
            }else{
                start = (unsigned int) cur->core.pos;
                strand = '-';
                tmpend = start - pre->core.isize;
            }
            end = min(cend, (unsigned int)tmpend);
            //remove dup first
            if (rmDup == 1){
                if (sprintf(key, "%s:%u:%u:%c", chr, start, end, strand) < 0)
                    errAbort("Mem ERROR");
                struct hashEl *hel = hashLookup(dup, key);
                if (hel == NULL) {
                    hashAddInt(dup, key, 1);
                } else {
                    continue;
                }
            }
            unique_reads++;
            //transfer coordinates
            int i, j;
            unsigned int qlen = end - start;
            struct binElement *hitList = NULL, *hit;
            struct hashEl *hel2 = hashLookup(hashRmsk, chr);
            if (hel2 != NULL) {
                struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
                hitList = binKeeperFind(bs2, start, end);
                if(hitList != NULL) {
                    for (hit = hitList; hit !=NULL; hit = hit->next) {
                        struct rmsk *ss = (struct rmsk *) hit->val;
                        if (filter == 0){
                            struct hashEl *hel3 = hashLookup(hashRep, ss->name);
                            if (hel3 != NULL){
                                struct rep *rs = (struct rep *) hel3->val;
                                rs->read_count++;
                                if (rs->length != 0){
                                    rstart = start - ss->start;
                                    rstart = (rstart < 0) ? 0 : rstart;
                                    rend = rstart + qlen;
                                    rend = (rend < ss->end) ? rend : ss->end;
                                    for (i = rstart; i < rend; i++) {
                                        j = i + ss->consensus_start;
                                        if (j >= ss->consensus_end) {
                                            break;
                                        }
                                        if (j >= rs->length) {
                                            break;
                                        }
                                        (rs->bp_total)[j]++;
                                    }
                                }
                            }
                            //fill hashFam
                            struct hashEl *hel4 = hashLookup(hashFam, ss->fname);
                            if (hel4 != NULL) {
                                struct repfam *fs = (struct repfam *) hel4->val;
                                fs->read_count++;
                            }
                            //fill hashCla
                            struct hashEl *hel5 = hashLookup(hashCla, ss->cname);
                            if (hel5 != NULL) {
                                struct repcla *cs = (struct repcla *) hel5->val;
                                cs->read_count++;
                            }
                        } else {
                            slNameAddHead(&(ss->sl), bam1_qname(cur));
                        }
                        break;
                    }
                    repeat_reads++;
                    slFreeList(hitList);
                }
            }
        }
    }
    fprintf(stderr, "\r* Processed reads pair: %llu\n", reads_num);
    samclose(samfp);
    bam_destroy1(b[0]);
    bam_destroy1(b[1]);
    freeHash(&nochr);
    freeHash(&dup);
    cnt[0] = reads_num;
    cnt[1] = mapped_reads_num;
    cnt[2] = reads_used;
    cnt[3] = unique_reads;
    cnt[4] = repeat_reads;
    return cnt;
}

void freermsk(struct rmsk *s){
    freeMem(s->name);
    freeMem(s->cname);
    freeMem(s->fname);
    free(s);
}

void rmsk2binKeeperHash(char *rmskfile, struct hash *chrHash, struct hash *repHash, struct hash **hashRmsk, struct hash **hashRep, struct hash **hashFam, struct hash **hashCla, int filterField, char *filterName){
    char strand;
    char *row[17];
    struct lineFile *repeat_stream = lineFileOpen2(rmskfile, TRUE);
    int repeat_num = 0;
    struct hash *hash1 = newHash(0); //hashRmsk
    struct hash *hash2 = newHash(0); //hashRep
    struct hash *hash3 = newHash(0); //hashFam
    struct hash *hash4 = newHash(0); //hashCla
    while( lineFileNextRow(repeat_stream, row, ArraySize(row))){
        if (filterField != 0){
            if (strcmp(filterName, row[filterField]) != 0)
                continue;
        }
        struct rmsk *s = malloc(sizeof(struct rmsk));
        repeat_num++;
        s->chr = cloneString(row[5]);
        strand = row[9][0];
        if ( strand == '+' ) {
            s->consensus_start = (unsigned int)strtol(row[13], NULL, 0);    
        } else {
            s->consensus_start = (unsigned int)strtol(row[15], NULL, 0);    
        }
        s->consensus_end = (unsigned int)strtol(row[14], NULL, 0);
        s->start = (unsigned int)strtol(row[6], NULL, 0);
        s->end = (unsigned int)strtol(row[7], NULL, 0);
        s->name = cloneString(row[10]);
        s->cname = cloneString(row[11]);
        s->fname = cloneString(row[12]);
        s->length = s->end - s->start;
        s->sl = NULL;
        s->sl_unique = NULL;
        s->cpgCount = 0;
        s->cpgTotalScore = 0;
        
        struct hashEl *hel = hashLookup(hash1, s->chr);
        if (hel != NULL) {
            struct binKeeper *bk = (struct binKeeper *) hel->val;
            binKeeperAdd(bk, s->start, s->end, s);
        } else {
            int size = hashIntValDefault(chrHash, s->chr, 0);
            if (size == 0) {
                freermsk(s);
                continue;
            }
            struct binKeeper *bk = binKeeperNew(0, size);
            binKeeperAdd(bk, s->start, s->end, s);
            hashAdd(hash1, s->chr, bk);
        }

        if (filterField != 0)
            continue;
        
        struct hashEl *hel2 = hashLookup(hash2, s->name);
        if (hel2 != NULL) {
            struct rep *rr = (struct rep *) hel2->val;
            rr->genome_count++;
            rr->total_length += s->length;
        } else {
            struct rep *ns = malloc(sizeof(struct rep));
            ns->name = cloneString(s->name);
            ns->cname = cloneString(s->cname);
            ns->fname = cloneString(s->fname);
            ns->genome_count = 1;
            ns->total_length = s->length;
            ns->read_count = 0;
            ns->read_count_unique = 0;
            ns->cpgCount = 0;
            ns->cpgTotalScore = 0;
            ns->length = hashIntValDefault(repHash, ns->name, 0);
            ns->bp_total = malloc(sizeof(unsigned int) * ns->length);
            ns->bp_total_unique = malloc(sizeof(unsigned int) * ns->length);
            ns->bp_total_plus = malloc(sizeof(unsigned int) * ns->length);
            ns->bp_total_unique_plus = malloc(sizeof(unsigned int) * ns->length);
            ns->bp_total_minus = malloc(sizeof(unsigned int) * ns->length);
            ns->bp_total_unique_minus = malloc(sizeof(unsigned int) * ns->length);
            ns->cpgScore = malloc(sizeof(double) * ns->length);
            int k;
            for (k = 0; k < ns->length; k++){ 
                (ns->bp_total)[k] = 0;
                (ns->bp_total_unique)[k] = 0;
                (ns->cpgScore)[k] = 0;
            }
            hashAdd(hash2, ns->name, ns);
        }
        //fill hash3
        struct hashEl *hel3 = hashLookup(hash3, s->fname);
        if (hel3 != NULL) {
            struct repfam *rr = (struct repfam *) hel3->val;
            rr->genome_count++;
            rr->total_length += s->length;
        } else {
            struct repfam *ns = malloc(sizeof(struct repfam));
            ns->fname = cloneString(s->fname);
            ns->cname = cloneString(s->cname);
            ns->genome_count = 1;
            ns->total_length = s->length;
            ns->read_count = 0;
            ns->read_count_unique = 0;
            ns->cpgCount = 0;
            ns->cpgTotalScore = 0;
            hashAdd(hash3, ns->fname, ns);
        }
        //fill hash4
        struct hashEl *hel4 = hashLookup(hash4, s->cname);
        if (hel4 != NULL) {
            struct repcla *rr = (struct repcla *) hel4->val;
            rr->genome_count++;
            rr->total_length += s->length;
        } else {
            struct repcla *ns = malloc(sizeof(struct repcla));
            ns->cname = cloneString(s->cname);
            ns->genome_count = 1;
            ns->total_length = s->length;
            ns->read_count = 0;
            ns->read_count_unique = 0;
            ns->cpgCount = 0;
            ns->cpgTotalScore = 0;
            hashAdd(hash4, ns->cname, ns);
        }
    }
    lineFileClose(&repeat_stream);
    if (filterField == 0){
        fprintf(stderr, "* Total %d repeats found.\n", repeat_num);
    } else{
        if ( repeat_num <= 0 )
            errAbort("* No repeats found related to [%s], typo? or specify wrong repName/Class/Family filter?", filterName);
        fprintf(stderr, "* Total %d repeats for [%s].\n", repeat_num, filterName);
    }
    *hashRmsk = hash1;
    *hashRep = hash2;
    *hashFam = hash3;
    *hashCla = hash4;
}

void writeFilterOut(struct hash *hash, char *out, int readlist, int threshold, char *subfam, unsigned long long int reads_num){ 
    FILE *out_stream;
    out_stream = mustOpen(out, "w");
    int j = 0;
    struct hashEl *he;
    struct hashCookie cookie = hashFirst(hash);
    if (readlist)
        fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#chr", "start", "end", "length", "repName", "repClass", "repFamily", "readsCount", "RPKM", "RPM", "readsList");
    else
        fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#chr", "start", "end", "length", "repName", "repClass", "repFamily", "readsCount", "RPKM", "RPM");
    while ( (he = hashNext(&cookie)) != NULL ) {
        struct binKeeper *bk = (struct binKeeper *) he->val;
        struct binKeeperCookie becookie = binKeeperFirst(bk);
        struct binElement *be;
        while( (be = binKeeperNext(&becookie)) != NULL ){
            struct rmsk *os = (struct rmsk *) (be->val);
            int count = slCount(os->sl);
            if (readlist) {
                if (count >= threshold){
                    j++;
                    slReverse(&(os->sl));
                    char *s = slNameListToString(os->sl, ',');
                    fprintf(out_stream, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%.3f\t%.3f\t%s\n", os->chr, os->start, os->end, os->length, os->name, os->cname, os->fname, count, cal_rpkm((unsigned long long int)count,(unsigned long long int)os->length, reads_num), cal_rpm((unsigned long long int)count, reads_num), s);
                    freeMem(s);
                }
            } else {
                if (count >= threshold){
                    j++;
                    fprintf(out_stream, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%.3f\t%.3f\n", os->chr, os->start, os->end, os->length, os->name, os->cname, os->fname, count, cal_rpkm((unsigned long long int)count, (unsigned long long int)os->length, reads_num), cal_rpm((unsigned long long int)count, reads_num));
                }
            }
            slFreeList(&(os->sl));
        }
        binKeeperFree(&bk);
    }
    fclose(out_stream);
    fprintf(stderr, "* Total %d [%s] TEs have at least %d reads mapped.\n", j, subfam, threshold);
}

void writeFilterOutMRE(struct hash *hash, char *out, char *subfam, double scoreThreshold){ 
    FILE *out_stream;
    out_stream = mustOpen(out, "w");
    int j = 0;
    struct hashEl *he;
    struct hashCookie cookie = hashFirst(hash);
        fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#chr", "start", "end", "length", "repName", "repClass", "repFamily", "covered_CpG_site", "total_CpG_score");
    while ( (he = hashNext(&cookie)) != NULL ) {
        struct binKeeper *bk = (struct binKeeper *) he->val;
        struct binKeeperCookie becookie = binKeeperFirst(bk);
        struct binElement *be;
        while( (be = binKeeperNext(&becookie)) != NULL ){
            struct rmsk *os = (struct rmsk *) (be->val);
            double score = os->cpgTotalScore;
            if (score > scoreThreshold){
                j++;
                fprintf(out_stream, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%.3f\n", os->chr, os->start, os->end, os->length, os->name, os->cname, os->fname, os->cpgCount, os->cpgTotalScore);
            }
        }
        binKeeperFree(&bk);
    }
    fclose(out_stream);
    fprintf(stderr, "* Total %d [%s] TEs have CpG score larger than %.3f.\n", j, subfam, scoreThreshold);
}

void sortBedfile(char *bedfile) {
    struct lineFile *lf = NULL;
    FILE *f = NULL;
    struct bedLine *blList = NULL, *bl;
    char *line;
    int lineSize;

    lf = lineFileOpen2(bedfile, TRUE);
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
