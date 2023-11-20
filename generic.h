#include "sam.h"
#include "from_kent.h"

#define ITERES_VERSION "0.3.3-r123"

//struct hold contens from rmsk line
struct rmsk {
    char *chr;
    unsigned int start, end, consensus_start, consensus_end, length;
    char *name, *fname, *cname;
    struct slName *sl;
    struct slName *sl_unique;
    unsigned int cpgCount;
    double cpgTotalScore;
};

//struct hold contents for repeat subfamily
struct rep {
    char *name, *fname, *cname;
    unsigned long long int read_count, read_count_unique, genome_count, total_length;
    unsigned int length;
    unsigned int *bp_total, *bp_total_unique;
    unsigned int *bp_total_plus, *bp_total_unique_plus;
    unsigned int *bp_total_minus, *bp_total_unique_minus;
    unsigned int cpgCount;
    double *cpgScore;
    double cpgTotalScore;
};

//struct hold contents for repeat family
struct repfam {
    char *fname, *cname;
    unsigned long long int read_count, read_count_unique, genome_count, total_length;
    unsigned int cpgCount;
    double cpgTotalScore;
};

//struct hold contents for repeat class
struct repcla {
    char *cname;
    unsigned long long int read_count, read_count_unique, genome_count, total_length;
    unsigned int cpgCount;
    double cpgTotalScore;
};



int stat_usage();
int main_stat(int argc, char *argv[]);
int filter_usage();
int main_filter (int argc, char *argv[]);
//int nearby_usage();
//int main_nearby(int argc, char *argv[]);
int cpgstat_usage();
int main_cpgstat(int argc, char *argv[]);
int cpgfilter_usage();
int main_cpgfilter (int argc, char *argv[]);


char *get_filename_without_ext(char *filename);
char *get_filename_ext(char *filename);
bool is_file(const char* path); 
bool is_dir(const char* path);
double cal_rpkm (unsigned long long int reads_count, unsigned long long int total_length, unsigned long long int mapped_reads_num);
double cal_rpm (unsigned long long int reads_count, unsigned long long int mapped_reads_num);
struct lineFile *lineFileOpen2(char *fileName, bool zTerm);
void writeReport(char *outfile, unsigned long long int *cnt, unsigned int mapQ, char *subfam);
void plotMappingStat(unsigned long long int *cnt, char *prefix);
void writeWigandStat(struct hash *hash, struct hash *hash1, struct hash *hash2, char *of1, char *of2, char *of3, char *of4, char *of5, unsigned long long int reads_num, unsigned long int reads_num_unique);
void writeWigandStatCage(struct hash *hash, struct hash *hash1, struct hash *hash2, char *of1, char *of2, char *of6, char *of3, char *of4, char *of5, char *of7, unsigned long long int reads_num, unsigned long int reads_num_unique);
void MREwriteWigandStat(struct hash *hash, struct hash *hash1, struct hash *hash2, char *of1, char *of2, char *of3, char *of4);
unsigned long long int *samFile2nodupRepbedFile(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr);
float getCov(unsigned int aStart, unsigned int aEnd, unsigned int start, unsigned int end);
int mapped2diffSubfam(struct hash *hashRmsk, char *subfam, int nm, char *ahstring, int qlen);
unsigned long long int *samFile2nodupRepbedFileNew(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, float minCoverage, int treat, char *outbed, char *outbed_unique, int diffSubfam);
unsigned long long int *samFile2nodupRepbedFileNewCage(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, float minCoverage, int treat, char *outbed, char *outbed_unique, int diffSubfam, int optcagewindow);
unsigned long long int *samFiles2nodupRepbedFileNew(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, float minCoverage, int treat, char *outbed, char *outbed_unique, int diffSubfam);
unsigned long long int *samFiles2nodupRepbedFileNewCage(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, float minCoverage, int treat, char *outbed, char *outbed_unique, int diffSubfam, int optcagewindow);
unsigned long long int *PEsamFile2nodupRepbedFile(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, unsigned int iSize);
void freermsk(struct rmsk *s);
void rmsk2binKeeperHash(char *rmskfile, struct hash *chrHash, struct hash *repHash, struct hash **hashRmsk, struct hash **hashRep, struct hash **hashFam, struct hash **hashCla, int filterField, char *filterName);
void writeFilterOut(struct hash *hash, char *out, int readlist, int threshold, char *subfam, unsigned long long int reads_num);
void sortBedfile(char *bedfile);
void writeReportDensity(char *outfile, unsigned long long int *cnt, unsigned int mapQ);
unsigned long long int *sam2bed(char *samfile, char *outbed, struct hash *chrHash, int isSam, unsigned int mapQ, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, int treat);
void cpgBedGraphOverlapRepeat(char *cpgBedGraphFile, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int filter);
void writeFilterOutMRE(struct hash *hash, char *out, char *subfam, double scoreThreshold); 
