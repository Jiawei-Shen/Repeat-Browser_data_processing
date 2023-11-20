#include "generic.h"

int cpgstat_usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "obtain CpG statistics for each repeat subfamily, family and class.\n\n");
    fprintf(stderr, "Usage:   iteres cpgstat [options] <chromosome size file> <repeat size file> <rmsk.txt> <CpG bedGraph file>\n\n");
    fprintf(stderr, "Options: -w       keep the wiggle file [off]\n");
    fprintf(stderr, "         -o       output prefix [basename of input without extension]\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    return 1;
}

/* main cpgstat function */
int main_cpgstat (int argc, char *argv[]) {
    
    char *output, *outWig, *outbigWig, *outStat, *outFam, *outCla;
    int optkeepWig = 0, c;
    char *optoutput = NULL;
    time_t start_time, end_time;
    struct hash *hashRmsk = newHash(0);
    struct hash *hashRep = newHash(0);
    struct hash *hashFam = newHash(0);
    struct hash *hashCla = newHash(0);
    start_time = time(NULL);
    while ((c = getopt(argc, argv, "wo:h?")) >= 0) {
        switch (c) {
            case 'w': optkeepWig = 1; break;
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return cpgstat_usage(); break;
            default: return 1;
        }
    }
    if (optind + 4 > argc)
        return cpgstat_usage();

    char *chr_size_file = argv[optind];
    char *rep_size_file = argv[optind+1];
    char *rmsk_file = argv[optind+2];
    char *bedgraph_file = argv[optind+3];
    
    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(bedgraph_file)));
    }
    
    if(asprintf(&outWig, "%s.CpGstat.wig", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outbigWig, "%s.CpGstat.bigWig", output) < 0)
        errAbort("Mem Error.\n");
    if (asprintf(&outStat, "%s.CpG.subfamily.stat", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outFam, "%s.CpG.family.stat", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outCla, "%s.CpG.class.stat", output) < 0)
        errAbort("Preparing output wrong");

    struct hash *chrHash = hashNameIntFile(chr_size_file);
    struct hash *repHash = hashNameIntFile(rep_size_file);
    
    fprintf(stderr, "* Start to parse the rmsk file\n");
    rmsk2binKeeperHash(rmsk_file, chrHash, repHash, &hashRmsk, &hashRep, &hashFam, &hashCla, 0, "ALL");
    
    //bedgraph file
    fprintf(stderr, "* Start to parse the bedGraph file\n");
    cpgBedGraphOverlapRepeat(bedgraph_file, hashRmsk, hashRep, hashFam, hashCla, 0);

    fprintf(stderr, "* Writing stats and Wig file\n");
    MREwriteWigandStat(hashRep, hashFam, hashCla, outStat, outWig, outFam, outCla);

    fprintf(stderr, "* Generating bigWig files\n");
    bigWigFileCreate(outWig, rep_size_file, 256, 1024, 0, 1, outbigWig);
    
    if (!optkeepWig){
        unlink(outWig);
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
    free(outStat);
    free(outFam);
    free(outCla);
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}

