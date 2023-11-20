#include "generic.h"

int cpgfilter_usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "obtain CpG statistics for each repeat locus.\n\n");
    fprintf(stderr, "Usage:   iteres cpgfilter [options] <chromosome size file> <repeat size file> <rmsk.txt> <CpG bedGraph file>\n\n");
    fprintf(stderr, "Options: -n       use repName (subfamily) as filter [null]\n");
    fprintf(stderr, "         -f       use repFamily as filter [null]\n");
    fprintf(stderr, "         -c       use repClass as filter [null]\n");
    fprintf(stderr, "         -t       only output repeats have more than [0] CpG score\n");
    fprintf(stderr, "         -o       output prefix [basename of input without extension]\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    return 1;
}

/* main cpgstat function */
int main_cpgfilter (int argc, char *argv[]) {
    
    char *output, *out, *subfam;
    int c, filterField = 0;
    double optthreshold = 0;
    char *optoutput = NULL, *optname = NULL, *optclass = NULL, *optfamily = NULL;
    time_t start_time, end_time;
    struct hash *hashRmsk = newHash(0);
    struct hash *hashRep = newHash(0);
    struct hash *hashFam = newHash(0);
    struct hash *hashCla = newHash(0);
    start_time = time(NULL);
    while ((c = getopt(argc, argv, "n:c:f:t:o:h?")) >= 0) {
        switch (c) {
            case 'n': optname = strdup(optarg); break;
            case 'c': optclass = strdup(optarg); break;
            case 'f': optfamily = strdup(optarg); break;
            case 't': optthreshold = strtod(optarg, NULL); break;
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return cpgfilter_usage(); break;
            default: return 1;
        }
    }
    if (optind + 4 > argc)
        return cpgfilter_usage();

    char *chr_size_file = argv[optind];
    char *rep_size_file = argv[optind+1];
    char *rmsk_file = argv[optind+2];
    char *bedgraph_file = argv[optind+3];
    
    if ( (optname && optclass) || (optname && optfamily) || (optclass && optfamily) || (optname && optclass && optfamily))
        errAbort("Please specify only one filter, either -n, -c or -f.");
    
    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(bedgraph_file)));
    }
    
    subfam = cloneString("ALL");
    if (optname) {
        optclass = NULL;
        optfamily = NULL;
        subfam = cloneString(optname);
        filterField = 10;
    }else if (optclass) {
        optname = NULL;
        optfamily = NULL;
        subfam = cloneString(optclass);
        filterField = 11;
    } else if (optfamily) {
        optname = NULL;
        optclass= NULL;
        subfam = cloneString(optfamily);
        filterField = 12;
    }
    if (sameString(subfam, "ALL")){
        fprintf(stderr, "* You didn't specify any filter, will output all repeats\n");
        filterField = 0;
    }
    

    struct hash *chrHash = hashNameIntFile(chr_size_file);
    struct hash *repHash = hashNameIntFile(rep_size_file);
    
    fprintf(stderr, "* Start to parse the rmsk file\n");
    rmsk2binKeeperHash(rmsk_file, chrHash, repHash, &hashRmsk, &hashRep, &hashFam, &hashCla, filterField, subfam);
    
    //bedgraph file
    fprintf(stderr, "* Start to parse the bedGraph file\n");
    cpgBedGraphOverlapRepeat(bedgraph_file, hashRmsk, hashRep, hashFam, hashCla, 1);

    fprintf(stderr, "* Preparing the output file\n");
    if (asprintf(&out, "%s_%s.CpG.loci", output, subfam) < 0)
        errAbort("Preparing output wrong");
    
    writeFilterOutMRE(hashRmsk, out, subfam, optthreshold);
    
    //cleaning
    hashFree(&chrHash);
    hashFree(&repHash);
    hashFree(&hashRmsk);
    hashFree(&hashRep);
    hashFree(&hashFam);
    hashFree(&hashCla);
    free(out);
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}

