#include "generic.h"
#include <stdbool.h>

int filter_usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "Obtain alignment statistics of individual loci of each repeat subfamily, family or class.\n\n");
    fprintf(stderr, "Usage:   iteres filter [options] <chromosome size file> <repeat size file> <rmsk.txt> <bam/sam alignment file>\n\n");
    fprintf(stderr, "Options: -S       input is SAM [off]\n");
    fprintf(stderr, "         -Q       mapping Quality threshold [10]\n");
    fprintf(stderr, "         -g       coverage threshold for overlapping [0.0001]\n");
    fprintf(stderr, "         -N       normalized by number of (0: non-redundant unique mapped reads, 1: unique reads, 2: mapped reads, 3: total reads) [0])\n");
    fprintf(stderr, "         -n       use repName (subfamily) as filter [null]\n");
    fprintf(stderr, "         -f       use repFamily as filter [null]\n");
    fprintf(stderr, "         -c       use repClass as filter [null]\n");
    fprintf(stderr, "         -t       only output repeats have more than [1] reads mapped\n");
    fprintf(stderr, "         -r       output the list of reads [off]\n");
    fprintf(stderr, "         -R       remove redundant reads [off]\n");
    fprintf(stderr, "         -T       treat 1 paired-end read as 2 single-end reads [off]\n");
    fprintf(stderr, "         -D       discard if only one end mapped in a paired end reads [off]\n");
    fprintf(stderr, "         -C       Add 'chr' string as prefix of reference sequence [off]\n");
    fprintf(stderr, "         -E       extend reads to represent fragment [150], specify 0 if want no extension\n");
    fprintf(stderr, "         -I       Insert length threshold [500]\n");
    fprintf(stderr, "         -o       output prefix [basename of input without extension]\n");
    fprintf(stderr, "         -W       input the length of the cage-seq window\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    return 1;
}

/* main filter function */
int main_filter(int argc, char *argv[]){
    char *output, *subfam, *out, *outReport;
    unsigned long long int *cnt;
    struct hash *hashRmsk = newHash(0);
    struct hash *hashRep = newHash(0);
    struct hash *hashFam = newHash(0);
    struct hash *hashCla = newHash(0);
    bool cageFlag = false;
    int optSam = 0, optthreshold = 1;
    char *optoutput = NULL, *optname = NULL, *optclass = NULL, *optfamily = NULL;
    unsigned int optreadlist = 0, optQual = 10, optisize = 500, optExt = 150;
    int filterField = 0, c, optDup = 0, optNorm = 0, optaddChr = 0, optDis =0, optTreat = 0, optcagewindow = 0;
    float optCov = 0.0001;

    time_t start_time, end_time;
    start_time = time(NULL);
    
    while ((c = getopt(argc, argv, "SQ:g:N:n:c:t:f:rRTDCE:I:o:h?W:")) >= 0) {
        switch (c) {
            case 'S': optSam = 1; break;
            case 'Q': optQual = (unsigned int)strtol(optarg, 0, 0); break;
            case 'g': optCov = atof(optarg); break;
            case 'N': optNorm = (unsigned int)strtol(optarg, 0, 0); break;
            case 't': optthreshold = (unsigned int)strtol(optarg, 0, 0); break;
            case 'r': optreadlist = 1; break;
            case 'R': optDup = 1; break;
            case 'T': optTreat = 1; break;
            case 'D': optDis = 1; break;
            case 'C': optaddChr = 1; break;
            case 'n': optname = strdup(optarg); break;
            case 'c': optclass = strdup(optarg); break;
            case 'f': optfamily = strdup(optarg); break;
            case 'E': optExt = (unsigned int)strtol(optarg, 0, 0); break;
            case 'I': optisize = (unsigned int)strtol(optarg, 0, 0); break;
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return filter_usage(); break;
            case 'W':  if (optarg == NULL) {
                    // Handle -W without parameters
                    optcagewindow = 50;
                    cageFlag = true;
                } else {
                    optcagewindow = (unsigned int)strtol(optarg, 0, 0);
                    cageFlag = true;
                } break;
            default: return 1;
        }
    }
    if (optind + 4 > argc)
        return filter_usage();

    char *chr_size_file = argv[optind];
    char *rep_size_file = argv[optind+1];
    char *rmsk_file = argv[optind+2];
    char *sam_file = argv[optind+3];
    
    if ( (optname && optclass) || (optname && optfamily) || (optclass && optfamily) || (optname && optclass && optfamily))
        errAbort("Please specify only one filter, either -n, -c or -f.");
    
    int nindex = 0;
    if (optNorm == 0){
        nindex = 7;
    } else if (optNorm == 1){
        nindex = 8;
    } else if (optNorm == 2) {
        nindex = 6;    
    } else if (optNorm == 3) {
        nindex = 4;    
    } else{
        errAbort("Wrong normalization method specified");
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
    
    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(sam_file)));
    }
    
    struct hash *chrHash = hashNameIntFile(chr_size_file);
    struct hash *repHash = hashNameIntFile(rep_size_file);
    
    fprintf(stderr, "* Start to parse the rmsk file\n");
    rmsk2binKeeperHash(rmsk_file, chrHash, repHash, &hashRmsk, &hashRep, &hashFam, &hashCla, filterField, subfam);
    
    //sam file
    fprintf(stderr, "* Start to parse the SAM/BAM file\n");
    //if (optPair){
    //    cnt = PEsamFile2nodupRepbedFile(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, optQual, 1, optDup, optaddChr, optisize);
    //} else {
    //    cnt = samFile2nodupRepbedFile(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, optQual, 1, optDup, optaddChr);
    //}
//    cnt = samFile2nodupRepbedFileNew(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, optQual, 1, optDup, optaddChr, optDis, optisize, optExt, optCov, optTreat, NULL, NULL, 0);

//    if(optcagewindow > 0)
//        cnt = samFile2nodupRepbedFileNewCage(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, optQual, 1, optDup, optaddChr, optDis, optisize, optExt, optCov, optTreat, NULL, NULL, 0, optcagewindow);

    if(cageFlag)
        cnt = samFile2nodupRepbedFileNewCage(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, optQual, 1, optDup, optaddChr, optDis, optisize, optExt, optCov, optTreat, NULL, NULL, 0, optcagewindow);
    else
        cnt = samFile2nodupRepbedFileNew(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, optQual, 1, optDup, optaddChr, optDis, optisize, optExt, optCov, optTreat, NULL, NULL, 0);

    fprintf(stderr, "* Preparing the output file\n");
    if (asprintf(&out, "%s_%s.iteres.loci", output, subfam) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outReport, "%s_%s.iteres.reportloci", output, subfam) < 0)
        errAbort("Preparing output wrong");
    
    writeFilterOut(hashRmsk, out, optreadlist, optthreshold, subfam, cnt[nindex]); 
    
    //write report file
    fprintf(stderr, "* Preparing report file\n");
    writeReport(outReport, cnt, optQual, subfam);
    
    hashFree(&chrHash);
    hashFree(&repHash);
    hashFree(&hashRmsk);
    hashFree(&hashRep);
    hashFree(&hashFam);
    hashFree(&hashCla);
    free(out);
    free(outReport);
    
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;    
}

