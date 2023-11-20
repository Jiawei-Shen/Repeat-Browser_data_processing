#include "generic.h"

int nearby_usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "Obtain nearby genes from locations listed in a bed file by querying UCSC database.\n\n");
    fprintf(stderr, "Usage:   iteres nearby [options] <bed file>\n\n");
    fprintf(stderr, "Options: -d       database to query [hg19]\n");
    fprintf(stderr, "         -n       output how many genes each direction [1]\n");
    fprintf(stderr, "         -o       output prefix [basename of input without extension]\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note: the bed file should contain at least 3 fields which were [chr] [start] [end]\n");
    fprintf(stderr, "      also you need to have an internet connection\n\n");
    return 1;
}

/* main nearby function */
int main_nearby(int argc, char *argv[]){
    struct lineFile *infile_stream;
    FILE *out_stream;
    char *output, *outfile, *line, *row[20];
    int c;
    char *optoutput = NULL;
    char *optdb = "hg19";
    unsigned int optthreshold = 1;
    time_t start_time, end_time;
    start_time = time(NULL);
    while ((c = getopt(argc, argv, "d:n:o:h?")) >= 0) {
        switch (c) {
            case 'n': optthreshold = (unsigned int)strtol(optarg, 0, 0); break;
            case 'd': optdb = strdup(optarg); break;
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return nearby_usage(); break;
            default: return 1;
        }
    }
    if (optind + 1 > argc)
        return nearby_usage();
    char *infile = argv[optind];
    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(infile)));
    }
    if(asprintf(&outfile, "%s.iteres.gene", output) < 0)
        errAbort("Mem Error.\n");
    out_stream = mustOpen(outfile, "w");
    infile_stream = lineFileOpen(infile, TRUE);
    fprintf(stderr, "* Making connection to UCSC database %s\n", optdb);
    struct sqlConnection *con = sqlConnectRemote("genome-mysql.cse.ucsc.edu", "genome", "", optdb);
    fprintf(stderr, "* Start to parse the input file\n");
    char *tchr, *query1, *query2;
    int tstart, tend;
    while( lineFileNextReal(infile_stream, &line)){
        int numFields = chopByWhite(line, row, ArraySize(row));
        if (numFields < 3)
            errAbort("file %s doesn't appear to be in bed format. At least 3 fields required, got %d", infile, numFields);
        tchr = cloneString(row[0]);
        tstart = (int)strtol(row[1], NULL, 0);
        tend = (int)strtol(row[2], NULL, 0);
        
        //upstream
        char **row1;
        if (asprintf(&query1, "select e.chrom,e.txStart,e.txEnd,e.alignID,j.geneSymbol FROM \
                   knownGene e, \
                   kgXref j \
                   WHERE e.alignID = j.kgID AND e.chrom='%s' AND e.txEnd < %d \
                   ORDER BY e.txEnd DESC limit %d", tchr, tstart, optthreshold) < 0)
            errAbort("Memory allocating for query string wrong");

        struct sqlResult *sr1 = sqlGetResult(con, query1);
        while ( (row1 = sqlNextRow(sr1)) != NULL ){
            fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\n", row1[0], row1[1], row1[2], row1[3], row1[4], "upstream");
        }
        sqlFreeResult(&sr1);
        
        //downstream
        char **row2;
        if (asprintf(&query2, "select e.chrom,e.txStart,e.txEnd,e.alignID,j.geneSymbol FROM \
                   knownGene e, \
                   kgXref j \
                   WHERE e.alignID = j.kgID AND e.chrom='%s' AND e.txStart > %d \
                   ORDER BY e.txStart ASC limit %d", tchr, tend, optthreshold) < 0)
            errAbort("Memory allocating for query string wrong");

        struct sqlResult *sr2 = sqlGetResult(con, query2);
        while ( (row2 = sqlNextRow(sr2)) != NULL ){
            fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\n", row2[0], row2[1], row2[2], row2[3], row2[4], "downstream");
        }
        sqlFreeResult(&sr2);
    }
    lineFileClose(&infile_stream);
    sqlDisconnect(&con);
    fclose(out_stream);
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}

