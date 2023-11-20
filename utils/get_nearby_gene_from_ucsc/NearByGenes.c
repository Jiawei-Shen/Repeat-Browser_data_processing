#include <argp.h>
#include "common.h"
#include "linefile.h"
#include "jksql.h"

const char *argp_program_version = "NearByGenes 0.1";

const char *argp_program_bug_address = "<http://wang.wustl.edu>";

/* definitions of structures*/

struct arguments {
    char *args[1];
    int upstream, downstream, threshold;
    char *output, *db;
};

/* definitions of functions */


static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = state->input;

    switch (key)
    {
        case 'd':
            arguments->db = arg;
            break;
        case 'u':
            arguments->upstream = 1;
            break;
        case 'r':
            arguments->downstream = 1;
            break;
        case 't':
            arguments->threshold = (int)strtol(arg, NULL, 0);
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

/* main funtion */
int main (int argc, char *argv[]) {
    
    struct lineFile *infile_stream;
    FILE *out_stream;
    char *output, *row[5], *db;
    struct arguments arguments;
    
    time_t start_time, end_time;

    start_time = time(NULL);

    static struct argp_option options[] = 
    {
        {"db", 'd', "database", 0, "database to query (default hg19)"},
        {"upstream", 'u', 0, 0, "output upstream genes (default: off)"},
        {"downstream", 'r', 0, 0, "output downstream genes (default: on)"},
        {"threshold", 't', "NUMTHRESHOLD", 0, "output how many genes (default: 1)"},
        {"output", 'o', "OUTPUTFILE", 0, "Name for output file (default: stdout)"},
        {0}
    };

    static char args_doc[] = "<location file, first 3 columns should be [chr] [start] [end]>";

    static char doc[] = "This program fetches nearby genes for locations.\n";

    static struct argp argp = {options, parse_opt, args_doc, doc};

    db = cloneString("hg19");

    /* set default parameters*/
    
    arguments.output = NULL;
    arguments.db = db;
    arguments.upstream = 0;
    arguments.downstream = 1;
    arguments.threshold = 1;
    
    /* parameter magic*/
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    
    char *infile = arguments.args[0];
    
    if(arguments.output) {
        output = arguments.output;
        out_stream = mustOpen(output, "w");
    } else {
        out_stream = stdout;
    }
    
    

    infile_stream = lineFileOpen(infile, TRUE);
    
    
    fprintf(stderr, "* Make connection to UCSC database %s...\n", arguments.db);
    
    struct sqlConnection *con = sqlConnectRemote("genome-mysql.cse.ucsc.edu", "genome", "", arguments.db);
    
    fprintf(stderr, "* Start to parse the input file ...\n");

    char *tchr, *query;
    int tstart, tend;
    while( lineFileNextRow(infile_stream, row, ArraySize(row))){
        
        tchr = cloneString(row[0]);
        tstart = (int)strtol(row[1], NULL, 0);
        tend = (int)strtol(row[2], NULL, 0);
    
        //upstream
        if(arguments.upstream){
            char **row;
            if (asprintf(&query, "select e.chrom,e.txStart,e.txEnd,e.alignID,j.geneSymbol FROM \
                                    knownGene e, \
                                    kgXref j \
                                    WHERE e.alignID = j.kgID AND e.chrom='%s' AND e.txEnd < %d \
                      ORDER BY e.txEnd DESC limit %d", tchr, tstart, arguments.threshold) < 0)
                errAbort("Memory allocating for query string wrong");

            struct sqlResult *sr = sqlGetResult(con, query);
            while ( (row = sqlNextRow(sr)) != NULL ){
                fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\n", row[0], row[1], row[2], row[3], row[4], "upstream");
            }
            sqlFreeResult(&sr);
        }
        
        //downstream
        if(arguments.downstream){
            char **row;
            if (asprintf(&query, "select e.chrom,e.txStart,e.txEnd,e.alignID,j.geneSymbol FROM \
                   knownGene e, \
                      kgXref j \
                      WHERE e.alignID = j.kgID AND e.chrom='%s' AND e.txStart > %d \
                      ORDER BY e.txStart ASC limit %d", tchr, tend, arguments.threshold) < 0)
                errAbort("Memory allocating for query string wrong");

            struct sqlResult *sr = sqlGetResult(con, query);
            while ( (row = sqlNextRow(sr)) != NULL ){
                fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\n", row[0], row[1], row[2], row[3], row[4], "downstream");
            }
            sqlFreeResult(&sr);
        }
        
    }
    lineFileClose(&infile_stream);
    sqlDisconnect(&con);

    fclose(out_stream);
    
    end_time = time(NULL);

    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    
    return 0;
}
