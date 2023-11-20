#include <argp.h>
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "fa.h"

const char *argp_program_version = "RevMaskFaExt 0.1";
const char *argp_program_bug_address = "<http://wang.wustl.edu>";

/* definitions of structures*/

struct arguments {
    char *args[2];
    int extend;
    char *output;
};

/* definitions of functions */

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = state->input;

    switch (key)
    {
        case 'e':
            arguments->extend = (int)strtol(arg, NULL, 0);
            break;
        case 'o':
            arguments->output = arg;
            break;
        case ARGP_KEY_ARG:
            if( state->arg_num >= 2)
            {
                argp_usage(state);
            }
            arguments->args[state->arg_num] = arg;
            break;
        case ARGP_KEY_END:
            if (state->arg_num < 2)
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
    struct lineFile *repeat_stream;
    char *chr;
    char *row[17];
    struct dnaSeq *seqList = NULL, *seq1List = NULL, *seq, *seq1;
    FILE *out_stream;
    struct arguments arguments;
    time_t start_time, end_time;
    struct hash *hash = hashNew(0);
    struct hash *hash1 = hashNew(0);
    start_time = time(NULL);

    static struct argp_option options[] = 
    {
        {"extend", 'e', "extend-length", 0, "extention length for both side"},
        {"output", 'o', "output-file", 0, "Name for output file (default: stdout)"},
        {0}
    };

    static char args_doc[] = "<rmsk.txt> <genome.fa>";

    static char doc[] = "This program mask out non-repeat region with N, except repeat regions and their flanking region.\n";

    static struct argp argp = {options, parse_opt, args_doc, doc};

    /* set default parameters*/
    
    arguments.output = NULL;
    arguments.extend = 0;
    
    /* parameter magic*/
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    fprintf(stderr, "* Start to read in sequence file ...\n");
    
    seqList = faReadAllMixed(arguments.args[1]);
    for (seq = seqList; seq != NULL; seq = seq->next){
        hashAdd(hash, seq->name, seq);
        //add N sequence
        //seq1 = cloneMem(seq, sizeof(seq));
        AllocVar(seq1);
        seq1->name = cloneString(seq->name);
        seq1->size = seq->size;
        seq1->dna = cloneMem(seq->dna, seq->size + 1);
        memset(seq1->dna, 'N', seq1->size);
        hashAdd(hash1, seq1->name, seq1);
        slAddHead(&seq1List, seq1);
    }
    slReverse(&seq1List);
    
    char *repeat_file = arguments.args[0];
    repeat_stream = lineFileOpen(repeat_file, TRUE);
    
    if(arguments.output) {
        out_stream = mustOpen(arguments.output, "w");
    } else {
        out_stream = stdout;
    }
    
    fprintf(stderr, "* Start to parse the rmsk file ...\n");
    
    //rmsk file
    while( lineFileNextRow(repeat_stream, row, ArraySize(row))){
        int repSize, seqSize, start, end, i;
        chr = cloneString(row[5]);
        start = (int)strtol(row[6], NULL, 0);
        end = (int)strtol(row[7], NULL, 0);
        start -= arguments.extend;
        end += arguments.extend;
        seq = hashFindVal(hash, chr);
        seq1 = hashFindVal(hash1, chr);
        seqSize = seq->size;
        if (start < 0 || end > seqSize || start > end) {
            if (start < 0) start = 0;
            if (end > seqSize) end = seqSize;
            if (end < start) end = start;
        }
        repSize = end - start;
        if (repSize > 0) {
            for ( i = start; i < end; ++i)
                seq1->dna[i] = seq->dna[i];
        }

    }
    lineFileClose(&repeat_stream);
    
    fprintf(stderr, "* Preparing the output file ...\n");
    for (seq = seq1List; seq != NULL; seq = seq->next){
        faWriteNext(out_stream, seq->name, seq->dna, seq->size);
    }
    //struct hashEl *hel;
    //struct hashCookie cookie = hashFirst(hash1);
    //while ( (hel = hashNext(&cookie)) != NULL ) {
    //    struct dnaSeq *s = hel->val;
    //    //fprintf(out_stream, ">%s\n%s\n", seq1->name, seq1->dna);
    //    faWriteNext(out_stream, s->name, s->dna, s->size);
    //}
    carefulClose(&out_stream);
    //freeHashAndVals(&hash); //bug, already free'ed
    hashFree(&hash);
    hashFree(&hash1);
    
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    
    return 0;
}
