#include <argp.h>
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "binRange.h"

const char *argp_program_version = "evaluateBed 0.1.1";

const char *argp_program_bug_address = "<http://wang.wustl.edu>";

/* definitions of structures*/

//struct hold contens from rmsk line
struct rmsk {
    char *chr;
    int start, end, length;
    char *name, *fname, *cname;
};

struct arguments {
    char *args[2];
    char *output;
};

/* definitions of functions */

void freeRmsk(struct rmsk *r) {
    free(r->chr);
    free(r->name);
    free(r->fname);
    free(r->cname);
    free(r);
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


static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = state->input;

    switch (key)
    {
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


struct hash * read_in_rmsk (char *repeat_file) {
    fprintf(stderr, "* Start to parse the rmsk file ...\n");
    struct hash *hash = hashNew(0);
    char *row[17], strand;
    struct lineFile *repeat_stream;
    repeat_stream = lineFileOpen(repeat_file, TRUE);
    while( lineFileNextRow(repeat_stream, row, ArraySize(row))){
        struct rmsk *s = malloc(sizeof(struct rmsk));
        s->chr = cloneString(row[5]);
        strand = row[9][0];
        s->start = (int)strtol(row[6], NULL, 0);
        s->end = (int)strtol(row[7], NULL, 0);
        s->name = cloneString(row[10]);
        s->cname = cloneString(row[11]);
        s->fname = cloneString(row[12]);
        s->length = s->end - s->start;
        
        struct hashEl *hel = hashLookup(hash, s->chr);
        if (hel != NULL) {
            struct binKeeper *bk = (struct binKeeper *) hel->val;
            binKeeperAdd(bk, s->start, s->end, s);
        } else {
            struct binKeeper *bk = binKeeperNew(0, 249250621);
            binKeeperAdd(bk, s->start, s->end, s);
            hashAdd(hash, s->chr, bk);
        }
    }
    lineFileClose(&repeat_stream);
    return hash;
}

/* main funtion */
int main (int argc, char *argv[]) {
    
    samfile_t *samfp;
    char *output, prn[200], *outerr, *outerr2;
    int mapped_reads_num, reads_num, mapped_to_loc, mapped_to_subfam;
    struct arguments arguments;
    struct hash *hash = hashNew(0);
    
    time_t start_time, end_time;
    

    start_time = time(NULL);

    static struct argp_option options[] = 
    {
        {"output", 'o', "OUTPUTFILE", 0, "Name for output file (default: the basename of input file.report)"},
        {0}
    };

    static char args_doc[] = "<rmsk.txt> <sam/bed alignment file>";

    static char doc[] = "This program evaluate sam/bed alignment file.\n";

    static struct argp argp = {options, parse_opt, args_doc, doc};


    reads_num = mapped_reads_num = mapped_to_loc = mapped_to_subfam = 0;

    /* set default parameters*/
    
    arguments.output = NULL;
    
    /* parameter magic*/
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    
    hash = read_in_rmsk(arguments.args[0]);
    
    char *sam_file = arguments.args[1];
    
    if(arguments.output) {
        output = arguments.output;
    } else {
        //output = get_filename_without_ext(basename(sam_file));
        if (asprintf(&output, "%s.report", get_filename_without_ext(basename(sam_file))) < 0)
            errAbort("Preparing output wrong");
    }

    if (asprintf(&outerr, "%s.subfamerr", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outerr2, "%s.notaligntorepeat", output) < 0)
        errAbort("Preparing output wrong");
    
    FILE *outerrhandle;
    outerrhandle = mustOpen(outerr, "w");
    FILE *outerrhandle2;
    outerrhandle2 = mustOpen(outerr2, "w");
    
    if (sameWord( get_filename_ext(sam_file), "sam")) {
        if ( (samfp = samopen(sam_file, "r", 0)) == 0) {
            fprintf(stderr, "Fail to open SAM file %s\n", sam_file);
            return 1;
        }
    } else {
        if ( (samfp = samopen(sam_file, "rb", 0)) == 0) {
            fprintf(stderr, "Fail to open BAM file %s\n", sam_file);
            return 1;
        }
    }
       
    if (0) {
        struct hashEl *he;
        struct hashCookie cookie = hashFirst(hash);
        while ( (he = hashNext(&cookie)) != NULL ) {
            struct rmsk *os = (struct rmsk *) (he->val);
            printf("%s\t%d\n", os->chr, os->start);
        }
    }
    
    fprintf(stderr, "* Start to parse the SAM/BAM file ...\n");
    
    //sam file
    strcpy(prn, "empty");
    bam1_t *b = bam_init1();
    int readStart, readRealstart, t1count, t2count;
    while ( samread(samfp, b) >= 0) {
        struct sam *sa = fetch_sa(b, samfp);
        
        if ( sameString (sa->name, prn)) {
            freeSam(sa);
            continue;
        }
        reads_num++;
        strcpy(prn, sa->name);

        if ( sameString(sa->chr, "*") ){
            freeSam(sa);
            continue;
        }
        mapped_reads_num++;
        
        char *tmpreadname = cloneString(sa->name);

        t1count = chopByChar(sa->name, ':', NULL, 0);
        char **t1 = (char **)needMem((size_t) (t1count * sizeof(char *)));
        chopByChar(sa->name, ':', t1, t1count);
        char *readChr = cloneString(t1[0]);
        char *readRep = cloneString(t1[3]);
        readStart = (int)strtol(t1[1], NULL, 0);
        t2count = chopByChar(t1[4], '_', NULL, 0);
        char **t2 = (char **)needMem((size_t) (t2count * sizeof(char *)));
        chopByChar(t1[4], '_', t2, t2count);
        readRealstart = (int)strtol(t2[1], NULL, 0);

        if (arguments.rmskfile){
            if (sameWord(readChr, sa->chr) && (abs(readStart + readRealstart - sa->start) < 5)) {
                mapped_to_loc++;
                mapped_to_subfam++;
            }else{
                struct binElement *hitList = NULL, *hit;
                struct hashEl *hel2 = hashLookup(hash, sa->chr);
                if (hel2 != NULL) {
                    struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
                    hitList = binKeeperFind(bs2, sa->start, sa->end);
                    if (hitList != NULL ) {
                        for ( hit = hitList; hit != NULL ; hit = hit->next) {
                            struct rmsk *ss = (struct rmsk *) (hit->val);
                            if (sameWord(ss->name, readRep)){
                                mapped_to_subfam++;
                            } else {
                                fprintf(outerrhandle, "%s\t%s\n", tmpreadname, ss->name);
                            }
                            //freeRmsk(ss);
                            break;
                        }
                        slFreeList(&hitList);
                    } else {
                        fprintf(outerrhandle2, "%s\t%s\t%d\t%d\t%c\n", tmpreadname, sa->chr, sa->start, sa->end, sa->strand);
                    }
                }
            }
            freeMem(t1);
            freeMem(t2);
            freeMem(readChr);
            freeMem(readRep);
        } else {
            char *readChr1, *readRep1;
            int readStart1;
            int t3count = chopByChar(sa->chr, ':', NULL, 0);
            if (t3count < 4)
                errAbort("[Error], seems the sam/bam was aligned to normal reference, not repeat index. Please specify rmsk.txt file by -r option.");
            char **t3 = (char **)needMem((size_t) (t3count * sizeof(char *)));
            chopByChar(sa->chr, ':', t3, t3count);
            readChr1 = cloneString(t3[0]);
            readRep1 = cloneString(t3[3]);
            readStart1 = (int)strtol(t3[1], NULL, 0);
            if ( sameWord(readChr, readChr1) && (readStart == readStart1) && (abs(readRealstart - sa->start) < 5)){
                mapped_to_loc++;
                mapped_to_subfam++;
            } else {
                if (sameWord(readRep1, readRep)){
                    mapped_to_subfam++;
                } else{
                    fprintf(outerrhandle, "%s\t%s\n", tmpreadname, readRep1);
                }
            }
            freeMem(t3);
            freeMem(readChr1);
            freeMem(readRep1);
        }
        freeSam(sa);
        freeMem(tmpreadname);
    }
    bam_destroy1(b);
    samclose(samfp);
    
    fclose(outerrhandle);
    fclose(outerrhandle2);

    FILE *out_stream;
    char *out;

    if (asprintf(&out, "%s", output) < 0)
        errAbort("Preparing output wrong");
    
    out_stream = mustOpen(out, "w");
    
    fprintf(stderr, "* Preparing the output file ...\n");

    fprintf(out_stream, "* Total %d reads.\n", reads_num);
    fprintf(out_stream, "* Total %d reads mapped.\n", mapped_reads_num);
    fprintf(out_stream, "* Total %d reads mapped to right location.\n", mapped_to_loc);
    fprintf(out_stream, "* Total %d reads mapped to right subfamily.\n", mapped_to_subfam);
    
    fclose(out_stream);
    
    if (arguments.rmskfile)
        hashFree(&hash);
    if(!arguments.output)
        free(output);
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}
