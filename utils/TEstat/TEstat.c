#include <argp.h>
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "binRange.h"
#include "sam.h"

const char *argp_program_version = "TEstat 0.2";

const char *argp_program_bug_address = "<http://wang.wustl.edu>";

/* definitions of structures*/

//struct hold contens from rmsk line
struct rmsk {
    char *chr;
    int start, end, length;
    char *name, *fname, *cname;
    struct slName *sl;
};

//struct hold contens from sam line
struct sam {
    int start, end, length;
    char *chr;
    char *name;
    uint32_t qual:8;
    char strand;
};

struct arguments {
    char *args[2];
    int readlist, threshold;
    char *output, *name, *class, *family;
};

/* definitions of functions */

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

double cal_rpkm (unsigned int reads_count, unsigned int total_length, unsigned int mapped_reads_num) {
    return reads_count / (mapped_reads_num * 1e-9 * total_length);
}

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = state->input;

    switch (key)
    {
        case 'n':
            arguments->name = arg;
            break;
        case 'c':
            arguments->class = arg;
            break;
        case 'f':
            arguments->family = arg;
            break;
        case 'r':
            arguments->readlist = 1;
            break;
        case 't':
            arguments->threshold = (int)strtol(arg, NULL, 0);
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

struct sam * fetch_sa (const bam1_t *b, void *data){
    struct sam *s = malloc(sizeof(struct sam));
    samfile_t *fp = (samfile_t *) data;
    uint32_t *cigar = bam1_cigar(b);
    const bam1_core_t *c = &b->core;
    int i, l;
    s->name = cloneString(bam1_qname(b));
    s->chr = cloneString("*");
    s->start = 0;
    s->end = 0;
    s->length = 0;
    s->qual = 0;
    s->strand = '*';
    if (b->core.tid < 0) return s;
    for (i = l = 0; i < c->n_cigar; ++i) {
        int op = cigar[i]&0xf;
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
            l += cigar[i]>>4;
    }
    s->chr = cloneString(fp->header->target_name[c->tid]);
    s->start = c->pos;
    s->end = c->pos + l;
    s->length = l;
    s->qual = c->qual;
    s->strand = (c->flag&BAM_FREVERSE) ? '-' : '+';
    return s;
}

/* main funtion */
int main (int argc, char *argv[]) {
    
    struct lineFile *repeat_stream;
    samfile_t *samfp;
    char strand, *output, *prn, *subfam;
    unsigned int mapped_reads_num;
    int repeat_num, sam_in_repeat, reads_num;
    char *row[17];
    struct arguments arguments;
    
    time_t start_time, end_time;
    
    struct hash *hash = hashNew(0);

    start_time = time(NULL);

    static struct argp_option options[] = 
    {
        {"name", 'n', "repName", 0, "Use repName (subfamily) as filter"},
        {"class", 'c', "repClass", 0, "Use repClass as filter"},
        {"family", 'f', "repFamily", 0, "Use repFamily as filter"},
        {"readlist", 'r', 0, 0, "Output list of read names in last column (default: off)"},
        {"threshold", 't', "NUMTHRESHOLD", 0, "Only output TEs have more than t reads mapped (default: 1)"},
        {"output", 'o', "OUTPUTFILE", 0, "Name for output file (default: the basename of input file without extension.bed)"},
        {0}
    };

    static char args_doc[] = "<rmsk.txt> <sam/bam alignment file>";

    static char doc[] = "This program parses sam/bam alignment file, obtain statistics of reads mapped to repeats.\nUsing either repName, repClass or repFamily as filter, if you didn't specify any filter, will output reads mapping stat for all TEs.\n";

    static struct argp argp = {options, parse_opt, args_doc, doc};


    repeat_num = sam_in_repeat = reads_num = mapped_reads_num = 0;

    /* set default parameters*/
    
    arguments.output = NULL;
    arguments.name = NULL;
    arguments.class = NULL;
    arguments.family = NULL;
    arguments.readlist = 0;
    arguments.threshold = 1;
    
    /* parameter magic*/
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    if ( (arguments.name && arguments.class) || (arguments.name && arguments.family) || (arguments.class && arguments.family) || (arguments.name && arguments.class && arguments.family))
        errAbort("Please specify only one filter, either -n, -c or -f.");

    //if ( !(arguments.name || arguments.class || arguments.family) )
        
    subfam = cloneString("ALL");
    
    if (arguments.name) {
        arguments.class = NULL;
        arguments.family = NULL;
        subfam = cloneString(arguments.name);
    }else if (arguments.class) {
        arguments.name = NULL;
        arguments.family = NULL;
        subfam = cloneString(arguments.class);
    } else if (arguments.family) {
        arguments.name = NULL;
        arguments.class= NULL;
        subfam = cloneString(arguments.family);
    }

    if (sameString(subfam, "ALL"))
        fprintf(stderr, "* You didn't specify any filter, will output all repeats\n");
    
    //char *subfam = arguments.args[0];
    char *repeat_file = arguments.args[0];
    char *sam_file = arguments.args[1];
    
    if(arguments.output) {
        output = arguments.output;
    } else {
        //output = get_filename_without_ext(basename(sam_file));
        if (asprintf(&output, "%s.bed", get_filename_without_ext(basename(sam_file))) < 0)
            errAbort("Preparing output wrong");
    }

    repeat_stream = lineFileOpen(repeat_file, TRUE);
    
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
    
    fprintf(stderr, "* Start to parse the rmsk file ...\n");
    
    char *target = "ALL";
    //rmsk file
    while( lineFileNextRow(repeat_stream, row, ArraySize(row))){
        if (arguments.name)
            target = cloneString(row[10]);
        else if (arguments.class)
            target = cloneString(row[11]);
        else if (arguments.family)
            target = cloneString(row[12]);
        if( sameString(target, subfam) ){
            struct rmsk *s = malloc(sizeof(struct rmsk));
            repeat_num++;
            s->chr = cloneString(row[5]);
            strand = row[9][0];
            s->start = (int)strtol(row[6], NULL, 0);
            s->end = (int)strtol(row[7], NULL, 0);
            s->name = cloneString(row[10]);
            s->cname = cloneString(row[11]);
            s->fname = cloneString(row[12]);
            s->length = s->end - s->start;
            s->sl = NULL;
        
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
    }
    lineFileClose(&repeat_stream);
    //freeMem(target);   
       
    if (0) {
        struct hashEl *he;
        struct hashCookie cookie = hashFirst(hash);
        while ( (he = hashNext(&cookie)) != NULL ) {
            struct rmsk *os = (struct rmsk *) (he->val);
            printf("%s\t%d\n", os->chr, os->start);
        }
    }

    //assert(repeat_num > 0);
    if ( repeat_num <= 0 )
        errAbort("* No repeats found related to [%s], typo? or specify wrong repName/Class/Family filter?", subfam);
    
    fprintf(stderr, "* Total %d repeats for [%s].\n", repeat_num, subfam);
    
    fprintf(stderr, "* Start to parse the SAM/BAM file ...\n");
    
    //sam file
    prn = cloneString("empty");
    bam1_t *b = bam_init1();
    while ( samread(samfp, b) >= 0) {
        struct sam *sa = fetch_sa(b, samfp);
        
        if ( sameString (sa->name, prn)) {
            free(sa);
            continue;
        }
        reads_num++;
        prn = cloneString(sa->name);

        if ( sameString(sa->chr, "*") ){
            free(sa);
            continue;
        }
        mapped_reads_num++;

        struct binElement *hitList = NULL, *hit;
        struct hashEl *hel2 = hashLookup(hash, sa->chr);
        if (hel2 != NULL) {
            struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
            hitList = binKeeperFind(bs2, sa->start, sa->end);
            if (hitList != NULL ) {
                for ( hit = hitList; hit != NULL ; hit = hit->next) {
                    struct rmsk *ss = (struct rmsk *) (hit->val);
                    slNameAddHead(&(ss->sl), sa->name);
                    break;
                }
                sam_in_repeat++;
                slFreeList(&hitList);
            }
        }
        free(sa);
    }
    bam_destroy1(b);
    samclose(samfp);

    freeMem(prn);   
    
    fprintf(stderr, "* Total %d reads.\n", reads_num);
    fprintf(stderr, "* Total %u reads mapped.\n", mapped_reads_num);
    fprintf(stderr, "* Total %d reads mapped to [%s] TE.\n", sam_in_repeat, subfam);

    FILE *out_stream;
    char *out;

    if (asprintf(&out, "%s", output) < 0)
        errAbort("Preparing output wrong");
    
    out_stream = mustOpen(out, "w");
    
    fprintf(stderr, "* Preparing the output file ...\n");
    
    int j = 0;
    struct hashEl *he;
    struct hashCookie cookie = hashFirst(hash);
    
    if (arguments.readlist)
        fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#chr", "start", "end", "length", "repName", "repClass", "repFamily", "readsCount", "readsList");
    else
        fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#chr", "start", "end", "length", "repName", "repClass", "repFamily", "readsCount");

    while ( (he = hashNext(&cookie)) != NULL ) {
        struct binKeeper *bk = (struct binKeeper *) he->val;
        struct binKeeperCookie becookie = binKeeperFirst(bk);
        struct binElement *be;
        while( (be = binKeeperNext(&becookie)) != NULL ){
            struct rmsk *os = (struct rmsk *) (be->val);
            int count = slCount(os->sl);
            if (arguments.readlist) {
                if (count >= arguments.threshold){
                    j++;
                    slReverse(&(os->sl));
                    char *s = slNameListToString(os->sl, ',');
                    fprintf(out_stream, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%s\n", os->chr, os->start, os->end, os->length, os->name, os->cname, os->fname, count, s);
                    freeMem(s);
                }
            } else {
                if (count >= arguments.threshold){
                    j++;
                    fprintf(out_stream, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\n", os->chr, os->start, os->end, os->length, os->name, os->cname, os->fname, count);
                }
            }
            slFreeList(&(os->sl));
        }
        binKeeperFree(&bk);
    }

    fclose(out_stream);
    
    //freeHashAndVals(&hash); //bug, already free'ed
    hashFree(&hash);

    if(!arguments.output)
        free(output);
    
    fprintf(stderr, "* Total %d [%s] TEs have at least %d reads mapped.\n", j, subfam, arguments.threshold);
    
    end_time = time(NULL);

    //freeMem(subfam);

    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    
    return 0;
}
