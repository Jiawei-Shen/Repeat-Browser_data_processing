#include <argp.h>
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "binRange.h"
#include "bigWig.h"
#include "sam.h"
#define debug 0

const char *argp_program_version = "RAP 0.2";

const char *argp_program_bug_address = "<http://wang.wustl.edu>";

/* definitions of structures*/

//struct hold contens from rmsk line
struct rmsk {
    char *chr;
    int start, end, consensus_start, consensus_end, length;
    char *name, *fname, *cname;
};

//struct hold contens from sam line
struct sam {
    int start, end, length;
    char *chr;
    char *name;
    uint32_t qual:8;
    char strand;
};

//struct hold contents from chrom size file
struct chr {
    char *name;
    int length, count_total;
    int *bp_total;
};

//struct hold contents from repeat size file
struct rep {
    char *name, *fname, *cname;
    int length;
    unsigned int count_total;
    int genome_count;
    unsigned int total_length;
    int *bp_total;
};

//struct hold contents for repeat which has no size info
struct repns {
    char *name, *fname, *cname;
    int length, count_total, genome_count, total_length;
};

struct arguments {
    char *args[4];
    int keepWig;
    char *output_base;
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
        case 'w':
            arguments->keepWig = 1;
            break;
        case 'o':
            arguments->output_base = arg;
            break;
        case ARGP_KEY_ARG:
            if( state->arg_num >= 4)
            {
                argp_usage(state);
            }
            arguments->args[state->arg_num] = arg;
            break;
        case ARGP_KEY_END:
            if (state->arg_num < 4)
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
    
    struct lineFile *repeat_stream, *chr_size_stream, *rep_size_stream;
    FILE *no_stream;
    samfile_t *samfp;
    char strand, *output, *prn;
    unsigned int chr_count, rep_count, mapped_reads_num;
    int repeat_num, i, j, sam_in_repeat, k , m, reads_num;
    char *chrRow[2], *repRow[2], *row[17];
    struct arguments arguments;
    
    time_t start_time, end_time;
    
    struct hash *hashchr = hashNew(0);
    struct hash *hashrep = hashNew(0);
    struct hash *hashbk = hashNew(0);
    struct hash *hashns = hashNew(0);

    start_time = time(NULL);

    static struct argp_option options[] = 
    {
        {"keepwig", 'w', 0, 0, "Specify this option to keep the temporary wiggle files (default: false)"},
        {"output_base", 'o', "OUTPUTBASE", 0, "Specify the basename for output files (default: the basename of input file without extension)"},
        {0}
    };

    static char args_doc[] = "<chromosome_size_file> <repeat_size_file> <rmsk.txt> <sam alignment file>";

    static char doc[] = "This program parses sam alignment file, obtain statistics of reads mapped to repeats.\n";

    static struct argp argp = {options, parse_opt, args_doc, doc};


    repeat_num = sam_in_repeat = reads_num = mapped_reads_num = 0;

    /* set default parameters*/
    
    arguments.output_base = NULL;
    arguments.keepWig = 0;
    
    /* parameter magic*/
    argp_parse(&argp, argc, argv, 0, 0, &arguments);
    
    char *chr_size_file = arguments.args[0];
    char *rep_size_file = arguments.args[1];
    char *repeat_file = arguments.args[2];
    char *sam_file = arguments.args[3];
    
    if(arguments.output_base) {
        output = arguments.output_base;
    } else {
        output = get_filename_without_ext(basename(sam_file));
    }

    chr_size_stream = lineFileOpen(chr_size_file, TRUE);
    rep_size_stream = lineFileOpen(rep_size_file, TRUE);
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
    
    //read in chr size file
    while( lineFileRow(chr_size_stream, chrRow)){
        char *name = chrRow[0];
        int size = lineFileNeedNum(chr_size_stream, chrRow, 1);
        struct chr *sc;
        sc = malloc(sizeof(struct chr));
        sc->name = cloneString(name);
        sc->length = size;
        sc->bp_total = malloc(sizeof(int) * sc->length);
        for (k = 0; k < size; k++) {
            (sc->bp_total)[k] = 0;
        }
        sc->count_total = 0;
        hashAdd(hashchr, name, sc);
    
        struct binKeeper *bk = binKeeperNew(0, size);
        assert(size > 1);
        hashAdd(hashbk, name, bk);

    }
    lineFileClose(&chr_size_stream);
    
    chr_count = hashNumEntries(hashchr);
    fprintf (stderr, "* Total %u chromosomes found in chr size file.\n", chr_count);
    
    //read in rep size file
    while( lineFileRow(rep_size_stream, repRow)){
        char *name = repRow[0];
        int size = lineFileNeedNum(rep_size_stream, repRow, 1);
        struct rep *sr;
        sr = malloc(sizeof(struct rep));
        sr->name = cloneString(name);
        sr->cname = cloneString("null");
        sr->fname = cloneString( "null");
        sr->length = size;
        sr->bp_total = malloc(sizeof(int) * sr->length);
        for (k = 0; k < size; k++) {
            (sr->bp_total)[k] = 0;
        }
        sr->count_total = 0;
        sr->genome_count = 0;
        sr->total_length = 0;
        hashAdd(hashrep, name, sr);
    }
    lineFileClose(&rep_size_stream);

    rep_count = hashNumEntries(hashrep);
    fprintf (stderr, "* Total %u repeat sub-families found in repeat size file.\n", rep_count);

    //repeat file
    while( lineFileNextRow(repeat_stream, row, ArraySize(row))){
        struct rmsk *s = malloc(sizeof(struct rmsk));
        repeat_num++;
        s->chr = cloneString(row[5]);
        strand = row[9][0];
        if ( strand == '+' ) {
            s->consensus_start = (int)strtol(row[13], NULL, 0);    
        } else {
            s->consensus_start = (int)strtol(row[15], NULL, 0);    
        }
        s->consensus_end = (int)strtol(row[14], NULL, 0);
        s->start = (int)strtol(row[6], NULL, 0);
        s->end = (int)strtol(row[7], NULL, 0);
        s->name = cloneString(row[10]);
        s->cname = cloneString(row[11]);
        s->fname = cloneString(row[12]);
        s->length = s->end - s->start;
        
        struct hashEl *hel = hashLookup(hashbk, s->chr);
        if (hel != NULL) {
            struct binKeeper *bk = (struct binKeeper *) hel->val;
            binKeeperAdd(bk, s->start, s->end, s);
        }

        struct hashEl *hel2 = hashLookup(hashrep, s->name);
        if (hel2 != NULL) {
            struct rep *rr = (struct rep *) hel2->val;
            rr->genome_count++;
            rr->total_length += s->length;
            if (strcmp(rr->cname, "null") == 0)
                rr->cname = cloneString(s->cname);
            if (strcmp(rr->fname, "null") == 0)
                rr->fname = cloneString(s->fname);
        } else {
            struct hashEl *hel3 = hashLookup(hashns, s->name);
            if (hel3 == NULL) {
                struct repns *ns = malloc(sizeof(struct repns));
                ns->name = cloneString(s->name);
                ns->cname = cloneString(s->cname);
                ns->fname = cloneString(s->fname);
                ns->genome_count = 1;
                ns->total_length = s->length;
                ns->length = s->length;
                ns->count_total = 0;
                hashAdd(hashns, ns->name, ns);
            } else {
                struct repns *ns = (struct repns *) hel3->val;
                ns->genome_count++;
                ns->total_length += s->length;
            }
        }
    }
    lineFileClose(&repeat_stream);
    
    fprintf(stderr, "* Total %d repeats found.\n", repeat_num);
    
    fprintf(stderr, "* Start to parse the SAM/BAM file ...\n");
    
    //sam file
    int rstart, rend;
    prn = cloneString("empty");
    char *no_out;
    if (asprintf(&no_out, "%s.reads_not_in_repeats.bed", output)< 0)
        errAbort("Preparing reads not in repeats wrong");
    no_stream = mustOpen(no_out, "w");
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

        //handle genome
        struct hashEl *hel = hashLookup(hashchr, sa->chr);
        if (hel != NULL) {
            struct chr *cs = (struct chr *) hel->val;
            cs->count_total++;
            for (i = sa->start; ; i++){
                if (i >= sa->end) {
                    break;
                }
                if (i >= cs->length) {
                    break;
                }
                (cs->bp_total)[i]++;
            }
        }

        //handle repeat
        struct binElement *hitList = NULL, *hit;
        struct hashEl *hel2 = hashLookup(hashbk, sa->chr);
        if (hel2 != NULL) {
            struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
            hitList = binKeeperFind(bs2, sa->start, sa->end);
            if (hitList != NULL ) {
                for ( hit = hitList; hit != NULL ; hit = hit->next) {
                    struct rmsk *ss = (struct rmsk *) (hit->val);
                    struct hashEl *hel3 = hashLookup(hashrep, ss->name);
                    if (hel3 != NULL) {
                        struct rep *rs = (struct rep *) hel3->val;
                        rs->count_total++;
                        rstart = sa->start - ss->start;
                        rstart = (rstart < 0) ? 0 : rstart;
                        rend = rstart + sa->length;
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
                    } else {
                        struct hashEl *hel4 = hashLookup(hashns, ss->name);
                        if (hel4 != NULL) {
                            struct repns *ns = (struct repns *) hel4->val;
                            ns->count_total++;
                        }        
                    }
                    break;
                }
                sam_in_repeat++;
                slFreeList(hitList);
            } else {
                fprintf(no_stream, "%s\t%d\t%d\t%s\t%d\t%c\n", sa->chr, sa->start, sa->end, sa->name, sa->qual, sa->strand);
            }
        }
        free(sa);
    }
    bam_destroy1(b);
    fclose(no_stream);
    samclose(samfp);

    freeMem(prn);   
    
    fprintf(stderr, "* Total %d reads.\n", reads_num);
    fprintf(stderr, "* Total %u reads mapped.\n", mapped_reads_num);
    fprintf(stderr, "* Total %d reads found overlaped with repeats.\n", sam_in_repeat);

    FILE *chr_out_stream;
    FILE *rep_out_stream;
    char *chr_out;
    char *rep_out;
    char *chr_out_bw;
    char *rep_out_bw;
    FILE *chr_stat_stream;
    FILE *rep_stat_stream;
    FILE *repns_stat_stream;
    char *chr_stat;
    char *rep_stat;
    char *repns_stat;

    if (asprintf(&chr_out, "%s.genome.wig", output) < 0)
        errAbort("Preparing genome wig wrong");
    if (asprintf(&rep_out, "%s.repeat.wig", output) < 0)
        errAbort("Preparing repeat wig wrong");
    if (asprintf(&chr_out_bw, "%s.genome.bigWig", output) < 0)
        errAbort("Preparing genome bigwig wrong");
    if (asprintf(&rep_out_bw, "%s.repeat.bigWig", output) < 0)
        errAbort("Preparing repeat bigwig wrong");
    if (asprintf(&chr_stat, "%s.genome.stat", output)< 0)
        errAbort("Preparing genome stat wrong");
    if (asprintf(&rep_stat, "%s.repeat.stat", output)< 0)
        errAbort("Preparing repeat stat wrong");
    if (asprintf(&repns_stat, "%s.repeat_no_size.stat", output)< 0)
        errAbort("Preparing repeat no size stat wrong");
    
    chr_out_stream = mustOpen(chr_out, "w");
    rep_out_stream = mustOpen(rep_out, "w");
    chr_stat_stream = mustOpen(chr_stat, "w");
    rep_stat_stream = mustOpen(rep_stat, "w");
    repns_stat_stream = mustOpen(repns_stat, "w");
    
    fprintf(stderr, "* Preparing the output wig and stat files...\n");

    struct hashEl *helc;
    struct hashCookie cookiec = hashFirst(hashchr);
    while ( (helc = hashNext(&cookiec)) != NULL ) {
        int zeroFlag, headFlag = 1;
        struct chr *oc = (struct chr *) (helc->val);
        fprintf(chr_stat_stream, "%s\t%d\n", oc->name, oc->count_total);
        for (k = 0; k < oc->length; k++) {
            if ((oc->bp_total)[k] == 0){
                zeroFlag = 1;
                headFlag = 1;
                continue;
            } else {
                zeroFlag = 0;
                if (headFlag) {
                    fprintf(chr_out_stream, "fixedStep chrom=%s start=%d step=1 span=1\n", oc->name, k + 1);
                    headFlag = 0;
                }
            }
            if ( zeroFlag == 0){
                fprintf(chr_out_stream, "%d\n", (oc->bp_total)[k]);
            }
        }
    }
    
    struct hashEl *helr;
    struct hashCookie cookier = hashFirst(hashrep);
    fprintf(rep_stat_stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#subfamily", "class", "family", "reads_count", "total_length", "genome_count", "RPKM");
    while ( (helr = hashNext(&cookier)) != NULL ) {
        struct rep *or = (struct rep *) (helr->val);
        if (or->genome_count != 0) {
            fprintf(rep_stat_stream, "%s\t%s\t%s\t%u\t%u\t%d\t%.3f\n", or->name, or->cname, or->fname, or->count_total, or->total_length, or->genome_count, cal_rpkm(or->count_total, or->total_length, mapped_reads_num));
        }
        fprintf(rep_out_stream, "fixedStep chrom=%s start=1 step=1 span=1\n", or->name);
        for (m = 0; m < or->length; m++) {
            fprintf(rep_out_stream, "%d\n", (or->bp_total)[m]);
        }
    }
    
    struct hashEl *helrs;
    struct hashCookie cookiers = hashFirst(hashns);
    fprintf(repns_stat_stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#subfamily", "class", "family", "reads_count", "total_length", "genome_count", "RPKM");
    while ( (helrs = hashNext(&cookiers)) != NULL ) {
        struct repns *os = (struct repns *) (helrs->val);
        fprintf(repns_stat_stream, "%s\t%s\t%s\t%u\t%u\t%d\t%.3f\n", os->name, os->cname, os->fname, os->count_total, os->total_length, os->genome_count, cal_rpkm(os->count_total, os->total_length, mapped_reads_num));
    }

    fclose(rep_out_stream);
    fclose(chr_out_stream);
    fclose(rep_stat_stream);
    fclose(repns_stat_stream);
    fclose(chr_stat_stream);

    freeHashAndVals(&hashchr);
    freeHashAndVals(&hashrep);
    freeHashAndVals(&hashns);
    freeHashAndVals(&hashbk);
    
    hashFree(&hashchr);
    hashFree(&hashrep);
    hashFree(&hashns);
    hashFree(&hashbk);
    
    fprintf(stderr, "* Generating bigWig files...\n");
    bigWigFileCreate(chr_out, chr_size_file, 256, 1024, 0, 1, chr_out_bw);
    bigWigFileCreate(rep_out, rep_size_file, 256, 1024, 0, 1, rep_out_bw);

    if (!arguments.keepWig) {
        unlink(chr_out);
        unlink(rep_out);
    }

    if(!arguments.output_base)
        free(output);
    
    end_time = time(NULL);

    fprintf(stderr, "* Done, time used %.1f seconds.\n", difftime(end_time, start_time));
    
    return 0;
}
