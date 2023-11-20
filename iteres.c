#include "generic.h"

static int usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: iteres (repeat analysis utils from Wang lab)\n");
    fprintf(stderr, "Version: %s\n\n", ITERES_VERSION);
    fprintf(stderr, "Usage:   iteres <command> [options]\n\n");
    fprintf(stderr, "Command: stat        get repeat alignment statistics\n");
    fprintf(stderr, "         filter      filter alignment statistic on repName/repFamily/repClass\n");
    //fprintf(stderr, "         nearby      fetch nearby genes for locations from bed file\n");
    fprintf(stderr, "         cpgstat     generate CpG density from MRE-Seq data for repeats\n");
    fprintf(stderr, "         cpgfilter   filter CpG statistic on repName/repFamily/repClass\n");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc < 2) return usage();
    if (strcmp(argv[1], "stat") == 0) return main_stat(argc-1, argv+1);
    else if (strcmp(argv[1], "filter") == 0) return main_filter(argc-1, argv+1);
    //else if (strcmp(argv[1], "nearby") == 0) return main_nearby(argc-1, argv+1);
    else if (strcmp(argv[1], "cpgstat") == 0) return main_cpgstat(argc-1, argv+1);
    else if (strcmp(argv[1], "cpgfilter") == 0) return main_cpgfilter(argc-1, argv+1);
    else {
        fprintf(stderr, "[iteres] unrecognized command '%s'\n", argv[1]);
        return 1;
    }
    return 0;
}
