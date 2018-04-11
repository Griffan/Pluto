#include <iostream>
#include <cstring>
#include "SinglePhasing/thunderVCF/PhasingMain.h"
using namespace std;

int usage()
{
    fprintf(stderr,"Usage:\n");
    fprintf(stderr,"\tPluto index\n");
    fprintf(stderr,"\tPluto phase\n");
    fprintf(stderr,"\tPluto iterative\n");
}

int main(int argc, char ** argv) {

    if (argc < 2) return usage();
    if (strcmp(argv[1], "iterative") == 0) return PhasingMain(argc - 1, argv + 1);
    else if (strcmp(argv[1], "index") == 0) return BuildGraph(argc - 1, argv + 1);
    else if (strcmp(argv[1], "phase") == 0) return PhaseByRefGraph(argc - 1, argv + 1);
    else {
        fprintf(stderr,"unrecognized command '%s'\n", argv[1]);
        return 1;
    }
}
