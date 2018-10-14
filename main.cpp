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
    return 0;
}

int main(int argc, char ** argv) {

    if (argc < 2) return usage();
    if (strcmp(argv[1], "iterative") == 0) return PhaseIntersect::PhasingMain(argc - 1, argv + 1);
    else if (strcmp(argv[1], "index") == 0) return BuildGraph::BuildGraph(argc - 1, argv + 1);
    else if (strcmp(argv[1], "phase") == 0) return ReadGraph::PhaseByRefGraph(argc - 1, argv + 1);
    else {
        fprintf(stderr, "unrecognized command '%s'\n", argv[1]);
        return 1;
    }
//
//    PhaseIntersect::PhasingMain(argc - 1, argv + 1);
////    BuildGraph::BuildGraph(argc - 1, argv + 1);
////    ReadGraph::PhaseByRefGraph(argc - 1, argv + 1);

}
