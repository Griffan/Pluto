//
// Created by Fan Zhang on 8/18/15.
//

#ifndef PLUTO_PHASINGMAIN_H
#define PLUTO_PHASINGMAIN_H
namespace BuildGraph{
int BuildGraph(int argc, char **argv) ;
}
namespace ReadGraph{
int PhaseByRefGraph(int argc, char **argv);
}
namespace PhaseIntersect {
    int PhasingMain(int argc, char **argv);
}
#endif //PLUTO_PHASINGMAIN_H
