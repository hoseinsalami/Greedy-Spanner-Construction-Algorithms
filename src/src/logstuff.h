#include "structures.h"
#include "memory.h"

#ifndef COMMANDLINE
#define STARTLOG
#define ENDLOG
#endif

#ifdef COMMANDLINE
    #define STARTLOG \
        Timing::reset(); \
        T->Log("Staring new spanner calculation"); \
        T->Log("using seed: ",T->seed);

    #define ENDLOG \
        T->Log("#edges: ",T->edgeCount); \
        T->LogDouble("longest edge: ",T->longestEdge); \
        double vm2, rss2; \
        process_mem_usage(vm2, rss2); \
        T->LogDouble("virtual memory size: ",vm2); \
        T->LogDouble("Resident memory set: ",rss2);
#endif
