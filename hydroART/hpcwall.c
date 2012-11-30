#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void hpcwall (double *retval)
{
    static long zsec = 0;
    static long zusec = 0;
    double esec;
    struct timeval tp;
    struct timezone tzp;

    gettimeofday (&tp, &tzp);

    if ( zsec == 0 ) zsec = tp.tv_sec;
    if ( zusec == 0 ) zusec = tp.tv_usec;

    *retval = (tp.tv_sec - zsec) + (tp.tv_usec - zusec) * 0.000001;
}

/* for different function name convention */
void hpcwall_(double *retval) { hpcwall(retval); }
