/*  Contains functions for extraneous functionality to the associated
 *  exercise. e.g. Printing PE placement, etc.
 *
 *  Controlled via macro definitions:.
 *    - C_SERIAL_PRACTICAL
 *    - F_SERIAL_PRACTICAL
 *    - C_MPI_PRACTICAL
 *    - F_MPI_PRACTICAL
 *    - C_OPENMP_PRACTICAL
 *    - F_OPENMP_PRACTICAL
 *    - C_HYBRID_PRACTICAL
 *    - F_HYBRID_PRACTICAL
 *    - C_OPENSHMEM_PRACTICAL
 *    - F_OPENSHMEM_PRACTICAL
 *
 *  Dominic Sloan-Murphy, EPCC, May 2014
 */
 
#define _GNU_SOURCE
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sched.h>

#if defined(C_HYBRID_PRACTICAL) || defined(F_HYBRID_PRACTICAL) || defined(C_MPI_PRACTICAL) || defined(F_MPI_PRACTICAL) || defined(C_OPENSHMEM_PRACTICAL) || defined(F_OPENSHMEM_PRACTICAL)
#include <mpi.h>
#endif /* C_HYBRID_PRACTICAL || F_HYBRID_PRACTICAL || C_MPI_PRACTICAL || F_MPI_PRACTICAL || C_OPENSHMEM_PRACTICAL || F_OPENSHMEM_PRACTICAL */

#if defined(C_OPENMP_PRACTICAL) || defined(F_OPENMP_PRACTICAL) || defined(C_HYBRID_PRACTICAL) || defined(F_HYBRID_PRACTICAL)
#include <omp.h>
#endif /* C_OPENMP_PRACTICAL || F_OPENMP_PRACTICAL || C_HYBRID_PRACTICAL || F_HYBRID_PRACTICAL */

/* Add a trailing underscore to functions so they can be called from Fortran. */
#if defined(F_OPENMP_PRACTICAL) || defined(F_HYBRID_PRACTICAL) || defined(F_MPI_PRACTICAL) || defined(F_OPENSHMEM_PRACTICAL) || defined(F_SERIAL_PRACTICAL)
#define printlocation printlocation_
#define wtime wtime_
#endif /* F_OPENMP_PRACTICAL || F_HYBRID_PRACTICAL || F_MPI_PRACTICAL || F_OPENSHMEM_PRACTICAL || F_SERIAL_PRACTICAL */

/* Borrowed from util-linux-2.13-pre7/schedutils/taskset.c */
static char *cpuset_to_cstr(cpu_set_t *mask, char *str)
{
    char *ptr = str;
    int i, j, entry_made = 0;
    for (i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, mask)) {
        int run = 0;
        entry_made = 1;
    for (j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, mask)) run++;
        else break;
        }
        if (!run)
        sprintf(ptr, "%d,", i);
        else if (run == 1) {
        sprintf(ptr, "%d,%d,", i, i + 1);
        i++;
        } else {
        sprintf(ptr, "%d-%d,", i, i + run);
        i += run;
        }
        while (*ptr != 0) ptr++;
    }
    }
    ptr -= entry_made;
    *ptr = 0;
    return(str);
}

void printlocation()
{

#if defined(C_OPENSHMEM_PRACTICAL) || defined(F_OPENSHMEM_PRACTICAL)
    MPI_Init(NULL, NULL);
#endif /* C_OPENSHMEM_PRACTICAL || F_OPENSHMEM_PRACTICAL */
      
#if defined(C_HYBRID_PRACTICAL) || defined(F_HYBRID_PRACTICAL) || defined(C_MPI_PRACTICAL) || defined(F_MPI_PRACTICAL) || defined(C_OPENSHMEM_PRACTICAL) || defined(F_OPENSHMEM_PRACTICAL)
    int rank, namelen;
    char hnbuf[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);    
    memset(hnbuf, 0, sizeof(hnbuf));
    MPI_Get_processor_name(hnbuf, &namelen);
#endif /* C_HYBRID_PRACTICAL || F_HYBRID_PRACTICAL || C_MPI_PRACTICAL || F_MPI_PRACTICAL || C_OPENSHMEM_PRACTICAL || F_OPENSHMEM_PRACTICAL */    

#if defined(C_OPENMP_PRACTICAL) || defined(F_OPENMP_PRACTICAL) || defined(C_HYBRID_PRACTICAL) || defined(F_HYBRID_PRACTICAL)
    int thread;
    thread = omp_get_thread_num();
#endif /* C_OPENMP_PRACTICAL || F_OPENMP_PRACTICAL || C_HYBRID_PRACTICAL || F_HYBRID_PRACTICAL */

    cpu_set_t coremask;
    char clbuf[7 * CPU_SETSIZE];
    memset(clbuf, 0, sizeof(clbuf));
    (void)sched_getaffinity(0, sizeof(coremask), &coremask);
    cpuset_to_cstr(&coremask, clbuf);

#if defined(C_HYBRID_PRACTICAL) || defined(F_HYBRID_PRACTICAL)
    printf("Rank %d / thread %d on core %s of node <%s>\n", rank, thread, clbuf, hnbuf); 
#elif defined(C_MPI_PRACTICAL) || defined(F_MPI_PRACTICAL) || defined(C_OPENSHMEM_PRACTICAL) || defined(F_OPENSHMEM_PRACTICAL)
    printf("Rank %d on core %s of node <%s>\n", rank, clbuf, hnbuf);
#elif defined(C_OPENMP_PRACTICAL) || defined(F_OPENMP_PRACTICAL)
    printf("Thread %d on core %s\n", thread, clbuf);
#else /* if defined(C_SERIAL_PRACTICAL) || defined(F_SERIAL_PRACTICAL) */
    printf("Program on core %s\n",clbuf);
#endif /* C_HYBRID_PRACTICAL || F_HYBRID_PRACTICAL */

#if defined(C_OPENSHMEM_PRACTICAL) || defined(F_OPENSHMEM_PRACTICAL)
    MPI_Finalize();
#endif /* C_OPENSHMEM_PRACTICAL || F_OPENSHMEM_PRACTICAL */

}

#if    defined(C_OPENSHMEM_PRACTICAL) || defined(F_OPENSHMEM_PRACTICAL) \
|| defined(C_SERIAL_PRACTICAL)    || defined(F_SERIAL_PRACTICAL)

/*
 * The following wall-clock code taken from the OpenSHMEM standards document
 */

#include <sys/time.h>

/* wall-clock time */

double wtime(void)
{
  struct timeval tp;
  gettimeofday (&tp, NULL);
  return tp.tv_sec + tp.tv_usec/(double)1.0e6;
}

#endif /*    C_OPENSHMEM_PRACTICAL || F_OPENSHMEM_PRACTICAL 
          || C_SERIAL_PRACTICAL    || F_SERIAL_PRACTICAL    */
