#include <stdlib.h>        /* malloc */
#include <sys/time.h>      /* gettimeofday */
#include <sys/times.h>     /* times */
#include <unistd.h>        /* gettimeofday */
#include <stdio.h>
#include <string.h>        /* memset, strcmp (via STRMATCH) */
#include <assert.h>

#ifdef UNICOSMP
#include <intrinsics.h>    /* rtc */
#endif

#ifdef SPMD
#include "mpi.h"
#endif

#include "private.h"

/* Use hash table assist to linked list because it is faster */
#define HASH

static Timer **timers = 0;       /* linked list of timers */
static Timer **last = 0;         /* last element in list */
static int *max_depth;           /* maximum indentation level encountered */
static int *max_name_len;        /* max length of timer name */

typedef struct {
  unsigned int depth;            /* depth in calling tree */
  int padding[31];               /* padding is to mitigate false cache sharing */
} Nofalse; 
static Nofalse *current_depth;

static int nthreads    = -1;     /* num threads. Init to bad value */
static int maxthreads  = -1;     /* max threads (=nthreads for OMP). Init to bad value */
static int depthlimit  = 99999;  /* max depth for timers (99999 is effectively infinite) */
static bool initialized = false; /* GPTLinitialize has been called */

typedef struct {
  const Option option;           /* wall, cpu, etc. */
  const char *str;               /* descriptive string for printing */
  bool enabled;                  /* flag */
} Settings;

/* Options, print strings, and default enable flags */

static Settings cpustats =      {GPTLcpu,      "Usr       sys       usr+sys   ", false};
static Settings wallstats =     {GPTLwall,     "   Wallclock    max       min          ", true };
static Settings overheadstats = {GPTLoverhead, "Overhead  ",                     true };

static const int tablesize = 128*MAX_CHARS;  /* 128 is size of ASCII char set */
static Hashentry **hashtable;    /* table of entries hashed by sum of chars */
static unsigned int *novfl;      /* microsecond overflow counter (only when DIAG set) */
static long ticks_per_sec;       /* clock ticks per second */

#ifdef UNICOSMP
static long long ticks_per_secI = -1;
#endif

/* Local function prototypes */

static void printstats (const Timer *, FILE *, const int, const bool);
static void add (Timer *, const Timer *);
static int get_cpustamp (long *, long *);
static Timer *getentry (const Hashentry *, const char *, int *);

/*
** GPTLsetoption: set option value to true or false.
**
** Input arguments:
**   option: option to be set
**   val:    value to which option should be set (for boolean values, 
**           nonzero=true, zero=false)
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLsetoption (const int option,  /* option */
		   const int val)     /* whether to enable */
{

#ifdef DISABLE_TIMERS
  return 0;
#endif

  if (initialized)
    return GPTLerror ("GPTLsetoption: must be called BEFORE GPTLinitialize\n");

  if (option == GPTLabort_on_error) {
    GPTLset_abort_on_error ((bool) val);
    printf ("GPTLsetoption: setting abort on error flag to %d\n", val);
    return 0;
  }

  switch (option) {
  case GPTLcpu:      
    cpustats.enabled = (bool) val; 
    printf ("GPTLsetoption: set cpustats to %d\n", val);
    return 0;
  case GPTLwall:     
    wallstats.enabled = (bool) val; 
    printf ("GPTLsetoption: set wallstats to %d\n", val);
    return 0;
  case GPTLoverhead: 
    overheadstats.enabled = (bool) val; 
    printf ("GPTLsetoption: set overheadstats to %d\n", val);
    return 0;
  case GPTLdepthlimit: 
    depthlimit = val; 
    printf ("GPTLsetoption: set depthlimit to %d\n", val);
    return 0;
  default:
    break;
  }

#ifdef HAVE_PAPI
  if (GPTL_PAPIsetoption (option, val) == 0)
    return 0;
#endif
  return GPTLerror ("GPTLsetoption: option %d not found\n", option);
}

/*
** GPTLinitialize (): Initialization routine must be called from single-threaded
**   region before any other timing routines may be called.  The need for this
**   routine could be eliminated if not targetting timing library for threaded
**   capability. 
**
** return value: 0 (success) or GPTLerror (failure)
*/

int GPTLinitialize (void)
{
  int i;          /* loop index */
  int t;          /* thread index */
#ifdef UNICOSMP
  extern long long rtc_rate_();
#endif

#ifdef DISABLE_TIMERS
  return 0;
#endif

  if (initialized)
    return GPTLerror ("GPTLinitialize: has already been called\n");

  if (threadinit (&nthreads, &maxthreads) < 0)
    return GPTLerror ("GPTLinitialize: bad return from threadinit\n");

  if (get_thread_num (&nthreads, &maxthreads) > 0) 
    return GPTLerror ("GPTLinitialize: must only be called by master thread\n");

#ifdef UNICOSMP
  ticks_per_secI = rtc_rate_();
#else
#ifndef SPMD
  if ((ticks_per_sec = sysconf (_SC_CLK_TCK)) == -1)
    return GPTLerror ("GPTLinitialize: sysconf (_SC_CLK_TCK) failed\n");
#endif
#endif

  /* Allocate space for global arrays */

  timers        = (Timer **)     GPTLallocate (maxthreads * sizeof (Timer *));
  last          = (Timer **)     GPTLallocate (maxthreads * sizeof (Timer *));
  current_depth = (Nofalse *)    GPTLallocate (maxthreads * sizeof (Nofalse));
  max_depth     = (int *)        GPTLallocate (maxthreads * sizeof (int));
  max_name_len  = (int *)        GPTLallocate (maxthreads * sizeof (int));
  hashtable     = (Hashentry **) GPTLallocate (maxthreads * sizeof (Hashentry *));
#ifdef DIAG
  novfl         = (unsigned int *) GPTLallocate (maxthreads * sizeof (unsigned int));
#endif

  /* Initialize array values */

  for (t = 0; t < maxthreads; t++) {
    timers[t] = 0;
    current_depth[t].depth = 0;
    max_depth[t]     = 0;
    max_name_len[t]  = 0;
    hashtable[t] = (Hashentry *) GPTLallocate (tablesize * sizeof (Hashentry));
#ifdef DIAG
    novfl[t] = 0;
#endif
    for (i = 0; i < tablesize; i++) {
      hashtable[t][i].nument = 0;
      hashtable[t][i].entries = 0;
    }
  }

#ifdef HAVE_PAPI
  if (GPTL_PAPIinitialize (maxthreads) < 0)
    return GPTLerror ("GPTLinitialize: GPTL_PAPIinitialize failure\n");
#endif

  initialized = true;
  return 0;
}

/*
** GPTLfinalize (): Finalization routine must be called from single-threaded
**   region. Free all malloc'd space
**
** return value: 0 (success) or GPTLerror (failure)
*/

int GPTLfinalize (void)
{
  int t;                /* thread index */
  Timer *ptr, *ptrnext; /* ll indices */

#ifdef DISABLE_TIMERS
  return 0;
#endif

  if ( ! initialized)
    return GPTLerror ("GPTLfinalize: GPTLinitialize() has not been called\n");

  if (get_thread_num (&nthreads, &maxthreads) > 0) 
    return GPTLerror ("GPTLfinalize: must only be called by master thread\n");

  for (t = 0; t < maxthreads; ++t) {
    free (hashtable[t]);
    for (ptr = timers[t]; ptr; ptr = ptrnext) {
      ptrnext = ptr->next;
      free (ptr);
    }
  }

  free (timers);
  free (current_depth);
  free (max_depth);
  free (max_name_len);
  free (hashtable);
#ifdef DIAG
  free (novfl);
#endif

  threadfinalize ();
#ifdef HAVE_PAPI
  GPTL_PAPIfinalize (maxthreads);
#endif
  initialized = false;
  return 0;
}

/*
** GPTLstart: start a timer
**
** Input arguments:
**   name: timer name
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLstart (const char *name)       /* timer name */
{
#ifdef SPMD
  double wtime1, wtime2;        /* returned from MPI_Wtime() */
#else
  struct timeval tp1, tp2;      /* argument returned from gettimeofday */
#endif

  Timer *ptr;                   /* linked list pointer */
  Timer **eptr;                 /* for realloc */

  int nchars;                   /* number of characters in timer */
  int t;                        /* thread index (of this thread) */
  int indx;                     /* hash table index */
  int nument;                   /* number of entries for a hash collision */

#ifdef UNICOSMP
  long long nticks1, nticks2;        /* returned from rtc() */

#ifndef SSP
  if (__streaming() == 0) return 0;  /* timers don't work in this situation so disable */
#endif
#endif

#ifdef DISABLE_TIMERS
  return 0;
#endif

  if ((t = get_thread_num (&nthreads, &maxthreads)) < 0)
    return GPTLerror ("GPTLstart\n");

  /*
  ** If current depth exceeds a user-specified limit for print, just
  ** increment and return
  */

  if (current_depth[t].depth >= depthlimit) {
    ++current_depth[t].depth;
    return 0;
  }

  /* 1st calls to overheadstart and gettimeofday are solely for overhead timing */

  if (overheadstats.enabled) {
#ifdef HAVE_PAPI
    (void) GPTL_PAPIoverheadstart (t);
#endif

    if (wallstats.enabled)
#ifdef UNICOSMP
      nticks1 = _rtc();
#else
#ifdef SPMD
      wtime1 = MPI_Wtime();
#else
      gettimeofday (&tp1, 0);
#endif
#endif
  }

  if ( ! initialized)
    return GPTLerror ("GPTLstart: GPTLinitialize has not been called\n");

  /* Look for the requested timer in the current list. */

#ifdef HASH
  ptr = getentry (hashtable[t], name, &indx);
  assert (indx < tablesize);
#else
  for (ptr = timers[t]; ptr && ! STRMATCH (name, ptr->name); ptr = ptr->next);
#endif

  if (ptr && ptr->onflg)
    return GPTLerror ("GPTLstart thread %d: timer %s was already on: "
		      "not restarting.\n", t, ptr->name);

  ++current_depth[t].depth;
  if (current_depth[t].depth > max_depth[t])
    max_depth[t] = current_depth[t].depth;

  if (ptr) {

    /*
    ** Reset indentation level to ambiguous value if inconsistent with
    ** current value. This will likely happen when the thing being timed is
    ** called from more than 1 branch in the call tree.
    */
  
    if (ptr->depth != current_depth[t].depth)
      ptr->depth = 0;

  } else {

    /* Add a new entry and initialize */

    ptr = (Timer *) GPTLallocate (sizeof (Timer));
    memset (ptr, 0, sizeof (Timer));

    /* Truncate input name if longer than MAX_CHARS characters  */

    nchars = MIN (strlen (name), MAX_CHARS);
    max_name_len[t] = MAX (nchars, max_name_len[t]);

    strncpy (ptr->name, name, nchars);
    ptr->name[nchars] = '\0';
    ptr->depth = current_depth[t].depth;

    if (timers[t])
      last[t]->next = ptr;
    else
      timers[t] = ptr;

    last[t] = ptr;
#ifdef HASH
    ++hashtable[t][indx].nument;
    nument = hashtable[t][indx].nument;

    eptr = realloc (hashtable[t][indx].entries, nument * sizeof (Timer *));
    if ( ! eptr)
      return GPTLerror ("GPTLstart: realloc error\n");

    hashtable[t][indx].entries = eptr;						 
    hashtable[t][indx].entries[nument-1] = ptr;
#endif
  }

  ptr->onflg = true;

  if (cpustats.enabled && get_cpustamp (&ptr->cpu.last_utime, &ptr->cpu.last_stime) < 0)
    return GPTLerror ("GPTLstart: get_cpustamp error");
  
  /*
  ** The 2nd system timer call is used both for overhead estimation and
  ** the input timer.  Under UNICOSMP, "ticks" are recorded by rtc instead 
  ** of sec as returned by MPI_Wtime or sec and usec as returned by gettimeofday(). 
  */
  
  if (wallstats.enabled) {
#ifdef UNICOSMP
    nticks2               = _rtc ();
    ptr->wall.last_nticks = nticks2; /* just store rtc to be used and converted later */
    if (overheadstats.enabled)
      ptr->wall.overhead += (float) (nticks2 - nticks1) / (float) ticks_per_secI;
#else
#ifdef SPMD
    wtime2               = MPI_Wtime();
    ptr->wall.last_sec   = wtime2; 
    if (overheadstats.enabled)
      ptr->wall.overhead += wtime2 - wtime1;
#else
    gettimeofday (&tp2, 0);
    ptr->wall.last_sec  = tp2.tv_sec;
    ptr->wall.last_usec = tp2.tv_usec;
    if (overheadstats.enabled)
      ptr->wall.overhead +=       (tp2.tv_sec  - tp1.tv_sec) + 
	                    1.e-6*(tp2.tv_usec - tp1.tv_usec);
#endif
#endif
  }

#ifdef HAVE_PAPI
  if (GPTL_PAPIstart (t, &ptr->aux) < 0)
    return GPTLerror ("GPTLstart: error from GPTL_PAPIstart\n");

  if (overheadstats.enabled)
    (void) GPTL_PAPIoverheadstop (t, &ptr->aux);
#endif

  return (0);
}

/*
** GPTLstop: stop a timer
**
** Input arguments:
**   name: timer name
**
** Return value: 0 (success) or -1 (failure)
*/

int GPTLstop (const char *name) /* timer name */
{
#ifdef SPMD
  double wtime1, wtime2;    /* returned from MPI_Wtime() */
#else
  struct timeval tp1, tp2;  /* argument to gettimeofday() */
  long delta_wtime_sec;     /* wallclock sec change fm GPTLstart() to GPTLstop() */    
  long delta_wtime_usec;    /* wallclock usec change fm GPTLstart() to GPTLstop() */
#endif
  float delta_wtime;        /* floating point wallclock change */

  Timer *ptr;               /* linked list pointer */

  int t;                    /* thread number for this process */
  int indx;                 /* index into hash table */

  long usr;                 /* user time (returned from get_cpustamp) */
  long sys;                 /* system time (returned from get_cpustamp) */

#ifdef UNICOSMP
  long long nticks1, nticks2;        /* returned from rtc() */
  long long delta_nticks;            /* change in tick count */

#ifndef SSP
  if (__streaming() == 0) return 0;  /* timers don't work in this situation so disable */
#endif
#endif

#ifdef DISABLE_TIMERS
  return 0;
#endif

  if ((t = get_thread_num (&nthreads, &maxthreads)) < 0)
    return GPTLerror ("GPTLstop\n");

  /*
  ** If current depth exceeds a user-specified limit for print, just
  ** decrement and return
  */

  if (current_depth[t].depth > depthlimit) {
    --current_depth[t].depth;
    return 0;
  }

#ifdef HAVE_PAPI
  if (overheadstats.enabled)
    (void) GPTL_PAPIoverheadstart (t);
#endif

  /*
  ** The 1st system timer call is used both for overhead estimation and
  ** the input timer
  */
    
  if (wallstats.enabled)
#ifdef UNICOSMP
    nticks1 = _rtc ();
#else
#ifdef SPMD
    wtime1 = MPI_Wtime();
#else
    gettimeofday (&tp1, 0);
#endif
#endif

  if (cpustats.enabled && get_cpustamp (&usr, &sys) < 0)
    return GPTLerror (0);

  if ( ! initialized)
    return GPTLerror ("GPTLstop: GPTLinitialize has not been called\n");

#ifdef HASH
  ptr = getentry (hashtable[t], name, &indx);
#else
  for (ptr = timers[t]; ptr && ! STRMATCH (name, ptr->name); ptr = ptr->next);
#endif

  if ( ! ptr) 
    return GPTLerror ("GPTLstop: timer for %s had not been started.\n", name);

  if ( ! ptr->onflg )
    return GPTLerror ("GPTLstop: timer %s was already off.\n",ptr->name);

#ifdef HAVE_PAPI
  if (GPTL_PAPIstop (t, &ptr->aux) < 0)
    return GPTLerror ("GPTLstart: error from GPTL_PAPIstop\n");
#endif

  --current_depth[t].depth;

  ptr->onflg = false;
  ptr->count++;

  if (wallstats.enabled) {

    /*
    ** Under UNICOSMP, rtc() is being used instead of gettimeofday().
    ** Under SPMD, MPI_Wtime is being used instead of gettimeofday().
    ** In both cases the usec component is not needed.
    */

#ifdef UNICOSMP
    delta_nticks            = nticks1 - ptr->wall.last_nticks;
    delta_wtime             = (float) delta_nticks / (float) ticks_per_secI;
    ptr->wall.accum_nticks += delta_nticks;
#else
#ifdef SPMD
    delta_wtime             = wtime1 - ptr->wall.last_sec;
    ptr->wall.accum_sec    += delta_wtime;
    ptr->wall.accum_usec    = 0;
#else
    delta_wtime_sec         = tp1.tv_sec  - ptr->wall.last_sec;
    delta_wtime_usec        = tp1.tv_usec - ptr->wall.last_usec;
    delta_wtime             = delta_wtime_sec + 1.e-6*delta_wtime_usec;
    ptr->wall.accum_sec    += delta_wtime_sec;
    ptr->wall.accum_usec   += delta_wtime_usec;
#endif
#endif

    if (ptr->count == 1) {
      ptr->wall.max = delta_wtime;
      ptr->wall.min = delta_wtime;
    } else {
      ptr->wall.max = MAX (ptr->wall.max, delta_wtime);
      ptr->wall.min = MIN (ptr->wall.min, delta_wtime);
    }

    /*
    ** Adjust accumulated wallclock values to guard against overflow in the
    ** microsecond accumulator.
    */

    if (ptr->wall.accum_usec > 10000000) {
      ptr->wall.accum_sec  += 10;
      ptr->wall.accum_usec -= 10000000;
#ifdef DIAG
      ++novfl[t];
#endif
    } else if (ptr->wall.accum_usec < -10000000) {
      ptr->wall.accum_sec  -= 10;
      ptr->wall.accum_usec += 10000000;
#ifdef DIAG
      ++novfl[t];
#endif
    }

    /* 2nd system timer call is solely for overhead timing */

    if (overheadstats.enabled) {
#ifdef UNICOSMP
      nticks2 = _rtc ();
      ptr->wall.overhead += (float) (nticks2 - nticks1) / (float) ticks_per_secI;
#else
#ifdef SPMD
      wtime2 = MPI_Wtime();
      ptr->wall.overhead += (wtime2 - wtime1);
#else
      gettimeofday (&tp2, 0);
      ptr->wall.overhead +=       (tp2.tv_sec  - tp1.tv_sec) + 
	                    1.e-6*(tp2.tv_usec - tp1.tv_usec);
#endif
#endif
    }
  }

  if (cpustats.enabled) {
    ptr->cpu.accum_utime += usr - ptr->cpu.last_utime;
    ptr->cpu.accum_stime += sys - ptr->cpu.last_stime;
    ptr->cpu.last_utime   = usr;
    ptr->cpu.last_stime   = sys;
  }

#ifdef HAVE_PAPI
  if (overheadstats.enabled)
    (void) GPTL_PAPIoverheadstop (t, &ptr->aux);
#endif

  return 0;
}

/*
** GPTLstamp: Compute timestamp of usr, sys, and wallclock time (seconds)
**
** Output arguments:
**   wall: wallclock
**   usr:  user time
**   sys:  system time
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLstamp (double *wall, double *usr, double *sys)
{
  struct timeval tp;         /* argument to gettimeofday */
  struct tms buf;            /* argument to times */

#ifdef DISABLE_TIMERS
  return 0;
#endif

  if ( ! initialized)
    return GPTLerror ("GPTLstamp: GPTLinitialize has not been called\n");

  *usr = 0;
  *sys = 0;

  if (times (&buf) == -1)
    return GPTLerror ("GPTLstamp: times() failed. Results bogus\n");

  *usr = buf.tms_utime / (double) ticks_per_sec;
  *sys = buf.tms_stime / (double) ticks_per_sec;

  gettimeofday (&tp, 0);
  *wall = tp.tv_sec + 1.e-6*tp.tv_usec;

  return 0;
}

/*
** GPTLreset: reset all known timers to 0
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLreset (void)
{
  int t;             /* index over threads */
  Timer *ptr;        /* linked list index */

#ifdef DISABLE_TIMERS
  return 0;
#endif

  if ( ! initialized)
    return GPTLerror ("GPTLreset: GPTLinitialize has not been called\n");

  /* Only allow the master thread to reset timers */

  if (get_thread_num (&nthreads, &maxthreads) > 0) 
    return GPTLerror ("GPTLreset: must only be called by master thread\n");

  for (t = 0; t < nthreads; t++) {
    for (ptr = timers[t]; ptr; ptr = ptr->next) {
      ptr->onflg = false;
      ptr->count = 0;
      memset (&ptr->wall, 0, sizeof (ptr->wall));
      memset (&ptr->cpu, 0, sizeof (ptr->cpu));
#ifdef HAVE_PAPI
      memset (&ptr->aux, 0, sizeof (ptr->aux));
#endif
    }
  }
  printf ("GPTLreset: accumulators for all timers set to zero\n");
  return 0;
}

/* 
** GPTLpr: Print values of all timers
**
** Input arguments:
**   id: integer to append to string "timing."
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLpr (const int id)   /* output file will be named "timing.<id>" */
{
  FILE *fp;                /* file handle to write to */
  Timer *ptr;              /* walk through master thread linked list */
  Timer *tptr;             /* walk through slave threads linked lists */
  Timer sumstats;          /* sum of same timer stats over threads */
  int i, ii, n, t;         /* indices */
  char outfile[12];        /* name of output file: timing.xxxx */
  float *sum;              /* sum of overhead values (per thread) */
  float osum;              /* sum of overhead over threads */
  bool found;              /* jump out of loop when name found */
  bool foundany;           /* whether summation print necessary */
  bool first;              /* flag 1st time entry found */

#ifdef DISABLE_TIMERS
  return 0;
#endif

  if ( ! initialized)
    return GPTLerror ("GPTLpr: GPTLinitialize() has not been called\n");

  if (get_thread_num (&nthreads, &maxthreads) > 0) 
    return GPTLerror ("GPTLpr: must only be called by master thread\n");

  if (id < 0 || id > 9999)
    return GPTLerror ("GPTLpr: bad id=%d for output file. Must be >= 0 and < 10000\n", id);

  sprintf (outfile, "timing.%d", id);

  if ( ! (fp = fopen (outfile, "w")))
    fp = stderr;

  sum = GPTLallocate (nthreads * sizeof (float));

  for (t = 0; t < nthreads; ++t) {
    if (t > 0)
      fprintf (fp, "\n");
    fprintf (fp, "Stats for thread %d:\n", t);

    for (n = 0; n < max_depth[t]; ++n)    /* max indent level (depth starts at 1) */
      fprintf (fp, "  ");
    for (n = 0; n < max_name_len[t]; ++n) /* longest timer name */
      fprintf (fp, " ");

    fprintf (fp, "Called   ");

    /* Print strings for enabled timer types */

    if (cpustats.enabled)
      fprintf (fp, "%s", cpustats.str);
    if (wallstats.enabled) {
      fprintf (fp, "%s", wallstats.str);
      if (overheadstats.enabled)
	fprintf (fp, "%s", overheadstats.str);
    }

#ifdef HAVE_PAPI
    GPTL_PAPIprstr (fp);
#endif

    fprintf (fp, "\n");        /* Done with titles, go to next line */

    for (ptr = timers[t]; ptr; ptr = ptr->next)
      printstats (ptr, fp, t, true);

    /* Sum of overhead across timers is meaningful */

    if (wallstats.enabled && overheadstats.enabled) {
      sum[t] = 0;
      for (ptr = timers[t]; ptr; ptr = ptr->next)
	sum[t] += ptr->wall.overhead;
      fprintf (fp, "Overhead sum = %9.3f wallclock seconds\n", sum[t]);
    }
  }

  /* Print per-name stats for all threads */

  if (nthreads > 1) {
    fprintf (fp, "\nSame stats sorted by timer for threaded regions:\n");
    fprintf (fp, "Thd ");

    for (n = 0; n < max_name_len[0]; ++n) /* longest timer name */
      fprintf (fp, " ");

    fprintf (fp, "Called   ");

    if (cpustats.enabled)
      fprintf (fp, "%s", cpustats.str);
    if (wallstats.enabled) {
      fprintf (fp, "%s", wallstats.str);
      if (overheadstats.enabled)
	fprintf (fp, "%s", overheadstats.str);
    }

#ifdef HAVE_PAPI
    GPTL_PAPIprstr (fp);
#endif

    fprintf (fp, "\n");

    for (ptr = timers[0]; ptr; ptr = ptr->next) {
      
      /* 
      ** To print sum stats, first create a new timer then copy thread 0
      ** stats into it. then sum using "add", and finally print.
      */

      foundany = false;
      first = true;
      sumstats = *ptr;
      for (t = 1; t < nthreads; ++t) {
	found = false;
	for (tptr = timers[t]; tptr && ! found; tptr = tptr->next) {
	  if (STRMATCH (ptr->name, tptr->name)) {

	    /* Only print thread 0 when this timer found for other threads */

	    if (first) {
	      first = false;
	      fprintf (fp, "%3.3d ", 0);
	      printstats (ptr, fp, 0, false);
	    }

	    found = true;
	    foundany = true;
	    fprintf (fp, "%3.3d ", t);
	    printstats (tptr, fp, 0, false);
	    add (&sumstats, tptr);
	  }
	}
      }

      if (foundany) {
	fprintf (fp, "SUM ");
	printstats (&sumstats, fp, 0, false);
	fprintf (fp, "\n");
      }
    }

    /* Repeat overhead print in loop over threads */

    if (wallstats.enabled && overheadstats.enabled) {
      osum = 0.;
      for (t = 0; t < nthreads; ++t) {
	fprintf (fp, "OVERHEAD.%3.3d (wallclock seconds) = %9.3f\n", t, sum[t]);
	osum += sum[t];
      }
      fprintf (fp, "OVERHEAD.SUM (wallclock seconds) = %9.3f\n", osum);
    }
  }

#ifdef DIAG
  fprintf (fp, "\n");
  for (t = 0; t < nthreads; ++t) 
    fprintf (fp, "novfl[%d]=%d\n", t, novfl[t]);
#endif

  /* Print hash table stats */

#ifdef HASH  
  for (t = 0; t < nthreads; t++) {
    first = true;
    for (i = 0; i < tablesize; i++) {
      int nument = hashtable[t][i].nument;
      if (nument > 1) {
	if (first) {
	  first = false;
	  fprintf (fp, "\nthread %d had some hash collisions:\n", t);
	}
	fprintf (fp, "hashtable[%d][%d] had %d entries:", t, i, nument);
	for (ii = 0; ii < nument; ii++)
	  fprintf (fp, " %s", hashtable[t][i].entries[ii]->name);
	fprintf (fp, "\n");
      }
    }
  }
#endif

  fclose (fp);
  free (sum);
  return 0;
}

/* 
** printstats: print a single timer
**
** Input arguments:
**   timer:    timer for which to print stats
**   fp:       file descriptor to write to
**   t:        thread number
**   doindent: whether to indent
*/

static void printstats (const Timer *timer,     /* timer to print */
			FILE *fp,               /* file descriptor to write to */
			const int t,            /* thread number */
			const bool doindent)    /* whether indenting will be done */
{
  int i;               /* index */
  int indent;          /* index for indenting */
  int extraspace;      /* for padding to length of longest name */
  long ticks_per_sec;  /* returned from sysconf */
  float usr;           /* user time */
  float sys;           /* system time */
  float usrsys;        /* usr + sys */
  float elapse;        /* elapsed time */

  if ((ticks_per_sec = sysconf (_SC_CLK_TCK)) == -1)
    (void) GPTLerror ("printstats: token _SC_CLK_TCK is not defined\n");

  /* Indent to depth of this timer */

  if (doindent)
    for (indent = 0; indent < timer->depth; ++indent)  /* depth starts at 1 */
      fprintf (fp, "  ");

  fprintf (fp, "%s", timer->name);

  /* Pad to length of longest name */

  extraspace = max_name_len[t] - strlen (timer->name);
  for (i = 0; i < extraspace; ++i)
    fprintf (fp, " ");

  /* Pad to max indent level */

  if (doindent)
    for (indent = timer->depth; indent < max_depth[t]; ++indent)
      fprintf (fp, "  ");

  fprintf (fp, "%8ld ", timer->count);

  if (cpustats.enabled) {
    usr = timer->cpu.accum_utime / (float) ticks_per_sec;
    sys = timer->cpu.accum_stime / (float) ticks_per_sec;
    usrsys = usr + sys;
    fprintf (fp, "%9.3f %9.3f %9.3f ", usr, sys, usrsys);
  }

  if (wallstats.enabled) {
#ifdef UNICOSMP
    elapse = timer->wall.accum_nticks / (float) ticks_per_secI;
#else
#ifdef SPMD
    elapse = timer->wall.accum_sec;
#else
    elapse = timer->wall.accum_sec + 1.e-6*timer->wall.accum_usec;
#endif
#endif
    fprintf (fp, "%12.6f %12.6f %12.6f ", elapse, timer->wall.max, timer->wall.min);
    if (overheadstats.enabled)
      fprintf (fp, "%12.6f ", timer->wall.overhead);
  }

#ifdef HAVE_PAPI
  GPTL_PAPIpr (fp, &timer->aux);
#endif

  fprintf (fp, "\n");
}

/* 
** add: add the contents of tin to tout
**
** Input arguments:
**   tin:  input timer
** Input/output arguments:
**   tout: output timer summed into
*/

static void add (Timer *tout,   
		 const Timer *tin)
{
  if (wallstats.enabled) {
    tout->count           += tin->count;
    tout->wall.accum_sec  += tin->wall.accum_sec;
    tout->wall.accum_usec += tin->wall.accum_usec;
    if (tout->wall.accum_usec > 10000000) {
      tout->wall.accum_sec  += 10;
      tout->wall.accum_usec -= 10000000;
    } else if (tout->wall.accum_usec < -10000000) {
      tout->wall.accum_sec  -= 10;
      tout->wall.accum_usec += 10000000;
    }
      
    tout->wall.max = MAX (tout->wall.max, tin->wall.max);
    tout->wall.min = MIN (tout->wall.min, tin->wall.min);
    if (overheadstats.enabled)
      tout->wall.overhead += tin->wall.overhead;
  }

  if (cpustats.enabled) {
    tout->cpu.accum_utime += tin->cpu.accum_utime;
    tout->cpu.accum_stime += tin->cpu.accum_stime;
  }
#ifdef HAVE_PAPI
  GPTL_PAPIadd (&tout->aux, &tin->aux);
#endif
}

/*
** get_cpustamp: Invoke the proper system timer and return stats.
**
** Output arguments:
**   usr: user time
**   sys: system time
**
** Return value: 0 (success)
*/

static int get_cpustamp (long *usr, long *sys)
{
  struct tms buf;

  (void) times (&buf);
  *usr = buf.tms_utime;
  *sys = buf.tms_stime;

  return 0;
}

/*
** getentry: find the entry in the hash table and return a pointer to it.
**
** Input args:
**   hashtable: the hashtable (array)
**   name:      string to be hashed on (specifically, summed)
** Output args:
**   indx:      hashtable index
**
** Return value: pointer to the entry, or NULL if not found
*/

static Timer *getentry (const Hashentry *hashtable, /* hash table */
			const char *name,           /* name to hash */
			int *indx)                  /* hash index */
{
  int i;                 /* loop index */
  const char *c = name;  /* pointer to elements of "name" */

  /* Generate the hash value by summing values of the chars in "name" */

  for (*indx = 0; *c; c++)
    *indx += *c;

  /* 
  ** If nument exceeds 1 there was a hash collision and we must search
  ** linearly through an array for a match
  */

  for (i = 0; i < hashtable[*indx].nument; i++)
    if (STRMATCH (name, hashtable[*indx].entries[i]->name))
      return hashtable[*indx].entries[i];

  return 0;
}
