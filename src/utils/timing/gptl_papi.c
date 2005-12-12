#if ( defined HAVE_PAPI )
#include <papi.h>
#endif
#include <stdlib.h>

#if ( defined THREADED_OMP )
#include <omp.h>
#elif ( defined THREADED_PTHREADS )
#include <pthread.h>
#endif

#include "private.h"

typedef struct {
  int counter;      /* PAPI counter */
  char *prstr;      /* print string for output timers (16 chars) */
  char *str;        /* descriptive print string (more descriptive than prstr) */
} Papientry;

#if ( defined HAVE_PAPI )

/* Mapping of PAPI counters to short and long printed strings */

static Papientry papitable [] = {
  {PAPI_L1_DCM, "L1 Dcache miss  ", "Level 1 data cache misses"},
  {PAPI_L1_ICM, "L1 Icache miss  ", "Level 1 instruction cache misses"},
  {PAPI_L2_DCM, "L2 Dcache miss  ", "Level 2 data cache misses"},             
  {PAPI_L2_ICM, "L2 Icache miss  ", "Level 2 instruction cache misses"},
  {PAPI_L3_DCM, "L3 Dcache miss  ", "Level 3 data cache misses"},             
  {PAPI_L3_ICM, "L3 Icache miss  ", "Level 3 instruction cache misses"},       
  {PAPI_L1_TCM, "L1 cache miss   ", "Level 1 total cache misses"},             
  {PAPI_L2_TCM, "L2 cache miss   ", "Level 2 total cache misses"},
  {PAPI_L3_TCM, "L3 cache miss   ", "Level 3 total cache misses"},
  {PAPI_CA_SNP, "Snoops          ", "Snoops          "},
  {PAPI_CA_SHR, "PAPI_CA_SHR     ", "Request for shared cache line (SMP)"},
  {PAPI_CA_CLN, "PAPI_CA_CLN     ", "Request for clean cache line (SMP)"},
  {PAPI_CA_INV, "PAPI_CA_INV     ", "Request for cache line Invalidation (SMP)"},
  {PAPI_CA_ITV, "PAPI_CA_ITV     ", "Request for cache line Intervention (SMP)"},
  {PAPI_L3_LDM, "L3 load misses  ", "Level 3 load misses"},
  {PAPI_L3_STM, "L3 store misses ", "Level 3 store misses"},
  {PAPI_BRU_IDL,"PAPI_BRU_IDL    ", "Cycles branch units are idle"},
  {PAPI_FXU_IDL,"PAPI_FXU_IDL    ", "Cycles integer units are idle"},
  {PAPI_FPU_IDL,"PAPI_FPU_IDL    ", "Cycles floating point units are idle"},
  {PAPI_LSU_IDL,"PAPI_LSU_IDL    ", "Cycles load/store units are idle"},
  {PAPI_TLB_DM, "Data TLB misses ", "Data translation lookaside buffer misses"},
  {PAPI_TLB_IM, "Inst TLB misses ", "Instr translation lookaside buffer misses"},
  {PAPI_TLB_TL, "Tot TLB misses  ", "Total translation lookaside buffer misses"},
  {PAPI_L1_LDM, "L1 load misses  ", "Level 1 load misses"},
  {PAPI_L1_STM, "L1 store misses ", "Level 1 store misses"},
  {PAPI_L2_LDM, "L2 load misses  ", "Level 2 load misses"},
  {PAPI_L2_STM, "L2 store misses ", "Level 2 store misses"},
  {PAPI_BTAC_M, "BTAC miss       ", "BTAC miss"},
  {PAPI_PRF_DM, "PAPI_PRF_DM     ", "Prefetch data instruction caused a miss"},
  {PAPI_L3_DCH, "L3 DCache Hit   ", "Level 3 Data Cache Hit"},
  {PAPI_TLB_SD, "PAPI_TLB_SD     ", "Xlation lookaside buffer shootdowns (SMP)"},
  {PAPI_CSR_FAL,"PAPI_CSR_FAL    ", "Failed store conditional instructions"},
  {PAPI_CSR_SUC,"PAPI_CSR_SUC    ", "Successful store conditional instructions"},
  {PAPI_CSR_TOT,"PAPI_CSR_TOT    ", "Total store conditional instructions"},
  {PAPI_MEM_SCY,"Cyc Stalled Mem ", "Cycles Stalled Waiting for Memory Access"},
  {PAPI_MEM_RCY,"Cyc Stalled MemR", "Cycles Stalled Waiting for Memory Read"},
  {PAPI_MEM_WCY,"Cyc Stalled MemW", "Cycles Stalled Waiting for Memory Write"},
  {PAPI_STL_ICY,"Cyc no InstrIss ", "Cycles with No Instruction Issue"},
  {PAPI_FUL_ICY,"Cyc Max InstrIss", "Cycles with Maximum Instruction Issue"},
  {PAPI_STL_CCY,"Cyc No InstrComp", "Cycles with No Instruction Completion"},
  {PAPI_FUL_CCY,"Cyc Max InstComp", "Cycles with Maximum Instruction Completion"},
  {PAPI_HW_INT, "HW interrupts   ", "Hardware interrupts"},
  {PAPI_BR_UCN, "Uncond br instr ", "Unconditional branch instructions executed"},
  {PAPI_BR_CN,  "Cond br instr ex", "Conditional branch instructions executed"},
  {PAPI_BR_TKN, "Cond br instr tk", "Conditional branch instructions taken"},
  {PAPI_BR_NTK, "Cond br instrNtk", "Conditional branch instructions not taken"},
  {PAPI_BR_MSP, "Cond br instrMPR", "Conditional branch instructions mispred"},
  {PAPI_BR_PRC, "Cond br instrCPR", "Conditional branch instructions corr. pred"},
  {PAPI_FMA_INS,"FMA instr comp  ", "FMA instructions completed"},
  {PAPI_TOT_IIS,"Total instr iss ", "Total instructions issued"},
  {PAPI_TOT_INS,"Total instr ex  ", "Total instructions executed"},
  {PAPI_INT_INS,"Int instr ex    ", "Integer instructions executed"},
  {PAPI_FP_INS, "FP instr ex     ", "Floating point instructions executed"},
  {PAPI_LD_INS, "Load instr ex   ", "Load instructions executed"},
  {PAPI_SR_INS, "Store instr ex  ", "Store instructions executed"},
  {PAPI_BR_INS, "br instr ex     ", "Total branch instructions executed"},
  {PAPI_VEC_INS,"Vec/SIMD instrEx", "Vector/SIMD instructions executed"},
  {PAPI_RES_STL,"Cyc proc stalled", "Cycles processor is stalled on resource"},
  {PAPI_FP_STAL,"Cyc any FP stall", "Cycles any FP units are stalled"},
  {PAPI_TOT_CYC,"Total cycles    ", "Total cycles"},
  {PAPI_LST_INS,"Tot L/S inst ex ", "Total load/store inst. executed"},
  {PAPI_SYC_INS,"Sync. inst. ex  ", "Sync. inst. executed"},
  {PAPI_L1_DCH, "L1 D Cache Hit  ", "L1 D Cache Hit"},
  {PAPI_L2_DCH, "L2 D Cache Hit  ", "L2 D Cache Hit"},
  {PAPI_L1_DCA, "L1 D Cache Acc  ", "L1 D Cache Access"},
  {PAPI_L2_DCA, "L2 D Cache Acc  ", "L2 D Cache Access"},
  {PAPI_L3_DCA, "L3 D Cache Acc  ", "L3 D Cache Access"},
  {PAPI_L1_DCR, "L1 D Cache Read ", "L1 D Cache Read"},
  {PAPI_L2_DCR, "L2 D Cache Read ", "L2 D Cache Read"},
  {PAPI_L3_DCR, "L3 D Cache Read ", "L3 D Cache Read"},
  {PAPI_L1_DCW, "L1 D Cache Write", "L1 D Cache Write"},
  {PAPI_L2_DCW, "L2 D Cache Write", "L2 D Cache Write"},
  {PAPI_L3_DCW, "L3 D Cache Write", "L3 D Cache Write"},
  {PAPI_L1_ICH, "L1 I cache hits ", "L1 instruction cache hits"},
  {PAPI_L2_ICH, "L2 I cache hits ", "L2 instruction cache hits"},
  {PAPI_L3_ICH, "L3 I cache hits ", "L3 instruction cache hits"},
  {PAPI_L1_ICA, "L1 I cache acc  ", "L1 instruction cache accesses"},
  {PAPI_L2_ICA, "L2 I cache acc  ", "L2 instruction cache accesses"},
  {PAPI_L3_ICA, "L3 I cache acc  ", "L3 instruction cache accesses"},
  {PAPI_L1_ICR, "L1 I cache reads", "L1 instruction cache reads"},
  {PAPI_L2_ICR, "L2 I cache reads", "L2 instruction cache reads"},
  {PAPI_L3_ICR, "L3 I cache reads", "L3 instruction cache reads"},
  {PAPI_L1_ICW, "L1 I cache write", "L1 instruction cache writes"},
  {PAPI_L2_ICW, "L2 I cache write", "L2 instruction cache writes"},
  {PAPI_L3_ICW, "L3 I cache write", "L3 instruction cache writes"},
  {PAPI_L1_TCH, "L1 cache hits   ", "L1 total cache hits"},
  {PAPI_L2_TCH, "L2 cache hits   ", "L2 total cache hits"},
  {PAPI_L3_TCH, "L3 cache hits   ", "L3 total cache hits"},
  {PAPI_L1_TCA, "L1 cache access ", "L1 total cache accesses"},
  {PAPI_L2_TCA, "L2 cache access ", "L2 total cache accesses"},
  {PAPI_L3_TCA, "L3 cache access ", "L3 total cache accesses"},
  {PAPI_L1_TCR, "L1 cache reads  ", "L1 total cache reads"},
  {PAPI_L2_TCR, "L2 cache reads  ", "L2 total cache reads"},
  {PAPI_L3_TCR, "L3 cache reads  ", "L3 total cache reads"},
  {PAPI_L1_TCW, "L1 cache writes ", "L1 total cache writes"},
  {PAPI_L2_TCW, "L2 cache writes ", "L2 total cache writes"},
  {PAPI_L3_TCW, "L3 cache writes ", "L3 total cache writes"},
  {PAPI_FML_INS,"FM ins          ", "FM ins"},
  {PAPI_FAD_INS,"FA ins          ", "FA ins"},
  {PAPI_FDV_INS,"FD ins          ", "FD ins"},
  {PAPI_FSQ_INS,"FSq ins         ", "FSq ins"},
  {PAPI_FNV_INS,"Finv ins        ", "Finv ins"},
  {PAPI_FP_OPS, "FP ops executed ", "Floating point operations executed"}};

static const int nentries = sizeof (papitable) / sizeof (Papientry);
static Papientry eventlist[MAX_AUX];     /* list of PAPI events to be counted */
static Papientry propeventlist[MAX_AUX]; /* list of PAPI events hoped to be counted */
static int nevents = 0;                  /* number of events: initialize to 0 */ 
static int nprop = 0;                    /* number of hoped events: initialize to 0 */ 
static int *EventSet;                    /* list of events to be counted by PAPI */
static long_long **papicounters;         /* counters return from PAPI */
static char papiname[PAPI_MAX_STR_LEN];  /* returned from PAPI_event_code_to_name */
static const int BADCOUNT = -999999;     /* Set counters to this when they are bad */
static int GPTLoverheadindx = -1;         /* index into counters array */
static long_long *lastoverhead;          /* needed because aux not available for overhead */

/* Function prototypes */

static int create_and_start_events (const int);

/*
** GPTL_PAPIsetoption: enable or disable PAPI event defined by "counter". Called 
**   from GPTLsetoption.  Since all events are off by default, val=false degenerates
**   to a no-op.  Coded this way to be consistent with the rest of GPTL
**
** Input args: 
**   counter: PAPI counter
**   val:     true or false for enable or disable
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTL_PAPIsetoption (const int counter,  /* PAPI counter */
			const int val)      /* true or false for enable or disable */
{
  int n;  /* loop index */
  
  /* Just return if the flag says disable an option, because default is off */

  if ( ! val)
    return 0;

  /*
  ** Loop through table looking for counter. If found, add the entry to the
  ** list of "proposed events".  Won't know till init time whether the event
  ** is available on this arch.
  */

  for (n = 0; n < nentries; n++)
    if (counter == papitable[n].counter) {
      if (nprop+1 > MAX_AUX) {
	return GPTLerror ("GPTL_PAPIsetoption: Event %s is too many\n", papitable[n].str);
      } else {
	propeventlist[nprop].counter = counter;
	propeventlist[nprop].prstr   = papitable[n].prstr;
	propeventlist[nprop].str     = papitable[n].str;
	printf ("GPTL_PAPIsetoption: will attempt to enable event %s\n", 
		propeventlist[nprop].str);
	++nprop;
      }
      return 0;
    }

  return GPTLerror ("GPTL_PAPIsetoption: counter %d does not exist\n", counter);
}

/*
** GPTL_PAPIinitialize(): Initialize the PAPI interface. Called from GPTLinitialize.
**   PAPI_library_init must be called before any other PAPI routines.  
**   PAPI_thread_init is called subsequently if threading is enabled.
**   Finally, allocate space for PAPI counters and start them.
**
** Input args: 
**   maxthreads: number of threads
**
** Return value: 0 (success) or GPTLerror or -1 (failure)
*/
 
int GPTL_PAPIinitialize (const int maxthreads)  /* number of threads */
{
  int ret;       /* return code */
  int n;         /* loop index */
  int counter;   /* PAPI counter */
  int t;         /* thread index */
  int *rc;       /* array of return codes from create_and_start_events */
  bool badret;   /* true if any bad return codes were found */

  /* PAPI_library_init needs to be called before ANY other PAPI routine */

  if ((ret = PAPI_library_init (PAPI_VER_CURRENT)) != PAPI_VER_CURRENT)
    return GPTLerror ("GPTL_PAPIinitialize: PAPI_library_init failure:%s\n",
		      PAPI_strerror (ret));

  /* PAPI_thread_init needs to be called if threading enabled */

#if ( defined THREADED_OMP )
  if (PAPI_thread_init ((unsigned long (*)(void)) (omp_get_thread_num)) != PAPI_OK)
    return GPTLerror ("GPTL_PAPIinitialize: PAPI_thread_init failure\n");
#elif ( defined THREADED_PTHREADS )
  if (PAPI_thread_init ((unsigned long (*)(void)) (pthread_self)) != PAPI_OK)
    return GPTLerror ("GPTL_PAPIinitialize: PAPI_thread_init failure\n");
#endif

  /* allocate and initialize static local space */

  EventSet     = (int *)        GPTLallocate (maxthreads * sizeof (int));
  papicounters = (long_long **) GPTLallocate (maxthreads * sizeof (long_long *));
  lastoverhead = (long_long *)  GPTLallocate (maxthreads * sizeof (long_long));

  for (t = 0; t < maxthreads; t++) {
    EventSet[t] = PAPI_NULL;
    papicounters[t] = (long_long *) GPTLallocate (MAX_AUX * sizeof (long_long));
    lastoverhead[t] = -1;
  }

  /* 
  ** Loop over events set by earlier calls to GPTL_PAPIsetoption. For the
  ** events which can be counted on this architecture, fill in the values
  ** (array "eventlist")
  */

  for (n = 0; n < nprop; n++) {
    counter = propeventlist[n].counter;
    if (PAPI_query_event (counter) != PAPI_OK) {
      (void) PAPI_event_code_to_name (counter, papiname);
      printf ("GPTL_PAPIinitialize: event %s not available on this arch\n", papiname);
    } else {
      if (nevents+1 > MAX_AUX) {
	(void) PAPI_event_code_to_name (counter, papiname);
	printf ("GPTL_PAPIinitialize: Event %s is too many\n", papiname);
      } else {
	if (counter == PAPI_TOT_CYC)
	  GPTLoverheadindx = nevents;

	eventlist[nevents].counter = counter;
	eventlist[nevents].prstr   = propeventlist[n].prstr;
	eventlist[nevents].str     = propeventlist[n].str;
	printf ("GPTL_PAPIinitialize: event %s enabled\n", eventlist[nevents].str);
	++nevents;
      }
    }
  }

  /* Event starting apparently must be within a threaded loop. */

  if (nevents > 0) {
    rc = (int *) GPTLallocate (maxthreads * sizeof (int));

#pragma omp parallel for private (t)

    for (t = 0; t < maxthreads; t++)
      rc[t] = create_and_start_events (t);
  
    badret = false;
    for (t = 0; t < maxthreads; t++)
      if (rc[t] < 0)
	badret = true;
    
    free (rc);

    if (badret)
      return -1;
  }

  return 0;
}

/*
** create_and_start_events: Create and start the PAPI eventset.  File-local.
**   Threaded routine to create the "event set" (PAPI terminology) and start
**   the counters. This is only done once, and is called from GPTL_PAPIinitialize 
** 
** Input args: 
**   t: thread number
**
** Return value: 0 (success) or GPTLerror (failure)
*/

static int create_and_start_events (const int t)  /* thread number */
{
  int ret;
  int n;

  /* Create the event set */

  if ((ret = PAPI_create_eventset (&EventSet[t])) != PAPI_OK)
    return GPTLerror ("GPTL_PAPIstart: failure creating eventset: %s\n", 
		      PAPI_strerror (ret));

  /* Add requested events to the event set */

  for (n = 0; n < nevents; n++) {
    if ((ret = PAPI_add_event (EventSet[t], eventlist[n].counter)) != PAPI_OK) {
      printf ("%s\n", PAPI_strerror (ret));
      return GPTLerror ("GPTL_PAPIstart: failure adding event: %s\n", 
			eventlist[n].str);
    }
  }

  /* Start the event set.  It will only be read from now on--never stopped */

  if ((ret = PAPI_start (EventSet[t])) != PAPI_OK)
    return GPTLerror ("%s\n", PAPI_strerror (ret));

  return 0;
}

/*
** GPTL_PAPIstart: Start the PAPI counters (actually they are just read).  
**   Called from GPTLstart.
**
** Input args:  
**   t: thread number
**
** Output args: 
**   aux: struct containing the counters
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTL_PAPIstart (const int t,          /* thread number */
		    Papistats *aux)       /* struct containing PAPI stats */
{
  int ret;  /* return code from PAPI lib calls */
  int n;    /* loop index */
  
  /* If no events are to be counted just return */

  if (nevents == 0)
    return 0;

  /* Read the counters */

  if ((ret = PAPI_read (EventSet[t], papicounters[t])) != PAPI_OK)
    return GPTLerror ("GPTL_PAPIstart: %s\n", PAPI_strerror (ret));

  /* 
  ** Store the counter values.  When GPTL_PAPIstop is called, the counters
  ** will again be read, and differenced with the values saved here.
  */

  for (n = 0; n < nevents; n++)
    aux->last[n] = papicounters[t][n];
  
  return 0;
}

/*
** GPTL_PAPIstop: Stop the PAPI counters (actually they are just read).  
**   Called from GPTLstop.
**
** Input args:
**   t: thread number
**
** Input/output args: 
**   aux: struct containing the counters
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTL_PAPIstop (const int t,         /* thread number */
		   Papistats *aux)      /* struct containing PAPI stats */
{
  int ret;          /* return code from PAPI lib calls */
  int n;            /* loop index */
  long_long delta;  /* change in counters from previous read */

  /* If no events are to be counted just return */

  if (nevents == 0)
    return 0;

  /* Read the counters */

  if ((ret = PAPI_read (EventSet[t], papicounters[t])) != PAPI_OK)
    return GPTLerror ("GPTL_PAPIstop: %s\n", PAPI_strerror (ret));
  
  /* 
  ** Accumulate the difference since timer start in aux.
  ** If negative accumulation has occurred (unfortunately this can and does
  ** happen, especially on AIX), store a flag value (BADCOUNT)
  */

  for (n = 0; n < nevents; n++) {
    delta = papicounters[t][n] - aux->last[n];
    if (delta < 0)
      aux->accum[n] = BADCOUNT;
    else if (aux->accum[n] != BADCOUNT)
      aux->accum[n] += delta;
  }
  return 0;
}

/* 
** GPTL_PAPIoverheadstart: Read the PAPI counters for overhead calcs (only
**   possible if total cycles are being counted). Called from GPTLstart and GPTLstop.
**   Have to set the static variable lastoverhead because a pointer to the correct 
**   aux timer is not yet available.
** 
** Input args: 
**   t: thread number
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTL_PAPIoverheadstart (const int t)  /* thread number */
{
  int ret;  /* return code from PAPI lib routine */

  /* If the overhead index hasn't been set we can't do anything */

  if (GPTLoverheadindx < 0)
    return -1;

  if ((ret = PAPI_read (EventSet[t], papicounters[t])) != PAPI_OK)
    return GPTLerror ("GPTL_PAPIoverheadstart: %s\n", PAPI_strerror (ret));

  lastoverhead[t] = papicounters[t][GPTLoverheadindx];

  return 0;
}

/*
** GPTL_PAPIoverheadstop: Read the PAPI counters and record overhead (only
**   possible if total cycles are being counted). Called from GPTLstart and GPTLstop.
**
** Input args:
**   t: thread number
**
** Input/output args:
**   aux: struct containing the overhead accumulator
**
** Return value: 0 (success) or GPTLerror (failure)
*/
    
int GPTL_PAPIoverheadstop (const int t,         /* thread number */               
			   Papistats *aux)      /* struct containing PAPI stats */
{
  int ret;         /* return code from PAPI_read */
  long_long diff;  /* difference in cycle count from when started */

  /* overheadindx <= 0 means cycle counting was not enabled */

  if (GPTLoverheadindx < 0)
    return -1;

  if ((ret = PAPI_read (EventSet[t], papicounters[t])) != PAPI_OK)
    return GPTLerror ("GPTL_PAPIoverheadstart: %s\n", PAPI_strerror (ret));

  /* Accumulate the overhead cycles.  Check for a negative increment */

  diff = papicounters[t][GPTLoverheadindx] - lastoverhead[t];
  if (diff < 0)
    aux->accum_cycles = BADCOUNT;
  else
    aux->accum_cycles += diff;

  return 0;
}

/*
** GPTL_PAPIprstr: Print the descriptive string for all enabled PAPI events.
**   Called from GPTLpr.
**
** Input args: 
**   fp: file descriptor
*/

void GPTL_PAPIprstr (FILE *fp)   /* file descriptor */
{
  int n;
  
  for (n = 0; n < nevents; n++)
    fprintf (fp, "%16s ", eventlist[n].prstr);

  if (lastoverhead[0] > -1)
    fprintf (fp, "Overhead (cycles)");
}

/*
** GPTL_PAPIpr: Print PAPI counter values for all enabled events. Called from
**   GPTLpr.
**
** Input args: 
**   fp: file descriptor
**   aux: struct containing the counters
*/

void GPTL_PAPIpr (FILE *fp,              /* file descriptor to write to */
		  const Papistats *aux)  /* stats to write */
{
  int n;
  
  for (n = 0; n < nevents; n++) {
    if (aux->accum[n] < 1000000)
      fprintf (fp, "%16ld ", (long) aux->accum[n]);
    else
      fprintf (fp, "%16.10e ", (double) aux->accum[n]);
  }

  /* The check on lastoverhead > -1 determines whether it was ever set */

  if (lastoverhead[0] > -1)
    if (aux->accum_cycles < 1000000)
      fprintf (fp, "%16ld ", (long) aux->accum_cycles);
    else
      fprintf (fp, "%16.10e ", (double) aux->accum_cycles);
}

/*
** GPTLPAPIprinttable:  Print table of PAPI native counters. Not all are
**   necessarily available on this architecture. This is the one routine 
**   in this file which is user-visible. No underscores in GPTLPAPIprinttable 
**   to avoid underscore weirdness of g77
*/

void GPTLPAPIprinttable ()
{
  int n;

  for (n = 0; n < nentries; n++)
    printf ("%d %s\n", papitable[n].counter, papitable[n].str);
}

/*
** GPTL_PAPIadd: Accumulate PAPI counters. Called from add.
**
** Input/Output args: 
**   auxout: auxout = auxout + auxin
**
** Input args:
**   auxin: counters to be summed into auxout
*/

void GPTL_PAPIadd (Papistats *auxout,      /* output struct */
		   const Papistats *auxin) /* input struct */
{
  int n;
  
  for (n = 0; n < nevents; n++)
    if (auxin->accum[n] == BADCOUNT || auxout->accum[n] == BADCOUNT)
      auxout->accum[n] = BADCOUNT;
    else
      auxout->accum[n] += auxin->accum[n];

  /* Overhead calcs */

  if (auxin->accum_cycles == BADCOUNT || auxout->accum_cycles == BADCOUNT)
    auxout->accum_cycles = BADCOUNT;
  else
    auxout->accum_cycles += auxin->accum_cycles;
}

/*
** GPTL_PAPIfinalize: finalization routine must be called from single-threaded
**   region. Free all malloc'd space
*/

void GPTL_PAPIfinalize (int maxthreads)
{
  int t;

  for (t = 0; t < maxthreads; t++)
    free (papicounters[t]);

  free (EventSet);
  free (papicounters);
  free (lastoverhead);
}
#endif
