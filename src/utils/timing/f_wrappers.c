/*
** $Id$
** 
** Fortran wrappers for timing library routines as called by CAM.
** When changing to the new interface (GPTLxxx), all t_ names below
** need to change to GPTL.
*/

#if ( defined T3D )
#include <fortran.h>
#endif

#include <cfort.h>
#include "private.h"

#if ( defined FORTRANCAPS )

#define gptlsetoption GPTLSETOPTION
#define gptlinitialize GPTLINITIALIZE
#define t_setoptionf T_SETOPTIONF
#define t_initializef T_INITIALIZEF
#define t_finalizef T_FINALIZEF
#define t_prf T_PRF
#define t_resetf T_RESETF
#define t_stampf T_STAMPF
#define t_startf T_STARTF
#define t_stopf T_STOPF

#elif ( defined FORTRANUNDERSCORE )

#define gptlsetoption gptlsetoption_
#define gptlinitialize gptlinitialize_
#define t_setoptionf t_setoptionf_
#define t_initializef t_initializef_
#define t_finalizef t_finalizef_
#define t_prf t_prf_
#define t_resetf t_resetf_
#define t_stampf t_stampf_
#define t_startf t_startf_
#define t_stopf t_stopf_

#elif ( defined FORTRANDOUBLEUNDERSCORE )

#define gptlsetoption gptlsetoption__
#define gptlinitialize gptlinitialize__
#define t_setoptionf t_setoptionf__
#define t_initializef t_initializef__
#define t_finalizef t_finalizef__
#define t_prf t_prf__
#define t_resetf t_resetf__
#define t_stampf t_stampf__
#define t_startf t_startf__
#define t_stopf t_stopf__

#endif

#if ( defined T3D )

int t_startf (_fcd);
int t_stopf (_fcd);

#else

int t_startf (char *, int);
int t_stopf (char *, int);

#endif

int gptlsetoption (int *option, int *val)
{
  return GPTLsetoption ( (Option) *option, *val);
}

int gptlinitialize ()
{
  return GPTLinitialize ();
}

int t_setoptionf (int *option, int *val)
{
  return GPTLsetoption ( (Option) *option, *val);
}

int t_initializef ()
{
  return GPTLinitialize ();
}

int t_finalizef ()
{
  return GPTLfinalize ();
}

int t_prf (int *procid)
{
  return GPTLpr (*procid);
}

void t_resetf ()
{
  GPTLreset();
  return;
}

int t_stampf (double *wall, double *usr, double *sys)
{
  return GPTLstamp (wall, usr, sys);
}

#if ( defined T3D )

int t_startf (_fcd name)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (_fcdlen (name), MAX_CHARS);
  strncpy (cname, _fcdtocp (name), numchars);
  cname[numchars] = '\0';
  return GPTLstart (cname);
}

int t_stopf (_fcd name)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (_fcdlen (name), MAX_CHARS);
  strncpy (cname, _fcdtocp (name), numchars);
  cname[numchars] = '\0';
  return GPTLstop (cname);
}

#else

int t_startf (char *name, int nc1)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc1, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';
  return GPTLstart (cname);
}

int t_stopf (char *name, int nc1)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc1, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';
  return GPTLstop (cname);
}

#endif
