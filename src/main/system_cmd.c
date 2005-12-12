#include <misc.h>
#include <cfort.h>
#include <stdlib.h>

#if ( defined FORTRANCAPS )
#define system_cmd SYSTEM_CMD
#elif ( defined FORTRANUNDERSCORE )
#define system_cmd system_cmd_
#elif ( defined FORTRANDOUBLEUNDERSCORE )
#define system_cmd system_cmd__
#endif

/*
//BOP
// !ROUTINE: system_cmd
// !INTERFACE:
int system_cmd(const char *text)
// !DESCRIPTION:
Wrapper to "C" system command.  Allows all platforms to have a consistent
interface that includes an error return code.
// !REVISION HISTORY:
//EOP
*/

int system_cmd(const char *text)
{
  int err;   /* Error return code */

  err = system(text);
  return( err );
}
