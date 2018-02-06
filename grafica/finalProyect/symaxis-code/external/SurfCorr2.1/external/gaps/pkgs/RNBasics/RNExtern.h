/* Include file for external stuff */



/* Standard dependency include files */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>



/* Machine dependent include files */

#if (RN_OS == RN_WINDOWSNT)
#   include <float.h>
#   include <windows.h>
#endif

#if ((RN_OS == RN_LINUX) || (RN_OS == RN_MACOSX) )
#   include <float.h>
#endif

#if ((RN_OS == RN_IRIX) || (RN_OS == RN_LINUX) || (RN_OS == RN_MACOSX))
#   include <sys/time.h>
#   include <sys/resource.h>
#endif

#if ((RN_OS == RN_IRIX) || (RN_OS == RN_LINUX) || (RN_OS == RN_WINDOWSNT))
#	include <malloc.h>
#endif


