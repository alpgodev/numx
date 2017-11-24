/*=======================================================================
	sysinfo-macosx.c
=========================================================================*/

/*
 * mac_addr_sys.c
 *
 * Return the MAC (ie, ethernet hardware) address by using system specific
 * calls.
 *
 * compile with: gcc -c -D "OS" mac_addr_sys.c
 * with "OS" is one of Linux, AIX, HPUX 
 */

#include <stdio.h>
#include <string.h>

/**
 * Return 0 if successful
 * Error code:
 *    1 : version date expired
 *    2 : license date expired
 *    3 : computer invalid
 *    4 : registry key not found 
 *    5 : decrypt failed
 */
int checkSysInfo() 
{	
	  return 0;
}	
	
	

