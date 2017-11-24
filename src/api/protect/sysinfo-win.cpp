/*=======================================================================

	sysinfo.cc					      version 1.0 
	
=========================================================================*/

#include <windows.h>
#include <wincon.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/timeb.h>
#include <sys/types.h>

typedef struct _ASTAT_
{
	ADAPTER_STATUS adapt;
	NAME_BUFFER    NameBuff [30];
} ASTAT, * PASTAT;

ASTAT Adapter;

int DecryptString(char *instr, char *decryptedstring);

extern "C" { int checkSysInfo(); };

#define MAC_ADDRESS_LENGTH 18
char Res[MAC_ADDRESS_LENGTH];

/* ethernet card - physical adress */
char* getMacAddress(void)
{
	NCB Ncb;
      	UCHAR uRetCode;
      	LANA_ENUM   lenum;

      	memset( &Ncb, 0, sizeof(Ncb) );
      	Ncb.ncb_command = NCBENUM;
      	Ncb.ncb_buffer = (UCHAR *)&lenum;
      	Ncb.ncb_length = sizeof(lenum);
      	uRetCode = Netbios( &Ncb );

	if ( lenum.length > 0 ) 
	{	
      		//for(i=0; i < lenum.length ;i++)
      		//{
          	memset( &Ncb, 0, sizeof(Ncb) );
          	Ncb.ncb_command = NCBRESET;
          	Ncb.ncb_lana_num = lenum.lana[0];

          	uRetCode = Netbios( &Ncb );
          	//printf( "The NCBRESET on LANA %d return code is: 0x%x \n",
                  //lenum.lana[i], uRetCode );

          	memset( &Ncb, 0, sizeof (Ncb) );
          	Ncb.ncb_command = NCBASTAT;
          	Ncb.ncb_lana_num = lenum.lana[0];

          	strcpy( (char *)Ncb.ncb_callname,  "*               " );
          	Ncb.ncb_buffer = (UCHAR *) &Adapter;
          	Ncb.ncb_length = sizeof(Adapter);
          	uRetCode = Netbios( &Ncb );
          	if ( uRetCode == 0 )
          	{

			sprintf( Res, "%02X-%02X-%02X-%02X-%02X-%02X", 
                  		Adapter.adapt.adapter_address[0],
                  		Adapter.adapt.adapter_address[1],
                  		Adapter.adapt.adapter_address[2],
                  		Adapter.adapt.adapter_address[3],
                  		Adapter.adapt.adapter_address[4],
                  		Adapter.adapt.adapter_address[5] );
#ifdef DEBUG
	 		printf( "local mac address = %s\n" , Res );
#endif
          	} else {
			sprintf( Res, "");
			//return (char*) 0;
		}
	} else {
		sprintf( Res, "");
		//return (char*) 0;
	}
	return Res;
}

time_t getSystemDate() 
{
	time_t t;
	time( &t );             
	return t;
}

/**
 * Return 0 if successful
 * Error code:
 *    1 : version date expired
 *    2 : license date expired
 *    3 : computer invalid (wrong physical adress)
 *    4 : registry key not found 
 *    5 : decrypt failed
 */
int checkSysInfo() 
{

    	HKEY hKey;
	char macAddress[200];
    	DWORD dwBufLen = 200;
	char date[200];
    	DWORD dwDateLen = 200;
    	//LONG lRet;
	char decaddr[251];
	char decdate[251];
	int decreturn;
	
	/* Expired date: a scalar with the number of seconds since Jan 1, 1970, 00:00 UTC (Unix Time Convention) */
	/* cf. http://www.onlineconversion.com/unix_time.htm */
	/* security: maximum expired date:  31-dec-2007 -> 1199052000 secondes */
	/* security: maximum expired date:  31-dec-2008 -> 1230674400 secondes */
	/* security: maximum expired date:  15-jan-2010 -> 1263513600 secondes */
	/* security: maximum expired date:  15-jan-2011 -> 1295049600 secondes */
	/* security: maximum expired date:   2-apr-2015 -> 1427977201 secondes */

	if ( 1427977201 < getSystemDate() ) return 1;

#ifndef NOCHECKMAC

	if( RegOpenKeyEx(HKEY_LOCAL_MACHINE,TEXT("SOFTWARE\\raisepartner\\norm"),
                     0, KEY_QUERY_VALUE, &hKey) != ERROR_SUCCESS) return 4;
	// Field1 = macAddress
    	if ( RegQueryValueEx(hKey, TEXT("field1"), NULL, NULL, (LPBYTE)macAddress, &dwBufLen) 
	 	!= ERROR_SUCCESS) return 4;	
		
#ifdef DEBUG
	printf( "macaddress in reg %s\n", macAddress);
#endif
	decreturn  = DecryptString( macAddress, decaddr );
	
	if ( decreturn == 1 ) return 5;

#ifdef DEBUG
	printf( "macaddress in reg %s\n", decaddr  );
#endif
		if ( strcmp( decaddr, getMacAddress() ) ) return 3;

      // Field2 = date
	if ( RegQueryValueEx(hKey, TEXT("field2"), NULL, NULL, (LPBYTE)date, &dwDateLen) 
	 	!= ERROR_SUCCESS) return 4;

	decreturn  = DecryptString( date, decdate);

	if ( decreturn == 1 ) return 5;



#ifdef DEBUG
		printf( "date in reg %s\n", decdate );
		printf( "current date %d\n", getSystemDate() );
#endif

    	if ( atol( decdate  ) < getSystemDate() ) return 2;

    	RegCloseKey(hKey);

#endif /* NOCHECKMAC */

	return 0;
}
