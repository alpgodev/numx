/*=======================================================================

	sysinfo-linux.c					      version 1.0 

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

#define Linux

#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#ifdef Linux
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <linux/if.h>
#endif

#ifdef HPUX
#include <netio.h>
#endif

#ifdef AIX
#include <sys/ndd_var.h>
#include <sys/kinfo.h>
#endif

#define IP_SIZE 128
#define OP_SIZE 200

/* typedef struct _ASTAT_
{
	ADAPTER_STATUS adapt;
	NAME_BUFFER    NameBuff [30];
} ASTAT, * PASTAT;

ASTAT Adapter; */

extern int decrypt (char *inbuff, char *outbuf);
extern unsigned char * base64_decode (unsigned char *bbuf, unsigned int *len);
extern int checkSysInfo();

#define MAC_ADDRESS_LENGTH 18
char Res[MAC_ADDRESS_LENGTH];

/* Reads system Macaddress for a banch of systems : provisional
AIX and HPUX, main function fir Linux Kernels */

long mac_addr_sys ( u_char *addr)
{
/* implementation for Linux */
#ifdef Linux
    struct ifreq ifr;
    struct ifreq *IFR;
    struct ifconf ifc;
    char buf[1024];
    int s, i;
    int ok = 0;

    s = socket(AF_INET, SOCK_DGRAM, 0);
    if (s==-1) {
        return -1;
    }

    ifc.ifc_len = sizeof(buf);
    ifc.ifc_buf = buf;
    ioctl(s, SIOCGIFCONF, &ifc);
 
    IFR = ifc.ifc_req;
    for (i = ifc.ifc_len / sizeof(struct ifreq); --i >= 0; IFR++) {

        strcpy(ifr.ifr_name, IFR->ifr_name);
        if (ioctl(s, SIOCGIFFLAGS, &ifr) == 0) {
            if (! (ifr.ifr_flags & IFF_LOOPBACK)) {
                if (ioctl(s, SIOCGIFHWADDR, &ifr) == 0) {
                    ok = 1;
                    break;
                }
            }
        }
    }

    close(s);
    if (ok) {
        bcopy( ifr.ifr_hwaddr.sa_data, addr, 6);
    }
    else {
        return -1;
    }
    return 0;
#endif

/* implementation for HP-UX */
#ifdef HPUX

#define LAN_DEV0 "/dev/lan0"

    int		fd;
    struct fis	iocnt_block;
    int		i;
    char	net_buf[sizeof(LAN_DEV0)+1];
    char	*p;

    (void)sprintf(net_buf, "%s", LAN_DEV0);
    p = net_buf + strlen(net_buf) - 1;

    /* 
     * Get 802.3 address from card by opening the driver and interrogating it.
     */
    for (i = 0; i < 10; i++, (*p)++) {
        if ((fd = open (net_buf, O_RDONLY)) != -1) {
			iocnt_block.reqtype = LOCAL_ADDRESS;
			ioctl (fd, NETSTAT, &iocnt_block);
			close (fd);

            if (iocnt_block.vtype == 6)
                break;
        }
    }

    if (fd == -1 || iocnt_block.vtype != 6) {
        return -1;
    }

	bcopy( &iocnt_block.value.s[0], addr, 6);
	return 0;

#endif /* HPUX */

/* implementation for AIX */
#ifdef AIX

    int size;
    struct kinfo_ndd *nddp;

    size = getkerninfo(KINFO_NDD, 0, 0, 0);
    if (size <= 0) {
        return -1;
    }
    nddp = (struct kinfo_ndd *)malloc(size);
          
    if (!nddp) {
        return -1;
    }
    if (getkerninfo(KINFO_NDD, nddp, &size, 0) < 0) {
        free(nddp);
        return -1;
    }
    bcopy(nddp->ndd_addr, addr, 6);
    free(nddp);
    return 0;
#endif

/* Not implemented platforms */
	return -1;
}

/***********************************************************************/
/*
 * Copies the  adresse Mac Address in the Res Field from the addr
  fileds : format with minus separator.
 */

char* getMacAddress(void)
{
    long stat;
    int i;
    u_char addr[6];

    stat = mac_addr_sys( addr);
    if (0 == stat) {

	sprintf( Res, "%02X-%02X-%02X-%02X-%02X-%02X", 
		 addr[0],
		 addr[1],
		 addr[2],
		 addr[3],
		 addr[4],
		 addr[5] );

        return (Res);

    }
    else {
      // fprintf( stderr, "can't get MAC address\n");
        return ((char *) 0);
    }
    return Res;
}

/**
 Function gets system date in seconds from the base
 date, 1 st of January 1970 (reference date)
*/
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
 *    3 : computer invalid
 *    4 : registry key not found 
 *    5 : decrypt failed
 */
int checkSysInfo() 
{
	char macAddress[200];
    	int dwBufLen = 200;
	char date[200];
    	int dwDateLen = 200;
      time_t sec;

	unsigned char outbuf[OP_SIZE];
	unsigned char decaddr[OP_SIZE];
	unsigned char decdate[OP_SIZE];
	char inbuff[IP_SIZE];
	int len;
      char * ptr2;

	FILE *file;
	int numbers[30];
	char str[200]; 
	
	// Additional Security Date (31 12 2004) 
	if ( 1291158000  < getSystemDate() ) return 1;

	file = fopen("/etc/raise.reg", "r");

	if(file==NULL) {
	  printf("Error: can't open file.\n");
	  return 4;
	};
	
	if (!feof(file)) {
	  fscanf(file, "%s", &str);
	} else {
	  return 1;
	}
	
	if (!feof(file)) {
	  fscanf(file, "%s", &macAddress);
	} else {
	  return 1;
	}	
	
	
#ifdef DEBUG
	printf( "macaddress in reg %s\n", macAddress);
      printf (" l'adresse mac %02X-%02X-%02X-%02X-%02X-%02X \n",  
              macAddress[0], macAddress[1], macAddress[3], 
              macAddress[4], macAddress[5], macAddress[5]);
#endif
#ifndef NOCHECKMAC
      ptr2 = base64_decode (macAddress, &len);

      decrypt(ptr2, decaddr);

      free(ptr2);
         
	if ( decaddr == 0 ) return 5;

#ifdef DEBUG
	printf( "macaddress in reg %s\n", decaddr  );
#endif

    	if ( strcmp( decaddr, getMacAddress() ) ) return 3;
#endif /* NOCHECKMAC */

      // Read the string field1 in the file of the registry

	if (!feof(file)) {
	  fscanf(file, "%s", &str);
	} else {
	  return 1;
	}

      // Read the date value encrypted and base64 encoded
	
	if (!feof(file)) {
	  fscanf(file, "%s", &date);
	} else {
	  return 1;
	}	

      ptr2 = base64_decode (date, &len);

      decrypt(ptr2, decdate);

      free(ptr2);

	// If the date looks ugly go out and return appropriate error
	
	if ( decdate == 0 ) return 5;
#ifdef DEBUG
	        sec = atoi (decdate);	
		printf( "date in reg %s  %s \n", decdate, ctime(&sec));
                sec = getSystemDate();
		printf( "current date %d  %s\n", getSystemDate(), ctime(&sec));
#endif
	
    	if ( atol( decdate  ) < getSystemDate() ) return 2;
 
	fclose(file);
	return 0;
}

