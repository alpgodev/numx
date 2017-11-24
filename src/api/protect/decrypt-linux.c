/*=======================================================================

	decrypt-linux.c					      version 1.0 

=========================================================================*/

#include <openssl/blowfish.h>
#include <openssl/evp.h>
#include <fcntl.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#define IP_SIZE 128
#define OP_SIZE 200

unsigned char key[16]="L5IeN_ELy=iq";
unsigned char iv[8]="k7^&%HF";

unsigned char *
base64_decode (unsigned char *bbuf, unsigned int *len)
{
  unsigned char *ret;
  unsigned int bin_len;

/* integer divide by 4 then multiply by 3, its binary so no NULL */
  bin_len = (((strlen (bbuf) + 3) / 4) * 3);
  ret = (unsigned char *) malloc (bin_len);
  *len = EVP_DecodeBlock (ret, bbuf, strlen (bbuf));
  return ret;
}


int
decrypt (char *inbuff, char *outbuf) 
{
	int olen, tlen, n;
	EVP_CIPHER_CTX ctx;
	EVP_CIPHER_CTX_init (&ctx);
	EVP_DecryptInit (&ctx, EVP_bf_cbc (), key, iv);


	if (EVP_DecryptUpdate (&ctx, outbuf, &olen, inbuff, 128) != 1)
	    {
		//printf ("error in decrypt update\n");
		return 0;
	    }
	
	if (EVP_DecryptFinal (&ctx, outbuf + olen, &tlen) != 1)
	    {
		//printf ("error in decrypt final\n");
		return 0;
	    }
	olen += tlen;

	//printf ("longeur de la decryption %i \n", olen);

	
	EVP_CIPHER_CTX_cleanup (&ctx);
	return 1;
}

