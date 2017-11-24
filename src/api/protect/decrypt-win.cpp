/*=======================================================================

     decrypt.cc                                               version 1.0

=======================================================================*/

#include "hex.h"
#include "default.h"
#include <stdio.h>


USING_NAMESPACE(CryptoPP)
	USING_NAMESPACE(std)
const int MAX_PHRASE_LENGTH=250;


extern "C" { int DecryptString(const char *ciphertext, char *decryptedstring); }
char passPhrase[15]="L5IÈN_ELy=Ôqzq";


int DecryptString(char *instr, char *decryptedstring)
{
try
	{
  string outstr;

  HexDecoder decryptor(new 
    DefaultDecryptorWithMAC(passPhrase, new StringSink(outstr)));
  decryptor.Put((byte*)instr, strlen(instr));
 
  //printf( "out= %s\n",outstr.c_str() ); 
  //if ( outstr.c_str() != NULL ) {
  	decryptor.MessageEnd();
  //}

    strcpy(decryptedstring, outstr.c_str());
	return 0 ;
   } catch(CryptoPP::Exception &e)
	{
		return 1;
	}
}

