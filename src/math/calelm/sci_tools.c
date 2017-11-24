/* Copyright (c) 1997 by Inria Lorraine.  All Rights Reserved */

/***
   NAME
     sci_tools
   PURPOSE
     
   NOTES
     
   HISTORY
     fleury - Dec 17, 1997: Created.
     $Log: sci_tools.c,v $
     Revision 1.3  2009-07-07 08:36:53  vernaz
     Y. Vernaz - add calemlm directory

     Revision 1.1  2007-07-06 09:49:15  thib
     *** empty log message ***

     Revision 1.1  2007/06/25 13:21:38  vernaz
     Y. Vernaz - add calelm/ sources (2)

     Revision 1.7  2005/10/22 18:53:10  cornet
     update memory management under Windows
     use HeapAlloc and VirtualAlloc (for scilab stack)
     Correction Bug 1576
     n=10000
     xbasc();
     plot2d([0,1],[0,n],0)
     xpols=[zeros(1,n); ones(2,n); zeros(1,n)];
     ypols=[2:n+1; 2:n+1; 1:n; 1:n];
     xfpolys(xpols,ypols,modulo((1:n),32))

     and for windows on PC with 256 mo
     stacksize(5000000*20)
     a=rand(9999,9999)

     Revision 1.6  2005/07/04 06:15:28  cornet
     correction compilation

     Revision 1.5  2005/07/03 18:16:10  cornet
     optimisation des MALLOC pour Windows ( A tester avec attention ) --> VirtualAlloc

     Revision 1.4  2005/07/01 07:08:08  cornet
     replace malloc, free, calloc & realloc by MALLOC,FREE,CALLOC & REALLOC defined in SCI/routines/sci_mem_alloc.h

     Revision 1.3  2001/10/16 08:23:38  chanceli
     change in includes

     Revision 1.2  2001/06/11 17:53:43  delebecq
     f772sci with 2 pointers

     Revision 1.1.1.1  2001/04/26 07:47:34  scilab
     Imported sources

     Revision 1.5  1998/03/27 12:20:22  fleury
     Version pvm OK.
     TODO: faire des tests de compil sur plateforme separee (POPC0
     TODO: commenter source (-;
     TODO: faire un peu de netoyage

     Revision 1.4  1998/03/17 11:49:33  fleury
     Broadcast OK.
     TODO: mettre les listes.
     TODO: faire qcq tests

     Revision 1.3  1998/03/13 13:57:07  fleury
     Version send/recv avec pack. A tester.
     TODO: ajouter les listes + BROADCAST
     TODO: faire un clean du dir et des fichiers...

     Revision 1.2  1998/01/06 13:23:48  fleury
     Use memcopy instead of fori

     Revision 1.1  1997/12/18 18:35:57  fleury
     Premier commit

     Ajout de fct permettant de tester le type d une variable scilab
     conversion de complex format scilab to complex format f77
     TODO:use memcpy
     TODO:use imatrix

***/
#include "../machine.h"

#ifdef __STDC__
#include <stdlib.h>
#else 
#include <malloc.h>
#endif 

#include <stdio.h>
#include <string.h>

#include "sci_tools.h"

#ifdef WIN32
 #include "../os_specific/win_mem_alloc.h" /* MALLOC */
#else
 #include "../os_specific/sci_mem_alloc.h" /* MALLOC */
#endif


#ifdef __STDC__
void 
C2F(ccomplexf)(int *n, double **ip, double *op)
#else
void 
C2F(ccomplexf)(n, ip, op)
  int *n;
  double **ip;
  double *op;
#endif 
{
  memcpy(op, *ip, *n * sizeof(double));

  /* int i */
  /*   for (i = *n; --i >= 0; ) { */
  /*     op[i] = (*ip)[i];*/		/* TODO: replace by memcpy */ 
  /*   } */
  
  SET_TYPE_COMPLEX(op);		        /* type is complex */
  SET_NB_ROW(op,  NB_ROW(op) / 2);	/* nb  row is halfed */

  FREE((char*) (*ip));
} /* ccomplexf */

#ifdef __STDC__
void 
(SciToF77)(double *ptr, int size, int lda)
#else
void 
SciToF77(ptr, size, lda)
  double *ptr;
  int size;
  int lda;
#endif 
{
  int i;
  double *tab;
  
  if ((tab = (double *) MALLOC(size * sizeof(double))) == NULL) {
    (void) fprintf(stderr, "SciToF77: Error malloc\n");
    return;
  }

  /* for (i = size; --i >= 0; ) { */
  /*     tab[i] = ptr[i]; */
  /*   } */

  memcpy(tab, ptr, size * sizeof(double));

  for (i = 0; i < size; ++i) {
    ptr[2*i] = tab[i];
    ptr[2*i+1] = ptr[lda+i];
  }

  FREE(tab);
} /* SciToF77 */


#ifdef __STDC__
void 
(F77ToSci)(double *ptr, int size, int lda)
#else
void 
F77ToSci(ptr, size, lda)
  double *ptr;
  int size;
  int lda;
#endif 
{
  int i;
  double *tab;
  
  if ((tab = (double *) MALLOC(size * sizeof(double))) == NULL) {
    (void) fprintf(stderr, "F77ToSci: Error malloc\n");
    return;
  }
  
  for (i = 0; i < size; ++i) {
    tab[i] = ptr[2*i+1];
    ptr[i] = ptr[2*i];
  }

  memcpy(ptr + lda, tab, size * sizeof(double));

  /*   for (i = size; --i >= 0; ) { */
  /*     ptr[lda+i] = tab[i]; */
  /*   } */

  FREE(tab);
} /* F77ToSci */


/* double2z and z2double : same as above with two pointers dest and src 
   double2z ptr = src, ptr77z = dest (z format)     
   z2double ptr = src (z format) , ptrsci = dest */  

#ifdef __STDC__
void 
(double2z)(double *ptr,  double *ptr77z, int size, int lda)
#else
void 
double2z(ptr, ptr77z, size, lda)
  double *ptr; double *ptr77z; 
  int size;
  int lda;
#endif 
{
  int i;
  double *tab;
  
  if ((tab = (double *) MALLOC(size * sizeof(double))) == NULL) {
    (void) fprintf(stderr, "Double2z: Error malloc\n");
    return;
  }

  memcpy(tab, ptr, size * sizeof(double));

  for (i = 0; i < size; ++i) {
    ptr77z[2*i] = tab[i];
    ptr77z[2*i+1] = ptr[lda+i];
  }

  FREE(tab);
} 


#ifdef __STDC__
void 
(z2double)(double *ptrz, double *ptrsci, int size, int lda)
#else
void 
z2double(ptrz, ptrsci, size, lda)
  double *ptrz; double *ptrsci;
  int size;
  int lda;
#endif 
{
  int i;
  double *tab;
  
  if ((tab = (double *) MALLOC(size * sizeof(double))) == NULL) {
    (void) fprintf(stderr, "z2double: Error malloc\n");
    return;
  }
  
  for (i = 0; i < size; ++i) {
    tab[i] = ptrz[2*i+1];
    ptrsci[i] = ptrz[2*i];
  }

  memcpy(ptrsci + lda, tab, size * sizeof(double));

  FREE(tab);
} 

