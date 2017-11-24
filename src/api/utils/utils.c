/*=======================================================================

     Java-C-Fortran Utilities                      

 =======================================================================*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef WINDOWS
#define DLLIMPORT __declspec (dllexport)
#endif /*WINDOWS */
int  lowestindex(int, int, int);
int  next(int, int, int);
int checkSysInfo();

void DisplayMatrixDouble( char *str, int ligne, int colonne, double *matptr)
{
  int i,j;
  printf("\n matrice %s :  \n", str);
  for (i=0; i<ligne; i++) {
    for (j=0; j<colonne; j++) {
      printf(" %5.8f",matptr[j+(i*colonne)]);
    }
    printf("\n");
  } 
}
void DisplayMatrixDoubleInFile( char *str, int ligne, int colonne, double *matptr, FILE *f)
{
  int i,j;
  fprintf(f, "\n matrice %s :  \n", str);
  for (i=0; i<ligne; i++) {
    for (j=0; j<colonne; j++) {
      fprintf(f, " %5.20f",matptr[j+(i*colonne)]);
    }
    fprintf(f, "\n");
  } 
}
void DisplayMatrixInt( char *str, int ligne, int colonne,int *matptr)
{
  int i,j;
  printf("\n matrice %s :  \n", str);
  for (i=0; i<ligne; i++) {
    for (j=0; j<colonne; j++) {
      printf(" %i",matptr[j+(i*colonne)]);
    }
    printf("\n");
  } 
}
void DisplayMatrixIntInFile( char *str, int ligne, int colonne,int *matptr, FILE *f)
{
  int i,j;
  fprintf(f, "\n matrice %s :  \n", str);
  for (i=0; i<ligne; i++) {
    for (j=0; j<colonne; j++) {
      fprintf(f, " %i",matptr[j+(i*colonne)]);
    }
    fprintf(f, "\n");
  } 
}
void DisplayVectorDouble( char *str, int ligne, double *matptr)
{
  int i;
  printf("\n vecteur %s : ", str);
  for (i=0; i<ligne; i++) {
    printf(" %5.8f",matptr[i]);
  }
}
void DisplayVectorDoubleInFile( char *str, int ligne, double *matptr, FILE *f)
{
  int i;
  fprintf(f, "\n vecteur %s : ", str);
  for (i=0; i<ligne; i++) {
    fprintf(f, " %5.20f",matptr[i]);
  }
}
void DisplayVectorInt( char *str, int ligne, int *matptr)
{
  int i;
  printf("\n vecteur %s : ", str);
  for (i=0; i<ligne; i++) {
    printf(" %i",matptr[i]);
  }
}
void DisplayVectorIntInFile( char *str, int ligne, int *matptr, FILE *f)
{
  int i;
  fprintf(f, "\n vecteur %s : ", str);
  for (i=0; i<ligne; i++) {
    fprintf(f, " %i",matptr[i]);
  }
}

void transposeDoubleSquareMatrix(double * B, int ligne, int colonne)
{

   int i, j;
   double Btemp;
   for (i=0; i<ligne; i++)
   {
     for (j=i+1; j<colonne; j++)
     {
       Btemp = B[i+(j*colonne)]; //B[j][i]; 
       B[i+(j*colonne)]= B[j+(i*colonne)];//     B[j][i] = B[i][j];
       B[j+(i*colonne)]=Btemp;    // B[i][j] = Btemp;
     }
   }

}
void transposeIntSquareMatrix(int * B, int ligne, int colonne)
{

   int i, j;
   int Btemp;
   for (i=0; i<ligne; i++)
   {
     for (j=i+1; j<colonne; j++)
     {
       Btemp = B[i+(j*colonne)]; //B[j][i]; 
       B[i+(j*colonne)]= B[j+(i*colonne)];//     B[j][i] = B[i][j];
       B[j+(i*colonne)]=Btemp;    // B[i][j] = Btemp;
     }
   }

}

//The algorithm:
//Iterate over the elements of a.  For each index x, determine whether
//x is the lowest-indexed element of the cycle that includes it by
//iterating over the positions in the cycle.  If it is, then iterate
//over the cycle again, this time moving data as you go.  Continue
//with the next element of a.

//C code:

void transposeDoubleMatrixOld(double *a, int m, int n)
     //    double *a;			/*  data = int, float, ...  */
     // int m, n;			/*  M and N from above  */
{
  int x, y;			/*  array indexes  */
  double t;			/*  temp storage for array values  */ 
  //bool lowestindex();
  //int next();
  
  for (x = 0; x < m * n; x++)
    if (lowestindex(x, m, n)!=1)
      for (y = next(x, m, n); y != x; y = next(y, m, n))
	{		/*  swap a[x] and a[y]  */
	  t = a[x];
	  a[x] = a[y];
	  a[y] = t;
	}
}

void transposeIntMatrixOld(int *a, int m, int n)
     //    double *a;			/*  data = int, float, ...  */
     // int m, n;			/*  M and N from above  */
{
  int x, y;			/*  array indexes  */
  int t;			/*  temp storage for array values  */ 
  //bool lowestindex();
  //int next();
  
  for (x = 0; x < m * n; x++)
    if (lowestindex(x, m, n)!=1)
      for (y = next(x, m, n); y != x; y = next(y, m, n))
	{		/*  swap a[x] and a[y]  */
	  t = a[x];
	  a[x] = a[y];
	  a[y] = t;
	}
}


int lowestindex(int x, int m,int  n)
{
  int y;
  //int next();
  
  for (y = next(x, m, n); y != x; y = next(y, m, n))
    if (y < x)
      return (0);
  return (1);
}


int next(int x, int m, int n)
{
  return ((x % m) * n + (x / m));
}

void transposeDoubleMatrix(double *a, int m, int n)
{
  int i, j;			
  double *b;

  b = malloc(m*n*sizeof(double));
  for (i = 0; i < m ; i++)
      for (j = 0; j < n ; j++)
	{ 
	  b[(j*m)+i] = a[(i*n)+j];
	}
  memcpy(a, b, m*n*sizeof(double));
  free(b);
}

void transposeIntMatrix(int *a, int m, int n)
{
  int i, j;			
  int *b;

  b = malloc(m*n*sizeof(int));
  for (i = 0; i < m ; i++)
      for (j = 0; j < n ; j++)
	{ 
	  b[(j*m)+i] = a[(i*n)+j];
	}
  memcpy(a, b, m*n*sizeof(int));
  free(b);
}
void mem_(int *a, int *m, int *n, int *b)
{
  int     *new;
  
  new = (int *)calloc(*m, *n);
  
  if (new == NULL)
    {
      printf("\n not enough memory \n");
      exit(1);
    }
  else
    *b = (int)(new - a) + 1;
  // printf("\n memd  \n");
  
}


void unmem_(int *a, int *m, int *n, int *b)
{
  // printf("\n (un) memd  \n");

  free(a + *b - 1);
}


void memd_(double *a, int *m, int *n, int *b)
{
    double     *new;

    new = (double *) calloc(*m, *n);
#ifdef DEBUG 
    printf("\n memd  alloue %i octets \n", (long ) *m * *n);
    printf("\n le pointeur est %i \n", (long ) new);
#endif /* DEBUG */
    if (new == NULL)
        {
        printf("\n not enough memory \n");
        exit(1);
        }
    else
        *b = (new - (double *) a) + 1;
#ifdef DEBUG 
    printf("\n l'offset b est %i \n", (int ) *b);
#endif /* DEBUG */

}

void unmemd_(double *a, int *m, int *n, int *b)
{
  free(a + *b - 1);
#ifdef DEBUG 
  printf("\n (un) memd libere %i \n", (long ) (a + *b -1));
#endif /* DEBUG */
  
}
#ifdef RELEASE_C
void checksysinfot_(int *retcode)
{
	*retcode = checkSysInfo();
	if ( *retcode != 0 ) {
		*retcode = 40 + *retcode;
	}
}
#endif

int getReturnCode(int info) {
	if (info < 0) {
		return 5;
	} else if (info > 0) {
		return 6;
	} else {
		return 0;
	}
}