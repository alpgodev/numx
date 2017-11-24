/*=======================================================================

     Java-C-Fortran utils headers

=======================================================================*/

/* C header functions are declared here */   
#ifdef WINDOWS 
#include <math.h>
#endif /*WINDOWS*/

extern void DisplayMatrixDouble( char *, int, int, double *);
extern void DisplayMatrixDoubleInFile( char *, int, int, double *, FILE *);
extern void DisplayMatrixInt( char *, int, int ,int *);
extern void DisplayMatrixIntInFile( char *, int, int ,int *, FILE *);
extern void DisplayVectorDouble( char *, int, double *);
extern void DisplayVectorDoubleInFile( char *, int, double *, FILE *);
extern void DisplayVectorInt( char *, int, int *);
extern void DisplayVectorIntInFile( char *, int, int *, FILE *);

extern void transposeDoubleMatrix(double * , int, int );
extern void transposeIntMatrix(int *, int, int);

extern void transposeDoubleMatrix(double *, int, int);
extern void transposeIntMatrixOld(int *, int, int);
extern void transposeDoubleMatrixOld(double *, int, int);
extern int  lowestindex(int, int, int);
extern int  next(int, int, int);
extern int checkSysInfo();
extern int getReturnCode(int);
