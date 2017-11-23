# numx
Library for Portfolio Optimization and Risk Analysis.

VERSION 1.0  :  October 2017

### Currently supported systems are:

- Windows
- Linux
- Mac OS X

### COMPILING NumX (C and JAVA APIs)

NumX can be installed with make. Configuration have to be set in the make.inc file.

Running ```./install.sh``` should do the job ;)

### INSTALLING NumX

1. Add “numx.jar" to your CLASSPATH

2. By default, the dynamic libraries are extracted from the jar file
   to the default temporary directory and loaded from there. If you
   don't want to do this, extract the dynamic library from numx.jar
   and copy it somewhere where it can be found.

   For Linux, use LD_LIBRARY_PATH, for Windows, PATH
 
### TESTING NumX

We recommend that you run the testing. 

```> ./test.sh```
