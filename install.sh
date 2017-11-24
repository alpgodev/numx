#!/bin/sh
echo ""
echo "***--- NumX v1.0 - Mac OS X installation ---***"
echo "               (c) 2014 - NumX                 "
echo ""
echo ""
echo "***--- Generating C/JAVA APIs, .h, javadoc ... ---***"
echo ""
#rm api/c/*/*.c api/jni/*/*.c # ne pas effacer mcvar.c
cd auto
ant clean
ant C_API
ant Java_API
cd ../src/api/include/jni
ant cleancompile -f buildnumx.xml
cd ../../..
echo ""
echo "***--- Generating Maths Libraries (Fotran/C) ---***"
echo "***--- output directory: numx/lib/math/"
echo ""
echo "***--- Compiling ALGLIN ... ---***"
echo ""
cd math/alglin
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling ANALYSIS ... ---***"
echo ""
cd ../analysis
make clean all || exit 1
make clean || exit 1
cd mcvar
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling CALELM ... ---***"
echo ""
cd ../../calelm
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling CLUSTER ... ---***"
echo ""
cd ../cluster
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling DATA ... ---***"
echo ""
cd ../missing_data
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling OPTIM ... ---***"
echo ""
cd ../optim
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling OPTIM PORT ... ---***"
echo ""
cd ../optim_portfolio
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling RND ... ---***"
echo ""
cd ../rnd
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling SIMUL ... ---***"
echo ""
cd ../simul
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling Maths Libraries (Fortran/C) OK ---***"
echo ""
echo ""
echo "***--- Generating API Libraries (C/C++) ---***"
echo "***--- output directory: numx/lib/api/"
cd ../../api/c
echo ""
echo "***--- Compiling Linear-Algebra Module (C/C++) ... ---***"
echo ""
cd linear-algebra
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling MODELING Module (C/C++) ... ---***"
echo ""
cd ../modeling
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling OPTIMIZATION Module (C/C++) ... ---***"
echo ""
cd ../optimization
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling REPORTING Module (C/C++) ... ---***"
echo ""
cd ../reporting
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling SIMULATION Module (C/C++) ... ---***"
echo ""
cd ../simulation
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- API Libraries (C/C++) for MacOSX OK ---***"
echo ""
echo "***--- Generating API Libraries (Java) ---***"
echo "***--- output directory: numx/lib/api/"
cd ../../jni
echo ""
echo "***--- Compiling Linear-Algebra Module (Java) ... ---***"
echo ""
cd linear-algebra
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling MODELING Module (Java) ... ---***"
echo ""
cd ../modeling
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling OPTIMIZATION Module (Java) ... ---***"
echo ""
cd ../optimization
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling REPORTING Module (Java) ... ---***"
echo ""
cd ../reporting
make clean all || exit 1
make clean || exit 1
echo ""
echo "***--- Compiling SIMULATION Module (Java) ... ---***"
echo ""
cd ../simulation
make clean all || exit 1
make clean || exit 1
cd ../../include/jni/
ant createjar -f buildnumx.xml
echo ""
echo "***--- API Libraries (Java) for MacOSX OK ---***"
echo ""
cd ../../../..
echo ""
echo "***--- BUILD SUCCESSFUL ---***"
echo ""
echo "**********************************************************"
echo "***                                                    ***"
echo "***  Files generated:                                  ***"
echo "***      numx/lib/api/        : *.a *.so *.dylib       ***"
echo "***      numx/lib/numx.jar    : JAR package            ***"
echo "***      numx/lib/api/javadoc : Javadoc for NumX       ***"
echo "***      numx/lib/api/include : C header               ***"
echo "***                                                    ***"
echo "***  Run test programs:                                ***"
echo "***      ./test.sh                                     ***"
echo "***                                                    ***"
echo "**********************************************************"
#done.


