C FILE: 'init_seed_linux.f'
C ---------------------------------------------------------------------
C
C       Seeding function for UNIX FORTRAN random number generators
C
C          The integer function SEED returns the product of the number
C     of minutes and seconds of the system clock, to serve as seed
C     for the random number generating functions: IRAND, RAND, and
C     DRAND.
C
C     Note: this function uses the UNIX FORTRAN library function
C           TIME to get the system time.
C
C     Command Syntax:  This function should be called as the
C                      argument to one of the system random
C                      number functions.  For example to seed
C                      the function RAND, the following command
C                      should be used:
C
C                              MUMBLE = RAND( SEED() )
C
C                      REMEMBER that SEED and RAND (or IRAND or DRAND)
C                      must be declared in your calling routine!
C
C
C   Variables used:
C        INIT = Logical flag which is true only on the first call to
C               SEED.
C
C        OLDSED = Previous value of the SEED function.
C
C ---------------------------------------------------------------------

      INTEGER FUNCTION INIT_SEED_LINUX()

      INTEGER TIME, IRAND, OLDSED
      LOGICAL INIT
      SAVE INIT, OLDSED
C    Initialize flag INIT to be true on the first run only.
      DATA INIT /.TRUE./

      IF (INIT) THEN
C        Call system function that returns time since 00:00 GMT 1/1/1970
C        and use that as an initial seed.
           OLDSED = TIME()
      ELSE
C         Else use old seed as to seed the next value.
           OLDSED = IRAND(OLDSED)
      ENDIF

C    Reset INIT to be false on any later calls to SEED
      INIT = .FALSE.

      INIT_SEED_LINUX = OLDSED
      RETURN
      END
