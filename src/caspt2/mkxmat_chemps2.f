************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
************************************************************************
      SUBROUTINE MKXMAT_CHEMPS2(TORB,XMAT)
      IMPLICIT NONE
#include "caspt2.fh"
* Build the full active-space transformation matrix from the separate
* symmetry and RAS blocks stored in TORB.
      REAL*8 TORB(*),XMAT(NASHT,NASHT)
      INTEGER ISYM,NI,NR1,NR2,NR3,NS,ITOSTA,ITOEND,ITO
      INTEGER ISTART,I,J,IR,JR

      IF (NASHT.LE.0) RETURN

      ITOEND=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NR1=NRAS1(ISYM)
        NR2=NRAS2(ISYM)
        NR3=NRAS3(ISYM)
        NS=NSSH(ISYM)
        ITOSTA=ITOEND+1
        ITOEND=ITOEND+NI**2+NR1**2+NR2**2+NR3**2+NS**2
        ITO=ITOSTA+NI**2

        ISTART=NAES(ISYM)
        DO JR=1,NR1
          J=ISTART+JR
          DO IR=1,NR1
            I=ISTART+IR
            XMAT(I,J)=TORB(ITO)
            ITO=ITO+1
          END DO
        END DO

        ISTART=NAES(ISYM)+NR1
        DO JR=1,NR2
          J=ISTART+JR
          DO IR=1,NR2
            I=ISTART+IR
            XMAT(I,J)=TORB(ITO)
            ITO=ITO+1
          END DO
        END DO

        ISTART=NAES(ISYM)+NR1+NR2
        DO JR=1,NR3
          J=ISTART+JR
          DO IR=1,NR3
            I=ISTART+IR
            XMAT(I,J)=TORB(ITO)
            ITO=ITO+1
          END DO
        END DO
      END DO

      END SUBROUTINE MKXMAT_CHEMPS2
