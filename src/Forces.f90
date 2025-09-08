!     Last change:  JE   24 Apr 2001   11:41 am
!***************     COPYRIGHT (C) 2001  JUSTIN E. KERWIN    ******************
      SUBROUTINE FORCES(NBLADE,MCP,ADVCO,WAKE,RV,RC,TANBC,UASTAR,UTSTAR,       &
                       VA,CHORD,CD,G,CT,CQ,CP,KT,KQ,EFFY,RHV,CTH,IHUB)
      IMPLICIT NONE
!--------------------- Declare the arguments -----------------------------------
      INTEGER, INTENT(IN) :: NBLADE,MCP,IHUB
      REAL, INTENT(IN) :: ADVCO,WAKE,RHV
      REAL, DIMENSION(:), INTENT(IN) :: RV,RC,TANBC,UASTAR,UTSTAR,VA,CHORD,CD,G
      REAL, INTENT(OUT) :: CT,CQ,CP,KT,KQ,EFFY,CTH

!--------------------- Declare the local variables -----------------------------
      REAL, PARAMETER :: PI=3.1415927E00, TWO=2.0E00, FOUR=4.0E00, EIGHT=8.0E00
      REAL :: DR,VSTAR,VTSTAR,VASTAR,VSTRSQ,DVISC,FKJ
      INTEGER :: M
      LOGICAL :: CD_LD

      CD_LD=.TRUE.  ! Default: Input CD interpreted as viscous drag coefficient
      IF(CD(1)>1.0) CD_LD=.FALSE. ! CD(1)>1 signals that input is L/D ----------
      CT=0.0
      CQ=0.0
      DO M=1,MCP
         DR=RV(M+1)-RV(M)
         VTSTAR=VA(M)/TANBC(M)+UTSTAR(M)
         VASTAR=VA(M)+UASTAR(M)
         VSTRSQ=VTSTAR**2+VASTAR**2
         VSTAR=SQRT(VSTRSQ)
         IF(CD_LD) THEN  ! Interpret CD as viscous drag coefficient, Cd---------
            DVISC=(VSTRSQ*CHORD(M)*CD(M))/(TWO*PI)
         ELSE            ! Interpret CD as the lift/drag ratio L/D -------------
            FKJ=VSTAR*G(M)
            DVISC=FKJ/CD(M)
         END IF
         CT=CT+(VTSTAR*G(M)-DVISC*VASTAR/VSTAR)*DR
         CQ=CQ+(VASTAR*G(M)+DVISC*VTSTAR/VSTAR)*RC(M)*DR
      END DO

      IF(IHUB/=0) THEN  ! Add hub vortex drag if hub image is present ----------
         CTH=0.5*(LOG(1.0/RHV)+3.0)*(REAL(NBLADE)*G(1))**2
      ELSE
         CTH=0.0
      END IF

      CT=CT*FOUR*REAL(NBLADE)-CTH
      CQ=CQ*TWO*REAL(NBLADE)
      CP=CQ*TWO*PI/ADVCO
      KT=CT*ADVCO**2*PI/EIGHT
      KQ=CQ*ADVCO**2*PI/EIGHT
      EFFY=CT*WAKE/CP

      RETURN
      END
