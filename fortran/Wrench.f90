!     Last change:  JEK  25 Apr 2001    8:48 pm
      SUBROUTINE WRENCH(NB,TANB,RC,RV,UA,UT)
      IMPLICIT NONE
!------------------Declare variables in argument list --------------------------
      INTEGER, INTENT(IN) :: NB
      DOUBLE PRECISION, INTENT(IN) :: TANB,RC,RV
      DOUBLE PRECISION, INTENT(OUT) :: UA,UT

!------------------Declare local variables -------------------------------------
      DOUBLE PRECISION :: C25=0.25D00, ONE=1.0D00, C15=1.5D00, TWO=2.0D00,     &
                          THREE=3.0D00, NINE=9.0D00, C24=24.0D00
      DOUBLE PRECISION :: BL,XG,ETA,H,XS,T,V,W,AE,U,R,XX,Y,Z,AF,AA,RATIO,AG,AB

      BL=DBLE(NB)

!-----Return infinite blade result if NB>20 JEK 9/19/98 ----------------
      IF(NB.GT.20) THEN
         IF(RC.GT.RV) THEN
            UA=0.0
            UT=BL*(RC-RV)/RC
         ELSE
            UA=-BL*(RC-RV)/(RV*TANB)
            UT=0.0
         END IF
         RETURN
      END IF
!-----End of infinite blade patch --------------------------------------

      XG=ONE/TANB
      ETA=RV/RC
      H=XG/ETA
      XS=ONE+H**2
      T=SQRT(XS)
      V=ONE+XG**2
      W=SQRT(V)
      AE=T-W
      U=EXP(AE)
      R=(((T-ONE)/H*(XG/(W-ONE)))*U)**BL
      XX=(ONE/(TWO*BL*XG))*((V/XS)**C25)
      Y=((NINE*XG**2)+TWO)/(V**C15)+((THREE*H**2-TWO)/(XS**C15))
      Z=ONE/(C24*BL)*Y
      IF(H.GE.XG) THEN
         AF=ONE+ONE/(R-ONE)
         AA=XX*(ONE/(R-ONE)-Z*LOG(AF))
         UA=TWO*BL**2*XG*H*(ONE-ETA)*AA
         UT=BL*(ONE-ETA)*(ONE+TWO*BL*XG*AA)
      ELSE
         IF(R.GT.1.0D-12) THEN
            RATIO=ONE/(ONE/R-ONE)
         ELSE
            RATIO=0.0
         END IF
         AG=ONE+RATIO
         AB=-XX*(RATIO+Z*LOG(AG))
         UA=BL*XG*(ONE-ONE/ETA)*(ONE-TWO*BL*XG*AB)
         UT=TWO*BL**2*XG*(ONE-ETA)*AB
      END IF
      RETURN
      END SUBROUTINE WRENCH

      REAL FUNCTION VOLWK(XR,XVA)
      USE DUCKMOD
      IMPLICIT NONE
      REAL :: YDX
      INTEGER :: NX,N
      REAL, DIMENSION(:), INTENT(IN) :: XR,XVA
      REAL, DIMENSION(:), ALLOCATABLE :: Y
      REAL, DIMENSION(:,:), ALLOCATABLE :: VWCUB
      NX=SIZE(XR)
      ALLOCATE ( Y(NX),VWCUB(NX-1,5) )
      Y(:)=XR(:)*XVA(:)
      CALL UGLYDK(0,0,XR,Y,0.0,0.0,VWCUB)
      CALL INTDK1(XR(1),XR(NX),YDX,VWCUB)
      VOLWK=2.0*YDX/(1.0-XR(1)**2)
      DEALLOCATE (Y,VWCUB)
      RETURN
      END FUNCTION VOLWK

