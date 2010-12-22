!     Last change:  JEK  30 Apr 2001    5:48 pm
!     *********** COPYRIGHT (C) 2001  JUSTIN E. KERWIN *************************
      USE PVLMOD
      USE DUCKMOD
!------------------------- Declare the Variables -------------------------------
      IMPLICIT NONE
      CHARACTER*36 :: FNAME,LABEL
      CHARACTER*72 :: TITLE
      INTEGER :: MT,NX,ITER,NBLADE,N,M,KTRY,IERR,IHUB
      REAL :: KT,KQ,DEL,HRR,RCWG,RM,DTANB,EDISK,ADVCO,CTDES,HR,HT,CRP,WAKE,    &
              CQ,CP,EFFY,HRF,CTH,RHV
      REAL, PARAMETER :: PI=3.1415927E00, TOL=0.000005, R2D=57.29578E00
      DOUBLE PRECISION :: TANBIW,RCW,RVW,UAIF,UTIF
      REAL, DIMENSION(:), ALLOCATABLE :: XR,XCHD,XCD,XVA,XVT,XRC,RV,TANBV,RC,  &
                                         TANBC,VAV,VTV,VAC,VTC,TANBIV,TANBIC,  &
                                         UAW,UTW,B,G,UASTAR,UTSTAR,T,CT,       &
                                         TANBXV,TANBXC,VBAV,VBAC,CD,CDC
      REAL, DIMENSION(:,:), ALLOCATABLE :: CHCUB,CDCUB,VACUB,VTCUB,UAHIF,      &
                                           UTHIF,A

!------------------ Start reading the input data -------------------------------
      WRITE(*,'(A)') ' ENTER INPUT FILE NAME.....    '
      READ(*,'(A)') FNAME
      OPEN(2,FILE=FNAME,STATUS='OLD',FORM='FORMATTED')
      READ(2,'(A)') TITLE  ! Title describing data file
      READ(2,*) MT   ! Number of vortex lattice panels
      READ(2,*) ITER ! Number of iterations to align wake
      READ(2,*) IHUB ! Hub image flag. IHUB=0 : No hub image, IHUB=1 : Image hub
      READ(2,*) RHV  ! Hub vortex radius/Hub radius. Only used if IHUB=1
      READ(2,*) NX   ! Number of radii used to specify the input data

!------------- Allocate all the arrays before reading rest of input-------------
      ALLOCATE ( XR(NX),XCHD(NX),XCD(NX),XVA(NX),XVT(NX),XRC(NX) )
      ALLOCATE ( CHCUB(NX-1,5),CDCUB(NX-1,5),VACUB(NX-1,5),VTCUB(NX-1,5) )
      ALLOCATE ( RV(MT+1),TANBV(MT+1),RC(MT),TANBC(MT),VAV(MT+1),VTV(MT+1),    &
                 VAC(MT),VTC(MT),TANBIV(MT+1),TANBIC(MT),UAW(MT+1),UTW(MT+1),  &
                 UAHIF(MT,MT),UTHIF(MT,MT),A(MT,MT),B(MT),G(MT),UASTAR(MT),    &
                 UTSTAR(MT),T(ITER),CT(ITER),TANBXV(MT+1),TANBXC(MT),          &
                 VBAV(MT+1),VBAC(MT),CD(MT),CDC(MT) )

!----------------All arrays allocated. read in rest of input data --------------
      READ(2,*) NBLADE     ! Number of blades
      READ(2,*) ADVCO      ! Advance coefficient based on ship speed
      READ(2,*) CTDES      ! Desires thrust coefficient CT (based on ship speed)
      READ(2,*) HR         ! Unloading ratio at hub
      READ(2,*) HT         ! Unloading ratio at tip
      READ(2,*) CRP        ! Tangential velocity cancellation factor
      READ(2,'(A)') LABEL  ! Alphanumeric label for output
!-------------------------------------------------------------------------------
!     XR=Input radii r/R, XCHD=Input chord length c/D, XCD=Input viscous drag
!     coefficient, Cd or Lift/Drag ratio, XVA,XVT=Input axial and tangential
!     velocities, Va/Vs, V_t/Vs
!-------------------------------------------------------------------------------
      Do N=1,NX
      READ(2,*)XR(N),XCHD(N),XCD(N),XVA(N),XVT(N)
!      READ(2,*) (XR(N),XCHD(N),XCD(N),XVA(N),XVT(N),N=1,NX)
      End do
      CLOSE(2)

!-----Compute volumetric mean inflow velocity ratio VA/VS ----------------------
      WAKE=VOLWK(XR,XVA)

!-----Spline chord over radius using square root stretched coordinates----------
      XRC(:)=1.0-SQRT(1.0-XR(:))
      CALL UGLYDK(0,0,XRC,XCHD,0.0,0.0,CHCUB)

!-----Spline Drag Coefficient Cd, Inflow Vx, Vt using radial coordinate directly
      CALL UGLYDK(0,0,XR,XCD,0.0,0.0,CDCUB)
      CALL UGLYDK(0,0,XR,XVA,0.0,0.0,VACUB)
      CALL UGLYDK(0,0,XR,XVT,0.0,0.0,VTCUB)

!-----Compute cosine spaced vortex radii and get Va,Vt,tanB,Vt*tanB/Va----------
      DEL=PI/(2.0*REAL(MT))
      HRR=0.5*(XR(NX)-XR(1))
      DO M=1,MT+1
         RV(M)=XR(1)+HRR*(1.0-COS(REAL(2*(M-1))*DEL))
         CALL EVALDK(RV(M),VAV(M),VACUB)
         CALL EVALDK(RV(M),VTV(M),VTCUB)
         TANBV(M)=VAV(M)/((PI*RV(M)/ADVCO)+VTV(M))
         VBAV(M)=VTV(M)*TANBV(M)/VAV(M)
      END DO

!-----Cosine spaced control point radii: Evaluate c/D,Va,Vt,tanB,Cd,Vt*tanB/Va -
      DO M=1,MT
         RC(M)=XR(1)+HRR*(1.0-COS(REAL(2*M-1)*DEL))
         RCWG=1.0-SQRT(1.0-RC(M))
         CALL EVALDK(RCWG,CDC(M),CHCUB)
         CALL EVALDK(RC(M),VAC(M),VACUB)
         CALL EVALDK(RC(M),VTC(M),VTCUB)
         TANBC(M)=VAC(M)/((PI*RC(M)/ADVCO)+VTC(M))
         CALL EVALDK(RC(M),CD(M),CDCUB)
         VBAC(M)=VTC(M)*TANBC(M)/VAC(M)
      END DO

!-----First estimate of tanBi based on 90 percent of actuator disk efficiency --
      EDISK=1.8/(1.0+SQRT(1.0+CTDES/WAKE**2))
      TANBXV(:)=TANBV(:)*SQRT(WAKE/(VAV(:)-VBAV(:)))/EDISK  ! Lerbs optimum-----
      TANBXC(:)=TANBC(:)*SQRT(WAKE/(VAC(:)-VBAC(:)))/EDISK

!-----Unload hub and tip as specified by input HR and HT -----------------------
      RM=0.5*(XR(1)+XR(NX))  ! Mid-radius. Unloading is quadratic, starting here
      DO M=1,MT+1
         IF(RV(M).LT.RM) THEN
            HRF=HR
         ELSE
            HRF=HT
         END IF
         DTANB=HRF*(TANBXV(M)-TANBV(M))*((RV(M)-RM)/(XR(1)-RM))**2
         TANBXV(M)=TANBXV(M)-DTANB
      END DO

      DO M=1,MT
         IF(RC(M).LT.RM) THEN
            HRF=HR
         ELSE
            HRF=HT
         END IF
         DTANB=HRF*(TANBXC(M)-TANBC(M))*((RC(M)-RM)/(XR(1)-RM))**2
         TANBXC(M)=TANBXC(M)-DTANB
      END DO

!-------------------------------------------------------------------------------
!     Iterations to scale tanBi to get desired value of thrust coefficient
!-------------------------------------------------------------------------------

      DO KTRY=1,ITER
         IF(KTRY.EQ.1) THEN
            T(KTRY)=1.0    ! T(KTRY) is the scale factor to apply to tanBi
         ELSE IF(KTRY.EQ.2) THEN
            T(KTRY)=1.0+(CTDES-CT(1))/(5.0*CTDES) ! Guess for second iteration
         ELSE IF(KTRY.GT.2) THEN
            T(KTRY)=T(KTRY-1)+(T(KTRY-1)-T(KTRY-2))*(CTDES-CT(KTRY-1))/        &
                    (CT(KTRY-1)-CT(KTRY-2))  ! Secant method for remaining iters
         END IF

            TANBIV(:)=T(KTRY)*TANBXV(:)    ! Scale tanBi at the vortex radii
            TANBIC(:)=T(KTRY)*TANBXC(:)    ! Scale tanbi at the control points

!-------------------------------------------------------------------------------
!        Compute axial and tangential horseshoe influence coefficients         !
!-------------------------------------------------------------------------------

         DO M=1,MT
            RCW=RC(M)
            DO N=1,MT+1
!--------------Induction of trailing vortices shed at RV(N)---------------------
               TANBIW=TANBIV(N)
               RVW=RV(N)
               CALL WRENCH(NBLADE,TANBIW,RCW,RVW,UAIF,UTIF)
               UAW(N)=-UAIF/(2.0*(RC(M)-RV(N)))
               UTIF=UTIF*CRP  ! Note if CRP=0, the tangential velocity is zero--
               UTW(N)=UTIF/(2.0*(RC(M)-RV(N)))
!--------------Induction of corresponding hub-image trailing vortices (if any)--
               IF(IHUB/=0) THEN
                  RVW=XR(1)**2/RV(N)
                  TANBIW=TANBIV(1)*RV(1)/RVW
                  CALL WRENCH(NBLADE,TANBIW,RCW,RVW,UAIF,UTIF)
                  UAW(N)=UAW(N)+UAIF/(2.0*(RC(M)-RVW))
                  UTIF=UTIF*CRP
                  UTW(N)=UTW(N)-UTIF/(2.0*(RC(M)-RVW))
               END IF
            END DO
!-----------Final step in building influence functions--------------------------
            DO N=1,MT
               UAHIF(M,N)=UAW(N+1)-UAW(N)
               UTHIF(M,N)=UTW(N+1)-UTW(N)
            END DO
         END DO

!-------------------------------------------------------------------------------
!        Solve simultaneous equations for circulation strengths G(M)           !
!-------------------------------------------------------------------------------
         DO M=1,MT
            B(M)=VAC(M)*((TANBIC(M)/TANBC(M))-1.0)      ! Right-hand side
            DO N=1,MT
               A(M,N)=UAHIF(M,N)-UTHIF(M,N)*TANBIC(M)   ! Coefficient matrix
            END DO
         END DO

         CALL SIMEQN(A,B,G,IERR)           ! Simultaneous equation solver
         IF(IERR/=0) EXIT                  ! Error return for singular matrix


!-------------------------------------------------------------------------------
!        Evaluate the induced velocities from the circulation GM)              !
!-------------------------------------------------------------------------------
         DO M=1,MT
            UASTAR(M)=0.0
            UTSTAR(M)=0.0
            DO N=1,MT
               UASTAR(M)=UASTAR(M)+G(N)*UAHIF(M,N)
               UTSTAR(M)=UTSTAR(M)+G(N)*UTHIF(M,N)
            END DO
         END DO

!-------------------------------------------------------------------------------
!        Compute the forces and test if Ct has converged to desired value      !
!-------------------------------------------------------------------------------
         CALL FORCES(NBLADE,MT,ADVCO,WAKE,RV,RC,TANBC,UASTAR,UTSTAR,VAC,       &
                  CDC,CD,G,CT(KTRY),CQ,CP,KT,KQ,EFFY,RHV,CTH,IHUB)
         WRITE(*,'(I5,'' CT='',F10.5,'' DESIRED VALUE='',F10.5)') KTRY,        &
               CT(KTRY),CTDES
         IF(ABS(CT(KTRY)-CTDES)<TOL) EXIT

      END DO

!-----Stop run if matrix is sigular---------------------------------------------
      IF(IERR/=0) THEN
         WRITE(*,'(A)') ' MATRIX SINGULAR. RUN TERMINATED..... '
         STOP
      ELSE

!-------------------------------------------------------------------------------
!                     Output results to Tecplot file                           !
!-------------------------------------------------------------------------------

         WRITE(*,'(//'' EFFICIENCY ='',F8.4)') EFFY
         WRITE(*,'('' Kt, Kq'',F8.4,F8.5)') KT,KQ
         WRITE(*,'('' HUB DRAG COEFFICIENT Cth='',F8.4)') CTH

         OPEN(1,FILE='APLOT.PLT',STATUS='UNKNOWN',FORM='FORMATTED')
         WRITE(1,'(A)') ' VARIABLES="R","G","VA","VT","UA","UT","BETA","BETAI",&
                          &"CDC","CD" '
         WRITE(1,'('' TEXT X=0.5, Y=0.50, T="    Ct='',F8.4,'' "'')')          &
               CT(KTRY)
         WRITE(1,'('' TEXT X=0.5, Y=0.46, T="    Cp='',F8.4,'' "'')')          &
               CP
         WRITE(1,'('' TEXT X=0.5, Y=0.42, T="    Kt='',F8.4,'' "'')')          &
               KT
         WRITE(1,'('' TEXT X=0.5, Y=0.38, T="    Kq='',F8.4,'' "'')')          &
               KQ
         WRITE(1,'('' TEXT X=0.5, Y=0.34, T=" Va/Vs='',F8.4,'' "'')')          &
              WAKE
         WRITE(1,'('' TEXT X=0.5, Y=0.30, T="     E='',F8.4,'' "'')')          &
              EFFY
         WRITE(1,'('' TEXT X=0.5, Y=0.26, T="'',A,'' " '')')                   &
              TITLE

         WRITE(1,'(F10.5,F10.6,4F10.5,2F10.3,2F10.5)') (RC(M),G(M),VAC(M),     &
               VTC(M),UASTAR(M),UTSTAR(M),R2D*ATAN(TANBC(M)),                  &
               R2D*ATAN(TANBIC(M)),CDC(M),CD(M),M=1,MT)
         CLOSE(1)
      END IF

      STOP
      END

