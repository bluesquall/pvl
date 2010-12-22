!     Last change:  JEK   2 Mar 99   11:44 am
      SUBROUTINE SIMEQN(A,B,X,IERR)
!-----This is a Fortran 90 version of Dave Greeley's FACTOR & SUBST (combined)--
      IMPLICIT NONE

!------------------- Declare variables in argument list ------------------------
      REAL, DIMENSION(:,:), INTENT(INOUT) :: A    ! Coefficient matrix
      REAL, DIMENSION(:), INTENT(INOUT)   :: B    ! Right hand side vector
      REAL, DIMENSION(:), INTENT(OUT)     :: X    ! Solution vector
      INTEGER, INTENT(OUT)                :: IERR ! Error flag

!-------------------- Declare local variables ----------------------------------
      REAL, DIMENSION(:), ALLOCATABLE    :: D       ! Row swapping storage
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPIVOT  ! Row swapping index
      INTEGER :: NEQ,I,J,NM1,K,KM1,KP1,IP,IPK,NP1MK
      REAL ::  ROWMAX,COLMAX,AWIKOV,SUMM,RATIO

!--------------------Allocate local arrays--------------------------------------
      NEQ=SIZE(B)
      ALLOCATE ( D(NEQ),IPIVOT(NEQ) )
      IERR=1

!-----Find |maximum| element in each row, and exit if a zero row is detected----
      IERR=1            ! Initialize error flag to 1 (denotes bad matrix)-------
      DO I=1,NEQ
         IPIVOT(I)=I
         ROWMAX=0.0
         DO J=1,NEQ
            ROWMAX=MAX(ROWMAX,ABS(A(I,J)))
         END DO
         IF(ROWMAX==0.0) RETURN     ! IERR=1 Matrix is singular ----------------
         D(I)=ROWMAX
      END DO

      NM1=NEQ-1
      IF(NM1>0) THEN                ! Otherwise special case of one equation----
         DO K=1,NM1
            J=K
            KP1=K+1
            IP=IPIVOT(K)
            COLMAX=ABS(A(IP,K))/D(IP)
            DO I=KP1,NEQ
               IP=IPIVOT(I)
               AWIKOV=ABS(A(IP,K))/D(IP)
               IF(AWIKOV>COLMAX) THEN
                  COLMAX=AWIKOV
                  J=I
               END IF
            END DO
            IF(COLMAX==0.0) RETURN  ! IERR=1 Matrix is singular ----------------
            IPK=IPIVOT(J)
            IPIVOT(J)=IPIVOT(K)
            IPIVOT(K)=IPK
            DO I=KP1,NEQ
               IP=IPIVOT(I)
               A(IP,K)=A(IP,K)/A(IPK,K)
               RATIO=-A(IP,K)
               DO J=KP1,NEQ
                  A(IP,J)=RATIO*A(IPK,J)+A(IP,J)
               END DO
            END DO
         END DO
         IF(A(IP,NEQ).EQ.0.) RETURN   ! IERR=1 Matrix is singular --------------
      END IF
      IERR=0                          ! Matrix survived singular tests ---------

!------------------Back substitute to obtain solution (X) ----------------------

      IF(NEQ==1) THEN         ! Special case of one equation again--------------
         X(1)=B(1)/A(1,1)
      ELSE

         IP=IPIVOT(1)
         X(1)=B(IP)
         DO K=2,NEQ
            IP=IPIVOT(K)
            KM1=K-1
            SUMM=0.0
            DO J=1,KM1
               SUMM=A(IP,J)*X(J)+SUMM
            END DO
            X(K)=B(IP)-SUMM
         END DO
         X(NEQ)=X(NEQ)/A(IP,NEQ)
         K=NEQ
         DO NP1MK=2,NEQ
            KP1=K
            K=K-1
            IP=IPIVOT(K)
            SUMM=0.0
            DO J=KP1,NEQ
               SUMM=A(IP,J)*X(J)+SUMM
            END DO
            X(K)=(X(K)-SUMM)/A(IP,K)
         END DO
      END IF

      DEALLOCATE (D,IPIVOT)
      RETURN
      END SUBROUTINE SIMEQN

