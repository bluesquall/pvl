!     Last change:  JEK  28 Mar 99    5:04 pm
      MODULE PVLMOD
         INTERFACE

         SUBROUTINE SIMEQN(A,B,X,IERR)
            REAL, DIMENSION(:,:), INTENT(IN) :: A    ! Coefficient matrix
            REAL, DIMENSION(:), INTENT(IN)   :: B    ! Right hand side vector
            REAL, DIMENSION(:), INTENT(OUT)  :: X    ! Solution vector
            INTEGER, INTENT(OUT)             :: IERR ! Error flag
         END SUBROUTINE SIMEQN

         SUBROUTINE WRENCH(NB,TANB,RC,RV,UA,UT)
            INTEGER, INTENT(IN) :: NB
            DOUBLE PRECISION, INTENT(IN) :: TANB,RC,RV
            DOUBLE PRECISION, INTENT(OUT) :: UA,UT
         END SUBROUTINE WRENCH

         SUBROUTINE FORCES(NBLADE,MCP,ADVCO,WAKE,RV,RC,TANBC,UASTAR,UTSTAR,    &
                           VA,CHORD,CD,G,CT,CQ,CP,KT,KQ,EFFY,RHV,CTH,IHUB)
            INTEGER, INTENT(IN) :: NBLADE,MCP,IHUB
            REAL, INTENT(IN) :: ADVCO,WAKE,RHV
            REAL, DIMENSION(:), INTENT(IN) :: RV,RC,TANBC,UASTAR,UTSTAR,VA,    &
                                              CHORD,CD,G
            REAL, INTENT(OUT) :: CT,CQ,CP,KT,KQ,EFFY,CTH
         END SUBROUTINE FORCES

         REAL FUNCTION VOLWK(XR,XVA)
            REAL, DIMENSION(:), INTENT(IN) :: XR,XVA
         END FUNCTION VOLWK

         END INTERFACE
      END MODULE PVLMOD
