
!**********************************************************************
!***********       SHOCK - VORTEX INTERACTION             ************!
!********   2D- DISCONTINUOUS GALERKIN METHOD              ***********!
!*****                   SATYVIR SINGH                         *******!
!**                        2017                                    ***!
!**********************************************************************

PROGRAM SATYA_SBI
      
USE FILE_INFO
USE GAS_INFO
USE TIME_INFO
USE INPUT_INFO
USE GRID_INFO
USE FLOWFIELD_INFO
USE DIMENSIONLESS_INFO
USE COEFFICIENTS
USE MACROSCOPIC_INFO

IMPLICIT NONE

! INTERNAL PARAMETERS AND ITERATION
INTEGER:: I,J,K,M
DOUBLE PRECISION :: ERROR_VAL(4)
DOUBLE PRECISION :: SIGMA_ERROR(4)
DOUBLE PRECISION :: SIGMA_ERRORII(4)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!**************************************************!
!*****          READ PROBLEM INFORMATION      *****!
!**************************************************!
CALL INPUT_FILE
CALL GAS_INFORMATION

!**************************************************!
!*****          MAKING DATA STRUCTURE         *****!
!**************************************************!
CALL PREPARATION


!**************************************************!
!*****          GRID GENERATION              *****!
!**************************************************!

CALL GRID_GENERATION

CALL INITIALIZE 
FILENAME="INITIAL_SOLUTION"
CALL OUTPUT (FILENAME)
        
!========================

CALL MASS_MATRIX
 
!======================       
      
ERROR_VAL(:) = 1E8

DO WHILE ( CURRENT_TIME .LT. FLOW_TIME_STAR .OR. MAXVAL(ERROR_VAL).GE. 1.e-9 )

CALL TIME_STEP

ITERATION = ITERATION +1
CURRENT_TIME=CURRENT_TIME+dt

!**************************************************!
!*****            SAVE OLD STEP DATA          *****!
!**************************************************!
      u0_OLD (:,:,:)= u0(:,:,:)
      ux_OLD (:,:,:)= ux(:,:,:)
      uy_OLD (:,:,:)= uy(:,:,:)
      uxx_OLD(:,:,:)= uxx(:,:,:)
      uxy_OLD(:,:,:)= uxy(:,:,:)
      uyy_OLD(:,:,:)= uyy(:,:,:)

!**************************************************!
!*****              RUNG KUTTA METHOD         *****!
!**************************************************!
!**** RUNG-KUTTA ZERO STEP      0 th STEP
      CALL CALRIGHT

      DO M=1,RGK_STEP
!**** RUNG-KUTTA                M th STEP
      DO  i=1,nx 
      DO  j=1,ny     
      DO  k=1,4 
       u0 (i,j,k)=u0_OLD (i,j,k)+dt*RGK_COEFF(M)*right0 (i,j,k)*rmass(1)
       ux (i,j,k)=ux_OLD (i,j,k)+dt*RGK_COEFF(M)*rightx (i,j,k)*rmass(2)
       uy (i,j,k)=uy_OLD (i,j,k)+dt*RGK_COEFF(M)*righty (i,j,k)*rmass(3)
       uxx(i,j,k)=uxx_OLD(i,j,k)+dt*RGK_COEFF(M)*rightxx(i,j,k)*rmass(4)
       uxy(i,j,k)=uxy_OLD(i,j,k)+dt*RGK_COEFF(M)*rightxy(i,j,k)*rmass(5)
       uyy(i,j,k)=uyy_OLD(i,j,k)+dt*RGK_COEFF(M)*rightyy(i,j,k)*rmass(6)
      END DO
      END DO 
      END DO

      CALL BOUNDARY
      CALL LIMITER
      CALL BOUNDARY
      END DO

!**************************************************!
!*****            ERROR CALCULATION           *****!
!**************************************************!
      SIGMA_ERROR(:) =0.0
      SIGMA_ERRORII(:)=0.0
      DO i=0,nx
      DO j=0,ny
      DO K=1,4
      SIGMA_ERROR(K)= SIGMA_ERROR(K)+dabs(u0_OLD(i,j,K)-u0(i,j,K))**2.
      SIGMA_ERRORII(K)= SIGMA_ERRORII(K)+ DABS(u0_OLD(i,j,K))
      END DO
      END DO
      END DO
      DO K=1,4
      ERROR_VAL (K)= DSQRT(SIGMA_ERROR(K))
      END DO

WRITE(*,'(A,I9,2(A,ES20.8),A,E14.6)') 'ITERATION:',ITERATION,' time:', &
CURRENT_TIME*REF_TIME,' time_star:', CURRENT_TIME,'  ERROR VALUE=', &
ERROR_VAL(1)/SIGMA_ERRORII(1)



IF(MOD(ITERATION,NPRINT).EQ.0.0)THEN
FILENAME="SBI-PROBLEM"
CALL OUTPUT(FILENAME)
END IF

END DO
	                         
FILENAME="FINAL_TIME_SOLUTION"
!**************************************************!
!*****           STEADY STATE RESULTS         *****!
!**************************************************!
CALL OUTPUT(FILENAME)
      
 END PROGRAM      
   
