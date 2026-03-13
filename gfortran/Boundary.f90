  
SUBROUTINE boundary     
USE GAS_INFO
USE TIME_INFO
USE INPUT_INFO
USE GRID_INFO
USE FLOWFIELD_INFO
USE COEFFICIENTS
USE DIMENSIONLESS_INFO
IMPLICIT NONE
INTEGER:: I,J,K                
!** DOUBLE PRECISION:: U,V,P,RHO


!*****************************************
!******* LEFT BOUNDARY
!*****************************************

IF(LEFT_BOUNDARY.EQ."INFLOW")THEN

DO j=0,ny+1   
DO k=1,4 
  u0(0,j,k)=u0(1,j,k) 
  ux(0,j,k)=0.0   
  uy(0,j,k)=0.0 
  uxx(0,j,k)=0.0
  uxy(0,j,k)=0.0
  uyy(0,j,k)=0.0 
END DO 
END DO  

END IF


!*****************************************
!******* RIGHT BOUNDARY
!***************************************** 

IF(RIGHT_BOUNDARY.EQ."OUTFLOW")THEN

DO j=0,ny+1 
  DO k=1,4   
   u0(nx+1,j,K) = u0(nx,j,K)   
   ux(nx+1,j,k)=0.0    
   uy(nx+1,j,k)=0.0 
   uxx(nx+1,j,k)=0.0 
   uxy(nx+1,j,k)=0.0 
   uyy(nx+1,j,k)=0.0
END DO
END DO

!******
ELSEIF(RIGHT_BOUNDARY.EQ."NO-SLIP")THEN
!******

DO j=0,ny+1 
  DO k=1,4   
   u0(nx+1,j,1) = u0(nx,j,1)   
   u0(nx+1,j,2) = 0.0
   u0(nx+1,j,3) = 0.0
   u0(nx+1,j,4) = u0(nx,j,4)   
   ux(nx+1,j,k)=0.0    
   uy(nx+1,j,k)=0.0 
   uxx(nx+1,j,k)=0.0 
   uxy(nx+1,j,k)=0.0 
   uyy(nx+1,j,k)=0.0
 END DO
END DO

!******
ELSEIF(RIGHT_BOUNDARY.EQ."REFLECTIVE")THEN
!******

DO j=0,ny+1 
  DO k=1,4   
   u0(nx+1,j,1) = u0(nx,j,1)   
   u0(nx+1,j,2) = 0.0              ! No penetration
   u0(nx+1,j,3) = u0(nx,j,3)       ! Reflective condition
   u0(nx+1,j,4) = u0(nx,j,4)   
   ux(nx+1,j,k)=0.0    
   uy(nx+1,j,k)=0.0 
   uxx(nx+1,j,k)=0.0 
   uxy(nx+1,j,k)=0.0 
   uyy(nx+1,j,k)=0.0
 END DO
END DO


!******
ELSEIF(RIGHT_BOUNDARY.EQ."SYMMETRY")THEN
!******
    
DO j=0,ny+1 
  DO k=1,4   
   u0(nx+1,j,1) = u0(nx,j,1)   
   u0(nx+1,j,2) = u0(nx,j,2)  
   u0(nx+1,j,3) = -u0(nx,j,3)  
   u0(nx+1,j,4) = u0(nx,j,4)   
   ux(nx+1,j,k)=0.0    
   uy(nx+1,j,k)=0.0 
   uxx(nx+1,j,k)=0.0 
   uxy(nx+1,j,k)=0.0 
   uyy(nx+1,j,k)=0.0
 END DO
END DO    
    
END IF

!*****************************************
!******* UPPER BOUNDARY
!*****************************************    

IF(UPPER_BOUNDARY.EQ."OUTFLOW")THEN

DO i=1,nx 
DO k=1,4
  u0(i,ny+1,k)=u0(i,ny,k)
  ux(i,ny+1,k)=0.0 
  uy(i,ny+1,k)=0.0 
  uxx(i,ny+1,k)=0.0 
  uxy(i,ny+1,k)=0.0    
  uyy(i,ny+1,k)=0.0 
END DO
END DO  

!***
ELSEIF(UPPER_BOUNDARY.EQ."NO-SLIP")THEN 
!*** 

DO i=1,nx 
 DO k=1,4
  u0(i,ny+1,1)=u0(i,ny,1)
    u0(i,ny+1,2)=0.0
    u0(i,ny+1,3)=0.0
    u0(i,ny+1,4)=u0(i,ny,4) 
  ux(i,ny+1,k)=0.0 
  uy(i,ny+1,k)=0.0 
  uxx(i,ny+1,k)=0.0 
  uxy(i,ny+1,k)=0.0    
  uyy(i,ny+1,k)=0.0 
END DO
END DO 

!***
ELSEIF(UPPER_BOUNDARY.EQ."REFLECTIVE")THEN 
!*** 

DO i=1,nx 
 DO k=1,4
  u0(i,ny+1,1)=u0(i,ny,1)
    u0(i,ny+1,2)=u0(i,ny,2)              ! Reflective condition
    u0(i,ny+1,3)=0.0                     ! No penetration
    u0(i,ny+1,4)=u0(i,ny,4) 
  ux(i,ny+1,k)=0.0 
  uy(i,ny+1,k)=0.0 
  uxx(i,ny+1,k)=0.0 
  uxy(i,ny+1,k)=0.0    
  uyy(i,ny+1,k)=0.0 
END DO
END DO 
 
!***
ELSEIF(UPPER_BOUNDARY.EQ."SYMMETRY")THEN 
!***    

 DO i=1,nx 
 DO k=1,4
  u0(i,ny+1,1)=u0(i,ny,1)
    u0(i,ny+1,2)=u0(i,ny,2)
    u0(i,ny+1,3)= -u0(i,ny,3)
    u0(i,ny+1,4)=u0(i,ny,4) 
  ux(i,ny+1,k)=0.0 
  uy(i,ny+1,k)=0.0 
  uxx(i,ny+1,k)=0.0 
  uxy(i,ny+1,k)=0.0    
  uyy(i,ny+1,k)=0.0 
END DO
END DO  
 
END IF 

!*****************************************
!******* BOTTOM BOUNDARY
!*****************************************    

IF(BOTTOM_BOUNDARY.EQ."OUTFLOW")THEN
    
DO i=1,nx
DO k=1,4
    u0(i,0,k)=u0(i,1,k) 
  ux(i,0,k)=0.0
  uy(i,0,k)=0.0
  uxx(i,0,k)=0.0 
  uxy(i,0,k)=0.0    
  uyy(i,0,k)=0.0 
END DO
END DO 

!***
ELSEIF(BOTTOM_BOUNDARY.EQ."NO-SLIP")THEN 
!*** 

DO i=1,nx
DO k=1,4
    u0(i,0,1)=u0(i,1,1) 
    u0(i,0,2)=0.0
    u0(i,0,3)=0.0
    u0(i,0,4)=u0(i,1,4) 
  ux(i,0,k)=0.0
  uy(i,0,k)=0.0
  uxx(i,0,k)=0.0 
  uxy(i,0,k)=0.0    
  uyy(i,0,k)=0.0 
END DO
END DO

!***
ELSEIF(BOTTOM_BOUNDARY.EQ."REFLECTIVE")THEN 
!*** 

DO i=1,nx
DO k=1,4
    u0(i,0,1)=u0(i,1,1) 
    u0(i,0,2)=u0(i,1,2)     ! Reflective condition
    u0(i,0,3)=0.0           ! No penetration
    u0(i,0,4)=u0(i,1,4) 
  ux(i,0,k)=0.0
  uy(i,0,k)=0.0
  uxx(i,0,k)=0.0 
  uxy(i,0,k)=0.0    
  uyy(i,0,k)=0.0 
END DO
END DO


!***
ELSEIF(BOTTOM_BOUNDARY.EQ."SYMMETRY")THEN 
!***  

DO i=1,nx
DO k=1,4
    u0(i,0,1)=u0(i,1,1) 
    u0(i,0,2)=u0(i,1,2) 
    u0(i,0,3)=-u0(i,1,3) 
    u0(i,0,4)=u0(i,1,4) 
  ux(i,0,k)=0.0
  uy(i,0,k)=0.0
  uxx(i,0,k)=0.0 
  uxy(i,0,k)=0.0    
  uyy(i,0,k)=0.0 
END DO
END DO

END IF

!*************************************************************
!*************************************************************
 
END SUBROUTINE