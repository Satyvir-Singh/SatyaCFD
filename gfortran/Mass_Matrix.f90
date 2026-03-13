      
!****************************************************************
!***************************************************************      
 
SUBROUTINE MASS_MATRIX
      
USE GRID_INFO
USE COEFFICIENTS
 
IMPLICIT NONE
DOUBLE PRECISION :: VOL, INV_VOL

VOL =DX*DY
INV_VOL = 1.0/VOL

RMASS(1) = 1.0*INV_VOL
RMASS(2) = 3.0*INV_VOL
RMASS(3) = 3.0*INV_VOL
RMASS(4) = (45.0/4.0)*INV_VOL
RMASS(5) = 9.0*INV_VOL
RMASS(6) = (45.0/4.0)*INV_VOL

END SUBROUTINE  
