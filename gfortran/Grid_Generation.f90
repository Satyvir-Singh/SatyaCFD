 
!**********************************************************************
!***********       SHOCK - VORTEX INTERACTION             ************!
!********   2D- DISCONTINUOUS GALERKIN METHOD              ***********!
!*****                   SATYVIR SINGH                         *******!
!**                        2017                                    ***!
!********************************************************************** 
 
SUBROUTINE GRID_GENERATION
         
USE GRID_INFO

IMPLICIT NONE
INTEGER :: I,J
                
DX = (xr - xl)/NX
DY = (yu - yb)/NY

DO I = 0, nx+1 
     X(I) = xl+I*DX
END DO 

DO J = 0, ny+1 
     Y(J) = yb+J*DY
END DO 

END SUBROUTINE