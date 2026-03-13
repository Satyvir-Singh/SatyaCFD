!**********************************************************************
!***********       SHOCK - VORTEX INTERACTION             ************!
!********   2D- DISCONTINUOUS GALERKIN METHOD              ***********!
!*****                   SATYVIR SINGH                         *******!
!**                        2017                                    ***!
!**********************************************************************  
 
 SUBROUTINE TIME_STEP
      
 USE GAS_INFO
 USE TIME_INFO
 USE GRID_INFO
 USE COEFFICIENTS
 USE DIMENSIONLESS_INFO

IMPLICIT NONE

      
  DOUBLE PRECISION :: rlamdamax
  DOUBLE PRECISION :: P,E,U,V,RHO
  DOUBLE PRECISION :: temp1,temp2,temp
  DOUBLE PRECISION :: MIN_DT
  INTEGER :: I,J

      rlamdamax=0. 
      do 101 i=1,nx 
        do 101 j=1,ny                 
          rho=u0(i,j,1) 
          u=u0(i,j,2)/u0(i,j,1)      
          v=u0(i,j,3)/u0(i,j,1)    
          e=u0(i,j,4) 
          p=(gama-1.)*gama*(e-0.5*rho*(u*u+v*v))
          temp1=dabs(u)+dsqrt(p/rho) 
          temp2=dabs(v)+dsqrt(p/rho) 
          temp=dsqrt(temp1*temp1+temp2*temp2) 
          if(temp.gt.rlamdamax) rlamdamax=temp 
  101 continue 
      dt=min(dx,dy)/rlamdamax*0.2 


!-----------------------------------------------------
 IF(CURRENT_TIME+dt .GT. FLOW_TIME_STAR )THEN
   dt = FLOW_TIME_STAR - CURRENT_TIME
 END IF

END SUBROUTINE  
