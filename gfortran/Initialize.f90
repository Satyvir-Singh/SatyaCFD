!**********************************************************************
!***********       SHOCK - VORTEX INTERACTION             ************!
!********   2D- DISCONTINUOUS GALERKIN METHOD              ***********!
!*****                   SATYVIR SINGH                         *******!
!**                        2017                                    ***!
!**********************************************************************

SUBROUTINE initialize
      
USE GAS_INFO
USE TIME_INFO
USE CONSTANTS
USE GRID_INFO
USE INPUT_INFO
USE FLOWFIELD_INFO
USE DIMENSIONLESS_INFO
USE COEFFICIENTS

IMPLICIT NONE
INTEGER:: I,J,K
DOUBLE PRECISION:: U,V,P,RHO
DOUBLE PRECISION:: origin_x,origin_y
DOUBLE PRECISION:: Inclined_angle, slope, x_distance
DOUBLE PRECISION:: origin_x1,origin_y1, origin_x2,origin_y2

!*********************** INPUT DATA
ITERATION = 0
CURRENT_TIME= 0.0

!===> SIMULATION REGIEM
MACH_REF= 1.0                            
REF_LENGTH = 40.0                        
  
REF_VELOCITY=MACH_REF*sqrt(gama*GASR*Ref_Temperature)       
REF_DENSITY = Ref_Pressure/ (GASR * Ref_Temperature)              

Reynolds = REF_DENSITY*REF_VELOCITY*REF_LENGTH/Viscosity 

Eckert=(gama-1.0)*MACH_REF*MACH_REF                  
Prandtl=Cp*Viscosity/Thermal_Conductivity                    
xndelta=Viscosity*REF_VELOCITY/Ref_Pressure           

REF_SHEAR_STRESS = Viscosity*REF_VELOCITY 
REF_HEAT_FLUX  = Thermal_Conductivity*Ref_Temperature 

REF_TIME       = REF_LENGTH / REF_VELOCITY
FLOW_TIME_STAR = FINAL_FLOW_TIME/REF_TIME

!===> REFERENCE VALUES
rhoref=1.25    
pref=1.0      

rho =0.0
p =0.0
u =0.0
v =0.0

!%%%%%%%%%%%%%%%%%%%%%%%
Inclined_angle = cotan(pi/12)
x_distance = 50*Inclined_angle
slope = (100-50)/(x_distance-25)

!slope = tan(5*pi/12)
slope = tan(pi/6)

!**************************************************************************
!**************************************************************************
!=============================================
! SVI PROBELM : MOVING SHOCK AND STATIONARY BUBBLE
!=============================================

origin_x1 =60.0
origin_x2 =25.0
origin_Y1 =55.0
origin_Y2 =55.0


DO i=0,nx+1 
DO j=0,ny+1   
    
IF(x(i) .LE. SHOCK_POSITION)THEN
RHO = rhoref*((gama+1)*SHOCK_MACH*SHOCK_MACH/(2+(gama-1)*SHOCK_MACH*SHOCK_MACH))    
P   = pref* (1.0 + 2*gama*(SHOCK_MACH*SHOCK_MACH-1.0)/(gama+1))    
U   = SHOCK_MACH*SQRT(gama*pref/rhoref)*(1.0 -rhoref/RHO)    
v=0.0     
!ELSEIF (  (50 + slope*(X(I) -25)) .GE. Y(J)  .AND.   (50 - slope*(X(I) -25)) .LE. Y(J) )THEN
ELSEIF (sqrt((x(i)-origin_x1)**2 + (Y(J)-origin_y1)**2)  .LE.  25.0)THEN 
rho=0.16
p=1.0
u=0.0
v=0.0 
ELSE
rho=1.25
p=1.0
u=0.0
v=0.0  
END IF 
    
!======> INITIALIZE      
u0(i,j,1)=rho 
u0(i,j,2)=rho*u 
u0(i,j,3)=rho*v 
u0(i,j,4)=p/(gama-1.)/gama+0.5*rho*(u*u+v*v) 

DO k=1,4 
ux(i,j,k)=0.0    
uy(i,j,k)=0.0 
uxx(i,j,k)=0.0 
uxy(i,j,k)=0.0 
uyy(i,j,k)=0.0 
END DO

END DO
END DO

END SUBROUTINE                       

