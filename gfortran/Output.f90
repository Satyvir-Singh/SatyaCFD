!**********************************************************************
!***********       SHOCK - VORTEX INTERACTION             ************!
!********   2D- DISCONTINUOUS GALERKIN METHOD              ***********!
!*****                   SATYVIR SINGH                         *******!
!**                        2017                                    ***!
!**********************************************************************


SUBROUTINE OUTPUT(FILENAME)   
USE CONSTANTS
USE INPUT_INFO
USE GAS_INFO
USE TIME_INFO
USE GRID_INFO
USE FLOWFIELD_INFO
USE DIMENSIONLESS_INFO
USE COEFFICIENTS
USE MACROSCOPIC_INFO

IMPLICIT NONE
INTEGER :: I,J
DOUBLE PRECISION :: MIU
DOUBLE PRECISION :: ZERO

DOUBLE PRECISION :: VIS_1_TERM, VIS_2_TERM, VIS_3_TERM, VIS_4_TERM, VIS_5_TERM, VIS_6_TERM 
DOUBLE PRECISION :: VIS_7_TERM, VIS_8_TERM, VIS_9_TERM, VIS_10_TERM, VIS_11_TERM, VIS_12_TERM
DOUBLE PRECISION::  XPUX, XPUY
DOUBLE PRECISION :: DISSIPATION_RATE_VAL, ENSTROPHY_VAL, CIRCULATION_VAL, ABS_CIRCULATION_VAL
DOUBLE PRECISION :: AVERAGE_VORTICITY, ABS_DILATATIONAL_VORTICITY, ABS_BAROCLINIC_VORTICITY, ABSLOUTE_CIRCULATION
DOUBLE PRECISION :: ABS_VISCOUS_VORTICITY, ABS_DELTA_VAL,R_VAL, HEAT_VAL, STRESS_VAL, AVERAGE_AREA
DOUBLE PRECISION :: ABS_CONVECTION_VORTICITY

DOUBLE PRECISION :: DIL_VORTICITY, BARO_VORTICITY, VIS_VORTICITY, NET_VORTCITY

DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: Rho, U, V, E, P, T
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: DRHO_DX, DRHOU_DX, DRHOV_DX,DRHOE_DX 
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: DRHO_DY, DRHOU_DY, DRHOV_DY,DRHOE_DY 
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: DU_DX,DU_DY, DV_DX,DV_DY, DP_DX,  DP_DY, DT_DX, DT_DY, dMIU_DX,dMIU_DY
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: PI_XX,PI_XY, PI_YY, Q_X, Q_Y, Delta, R_PARAMETER, STRESS_TENSOR, HEAT_FLUX
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: DIVERGENCE_U
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: KINETIC_ENERGY
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: Velocity

DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: ENSTROPHY, DISSIPATION_RATE, CONVECTION_VORTICITY
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: VORTICITY, DILATATIONAL_VORTICITY, BAROCLINIC_VORTICITY, TOTAL_VISCOUS_VORTCITIY 

DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: dVORTCITY_DX, dVORTCITY_DY, ddVORTCITY_DXX, ddVORTCITY_DYY
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: ddU_dXX, ddV_dXX, ddV_dYY, ddU_dxy, ddV_dxy, ddU_dYY
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: dDIVERGENCE_U_dx, dDIVERGENCE_U_dy
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: ddMIU_dxx, ddMIU_dyy, ddMIU_dxy

DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:):: S_xx, S_yy,  S_xy
DOUBLE PRECISION :: VISCOUS_VORTICITY_PRODUCTION
DOUBLE PRECISION :: BAROCLINIC_VORTICITY_PRODUCTION
DOUBLE PRECISION :: DILATATIONAL_VORTICITY_PRODUCTION
DOUBLE PRECISION :: CONVECTION_VORTICITY_PRODUCTION
DOUBLE PRECISION :: KINETIC_ENERGY_VAL

DOUBLE PRECISION :: pos_circulation
DOUBLE PRECISION :: neg_circulation
DOUBLE PRECISION :: total_circulation

DOUBLE PRECISION :: vorticity_max, vorticity_max_integrated
DOUBLE PRECISION :: vorticity_min, vorticity_min_integrated


!=============================================================
CHARACTER(LEN=20) FILENAME
CHARACTER(LEN=20) ZONENAME
CHARACTER(LEN=20) FILETYPE
CHARACTER(LEN=40) DIRECTORY
CHARACTER(LEN=40) DIRECTORY_3
CHARACTER(LEN=40) SUFFIX1

DIRECTORY ='RESULT/'
DIRECTORY_3 = 'Result_Profiles/'
FILETYPE=".plt"
!=============================================================
!=============================================================
ALLOCATE ( Rho(0:NX+1,0:NY+1))
ALLOCATE ( U(0:NX+1,0:NY+1))
ALLOCATE ( V(0:NX+1,0:NY+1))
ALLOCATE ( E(0:NX+1,0:NY+1))
ALLOCATE ( P(0:NX+1,0:NY+1))
ALLOCATE ( T(0:NX+1,0:NY+1))
Rho(:,:) =0.0
U(:,:) =0.0
V(:,:) =0.0
E(:,:) =0.0
P(:,:) =0.0
T(:,:) =0.0

ALLOCATE ( DRHO_DX(0:NX+1,0:NY+1))
ALLOCATE ( DRHOU_DX(0:NX+1,0:NY+1))
ALLOCATE ( DRHOV_DX(0:NX+1,0:NY+1))
ALLOCATE ( DRHOE_DX(0:NX+1,0:NY+1))
DRHO_DX(:,:) =0.0
DRHOU_DX(:,:) =0.0
DRHOV_DX(:,:) =0.0
DRHOE_DX(:,:) =0.0

ALLOCATE ( DRHO_DY(0:NX+1,0:NY+1))
ALLOCATE ( DRHOU_DY(0:NX+1,0:NY+1))
ALLOCATE ( DRHOV_DY(0:NX+1,0:NY+1))
ALLOCATE ( DRHOE_DY(0:NX+1,0:NY+1))
DRHO_DY(:,:) =0.0
DRHOU_DY(:,:) =0.0
DRHOV_DY(:,:) =0.0
DRHOE_DY(:,:) =0.0

ALLOCATE ( DU_DX(0:NX+1,0:NY+1))
ALLOCATE ( DU_DY(0:NX+1,0:NY+1))
ALLOCATE ( DV_DX(0:NX+1,0:NY+1))
ALLOCATE ( DV_DY(0:NX+1,0:NY+1))
ALLOCATE ( DP_DX(0:NX+1,0:NY+1))
ALLOCATE ( DP_DY(0:NX+1,0:NY+1))
ALLOCATE ( DT_DX(0:NX+1,0:NY+1))
ALLOCATE ( DT_DY(0:NX+1,0:NY+1))
ALLOCATE ( dMIU_DX(0:NX+1,0:NY+1))
ALLOCATE ( dMIU_DY(0:NX+1,0:NY+1))
DU_DX(:,:) =0.0
DU_DY(:,:) =0.0
DV_DX(:,:) =0.0
DV_DY(:,:) =0.0
DP_DX(:,:) =0.0
DP_DY(:,:) =0.0
DT_DX(:,:) =0.0
DT_DY(:,:) =0.0
dMIU_DX(:,:) =0.0
dMIU_DY(:,:) =0.0

ALLOCATE ( PI_XX(0:NX+1,0:NY+1))
ALLOCATE ( PI_XY(0:NX+1,0:NY+1))
ALLOCATE ( PI_YY(0:NX+1,0:NY+1))
ALLOCATE ( Q_X(0:NX+1,0:NY+1))
ALLOCATE ( Q_Y(0:NX+1,0:NY+1))
ALLOCATE ( Delta(0:NX+1,0:NY+1))
ALLOCATE ( R_PARAMETER(0:NX+1,0:NY+1))
ALLOCATE ( STRESS_TENSOR(0:NX+1,0:NY+1))
ALLOCATE ( HEAT_FLUX(0:NX+1,0:NY+1))
ALLOCATE ( DIVERGENCE_U(0:NX+1,0:NY+1))
PI_XX(:,:) =0.0
PI_XY(:,:) =0.0
PI_YY(:,:) =0.0
Q_X(:,:) =0.0
Q_Y(:,:) =0.0
Delta(:,:) =0.0
R_PARAMETER(:,:) =0.0
STRESS_TENSOR(:,:) =0.0
HEAT_FLUX(:,:) =0.0
DIVERGENCE_U(:,:) =0.0

ALLOCATE ( ENSTROPHY(0:NX+1,0:NY+1))
ALLOCATE ( DISSIPATION_RATE(0:NX+1,0:NY+1))
ALLOCATE ( VORTICITY(0:NX+1,0:NY+1))
ALLOCATE ( BAROCLINIC_VORTICITY(0:NX+1,0:NY+1))
ALLOCATE ( DILATATIONAL_VORTICITY(0:NX+1,0:NY+1))
ALLOCATE ( CONVECTION_VORTICITY(0:NX+1,0:NY+1))
ALLOCATE ( TOTAL_VISCOUS_VORTCITIY(0:NX+1,0:NY+1))
ALLOCATE ( dVORTCITY_DX(0:NX+1,0:NY+1))
ALLOCATE ( dVORTCITY_DY(0:NX+1,0:NY+1))
ALLOCATE ( ddVORTCITY_DXX(0:NX+1,0:NY+1))
ALLOCATE ( ddVORTCITY_DYY(0:NX+1,0:NY+1))
ENSTROPHY(:,:) =0.0
VORTICITY(:,:) =0.0
DISSIPATION_RATE(:,:) =0.0
BAROCLINIC_VORTICITY(:,:) =0.0
CONVECTION_VORTICITY(:,:) =0.0
DILATATIONAL_VORTICITY(:,:) =0.0
TOTAL_VISCOUS_VORTCITIY(:,:) =0.0
dVORTCITY_DX(:,:) =0.0
dVORTCITY_DY(:,:) =0.0
ddVORTCITY_DXX(:,:) =0.0
ddVORTCITY_DYY(:,:) =0.0

ALLOCATE ( ddU_dXX(0:NX+1,0:NY+1))
ALLOCATE ( ddV_dXX(0:NX+1,0:NY+1))
ALLOCATE ( ddV_dYY(0:NX+1,0:NY+1))
ALLOCATE ( ddU_dxy(0:NX+1,0:NY+1))
ALLOCATE ( ddV_dxy(0:NX+1,0:NY+1))
ALLOCATE ( ddU_dYY(0:NX+1,0:NY+1))
ALLOCATE ( dDIVERGENCE_U_dx(0:NX+1,0:NY+1))
ALLOCATE ( dDIVERGENCE_U_dy(0:NX+1,0:NY+1))
ALLOCATE ( ddMIU_dxx(0:NX+1,0:NY+1))
ALLOCATE ( ddMIU_dyy(0:NX+1,0:NY+1))
ALLOCATE ( ddMIU_dxy(0:NX+1,0:NY+1))
ddU_dXX(:,:) =0.0
ddV_dXX(:,:) =0.0
ddV_dYY(:,:) =0.0
ddU_dxy(:,:) =0.0
ddV_dxy(:,:) =0.0
ddU_dYY(:,:) =0.0
dDIVERGENCE_U_dx(:,:) =0.0
dDIVERGENCE_U_dy(:,:) =0.0
ddMIU_dxx(:,:) =0.0
ddMIU_dyy(:,:) =0.0
ddMIU_dxy(:,:) =0.0

ALLOCATE ( S_xx(0:nx+1,0:ny+1))
ALLOCATE ( S_yy(0:nx+1,0:ny+1))
ALLOCATE ( S_xy(0:nx+1,0:ny+1))
S_xx(:,:)=0.0
S_yy(:,:)=0.0
S_xy(:,:)=0.0

ALLOCATE ( KINETIC_ENERGY(0:nx+1,0:ny+1))
ALLOCATE ( Velocity(0:nx+1,0:ny+1))
KINETIC_ENERGY(:,:) =0.0
Velocity(:,:) =0.0


!*********************************************************************************
!*********************************************************************************
WRITE(SUFFIX1,"(F16.9)")CURRENT_TIME*REF_TIME
OPEN(3,FILE=TRIM(DIRECTORY)//TRIM(FILENAME)//" -"//TRIM(GOV_EQ_SWITCH)//" -"//TRIM(GAS_TYPE) // &
" - t= "//TRIM(SUFFIX1)//" "//TRIM(FILETYPE) )
WRITE(3,*)'TITLE = "TESTCASE: BUBBLE-SHOCK INTERACTION" ' 

WRITE(3,'(A)')'VARIABLES = "X","Y","Rho","U","V","P"'

WRITE(3,*)'ZONE T="Floor", I=',ny,' J=',nx,' F=POINT ' 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DO I=0,nx+1 
DO J=0,ny+1 

ZERO = 0.0
    
!=== Primtiive Variables    
Rho(I,J) = u0(I,J,1)    
U(I,J)   = u0(I,J,2)/u0(I,J,1)  
V(I,J)   = u0(I,J,3)/u0(I,J,1) 
E(I,J)   = u0(I,J,4)/u0(I,J,1) 
P(I,J)   = gama*(gama-1.0)*Rho(i,j)*(E(I,J) -0.5*(U(I,J)**2 + V(i,j)**2))   
T(I,J)   = P(I,J)/Rho(I,J)  

Velocity(i,j) = SQRT(U(i,j)*U(i,j) + V(i,j)*V(i,j))
KINETIC_ENERGY(i,j) = 0.5*Rho(I,J)*Velocity(i,j)*Velocity(i,j)

!==== Derivative Values    

DRHO_DX(I,J)   = xpux(dx,ux(i,j,1),uxx(i,j,1),uxy(i,j,1),ZERO,ZERO)      ! dRho_dx
DRHOU_DX(I,J)  = xpux(dx,ux(i,j,2),uxx(i,j,2),uxy(i,j,2),ZERO,ZERO)      ! dRhoU_dx

DRHOV_DX(I,J)  = xpux(dx,ux(i,j,3),uxx(i,j,3),uxy(i,j,3),ZERO,ZERO)      ! dRhoV_dx
DRHOE_DX(I,J)  = xpux(dx,ux(i,j,4),uxx(i,j,4),uxy(i,j,4),ZERO,ZERO)      ! dRhoE_dx

DRHO_DY(I,J)   = xpuy(dy,uy(i,j,1),uxy(i,j,1),uyy(i,j,1),ZERO,ZERO)      ! dRho_dy
DRHOU_DY(I,J)  = xpuy(dy,uy(i,j,2),uxy(i,j,2),uyy(i,j,2),ZERO,ZERO)      ! dRhoU_dy

DRHOV_DY(I,J)  = xpuy(dy,uy(i,j,3),uxy(i,j,3),uyy(i,j,3),ZERO,ZERO)      ! dRhoV_dy
DRHOE_DY(I,J)  = xpuy(dy,uy(i,j,4),uxy(i,j,4),uyy(i,j,4),ZERO,ZERO)      ! dRhoE_dy    
   
DU_DX(I,J)    = (1.0/Rho(I,J))*(DRHOU_DX(I,J) - U(I,J)*DRHO_DX(I,J))            ! du_dx
DU_DY(I,J)    = (1.0/Rho(I,J))*(DRHOU_DY(I,J) - U(I,J)*DRHO_DY(I,J))            ! du_dy 

DV_DX(I,J)    = (1.0/Rho(I,J))*(DRHOV_DX(I,J) - V(I,J)*DRHO_DX(I,J))            ! dv_dx 
DV_DY(I,J)    = (1.0/Rho(I,J))*(DRHOV_DY(I,J) - V(I,J)*DRHO_DY(I,J))            ! dv_dy 
    
DP_DX(I,J)    = gama*(gama-1.0)*(DRHOE_DX(I,J) - (U(I,J)*DRHOU_DX(I,J) + V(I,J)*DRHOV_DX(I,J)) + &
 0.5*(U(I,J)**2 + V(I,J)**2)*DRHO_DX(I,J))
DP_DY(I,J)    = gama*(gama-1.0)*(DRHOE_DY(I,J) - (U(I,J)*DRHOU_DY(I,J) + V(I,J)*DRHOV_DY(I,J)) + &
 0.5*(U(I,J)**2 + V(I,J)**2)*DRHO_DY(I,J))

DT_DX(I,J)    = (1.0/Rho(I,J))*gama*(gama-1.0)*(DRHOE_DX(I,J) - (U(I,J)*DRHOU_DX(I,J) + V(I,J)*DRHOV_DX(I,J)) + &
                0.5*(U(I,J)**2 + V(I,J)**2)*DRHO_DX(I,J) - E(I,J) - 0.5*(U(I,J)**2 + V(I,J)**2)*DRHO_DX(I,J))
DT_DY(I,J)    = (1.0/Rho(I,J))*gama*(gama-1.0)*(DRHOE_DY(I,J) - (U(I,J)*DRHOU_DY(I,J) + V(I,J)*DRHOV_DY(I,J)) + &
                0.5*(U(I,J)**2 + V(I,J)**2)*DRHO_DY(I,J) - E(I,J) - 0.5*(U(I,J)**2 + V(I,J)**2)*DRHO_DY(I,J))
 
dMIU_DX(I,J)  = VISCOSITY_INDEX*T(I,J)**(VISCOSITY_INDEX-1.0)*DT_DX(I,J)
dMIU_DY(I,J)  = VISCOSITY_INDEX*T(I,J)**(VISCOSITY_INDEX-1.0)*DT_DY(I,J)

!=== Non-conservative Variables    
PI_xx(I,J) = pixx(I,J,0,0)*REF_SHEAR_STRESS 
PI_xy(I,J) = pixy(I,J,0,0)*REF_SHEAR_STRESS 
PI_yy(I,J) = piyy(I,J,0,0)*REF_SHEAR_STRESS 
Q_x(I,J)   = qx(I,J,0,0)*REF_HEAT_FLUX 
Q_y(I,J)   = qy(I,J,0,0)*REF_HEAT_FLUX
Delta(I,J) = bx(I,J,0,0)*REF_SHEAR_STRESS 

!=== FULL TENSORS Non-conservative Variables
STRESS_TENSOR(I,J) = SQRT(PI_xx(i,j)*PI_xx(i,j) + 2.0*PI_xy(i,j)*PI_xy(i,j) + PI_yy(i,j)*PI_yy(i,j))

HEAT_FLUX(I,J) = SQRT(Q_x(i,j)*Q_x(i,j) + Q_y(i,j)*Q_y(i,j))

R_PARAMETER(i,j) = SQRT(PI_XX(I,J)*PI_XX(I,J) + 2.0*PI_XY(I,J)*PI_XY(I,J) + PI_YY(I,J)*PI_XX(I,J) + &
                   2.0*GAMA_BULK*Bulk_Ratio*Delta(i,j)*Delta(i,j) + Q_X(i,j)*Q_X(i,j) + Q_Y(i,j)*Q_Y(i,j))

!================================================
!***** Strain rate tensor
S_xx(I,J) = DU_DX(I,J)
S_yy(I,J) = DV_DY(I,J)
S_xy(I,J) = 0.5*(DU_DY(I,J) + DV_DX(I,J))   ! =  S_yx(I,J)
!================================================
! DISSIPATION RATE OR  KINETIC ENERGY
!================================================
DISSIPATION_RATE(i,j) = -((PI_xx(i,j) + Bulk_Ratio*Delta(i,j))*S_xx(I,J) + 2.0*PI_xy(i,j)*S_xy(I,J) + &
                         (PI_yy(i,j) + Bulk_Ratio*Delta(i,j))*S_yy(I,J))

!===========================================================================================
!     EXTRA TERMS FOR CALCUALTING VISCOUS VORTCITY EUQATUATIONS
!===========================================================================================

DIVERGENCE_U(I,J) = DU_DX(I,J) + DV_DY(I,J)

VORTICITY(I,J)   = DV_DX(I,J) - DU_DY(I,J)

ENSTROPHY(I,J)   = 0.5*VORTICITY(I,J)*VORTICITY(I,J)

DILATATIONAL_VORTICITY(I,J) = VORTICITY(I,J)*DIVERGENCE_U(I,J)

BAROCLINIC_VORTICITY(I,J) = (1.0/(Rho(I,J)**2))*(DRHO_DX(I,J)*DP_DY(I,J) - DRHO_DY(I,J)*DP_DX(I,J))

END DO
END DO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!*********************************************************************************
!*********************************************************************************
!                           single/double derivatvies
!*********************************************************************************
!*********************************************************************************
DO I=1,NX
DO J=1,NY

dVORTCITY_DX(I,J) = (VORTICITY(I+1,J) - VORTICITY(I-1,J))/(2.0*DX)      ! domega/dx
dVORTCITY_DY(I,J) = (VORTICITY(I,J+1) - VORTICITY(I,J-1))/(2.0*DY)      ! domega/dy

ddVORTCITY_DXX(I,J) = (dVORTCITY_DX(I+1,J) - dVORTCITY_DX(I-1,J))/(2.0*DX)  ! d^2omega/dx^2
ddVORTCITY_DYY(I,J) = (dVORTCITY_DY(I,J+1) - dVORTCITY_DY(I,J-1))/(2.0*DY)  ! d^2omega/dy^2   
    
ddU_dXX(I,J) = (DU_DX(I+1,J) - DU_DX(I-1,J))/(2.0*DX)   ! d^2u/dx^2
ddU_dYY(I,J) = (DU_DY(I,J+1) - DU_DY(I,J+1))/(2.0*DY)   ! d^2v/dy^2
ddV_dXX(I,J) = (DV_DX(I+1,J) - DV_DX(I-1,J))/(2.0*DX)   ! d^2u/dx^2
ddV_dYY(I,J) = (DV_DY(I,J+1) - DV_DY(I,J+1))/(2.0*DY)   ! d^2v/dy^2
ddU_dxy(I,J) = (DU_DX(I,J+1) - DU_DX(I,J+1))/(2.0*DY)   ! d^2u/dxdy
ddV_dxy(I,J) = (DV_DX(I,J+1) - DV_DX(I,J+1))/(2.0*DY)   ! d^2u/dxdy

dDIVERGENCE_U_dx(I,J) = (DIVERGENCE_U(I+1,J) - DIVERGENCE_U(I-1,J))/(2.0*DX)   ! d/dx (Nebla.u)
dDIVERGENCE_U_dy(I,J) = (DIVERGENCE_U(I,J+1) - DIVERGENCE_U(I,J-1))/(2.0*DY)   ! d/dy (Nebla.u)

ddMIU_dxx(I,J)   = (dMIU_DX(I+1,J) - dMIU_DX(I-1,J))/(2.0*DX)    ! d^2 miu/dx^2
ddMIU_dyy(I,J)   = (dMIU_DY(I,J+1) - dMIU_DY(I,J-1))/(2.0*DY)    ! d^2 miu/dy^2 
ddMIU_dxy(I,J)  = (dMIU_DX(I,J+1) - dMIU_DX(I,J-1))/(2.0*DY)    ! d^2 miu/dxdy

END DO
END DO

!*********************************************************************************
!              TOTAL VISCOUS VORTCITY
!*********************************************************************************         
DO i=1,nx
DO j=1,ny      

MIU  = T(I,J)**VISCOSITY_INDEX    
       
VIS_1_TERM = (MIU/Rho(I,J))*(ddVORTCITY_DXX(I,J) + ddVORTCITY_DYY(I,J) )
VIS_2_TERM = -(MIU/(Rho(I,J)**2))*(DRHO_DX(I,J)*ddV_dYY(I,J) -  DRHO_DY(I,J)*ddU_dXX(I,J))   
VIS_3_TERM = -(1.0/3.0 + Bulk_Ratio)*(MIU/(Rho(I,J)**2))*(DRHO_DX(I,J)*dDIVERGENCE_U_dy(I,J) - &
DRHO_DY(I,J)*dDIVERGENCE_U_dx(I,J))  
VIS_4_TERM =  (1.0/Rho(I,J))*(dMIU_DX(I,J)*dVORTCITY_DX(I,J) + dMIU_DY(I,J)*dVORTCITY_DY(I,J))
VIS_5_TERM =  0.0
VIS_6_TERM = -(VORTICITY(I,J)/Rho(I,J))*(ddMIU_dxx(I,J) + ddMIU_dyy(I,J) )
VIS_7_TERM = (1.0/Rho(I,J))*(dMIU_DX(I,J)*ddU_dYY(I,J)  - dMIU_DY(I,J)*ddU_dXX(I,J) )
VIS_8_TERM = (2.0/Rho(I,J))*(dMIU_DX(I,J)*ddV_dXX(I,J) + DV_DX(I,J)*ddMIU_dxx(I,J) + DV_DY(I,J)*ddMIU_dxy(I,J) + &
dMIU_DY(I,J)*ddV_dxy(I,J) - &
                 dMIU_DX(I,J)*ddU_dXY(I,J) - DU_DX(I,J)*ddMIU_dxy(I,J) - dMIU_DY(I,J)*ddU_dYY(I,J) - DU_DY(I,J)*ddMIU_dyy(I,J) )
VIS_9_TERM = (1.0/(Rho(I,J)**2))*(DRHO_DX(I,J)*DU_DX(I,J) + DRHO_DY(I,J)*DV_DX(I,J))*VORTICITY(I,J)
VIS_10_TERM = -(2.0/(Rho(I,J)**2))*(DRHO_DX(I,J)*dMIU_DY(I,J)*DV_DY(I,J) - DRHO_DY(I,J)*dMIU_DX(I,J)*DU_DX(I,J) )
VIS_11_TERM = (2.0/3.0 - Bulk_Ratio)*(DIVERGENCE_U(I,J)/(Rho(I,J)**2))*(DRHO_DX(I,J)*dMIU_DY(I,J) - DRHO_DY(I,J)*dMIU_DX(I,J) )
VIS_12_TERM =  (1.0/Rho(I,J))*(dMIU_DX(I,J)*dDIVERGENCE_U_dy(I,J) - dMIU_DY(I,J)*dDIVERGENCE_U_dx(I,J))

TOTAL_VISCOUS_VORTCITIY(I,J) =  VIS_1_TERM + VIS_2_TERM + VIS_3_TERM + VIS_4_TERM + VIS_5_TERM +  VIS_6_TERM + VIS_7_TERM + &
                                VIS_8_TERM + VIS_9_TERM + VIS_10_TERM + VIS_11_TERM + VIS_12_TERM


CONVECTION_VORTICITY(I,J) = -U(I,J)*dVORTCITY_DX(I,J) - V(I,J)*dVORTCITY_DY(I,J)

END DO
END DO

!*********************************************************************************
!*********************************************************************************
!                        FLOW FIELDS
!*********************************************************************************
!*********************************************************************************  
!==========================
DO i=1,nx
DO j=1,ny    
    
WRITE(3,444) X(I),Y(J),RHO(I,J),U(I,J),V(I,J),P(I,J)
END DO
END DO 

CLOSE(3)   


!==========================
!*********************************************************************************
!*********************************************************************************  

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  
!     SBI-PROFILES
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!*********************************************************************************
!  1. Dissipation Rate 
!*********************************************************************************
OPEN(11,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"Dissipation_Rate"//TRIM(FILETYPE))

DISSIPATION_RATE_VAL =0.0

DO I=1,NX
DO J=1,NY
    DISSIPATION_RATE_VAL = DISSIPATION_RATE_VAL + DISSIPATION_RATE(i,j)*DX*DY
END DO
END DO

WRITE(11,444) CURRENT_TIME*REF_TIME, DISSIPATION_RATE_VAL
            
!*********************************************************************************
!  2. Enstrophy
!*********************************************************************************
OPEN(12,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"Enstrophy"//TRIM(FILETYPE))

ENSTROPHY_VAL = 0.0

DO I=1,NX
DO J=1,NY
    ENSTROPHY_VAL = ENSTROPHY_VAL + ENSTROPHY(I,J)*DX*DY
END DO
END DO

WRITE(12,444) CURRENT_TIME*REF_TIME, ENSTROPHY_VAL


!*********************************************************************************
!  5. Average Vorticity
!*********************************************************************************
OPEN(13,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"Average_Vorticity"//TRIM(FILETYPE))

AVERAGE_AREA =0.0
AVERAGE_VORTICITY =0.0

DO I=1,NX
DO J=1,NY
    ABSLOUTE_CIRCULATION = ABSLOUTE_CIRCULATION + ABS(Vorticity(i,j))*DX*DY
    AVERAGE_AREA = AVERAGE_AREA + DX*DY
    AVERAGE_VORTICITY = ABSLOUTE_CIRCULATION/AVERAGE_AREA 
END DO
END DO

WRITE(13,444) CURRENT_TIME*REF_TIME, AVERAGE_VORTICITY

!*********************************************************************************
!  6. CONVECTION Vorticity PRODUCTION
!*********************************************************************************
OPEN(1400,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"Convection_Vorticity_Production"//TRIM(FILETYPE))

ABS_CONVECTION_VORTICITY = 0.0
AVERAGE_AREA =0.0

DO I=1,NX
DO J=1,NY
    ABS_CONVECTION_VORTICITY = ABS_CONVECTION_VORTICITY + ABS(CONVECTION_VORTICITY(I,J))*DX*DY
    AVERAGE_AREA = AVERAGE_AREA + DX*DY
    CONVECTION_VORTICITY_PRODUCTION = ABS_CONVECTION_VORTICITY/AVERAGE_AREA
END DO
END DO

WRITE(1400,444) CURRENT_TIME*REF_TIME, CONVECTION_VORTICITY_PRODUCTION

!*********************************************************************************



!*********************************************************************************
!  6. Dillatational Vorticity PRODUCTION
!*********************************************************************************
OPEN(14,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"Dillatational_Vorticity_Production"//TRIM(FILETYPE))

ABS_DILATATIONAL_VORTICITY = 0.0
AVERAGE_AREA =0.0

DO I=1,NX
DO J=1,NY
    ABS_DILATATIONAL_VORTICITY = ABS_DILATATIONAL_VORTICITY + ABS(DILATATIONAL_VORTICITY(I,J))*DX*DY
    AVERAGE_AREA = AVERAGE_AREA + DX*DY
    DILATATIONAL_VORTICITY_PRODUCTION = ABS_DILATATIONAL_VORTICITY/AVERAGE_AREA
END DO
END DO

WRITE(14,444) CURRENT_TIME*REF_TIME, DILATATIONAL_VORTICITY_PRODUCTION

!*********************************************************************************
!  7. Baroclinic Vorticity PRODUCTION
!*********************************************************************************
OPEN(15,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"Baroclinic_Vorticity_Production"//TRIM(FILETYPE))

ABS_BAROCLINIC_VORTICITY =0.0
AVERAGE_AREA =0.0

DO I=1,NX
DO J=1,NY
    ABS_BAROCLINIC_VORTICITY = ABS_BAROCLINIC_VORTICITY + ABS(BAROCLINIC_VORTICITY(I,J))*DX*DY
    AVERAGE_AREA = AVERAGE_AREA + DX*DY
    BAROCLINIC_VORTICITY_PRODUCTION = ABS_BAROCLINIC_VORTICITY/AVERAGE_AREA
END DO
END DO

WRITE(15,444) CURRENT_TIME*REF_TIME, BAROCLINIC_VORTICITY_PRODUCTION


!*********************************************************************************
!  8. Viscous Vorticity
!*********************************************************************************
OPEN(16,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"Viscous_Vorticity_Production"//TRIM(FILETYPE))

ABS_VISCOUS_VORTICITY =0.0
AVERAGE_AREA =0.0

DO I=1,NX
DO J=1,NY
    ABS_VISCOUS_VORTICITY = ABS_VISCOUS_VORTICITY + ABS(TOTAL_VISCOUS_VORTCITIY(I,J))*DX*DY
    AVERAGE_AREA = AVERAGE_AREA + DX*DY
    VISCOUS_VORTICITY_PRODUCTION = ABS_VISCOUS_VORTICITY/AVERAGE_AREA
END DO
END DO

WRITE(16,444) CURRENT_TIME*REF_TIME, VISCOUS_VORTICITY_PRODUCTION

!*********************************************************************************
!  9. Delta profile
!*********************************************************************************
OPEN(19,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"Delta_Value"//TRIM(FILETYPE))

ABS_DELTA_VAL =0.0

DO I=1,NX
DO J=1,NY
    ABS_DELTA_VAL = ABS_DELTA_VAL + abs(Delta(i,j))*DX*DY
END DO
END DO

WRITE(19,444) CURRENT_TIME*REF_TIME, ABS_DELTA_VAL

!*********************************************************************************
!*********************************************************************************

!*********************************************************************************
!  20. R-PARAMETER profile
!*********************************************************************************
OPEN(20,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"R_Parameter"//TRIM(FILETYPE))

R_VAL =0.0

DO I=1,NX
DO J=1,NY
    R_VAL = R_VAL + R_PARAMETER(I,J)*DX*DY
END DO
END DO

WRITE(20,444) CURRENT_TIME*REF_TIME, R_VAL

!*********************************************************************************
!*********************************************************************************

!*********************************************************************************
!  21. stress tensor profile
!*********************************************************************************
OPEN(21,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"Stress_Tensor_Value"//TRIM(FILETYPE))

STRESS_VAL =0.0

DO I=1,NX
DO J=1,NY
    STRESS_VAL = STRESS_VAL + STRESS_TENSOR(I,J)*DX*DY
END DO
END DO

WRITE(21,444) CURRENT_TIME*REF_TIME, STRESS_VAL

!*********************************************************************************
!*********************************************************************************

!*********************************************************************************
!  22. Heat flux tensor profile
!*********************************************************************************
OPEN(22,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"Heat_Flux_Value"//TRIM(FILETYPE))

HEAT_VAL =0.0

DO I=1,NX
DO J=1,NY
    HEAT_VAL = HEAT_VAL + HEAT_FLUX(I,J)*DX*DY
END DO
END DO

WRITE(22,444) CURRENT_TIME*REF_TIME, HEAT_VAL


!*********************************************************************************
!  2. Kinetic_Energy_
!*********************************************************************************
OPEN(10001,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"Kinetic_Energy"//TRIM(FILETYPE))

KINETIC_ENERGY_VAL =0.0

DO I=1,NX
DO J=1,NY
    KINETIC_ENERGY_VAL = KINETIC_ENERGY_VAL + KINETIC_ENERGY(i,j)*DX*DY
END DO
END DO

WRITE(10001,444) CURRENT_TIME*REF_TIME, KINETIC_ENERGY_VAL



!*********************************************************************************
!  100. Circulation
!*********************************************************************************
OPEN(10002,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"Pos_Neg_tot_circulation"//TRIM(FILETYPE))

pos_circulation = 0.0
neg_circulation = 0.0
total_circulation = 0.0

DO I=1,NX
DO J=1,NY
    if(Vorticity(i,j) > 0.0)then
      pos_circulation = pos_circulation + Vorticity(i,j) * dx * dy
    else
     neg_circulation = neg_circulation + Vorticity(i,j) * dx * dy
    end if
    total_circulation = total_circulation + Vorticity(i,j) * dx * dy
END DO
END DO

WRITE(10002,444) CURRENT_TIME*REF_TIME, pos_circulation, neg_circulation, total_circulation



!*********************************************************************************
!  101. Vorticity max and Vortcity min
!*********************************************************************************
OPEN(10003,FILE=TRIM(DIRECTORY_3)//TRIM(GAS_TYPE)// "-"//"Max_Min_Vortcity"//TRIM(FILETYPE))

vorticity_max = -1.0e30
vorticity_min = 1.0e30

vorticity_max_integrated = 0.0
vorticity_min_integrated = 0.0

DO I=1,NX
DO J=1,NY
    if (Vorticity(i,j) > vorticity_max) then
        vorticity_max = Vorticity(i,j)
      end if
    if (Vorticity(i,j) < vorticity_min) then
        vorticity_min = Vorticity(i,j)
    end if
   ! Integrate max and min vorticity values over the grid
    vorticity_max_integrated = vorticity_max_integrated + vorticity_max * dx * dy
    vorticity_min_integrated = vorticity_min_integrated + vorticity_min * dx * dy   
END DO
END DO

WRITE(10003,444) CURRENT_TIME*REF_TIME, vorticity_max_integrated, vorticity_min_integrated




!*********************************************************************************
!*********************************************************************************

444 FORMAT(2x,37e25.10)  	
      
END SUBROUTINE