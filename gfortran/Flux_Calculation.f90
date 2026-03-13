          
subroutine calfv_edge
USE GAS_INFO
USE GRID_INFO
USE INPUT_INFO
USE DIMENSIONLESS_INFO
USE COEFFICIENTS
      
IMPLICIT NONE

INTEGER :: I,J,K,M

DOUBLE PRECISION,ALLOCATABLE,DIMENSION (:) :: GP
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: uplus, umins
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: fplus, fmins
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: fvplus, fvmins

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: pixxr, piyyr, pixyr, qxr, qyr,bxr
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: pixxl, piyyl, pixyl, qxl, qyl,bxl
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: pixxu, piyyu, pixyu, qxu, qyu,bxu
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: pixxb, piyyb, pixyb, qxb, qyb,bxb

DOUBLE PRECISION:: ur, ul, uu, ub
DOUBLE PRECISION:: f1,f2,f3,f4,g1,g2,g3,g4
DOUBLE PRECISION:: f1v,f2v,f3v,f4v,g1v,g2v,g3v,g4v
DOUBLE PRECISION:: xpuxr,xpuxl,xpuxu,xpuxb,xpuyu,xpuyb,xpuyr,xpuyl
    
DOUBLE PRECISION:: RHO0,P0,SONIC,alfa
DOUBLE PRECISION :: RHO_VAL,U_VAL,V_VAL,E_VAL,P_VAL,T_VAL       
DOUBLE PRECISION :: D_RHO_DX,D_RHO_U_DX,D_RHO_V_DX,D_RHO_E_DX
DOUBLE PRECISION :: D_RHO_DY,D_RHO_U_DY,D_RHO_V_DY,D_RHO_E_DY
DOUBLE PRECISION :: DU_DX,DU_DY,DV_DX,DV_DY,DE_DX,DE_DY
DOUBLE PRECISION :: DT_DX,DT_DY
DOUBLE PRECISION :: MUI,KAPA

ALLOCATE ( GP(-1:1) )
ALLOCATE ( uplus(-1:1,4) )
ALLOCATE ( umins(-1:1,4) )
ALLOCATE ( fplus(-1:1,4) )
ALLOCATE ( fmins(-1:1,4) )
ALLOCATE ( fvplus(-1:1,4) )
ALLOCATE ( fvmins(-1:1,4) )

ALLOCATE ( pixxr(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( piyyr(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( pixyr(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( qxr(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( qyr(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( bxr(0:nx+1,0:ny+1,-1:1) )

ALLOCATE ( pixxl(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( piyyl(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( pixyl(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( qxl(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( qyl(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( bxl(0:nx+1,0:ny+1,-1:1) )

ALLOCATE ( pixxu(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( piyyu(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( pixyu(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( qxu(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( qyu(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( bxu(0:nx+1,0:ny+1,-1:1) )

ALLOCATE ( pixxb(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( piyyb(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( pixyb(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( qxb(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( qyb(0:nx+1,0:ny+1,-1:1) )
ALLOCATE ( bxb(0:nx+1,0:ny+1,-1:1) )

!======> RESET DATA
!=============================================================
!======> GAUSS POINTS 
 gp(-1)=-sqrt(0.6) 
 gp(0)=0. 
 gp(1)=sqrt(0.6)
	    
 do 1500 i=0,nx+1 
 do 1500 j=0,ny+1
 do m=-1,1

 !******** SIDE 1 (RIGHT)
 uplus(m,1)=ur(u0(i,j,1),ux(i,j,1),uy(i,j,1),uxx(i,j,1),uxy(i,j,1),uyy(i,j,1),gp(m))      ! RHO
 uplus(m,2)=ur(u0(i,j,2),ux(i,j,2),uy(i,j,2),uxx(i,j,2),uxy(i,j,2),uyy(i,j,2),gp(m))      ! RHO_U
 uplus(m,3)=ur(u0(i,j,3),ux(i,j,3),uy(i,j,3),uxx(i,j,3),uxy(i,j,3),uyy(i,j,3),gp(m))      ! RHO_V
 uplus(m,4)=ur(u0(i,j,4),ux(i,j,4),uy(i,j,4),uxx(i,j,4),uxy(i,j,4),uyy(i,j,4),gp(m))      ! RHO_E
  
RHO_VAL = uplus(m,1)
U_VAL   = uplus(m,2)/uplus(m,1)
V_VAL   = uplus(m,3)/uplus(m,1)
E_VAL   = uplus(m,4)/uplus(m,1)
P_VAL   = gama*(gama-1.0)*RHO_VAL*(E_VAL -0.5*(U_VAL**2 + V_VAL**2))
T_VAL   = P_VAL/RHO_VAL 
  
 D_RHO_DX  =xpuxr(dx,ux(i,j,1),uxx(i,j,1),uxy(i,j,1),gp(m))         !  D_RHO_DX
 D_RHO_U_DX=xpuxr(dx,ux(i,j,2),uxx(i,j,2),uxy(i,j,2),gp(m))         !  D_RHO_U_DX
 D_RHO_V_DX=xpuxr(dx,ux(i,j,3),uxx(i,j,3),uxy(i,j,3),gp(m))         !  D_RHO_V_DX
 D_RHO_E_DX=xpuxr(dx,ux(i,j,4),uxx(i,j,4),uxy(i,j,4),gp(m))         !  D_RHO_E_DX

 D_RHO_DY  =xpuyr(dy,uy(i,j,1),uxy(i,j,1),uyy(i,j,1),gp(m))         !  D_RHO_DY
 D_RHO_U_DY=xpuyr(dy,uy(i,j,2),uxy(i,j,2),uyy(i,j,2),gp(m))         !  D_RHO_U_DY
 D_RHO_V_DY=xpuyr(dy,uy(i,j,3),uxy(i,j,3),uyy(i,j,3),gp(m))         !  D_RHO_V_DY
 D_RHO_E_DY=xpuyr(dy,uy(i,j,4),uxy(i,j,4),uyy(i,j,4),gp(m))         !  D_RHO_E_DY

DU_DX = (1.0/RHO_VAL)*(D_RHO_U_DX - U_VAL*D_RHO_DX)
DU_DY = (1.0/RHO_VAL)*(D_RHO_U_DY - U_VAL*D_RHO_DY)

DV_DX = (1.0/RHO_VAL)*(D_RHO_V_DX - V_VAL*D_RHO_DX)
DV_DY = (1.0/RHO_VAL)*(D_RHO_V_DY - V_VAL*D_RHO_DY)

DE_DX = (1.0/RHO_VAL)*(D_RHO_E_DX - E_VAL*D_RHO_DX)
DE_DY = (1.0/RHO_VAL)*(D_RHO_E_DY - E_VAL*D_RHO_DY)

DT_DX = (gama-1.0)*gama*(DE_DX - U_VAL*DU_DX - V_VAL*DV_DX)
DT_DY = (gama-1.0)*gama*(DE_DY - U_VAL*DU_DY - V_VAL*DV_DY)

MUI  = T_VAL**VISCOSITY_INDEX
KAPA = T_VAL**VISCOSITY_INDEX 
 
IF(GOV_EQ_SWITCH.EQ."EULER")THEN
pixxr(i,j,m) =0.0
piyyr(i,j,m) =0.0
pixyr(i,j,m) =0.0
qxr(i,j,m) =0.0
qyr(i,j,m) =0.0
bxr(i,j,m) =0.0
ELSEIF(GOV_EQ_SWITCH.EQ."NSF")THEN
pixxr(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DU_DX - (2.0/3.0)*DV_DY))
piyyr(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DV_DY - (2.0/3.0)*DU_DX))
pixyr(i,j,m)= (1.0/Reynolds)*(-MUI*(DU_DY + DV_DX))
qxr(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DX)
qyr(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DY)
bxr(i,j,m)  =  0.0
ELSEIF(GOV_EQ_SWITCH.EQ."NF")THEN
pixxr(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DU_DX - (2.0/3.0)*DV_DY))
piyyr(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DV_DY - (2.0/3.0)*DU_DX))
pixyr(i,j,m)= (1.0/Reynolds)*(-MUI*(DU_DY + DV_DX))
qxr(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DX)
qyr(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DY)
bxr(i,j,m)  = (1.0/Reynolds)*(-Bulk_Ratio*MUI*(DU_DX+ DV_DY))     
END IF
 
!*************************** 
!***************************
!*****LEFT SIDE

umins(m,1)=ul(u0(i,j,1),ux(i,j,1),uy(i,j,1),uxx(i,j,1),uxy(i,j,1),uyy(i,j,1),gp(m))    ! RHO
umins(m,2)=ul(u0(i,j,2),ux(i,j,2),uy(i,j,2),uxx(i,j,2),uxy(i,j,2),uyy(i,j,2),gp(m))    !RHO_U
umins(m,3)=ul(u0(i,j,3),ux(i,j,3),uy(i,j,3),uxx(i,j,3),uxy(i,j,3),uyy(i,j,3),gp(m))    !RHO_V
umins(m,4)=ul(u0(i,j,4),ux(i,j,4),uy(i,j,4),uxx(i,j,4),uxy(i,j,4),uyy(i,j,4),gp(m))    !RHO_E

RHO_VAL = umins(m,1)
U_VAL   = umins(m,2)/umins(m,1)
V_VAL   = umins(m,3)/umins(m,1)
E_VAL   = umins(m,4)/umins(m,1)
P_VAL   = gama*(gama-1.0)*RHO_VAL*(E_VAL -0.5*(U_VAL**2 + V_VAL**2))
T_VAL   = P_VAL/RHO_VAL 

D_RHO_DX   = xpuxl(dx,ux(i,j,1),uxx(i,j,1),uxy(i,j,1),gp(m))
D_RHO_U_DX = xpuxl(dx,ux(i,j,2),uxx(i,j,2),uxy(i,j,2),gp(m))
D_RHO_V_DX = xpuxl(dx,ux(i,j,3),uxx(i,j,3),uxy(i,j,3),gp(m))
D_RHO_E_DX = xpuxl(dx,ux(i,j,4),uxx(i,j,4),uxy(i,j,4),gp(m))

D_RHO_DY   = xpuyl(dy,uy(i,j,1),uxy(i,j,1),uyy(i,j,1),gp(m))
D_RHO_U_DY = xpuyl(dy,uy(i,j,2),uxy(i,j,2),uyy(i,j,2),gp(m))
D_RHO_V_DY = xpuyl(dy,uy(i,j,3),uxy(i,j,3),uyy(i,j,3),gp(m))
D_RHO_E_DY = xpuyl(dy,uy(i,j,4),uxy(i,j,4),uyy(i,j,4),gp(m))

DU_DX = (1.0/RHO_VAL)*(D_RHO_U_DX - U_VAL*D_RHO_DX)
DU_DY = (1.0/RHO_VAL)*(D_RHO_U_DY - U_VAL*D_RHO_DY)

DV_DX = (1.0/RHO_VAL)*(D_RHO_V_DX - V_VAL*D_RHO_DX)
DV_DY = (1.0/RHO_VAL)*(D_RHO_V_DY - V_VAL*D_RHO_DY)

DE_DX = (1.0/RHO_VAL)*(D_RHO_E_DX - E_VAL*D_RHO_DX)
DE_DY = (1.0/RHO_VAL)*(D_RHO_E_DY - E_VAL*D_RHO_DY)

DT_DX = (gama-1.0)*gama*(DE_DX - U_VAL*DU_DX - V_VAL*DV_DX)
DT_DY = (gama-1.0)*gama*(DE_DY - U_VAL*DU_DY - V_VAL*DV_DY)

MUI  = T_VAL**VISCOSITY_INDEX
KAPA = T_VAL**VISCOSITY_INDEX 

IF(GOV_EQ_SWITCH.EQ."EULER")THEN
pixxl(i,j,m) =0.0
piyyl(i,j,m) =0.0
pixyl(i,j,m) =0.0
qxl(i,j,m) =0.0
qyl(i,j,m) =0.0
bxl(i,j,m) =0.0
ELSEIF(GOV_EQ_SWITCH.EQ."NSF")THEN
pixxl(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DU_DX - (2.0/3.0)*DV_DY))
piyyl(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DV_DY - (2.0/3.0)*DU_DX))
pixyl(i,j,m)= (1.0/Reynolds)*(-MUI*(DU_DY + DV_DX))
qxl(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DX)
qyl(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DY)
bxl(i,j,m)  =  0.0
ELSEIF(GOV_EQ_SWITCH.EQ."NF")THEN
pixxl(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DU_DX - (2.0/3.0)*DV_DY))
piyyl(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DV_DY - (2.0/3.0)*DU_DX))
pixyl(i,j,m)= (1.0/Reynolds)*(-MUI*(DU_DY + DV_DX))
qxl(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DX)
qyl(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DY)
bxl(i,j,m)  = (1.0/Reynolds)*(-Bulk_Ratio*MUI*(DU_DX+ DV_DY))     
END IF

!***************************************************************************************
!***************************************************************************************      
!***** UPPER SIDE
      
 uplus(m,1)=uu(u0(i,j,1),ux(i,j,1),uy(i,j,1),uxx(i,j,1),uxy(i,j,1),uyy(i,j,1),gp(m))   ! RHO
 uplus(m,2)=uu(u0(i,j,2),ux(i,j,2),uy(i,j,2),uxx(i,j,2),uxy(i,j,2),uyy(i,j,2),gp(m))   ! RHO_U
 uplus(m,3)=uu(u0(i,j,3),ux(i,j,3),uy(i,j,3),uxx(i,j,3),uxy(i,j,3),uyy(i,j,3),gp(m))   ! RHO_V
 uplus(m,4)=uu(u0(i,j,4),ux(i,j,4),uy(i,j,4),uxx(i,j,4),uxy(i,j,4),uyy(i,j,4),gp(m))   ! RHO_E

RHO_VAL = uplus(m,1)
U_VAL   = uplus(m,2)/uplus(m,1)
V_VAL   = uplus(m,3)/uplus(m,1)
E_VAL   = uplus(m,4)/uplus(m,1)
P_VAL   = gama*(gama-1.0)*RHO_VAL*(E_VAL -0.5*(U_VAL**2 + V_VAL**2))
T_VAL   = P_VAL/RHO_VAL 

 D_RHO_DX  =xpuxu(dx,ux(i,j,1),uxx(i,j,1),uxy(i,j,1),gp(m))    
 D_RHO_U_DX=xpuxu(dx,ux(i,j,2),uxx(i,j,2),uxy(i,j,2),gp(m))
 D_RHO_V_DX=xpuxu(dx,ux(i,j,3),uxx(i,j,3),uxy(i,j,3),gp(m))
 D_RHO_E_DX=xpuxu(dx,ux(i,j,4),uxx(i,j,4),uxy(i,j,4),gp(m))

 D_RHO_DY  =xpuyu(dy,uy(i,j,1),uxy(i,j,1),uyy(i,j,1),gp(m))
 D_RHO_U_DY=xpuyu(dy,uy(i,j,2),uxy(i,j,2),uyy(i,j,2),gp(m))
 D_RHO_V_DY=xpuyu(dy,uy(i,j,3),uxy(i,j,3),uyy(i,j,3),gp(m))
 D_RHO_E_DY=xpuyu(dy,uy(i,j,4),uxy(i,j,4),uyy(i,j,4),gp(m))
 	  
DU_DX = (1.0/RHO_VAL)*(D_RHO_U_DX - U_VAL*D_RHO_DX)
DU_DY = (1.0/RHO_VAL)*(D_RHO_U_DY - U_VAL*D_RHO_DY)

DV_DX = (1.0/RHO_VAL)*(D_RHO_V_DX - V_VAL*D_RHO_DX)
DV_DY = (1.0/RHO_VAL)*(D_RHO_V_DY - V_VAL*D_RHO_DY)

DE_DX = (1.0/RHO_VAL)*(D_RHO_E_DX - E_VAL*D_RHO_DX)
DE_DY = (1.0/RHO_VAL)*(D_RHO_E_DY - E_VAL*D_RHO_DY)

DT_DX = (gama-1.0)*gama*(DE_DX - U_VAL*DU_DX - V_VAL*DV_DX)
DT_DY = (gama-1.0)*gama*(DE_DY - U_VAL*DU_DY - V_VAL*DV_DY)

MUI  = T_VAL**VISCOSITY_INDEX
KAPA = T_VAL**VISCOSITY_INDEX  
        
IF(GOV_EQ_SWITCH.EQ."EULER")THEN
pixxu(i,j,m) =0.0
piyyu(i,j,m) =0.0
pixyu(i,j,m) =0.0
qxu(i,j,m) =0.0
qyu(i,j,m) =0.0
bxu(i,j,m) =0.0
ELSEIF(GOV_EQ_SWITCH.EQ."NSF")THEN
pixxu(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DU_DX - (2.0/3.0)*DV_DY))
piyyu(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DV_DY - (2.0/3.0)*DU_DX))
pixyu(i,j,m)= (1.0/Reynolds)*(-MUI*(DU_DY + DV_DX))
qxu(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DX)
qyu(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DY)
bxu(i,j,m)  =  0.0
ELSEIF(GOV_EQ_SWITCH.EQ."NF")THEN
pixxu(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DU_DX - (2.0/3.0)*DV_DY))
piyyu(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DV_DY - (2.0/3.0)*DU_DX))
pixyu(i,j,m)= (1.0/Reynolds)*(-MUI*(DU_DY + DV_DX))
qxu(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DX)
qyu(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DY)
bxu(i,j,m)  = (1.0/Reynolds)*(-Bulk_Ratio*MUI*(DU_DX+ DV_DY))     
END IF
!***************************************************************************************
!***************************************************************************************    
!**** BOTTOM SIDE

umins(m,1)=ub(u0(i,j,1),ux(i,j,1),uy(i,j,1),uxx(i,j,1),uxy(i,j,1),uyy(i,j,1),gp(m))
umins(m,2)=ub(u0(i,j,2),ux(i,j,2),uy(i,j,2),uxx(i,j,2),uxy(i,j,2),uyy(i,j,2),gp(m))
umins(m,3)=ub(u0(i,j,3),ux(i,j,3),uy(i,j,3),uxx(i,j,3),uxy(i,j,3),uyy(i,j,3),gp(m))
umins(m,4)=ub(u0(i,j,4),ux(i,j,4),uy(i,j,4),uxx(i,j,4),uxy(i,j,4),uyy(i,j,4),gp(m))

RHO_VAL = umins(m,1)
U_VAL   = umins(m,2)/umins(m,1)
V_VAL   = umins(m,3)/umins(m,1)
E_VAL   = umins(m,4)/umins(m,1)
P_VAL   = gama*(gama-1.0)*RHO_VAL*(E_VAL -0.5*(U_VAL**2 + V_VAL**2))
T_VAL   = P_VAL/RHO_VAL 

D_RHO_DX   = xpuxb(dx,ux(i,j,1),uxx(i,j,1),uxy(i,j,1),gp(m))
D_RHO_U_DX = xpuxb(dx,ux(i,j,2),uxx(i,j,2),uxy(i,j,2),gp(m))
D_RHO_V_DX = xpuxb(dx,ux(i,j,3),uxx(i,j,3),uxy(i,j,3),gp(m))
D_RHO_E_DX = xpuxb(dx,ux(i,j,4),uxx(i,j,4),uxy(i,j,4),gp(m))

D_RHO_DY   = xpuyb(dy,uy(i,j,1),uxy(i,j,1),uyy(i,j,1),gp(m))
D_RHO_U_DY = xpuyb(dy,uy(i,j,2),uxy(i,j,2),uyy(i,j,2),gp(m))
D_RHO_V_DY = xpuyb(dy,uy(i,j,3),uxy(i,j,3),uyy(i,j,3),gp(m))
D_RHO_E_DY = xpuyb(dy,uy(i,j,4),uxy(i,j,4),uyy(i,j,4),gp(m))
	
DU_DX = (1.0/RHO_VAL)*(D_RHO_U_DX - U_VAL*D_RHO_DX)
DU_DY = (1.0/RHO_VAL)*(D_RHO_U_DY - U_VAL*D_RHO_DY)

DV_DX = (1.0/RHO_VAL)*(D_RHO_V_DX - V_VAL*D_RHO_DX)
DV_DY = (1.0/RHO_VAL)*(D_RHO_V_DY - V_VAL*D_RHO_DY)

DE_DX = (1.0/RHO_VAL)*(D_RHO_E_DX - E_VAL*D_RHO_DX)
DE_DY = (1.0/RHO_VAL)*(D_RHO_E_DY - E_VAL*D_RHO_DY)

DT_DX = (gama-1.0)*gama*(DE_DX - U_VAL*DU_DX - V_VAL*DV_DX)
DT_DY = (gama-1.0)*gama*(DE_DY - U_VAL*DU_DY - V_VAL*DV_DY)

MUI  = T_VAL**VISCOSITY_INDEX
KAPA = T_VAL**VISCOSITY_INDEX  

IF(GOV_EQ_SWITCH.EQ."EULER")THEN
pixxb(i,j,m) =0.0
piyyb(i,j,m) =0.0
pixyb(i,j,m) =0.0
qxb(i,j,m) =0.0
qyb(i,j,m) =0.0
bxb(i,j,m) =0.0
ELSEIF(GOV_EQ_SWITCH.EQ."NSF")THEN
pixxb(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DU_DX - (2.0/3.0)*DV_DY))
piyyb(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DV_DY - (2.0/3.0)*DU_DX))
pixyb(i,j,m)= (1.0/Reynolds)*(-MUI*(DU_DY + DV_DX))
qxb(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DX)
qyb(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DY)
bxb(i,j,m)  =  0.0
ELSEIF(GOV_EQ_SWITCH.EQ."NF")THEN
pixxb(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DU_DX - (2.0/3.0)*DV_DY))
piyyb(i,j,m)= (1.0/Reynolds)*(-MUI*((4.0/3.0)*DV_DY - (2.0/3.0)*DU_DX))
pixyb(i,j,m)= (1.0/Reynolds)*(-MUI*(DU_DY + DV_DX))
qxb(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DX)
qyb(i,j,m)  = (1.0/(Reynolds*Eckert*Prandtl))*(-KAPA*DT_DY)
bxb(i,j,m)  = (1.0/Reynolds)*(-Bulk_Ratio*MUI*(DU_DX+ DV_DY))     
END IF

!***************************************************************************
!***************************************************************************
!****************************************************************

END DO
1500   continue      

       
 do 500 i=0,nx 
 do 500 j=0,ny 
         
   do 501 k=1,4 
   do 501 m=-1,1   
     uplus(m,k)=ur(u0(i+1,j,k),ux(i+1,j,k),uy(i+1,j,k),uxx(i+1,j,k),uxy(i+1,j,k),uyy(i+1,j,k),gp(m)) 
     umins(m,k)=ul(u0(i,j,k),ux(i,j,k),uy(i,j,k),uxx(i,j,k),uxy(i,j,k),uyy(i,j,k),gp(m)) 
  501     continue  
                                     
     do m=-1,1 
        fplus(m,1)=f1(uplus(m,1),uplus(m,2),uplus(m,3),uplus(m,4))  
        fplus(m,2)=f2(uplus(m,1),uplus(m,2),uplus(m,3),uplus(m,4),gama) 
        fplus(m,3)=f3(uplus(m,1),uplus(m,2),uplus(m,3),uplus(m,4))
        fplus(m,4)=f4(uplus(m,1),uplus(m,2),uplus(m,3),uplus(m,4),gama)
        fmins(m,1)=f1(umins(m,1),umins(m,2),umins(m,3),umins(m,4))
        fmins(m,2)=f2(umins(m,1),umins(m,2),umins(m,3),umins(m,4),gama)
        fmins(m,3)=f3(umins(m,1),umins(m,2),umins(m,3),umins(m,4))
        fmins(m,4)=f4(umins(m,1),umins(m,2),umins(m,3),umins(m,4),gama)
          enddo 
          rho0=u0(i,j,1) 
          p0=(gama-1.)*gama*(u0(i,j,4)-0.5*(u0(i,j,2)**2+u0(i,j,3)**2)/rho0) 
          sonic=dsqrt(p0/rho0) 
          alfa=dabs(u0(i,j,2)/u0(i,j,1))+sonic 


  do 505 m=-1,1           
      do 505 k=1,4 
      fv_vedge(i,j,m,k)=0.5*(fplus(m,k)+fmins(m,k)-alfa*(uplus(m,k)-umins(m,k))) 
  505     continue 
  
  !     calculate visous stress  !
     do m=-1,1
       fvplus(m,1)=f1v(pixxr(i+1,j,m),pixyr(i+1,j,m),piyyr(i+1,j,m),uplus(m,2)/uplus(m,1), &
       uplus(m,3)/uplus(m,1),qxr(i+1,j,m),qyr(i+1,j,m))
	   fvplus(m,2)=f2v(pixxr(i+1,j,m),pixyr(i+1,j,m),piyyr(i+1,j,m),uplus(m,2)/uplus(m,1), &
     uplus(m,3)/uplus(m,1),qxr(i+1,j,m),qyr(i+1,j,m),bxr(i+1,j,m),Bulk_Ratio)
	   fvplus(m,3)=f3v(pixxr(i+1,j,m),pixyr(i+1,j,m),piyyr(i+1,j,m),uplus(m,2)/uplus(m,1), &
     uplus(m,3)/uplus(m,1),qxr(i+1,j,m),qyr(i+1,j,m))
       fvplus(m,4)=f4v(pixxr(i+1,j,m),pixyr(i+1,j,m),piyyr(i+1,j,m),uplus(m,2)/uplus(m,1), &
       uplus(m,3)/uplus(m,1),qxr(i+1,j,m),qyr(i+1,j,m),bxr(i+1,j,m),Bulk_Ratio) 
          
       fvmins(m,1)=f1v(pixxl(i,j,m),pixyl(i,j,m),piyyl(i,j,m),umins(m,2)/umins(m,1), &
       umins(m,3)/umins(m,1),qxl(i,j,m),qyl(i,j,m))
	   fvmins(m,2)=f2v(pixxl(i,j,m),pixyl(i,j,m),piyyl(i,j,m),umins(m,2)/umins(m,1), &
     umins(m,3)/umins(m,1),qxl(i,j,m),qyl(i,j,m),bxl(i,j,m),Bulk_Ratio)     
	   fvmins(m,3)=f3v(pixxl(i,j,m),pixyl(i,j,m),piyyl(i,j,m),umins(m,2)/umins(m,1), &
     umins(m,3)/umins(m,1),qxl(i,j,m),qyl(i,j,m))
	   fvmins(m,4)=f4v(pixxl(i,j,m),pixyl(i,j,m),piyyl(i,j,m),umins(m,2)/umins(m,1), &
     umins(m,3)/umins(m,1),qxl(i,j,m),qyl(i,j,m),bxl(i,j,m),Bulk_Ratio)
	   enddo

	 do m=-1,1
	 do k=1,4
       fv_vedge(i,j,m,k)=fv_vedge(i,j,m,k)+0.5*(fvplus(m,k)+fvmins(m,k))
	 enddo
	 enddo 

!     calculate visous stress flux !
               
           
  do 511 k=1,4 
  do 511 m=-1,1   
     uplus(m,k)=uu(u0(i,j+1,k),ux(i,j+1,k),uy(i,j+1,k),uxx(i,j+1,k),uxy(i,j+1,k),uyy(i,j+1,k),gp(m))
     umins(m,k)=ub(u0(i,j,k),ux(i,j,k),uy(i,j,k),uxx(i,j,k),uxy(i,j,k),uyy(i,j,k),gp(m)) 
  511     continue                                     
  
  do m=-1,1 
        fplus(m,1)=g1(uplus(m,1),uplus(m,2),uplus(m,3),uplus(m,4))
        fplus(m,2)=g2(uplus(m,1),uplus(m,2),uplus(m,3),uplus(m,4))
        fplus(m,3)=g3(uplus(m,1),uplus(m,2),uplus(m,3),uplus(m,4),gama)
        fplus(m,4)=g4(uplus(m,1),uplus(m,2),uplus(m,3),uplus(m,4),gama)
        fmins(m,1)=g1(umins(m,1),umins(m,2),umins(m,3),umins(m,4))
        fmins(m,2)=g2(umins(m,1),umins(m,2),umins(m,3),umins(m,4))
        fmins(m,3)=g3(umins(m,1),umins(m,2),umins(m,3),umins(m,4),gama)
        fmins(m,4)=g4(umins(m,1),umins(m,2),umins(m,3),umins(m,4),gama)
   enddo 
          rho0=u0(i,j,1) 
          p0=(gama-1.)*gama*(u0(i,j,4)-0.5*(u0(i,j,2)**2+u0(i,j,3)**2)/rho0) 
          sonic=dsqrt(p0/rho0) 
          alfa=dabs(u0(i,j,3)/u0(i,j,1))+sonic 

          do 515 m=-1,1           
            do 515 k=1,4 
              fv_hedge(i,j,m,k)=0.5*(fplus(m,k)+fmins(m,k)-alfa*(uplus(m,k)-umins(m,k))) 
  515     continue 
  
  !     calculate visous stress  !
         do m=-1,1
       fvplus(m,1)=g1v(pixxu(i,j+1,m),pixyu(i,j+1,m),piyyu(i,j+1,m),uplus(m,2)/uplus(m,1), &
       uplus(m,3)/uplus(m,1),qxu(i,j+1,m),qyu(i,j+1,m))
	   fvplus(m,2)=g2v(pixxu(i,j+1,m),pixyu(i,j+1,m),piyyu(i,j+1,m),uplus(m,2)/uplus(m,1), &
     uplus(m,3)/uplus(m,1),qxu(i,j+1,m),qyu(i,j+1,m))
	   fvplus(m,3)=g3v(pixxu(i,j+1,m),pixyu(i,j+1,m),piyyu(i,j+1,m),uplus(m,2)/uplus(m,1), &
     uplus(m,3)/uplus(m,1),qxu(i,j+1,m),qyu(i,j+1,m),bxu(i,j+1,m),Bulk_Ratio)
	   fvplus(m,4)=g4v(pixxu(i,j+1,m),pixyu(i,j+1,m),piyyu(i,j+1,m),uplus(m,2)/uplus(m,1), &
     uplus(m,3)/uplus(m,1),qxu(i+1,j,m),qyu(i+1,j,m),bxu(i,j+1,m),Bulk_Ratio)

	   fvmins(m,1)=g1v(pixxb(i,j,m),pixyb(i,j,m),piyyb(i,j,m),umins(m,2)/umins(m,1), &
     umins(m,3)/umins(m,1),qxb(i,j,m),qyb(i,j,m))
	   fvmins(m,2)=g2v(pixxb(i,j,m),pixyb(i,j,m),piyyb(i,j,m),umins(m,2)/umins(m,1), &
     umins(m,3)/umins(m,1),qxb(i,j,m),qyb(i,j,m))
	   fvmins(m,3)=g3v(pixxb(i,j,m),pixyb(i,j,m),piyyb(i,j,m),umins(m,2)/umins(m,1), &
     umins(m,3)/umins(m,1),qxb(i,j,m),qyb(i,j,m),bxb(i,j,m),Bulk_Ratio)      
	   fvmins(m,4)=g4v(pixxb(i,j,m),pixyb(i,j,m),piyyb(i,j,m),umins(m,2)/umins(m,1), &
     umins(m,3)/umins(m,1),qxb(i,j,m),qyb(i,j,m),bxb(i,j,m),Bulk_Ratio)

	   enddo
  
  do m=-1,1
  do k=1,4

  fv_hedge(i,j,m,k)=fv_hedge(i,j,m,k)+0.5*(fvplus(m,k)+fvmins(m,k))

 enddo
  enddo
   
   
500 continue    
 
    
DEALLOCATE(GP)
DEALLOCATE(uplus) 
DEALLOCATE(umins) 
DEALLOCATE(fplus) 
DEALLOCATE(fmins) 
DEALLOCATE(fvplus) 
DEALLOCATE(fvmins) 

DEALLOCATE(pixxr) 
DEALLOCATE(piyyr) 
DEALLOCATE(pixyr) 
DEALLOCATE(qxr)
DEALLOCATE(qyr)
DEALLOCATE(bxr)

DEALLOCATE(pixxl) 
DEALLOCATE(piyyl) 
DEALLOCATE(pixyl) 
DEALLOCATE(qxl)
DEALLOCATE(qyl)
DEALLOCATE(bxl)

DEALLOCATE(pixxu) 
DEALLOCATE(piyyu) 
DEALLOCATE(pixyu) 
DEALLOCATE(qxu)
DEALLOCATE(qyu)
DEALLOCATE(bxu)

DEALLOCATE(pixxb) 
DEALLOCATE(piyyb) 
DEALLOCATE(pixyb) 
DEALLOCATE(qxb)
DEALLOCATE(qyb)
DEALLOCATE(bxb)
       
    
END SUBROUTINE  
      
