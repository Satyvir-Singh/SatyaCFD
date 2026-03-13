     
 
SUBROUTINE calright
USE GAS_INFO
USE GRID_INFO
USE INPUT_INFO
USE DIMENSIONLESS_INFO
USE COEFFICIENTS
USE MACROSCOPIC_INFO
IMPLICIT NONE 
INTEGER :: I,J,K
INTEGER :: M,N
INTEGER :: i1,i2

DOUBLE PRECISION:: uv(-1:1,-1:1,4)
DOUBLE PRECISION:: gp(-1:1)
DOUBLE PRECISION:: rint(0:5,4)
DOUBLE PRECISION:: flux(0:5,4) 
DOUBLE PRECISION:: fluxr(0:5,4)
DOUBLE PRECISION:: fluxl(0:5,4)
DOUBLE PRECISION:: fluxu(0:5,4)
DOUBLE PRECISION:: fluxb(0:5,4) 
DOUBLE PRECISION:: fv(-1:1,-1:1)

DOUBLE PRECISION :: xpux,xpuy
DOUBLE PRECISION :: F1,F2,F3,F4,G1,G2,G3,G4
      
DOUBLE PRECISION :: f1_int,f1x_int,f1y_int
DOUBLE PRECISION :: f2_int,f2x_int,f2y_int
DOUBLE PRECISION :: f3_int,f3x_int,f3y_int
DOUBLE PRECISION :: f4_int,f4x_int,f4y_int

DOUBLE PRECISION :: G1_int,G1x_int,G1y_int
DOUBLE PRECISION :: G2_int,G2x_int,G2y_int
DOUBLE PRECISION :: G3_int,G3x_int,G3y_int
DOUBLE PRECISION :: G4_int,G4x_int,G4y_int
DOUBLE PRECISION :: rint_area, rint_edge
DOUBLE PRECISION :: w_1,w0,w1

DOUBLE PRECISION :: RHO_VAL,U_VAL,V_VAL,E_VAL,P_VAL,T_VAL       
DOUBLE PRECISION :: D_RHO_DX,D_RHO_U_DX,D_RHO_V_DX,D_RHO_E_DX
DOUBLE PRECISION :: D_RHO_DY,D_RHO_U_DY,D_RHO_V_DY,D_RHO_E_DY
DOUBLE PRECISION :: DU_DX,DU_DY,DV_DX,DV_DY,DE_DX,DE_DY
DOUBLE PRECISION :: DT_DX,DT_DY
DOUBLE PRECISION :: MUI,KAPA

!===========================================================
!--------> gauss points          
 gp(-1)=-sqrt(0.6) 
 gp(0)=0. 
 gp(1)=sqrt(0.6)

!--------> flux calculation
call calfv_edge
 
!-------->  calculate right hand side value
  DO 201 i=1,nx 
  DO 201 j=1,ny     
 
 !*** U VALUE AT GAUSS POINT IS CALCULATED
  DO k=1,4  
  DO i1=-1,1  
  DO i2=-1,1 
     uv(i1,i2,k)= u0(i,j,k) +ux(i,j,k)*gp(i1) +uy(i,j,k)*gp(i2) + &
                  uxx(i,j,k)*(gp(i1)*gp(i1)-1.0/3.0) +uxy(i,j,k)*gp(i1)*gp(i2) +uyy(i,j,k)*(gp(i2)*gp(i2)-1.0/3.0)
  END DO      
  END DO
  END DO

!*** 
DO i1=-1,1
DO i2=-1,1

RHO_VAL = uv(i1,i2,1) 
U_VAL   = uv(i1,i2,2)/uv(i1,i2,1)  
V_VAL   = uv(i1,i2,3)/uv(i1,i2,1)  
E_VAL   = uv(i1,i2,4)/uv(i1,i2,1)
P_VAL   = gama*(gama-1.0)*RHO_VAL*(E_VAL -0.5*(U_VAL**2 + V_VAL**2))
T_VAL   = P_VAL/RHO_VAL
  
D_RHO_DX   = xpux(dx,ux(i,j,1),uxx(i,j,1),uxy(i,j,1),gp(i1),gp(i2))     
D_RHO_U_DX = xpux(dx,ux(i,j,2),uxx(i,j,2),uxy(i,j,2),gp(i1),gp(i2))
D_RHO_V_DX = xpux(dx,ux(i,j,3),uxx(i,j,3),uxy(i,j,3),gp(i1),gp(i2))
D_RHO_E_DX = xpux(dx,ux(i,j,4),uxx(i,j,4),uxy(i,j,4),gp(i1),gp(i2))

D_RHO_DY   = xpuy(dy,uy(i,j,1),uxy(i,j,1),uyy(i,j,1),gp(i1),gp(i2))
D_RHO_U_DY = xpuy(dy,uy(i,j,2),uxy(i,j,2),uyy(i,j,2),gp(i1),gp(i2))
D_RHO_V_DY = xpuy(dy,uy(i,j,3),uxy(i,j,3),uyy(i,j,3),gp(i1),gp(i2))
D_RHO_E_DY = xpuy(dy,uy(i,j,4),uxy(i,j,4),uyy(i,j,4),gp(i1),gp(i2))
      
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
pixx(i,j,i1,i2) =0.0
piyy(i,j,i1,i2) =0.0
pixy(i,j,i1,i2) =0.0
qx(i,j,i1,i2) =0.0
qy(i,j,i1,i2) =0.0
bx(i,j,i1,i2) =0.0
ELSEIF(GOV_EQ_SWITCH.EQ."NSF")THEN
pixx(i,j,i1,i2)= -MUI*((4.0/3.0)*DU_DX - (2.0/3.0)*DV_DY)
piyy(i,j,i1,i2)= -MUI*((4.0/3.0)*DV_DY - (2.0/3.0)*DU_DX)
pixy(i,j,i1,i2)= -MUI*(DU_DY + DV_DX)
qx(i,j,i1,i2)  = -KAPA*DT_DX
qy(i,j,i1,i2)  = -KAPA*DT_DY
bx(i,j,i1,i2)  =  0.0
ELSEIF(GOV_EQ_SWITCH.EQ."NF")THEN
pixx(i,j,i1,i2)= -MUI*((4.0/3.0)*DU_DX - (2.0/3.0)*DV_DY)
piyy(i,j,i1,i2)= -MUI*((4.0/3.0)*DV_DY - (2.0/3.0)*DU_DX)
pixy(i,j,i1,i2)= -MUI*(DU_DY + DV_DX)
qx(i,j,i1,i2)  = -KAPA*DT_DX
qy(i,j,i1,i2)  = -KAPA*DT_DY
bx(i,j,i1,i2)  = -Bulk_Ratio*MUI*(DU_DX+DV_DY)
END IF

END DO
END DO
           
 do 251 m=-1,1 
 do 251 n=-1,1 
   fv(m,n)=f1(uv(m,n,1),uv(m,n,2),uv(m,n,3),uv(m,n,4)) + 0.0 
  251    continue
      
  f1_int=rint_area(fv(-1,-1),fv(-1,1),fv(1,-1),fv(1,1),fv(0,-1),fv(0,1),fv(-1,0),fv(1,0),fv(0,0))   
  f1x_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(-1),fv(1,-1)*gp(1),fv(1,1)*gp(1),fv(0,-1)*gp(0),  &
                           fv(0,1)*gp(0),fv(-1,0)*gp(-1),fv(1,0)*gp(1),fv(0,0)*gp(0))
  f1y_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(1),fv(1,-1)*gp(-1),fv(1,1)*gp(1),fv(0,-1)*gp(-1), &
                  fv(0,1)*gp(1),fv(-1,0)*gp(0),fv(1,0)*gp(0),fv(0,0)*gp(0)) 
 
  do 252 m=-1,1 
  do 252 n=-1,1 
       fv(m,n)=f2(uv(m,n,1),uv(m,n,2),uv(m,n,3),uv(m,n,4),gama)+(1.0/Reynolds)*(pixx(i,j,m,n)+Bulk_Ratio*bx(i,j,m,n))
  252   continue

  f2_int=rint_area(fv(-1,-1),fv(-1,1),fv(1,-1),fv(1,1),fv(0,-1),fv(0,1),fv(-1,0),fv(1,0),fv(0,0))         
  f2x_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(-1),fv(1,-1)*gp(1),fv(1,1)*gp(1),fv(0,-1)*gp(0), &
            fv(0,1)*gp(0),fv(-1,0)*gp(-1),fv(1,0)*gp(1),fv(0,0)*gp(0)) 
  f2y_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(1),fv(1,-1)*gp(-1),fv(1,1)*gp(1),fv(0,-1)*gp(-1), &
            fv(0,1)*gp(1),fv(-1,0)*gp(0),fv(1,0)*gp(0),fv(0,0)*gp(0)) 
 
 do 253 m=-1,1 
 do 253 n=-1,1 
       fv(m,n)=f3(uv(m,n,1),uv(m,n,2),uv(m,n,3),uv(m,n,4))+(1.0/Reynolds)*pixy(i,j,m,n)          
  253     continue 

  f3_int=rint_area(fv(-1,-1),fv(-1,1),fv(1,-1),fv(1,1),fv(0,-1),fv(0,1),fv(-1,0),fv(1,0),fv(0,0))         
  f3x_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(-1),fv(1,-1)*gp(1),fv(1,1)*gp(1),fv(0,-1)*gp(0), &
                          fv(0,1)*gp(0),fv(-1,0)*gp(-1),fv(1,0)*gp(1),fv(0,0)*gp(0)) 
  f3y_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(1),fv(1,-1)*gp(-1),fv(1,1)*gp(1),fv(0,-1)*gp(-1), &
                   fv(0,1)*gp(1),fv(-1,0)*gp(0),fv(1,0)*gp(0),fv(0,0)*gp(0))
           
 do 254 m=-1,1 
 do 254 n=-1,1 
   fv(m,n)=f4(uv(m,n,1),uv(m,n,2),uv(m,n,3),uv(m,n,4),gama)+ &
           (1.0/Reynolds)*((pixx(i,j,m,n)+Bulk_Ratio*bx(i,j,m,n))*uv(m,n,2)/uv(m,n,1) &
           + pixy(i,j,m,n)*uv(m,n,3)/uv(m,n,1) + qx(i,j,m,n)/(Eckert*Prandtl) )
  254     continue 
         
 f4_int=rint_area(fv(-1,-1),fv(-1,1),fv(1,-1),fv(1,1),fv(0,-1),fv(0,1),fv(-1,0),fv(1,0),fv(0,0))        
 f4x_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(-1),fv(1,-1)*gp(1),fv(1,1)*gp(1),fv(0,-1)*gp(0), &
          fv(0,1)*gp(0),fv(-1,0)*gp(-1),fv(1,0)*gp(1),fv(0,0)*gp(0)) 
 f4y_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(1),fv(1,-1)*gp(-1),fv(1,1)*gp(1),fv(0,-1)*gp(-1), &
                  fv(0,1)*gp(1),fv(-1,0)*gp(0),fv(1,0)*gp(0),fv(0,0)*gp(0)) 
      
 do 261 m=-1,1 
 do 261 n=-1,1 
 fv(m,n)=g1(uv(m,n,1),uv(m,n,2),uv(m,n,3),uv(m,n,4))
  261     continue        
 g1_int=rint_area(fv(-1,-1),fv(-1,1),fv(1,-1),fv(1,1),fv(0,-1),fv(0,1),fv(-1,0),fv(1,0),fv(0,0))         
 g1x_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(-1),fv(1,-1)*gp(1),fv(1,1)*gp(1),fv(0,-1)*gp(0),  &
                  fv(0,1)*gp(0),fv(-1,0)*gp(-1),fv(1,0)*gp(1),fv(0,0)*gp(0)) 
 g1y_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(1),fv(1,-1)*gp(-1),fv(1,1)*gp(1),fv(0,-1)*gp(-1), &
                  fv(0,1)*gp(1),fv(-1,0)*gp(0),fv(1,0)*gp(0),fv(0,0)*gp(0)) 
                 
 do 262 m=-1,1 
 do 262 n=-1,1 
       fv(m,n)=g2(uv(m,n,1),uv(m,n,2),uv(m,n,3),uv(m,n,4))+(1.0/Reynolds)*pixy(i,j,m,n)                         
  262     continue        
 g2_int=rint_area(fv(-1,-1),fv(-1,1),fv(1,-1),fv(1,1),fv(0,-1),fv(0,1),fv(-1,0),fv(1,0),fv(0,0))         
 g2x_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(-1),fv(1,-1)*gp(1),fv(1,1)*gp(1),fv(0,-1)*gp(0), &
                  fv(0,1)*gp(0),fv(-1,0)*gp(-1),fv(1,0)*gp(1),fv(0,0)*gp(0)) 
 g2y_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(1),fv(1,-1)*gp(-1),fv(1,1)*gp(1),fv(0,-1)*gp(-1), & 
            fv(0,1)*gp(1),fv(-1,0)*gp(0),fv(1,0)*gp(0),fv(0,0)*gp(0))
                 
 do 263 m=-1,1 
 do 263 n=-1,1 
      fv(m,n)=g3(uv(m,n,1),uv(m,n,2),uv(m,n,3),uv(m,n,4),gama)+(1.0/Reynolds)*(piyy(i,j,m,n)+Bulk_Ratio*bx(i,j,m,n))
  263     continue        
 g3_int=rint_area(fv(-1,-1),fv(-1,1),fv(1,-1),fv(1,1),fv(0,-1),fv(0,1),fv(-1,0),fv(1,0),fv(0,0))        
 g3x_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(-1),fv(1,-1)*gp(1),fv(1,1)*gp(1),fv(0,-1)*gp(0), &
                           fv(0,1)*gp(0),fv(-1,0)*gp(-1),fv(1,0)*gp(1),fv(0,0)*gp(0)) 
 g3y_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(1),fv(1,-1)*gp(-1),fv(1,1)*gp(1),fv(0,-1)*gp(-1), & 
                  fv(0,1)*gp(1),fv(-1,0)*gp(0),fv(1,0)*gp(0),fv(0,0)*gp(0)) 
 
 do 264 m=-1,1 
 do 264 n=-1,1 
     fv(m,n)=g4(uv(m,n,1),uv(m,n,2),uv(m,n,3),uv(m,n,4),gama)+ &
     (1.0/Reynolds)*( pixy(i,j,m,n)*uv(m,n,2)/uv(m,n,1) + (piyy(i,j,m,n)+Bulk_Ratio*bx(i,j,m,n))*uv(m,n,3)/uv(m,n,1) &
     + qy(i,j,m,m)/(Eckert*Prandtl))   
  264     continue  
  
  g4_int=rint_area(fv(-1,-1),fv(-1,1),fv(1,-1),fv(1,1),fv(0,-1),fv(0,1),fv(-1,0),fv(1,0),fv(0,0))         
  g4x_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(-1), fv(1,-1)*gp(1),fv(1,1)*gp(1),fv(0,-1)*gp(0), &
                  fv(0,1)*gp(0),fv(-1,0)*gp(-1),fv(1,0)*gp(1),fv(0,0)*gp(0)) 
  g4y_int=rint_area(fv(-1,-1)*gp(-1),fv(-1,1)*gp(1),fv(1,-1)*gp(-1),fv(1,1)*gp(1),fv(0,-1)*gp(-1), &
                  fv(0,1)*gp(1),fv(-1,0)*gp(0),fv(1,0)*gp(0),fv(0,0)*gp(0))
                 
 do k=1,4 
   rint(0,k)=0. 
 enddo  

 rint(1,1)=0.5*dy*f1_int+0 
 rint(2,1)=0+0.5*dx*g1_int 
 rint(3,1)=dy*f1x_int+0 
 rint(4,1)=dy/2.*f1y_int+dx/2.*g1x_int 
 rint(5,1)=0+dx*g1y_int 
           
 rint(1,2)=0.5*dy*f2_int+0 
 rint(2,2)=0+0.5*dx*g2_int 
 rint(3,2)=dy*f2x_int+0 
 rint(4,2)=dy/2.*f2y_int+dx/2.*g2x_int 
 rint(5,2)=0+dx*g2y_int 
           
 rint(1,3)=0.5*dy*f3_int+0 
 rint(2,3)=0+0.5*dx*g3_int 
 rint(3,3)=dy*f3x_int+0 
 rint(4,3)=dy/2.*f3y_int+dx/2.*g3x_int 
 rint(5,3)=0+dx*g3y_int 
           
 rint(1,4)=0.5*dy*f4_int+0 
 rint(2,4)=0+0.5*dx*g4_int 
 rint(3,4)=dy*f4x_int+0 
 rint(4,4)=dy/2.*f4y_int+dx/2.*g4x_int 
 rint(5,4)=0+dx*g4y_int 


          do k=1,4 
             fluxr(0,k)=rint_edge(fv_vedge(i,j,-1,k),fv_vedge(i,j,0,k),fv_vedge(i,j,1,k))                    
          enddo 
          do k=1,4 
             fluxr(1,k)=rint_edge(fv_vedge(i,j,-1,k),fv_vedge(i,j,0,k),fv_vedge(i,j,1,k))                    
          enddo 
          do k=1,4 
             fluxr(2,k)=rint_edge(fv_vedge(i,j,-1,k)*gp(-1),fv_vedge(i,j,0,k)*gp(0),fv_vedge(i,j,1,k)*gp(1))                    
          enddo 
          do k=1,4 
             fluxr(3,k)=r23*rint_edge(fv_vedge(i,j,-1,k),fv_vedge(i,j,0,k),fv_vedge(i,j,1,k))                    
          enddo 
          do k=1,4 
             fluxr(4,k)=rint_edge(fv_vedge(i,j,-1,k)*gp(-1),fv_vedge(i,j,0,k)*gp(0),fv_vedge(i,j,1,k)*gp(1))                    
          enddo 
          w_1=gp(-1)*gp(-1)-r13 
          w0 =gp( 0)*gp( 0)-r13 
          w1 =gp( 1)*gp( 1)-r13 
          do k=1,4 
             fluxr(5,k)=rint_edge(fv_vedge(i,j,-1,k)*w_1,fv_vedge(i,j,0,k)*w0,fv_vedge(i,j,1,k)*w1)                    
          enddo 
         
          do k=1,4 
             fluxl(0,k)=rint_edge(fv_vedge(i-1,j,-1,k),fv_vedge(i-1,j,0,k),fv_vedge(i-1,j,1,k))               !!!!!!!*dy/2                                 
          enddo 
          do k=1,4 
             fluxl(1,k)=-rint_edge(fv_vedge(i-1,j,-1,k),fv_vedge(i-1,j,0,k),fv_vedge(i-1,j,1,k))                    
          enddo 
          do k=1,4 
             fluxl(2,k)=rint_edge(fv_vedge(i-1,j,-1,k)*gp(-1),fv_vedge(i-1,j,0,k)*gp(0),fv_vedge(i-1,j,1,k)*gp(1))                    
          enddo 
          do k=1,4 
             fluxl(3,k)=r23*rint_edge(fv_vedge(i-1,j,-1,k),fv_vedge(i-1,j,0,k),fv_vedge(i-1,j,1,k))                    
          enddo 
          do k=1,4 
             fluxl(4,k)=-rint_edge(fv_vedge(i-1,j,-1,k)*gp(-1),fv_vedge(i-1,j,0,k)*gp(0),fv_vedge(i-1,j,1,k)*gp(1))                    
          enddo 
          w_1=gp(-1)*gp(-1)-r13 
          w0 =gp( 0)*gp( 0)-r13 
          w1 =gp( 1)*gp( 1)-r13 
          do k=1,4 
             fluxl(5,k)=rint_edge(fv_vedge(i,j,-1,k)*w_1,fv_vedge(i,j,0,k)*w0,fv_vedge(i,j,1,k)*w1)                    
          enddo 
           
          do k=1,4 
             fluxu(0,k)=rint_edge(fv_hedge(i,j,-1,k),fv_hedge(i,j,0,k), fv_hedge(i,j,1,k))                    
          enddo 
          do k=1,4 
             fluxu(1,k)=rint_edge(fv_hedge(i,j,-1,k)*gp(-1),fv_hedge(i,j,0,k)*gp(0),fv_hedge(i,j,1,k)*gp(1)) 
          enddo 
          do k=1,4 
             fluxu(2,k)=rint_edge(fv_hedge(i,j,-1,k),fv_hedge(i,j,0,k),fv_hedge(i,j,1,k)) 
          enddo      
          w_1=gp(-1)*gp(-1)-r13 
          w0 =gp( 0)*gp( 0)-r13 
          w1 =gp( 1)*gp( 1)-r13 
          do k=1,4 
             fluxu(3,k)=rint_edge(fv_hedge(i,j,-1,k)*w_1,fv_hedge(i,j,0,k)*w0,fv_hedge(i,j,1,k)*w1)                    
          enddo 
          do k=1,4 
             fluxu(4,k)=rint_edge(fv_hedge(i,j,-1,k)*gp(-1),fv_hedge(i,j,0,k)*gp(0),fv_hedge(i,j,1,k)*gp(1))                    
          enddo 
           
          do k=1,4 
             fluxu(5,k)=r23*rint_edge(fv_hedge(i,j,-1,k),fv_hedge(i,j,0,k),fv_hedge(i,j,1,k))                    
          enddo 
 
!         bottom edge i,j-1/2 
          do k=1,4 
             fluxb(0,k)=rint_edge(fv_hedge(i,j-1,-1,k),fv_hedge(i,j-1,0,k),fv_hedge(i,j-1,1,k))              !!!!!!!*dx/2                          
          enddo 
          do k=1,4 
             fluxb(1,k)=rint_edge(fv_hedge(i,j-1,-1,k)*gp(-1),fv_hedge(i,j-1,0,k)*gp(0),fv_hedge(i,j-1,1,k)*gp(1)) 
          enddo 
          do k=1,4 
             fluxb(2,k)=-rint_edge(fv_hedge(i,j-1,-1,k),fv_hedge(i,j-1,0,k),fv_hedge(i,j-1,1,k)) 
          enddo      
          w_1=gp(-1)*gp(-1)-r13 
          w0 =gp( 0)*gp( 0)-r13 
          w1 =gp( 1)*gp( 1)-r13 
          do k=1,4 
             fluxb(3,k)=rint_edge(fv_hedge(i,j-1,-1,k)*w_1,fv_hedge(i,j-1,0,k)*w0,fv_hedge(i,j-1,1,k)*w1)
          enddo 
          do k=1,4 
             fluxb(4,k)=-rint_edge(fv_hedge(i,j-1,-1,k)*gp(-1),fv_hedge(i,j-1,0,k)*gp(0),fv_hedge(i,j-1,1,k)*gp(1))
          enddo 
           
          do k=1,4 
             fluxb(5,k)=r23*rint_edge(fv_hedge(i,j-1,-1,k),fv_hedge(i,j-1,0,k),fv_hedge(i,j-1,1,k))
          enddo 
           
          do 280 m=0,5 
            do 280 k=1,4 
              flux(m,k)=(-fluxr(m,k)+fluxl(m,k))*0.5*dy+(-fluxu(m,k)+fluxb(m,k))*0.5*dx 
  280     continue 
!c    *******************     the above is flux     ***************** 

          do k=1,4 
             right0(i,j,k)=rint(0,k)+flux(0,k) 
             rightx(i,j,k)=rint(1,k)+flux(1,k) 
             righty(i,j,k)=rint(2,k)+flux(2,k) 
             rightxx(i,j,k)=rint(3,k)+flux(3,k) 
             rightxy(i,j,k)=rint(4,k)+flux(4,k) 
             rightyy(i,j,k)=rint(5,k)+flux(5,k) 
          enddo 
  201 continue 


END SUBROUTINE     


  