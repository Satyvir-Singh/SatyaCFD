      subroutine LIMITER
       
       USE GAS_INFO
       USE GRID_INFO
       USE COEFFICIENTS	

      IMPLICIT NONE
      
      DOUBLE PRECISION:: s(4,4)
      DOUBLE PRECISION:: sv(4,4)               !left eignvector s! sAs^1=diag(u,u,u+c,u-c) 
      DOUBLE PRECISION:: vx(4)
      DOUBLE PRECISION:: vy(4)
      DOUBLE PRECISION:: dvplus(4)
      DOUBLE PRECISION:: dvminus(4)
      DOUBLE PRECISION:: tux(4)
      DOUBLE PRECISION:: tuy(4) 
      INTEGER:: I,J,K,M
      
      DOUBLE PRECISION:: RHO,U,V,P,C
      DOUBLE PRECISION:: ALFA,RMINMOD,ERROR

      do 301 i=1,nx                     
        do 301 j=1,ny 
          rho=u0(i,j,1) 
          u=u0(i,j,2)/u0(i,j,1) 
          v=u0(i,j,3)/u0(i,j,1) 
          p=(gama-1.)*gama*(u0(i,j,4)-0.5*rho*(u*u+v*v))
		if(p.lt.0.0)then 
	    write(*,*)'1227',i,j,p,gama,u0(i,j,4),rho,u,v
	    stop
	    endif

          c=dsqrt(p/rho) 
          alfa=(gama-1.)/2./c/c 
          s(1,1)=1.-alfa*(u*u+v*v) 
          s(1,2)=2.*alfa*u 
          s(1,3)=2.*alfa*v 
          s(1,4)=-2.*alfa 
          s(2,1)=-v 
          s(2,2)=0. 
          s(2,3)=1. 
          s(2,4)=0. 
          s(3,1)=0.5*(-u/c+alfa*(u*u+v*v)) 
          s(3,2)=1./2./c-alfa*u 
          s(3,3)=-alfa*v 
          s(3,4)=alfa 
          s(4,1)=0.5*(u/c+alfa*(u*u+v*v)) 
          s(4,2)=-1./2./c-alfa*u 
          s(4,3)=-alfa*v 
          s(4,4)=alfa 
          sv(1,1)=1. 
          sv(1,2)=0. 
          sv(1,3)=1. 
          sv(1,4)=1. 
          sv(2,1)=u 
          sv(2,2)=0. 
          sv(2,3)=u+c 
          sv(2,4)=u-c 
          sv(3,1)=v 
          sv(3,2)=1. 
          sv(3,3)=v 
          sv(3,4)=v 
          sv(4,1)=0.5*(u*u+v*v) 
          sv(4,2)=v 
          sv(4,3)=0.5*(u*u+v*v)+c*u+c*c/(gama-1.) 
          sv(4,4)=0.5*(u*u+v*v)-c*u+c*c/(gama-1.) 
           
          do 310 k=1,4 
            vx(k)=0. 
            dvplus(k)=0. 
            dvminus(k)=0. 
            do 310 m=1,4 
              vx(k)=vx(k)+s(k,m)*ux(i,j,m) 
              dvplus(k)=dvplus(k)+s(k,m)*(u0(i+1,j,m)-u0(i,j,m)) 
              dvminus(k)=dvminus(k)+s(k,m)*(u0(i,j,m)-u0(i-1,j,m)) 
  310     continue 
          do k=1,4 
            vx(k)=rminmod(vx(k),dvplus(k),dvminus(k),dx) 
          enddo 
          do 311 k=1,4 
            tux(k)=0. 
            do 311 m=1,4 
              tux(k)=tux(k)+sv(k,m)*vx(m) 
  311     continue 
           
           
          s(1,1)=1.-alfa*(u*u+v*v) 
          s(1,2)=2.*alfa*u 
          s(1,3)=2.*alfa*v 
          s(1,4)=-2.*alfa 
          s(2,1)=-u 
          s(2,2)=1. 
          s(2,3)=0. 
          s(2,4)=0. 
          s(3,1)=0.5*(-v/c+alfa*(u*u+v*v)) 
          s(3,2)=-alfa*u 
          s(3,3)=1./2./c-alfa*v 
          s(3,4)=alfa 
          s(4,1)=0.5*(v/c+alfa*(u*u+v*v)) 
          s(4,2)=-alfa*u 
          s(4,3)=-1./2./c-alfa*v 
          s(4,4)=alfa 
          sv(1,1)=1. 
          sv(1,2)=0. 
          sv(1,3)=1. 
          sv(1,4)=1. 
          sv(2,1)=u 
          sv(2,2)=1. 
          sv(2,3)=u 
          sv(2,4)=u 
          sv(3,1)=v 
          sv(3,2)=0. 
          sv(3,3)=v+c 
          sv(3,4)=v-c 
          sv(4,1)=0.5*(u*u+v*v) 
          sv(4,2)=u 
          sv(4,3)=0.5*(u*u+v*v)+c*v+c*c/(gama-1.) 
          sv(4,4)=0.5*(u*u+v*v)-c*v+c*c/(gama-1.) 
           
          do 320 k=1,4 
            vy(k)=0. 
            dvplus(k)=0. 
            dvminus(k)=0. 
            do 320 m=1,4 
              vy(k)=vy(k)+s(k,m)*uy(i,j,m) 
              dvplus(k)=dvplus(k)+s(k,m)*(u0(i,j+1,m)-u0(i,j,m)) 
              dvminus(k)=dvminus(k)+s(k,m)*(u0(i,j,m)-u0(i,j-1,m)) 
  320     continue 
          do k=1,4 
            vy(k)=rminmod(vy(k),dvplus(k),dvminus(k),dy) 
          enddo 
          do 321 k=1,4 
            tuy(k)=0. 
            do 321 m=1,4 
              tuy(k)=tuy(k)+sv(k,m)*vy(m) 
  321     continue                               
          error=0. 
          do k=1,4 
            error=error+dabs(tux(k)-ux(i,j,k)) 
          enddo 
          do k=1,4 
            error=error+dabs(tuy(k)-uy(i,j,k)) 
          enddo 
          if(error.gt.1.e-6) then 
            do k=1,4 
              ux(i,j,k)=tux(k) 
              uy(i,j,k)=tuy(k) 
              uxx(i,j,k)=0. 
              uxy(i,j,k)=0. 
              uyy(i,j,k)=0. 
            enddo 
          endif      
  301 continue 
      end                           
            
           
!!   minmod function      minmod(x,y,z) 

      function rminmod(x,y,z,dx)
	implicit real*8 (a-h,o-z)  
      rm=0. 
      if(dabs(x).le.rm*dx*dx) then 
        rminmod=x 
      else if(x.gt.0.and.y.gt.0.and.z.gt.0) then 
        rminmod=min(x,y,z) 
      else if(x.lt.0.and.y.lt.0.and.z.lt.0) then 
        rminmod=max(x,y,z) 
      else 
        rminmod=0.0 
      endif         
      end  