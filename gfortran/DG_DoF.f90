  
!****************************************************************
!*******INVISICD FLUX FUNCTIONS IN X-DIRECTIONS
!****************************************************************  

      FUNCTION ur(u0,u1,u2,u3,u4,u5,gp)
      IMPLICIT NONE
      DOUBLE PRECISION :: ur
      DOUBLE PRECISION :: u0,u1,u2,u3,u4,u5,gp
      DOUBLE PRECISION :: R13=1./3.
      DOUBLE PRECISION :: R23=2./3.

      ur=u0-u1+r23*u3-r13*u5+(u2-u4)*gp+u5*gp*gp  !from right|<--- 
      RETURN 
      END FUNCTION

      FUNCTION ul(u0,u1,u2,u3,u4,u5,gp)
      IMPLICIT NONE
      DOUBLE PRECISION :: ul
      DOUBLE PRECISION :: u0,u1,u2,u3,u4,u5,gp
      DOUBLE PRECISION :: R13=1./3.
      DOUBLE PRECISION :: R23=2./3.

      ul=u0+u1+r23*u3-r13*u5+(u2+u4)*gp+u5*gp*gp  !from left  --->| 
      RETURN 
      END FUNCTION

      FUNCTION uu(u0,u1,u2,u3,u4,u5,gp)
      IMPLICIT NONE
      DOUBLE PRECISION :: uu
      DOUBLE PRECISION :: u0,u1,u2,u3,u4,u5,gp
      DOUBLE PRECISION :: R13=1./3.
      DOUBLE PRECISION :: R23=2./3.

      uu=u0-u2-r13*u3+r23*u5+(u1-u4)*gp+u3*gp*gp  !from up       
      RETURN 
      END FUNCTION

      FUNCTION ub(u0,u1,u2,u3,u4,u5,gp)
      IMPLICIT NONE
      DOUBLE PRECISION :: ub
      DOUBLE PRECISION :: u0,u1,u2,u3,u4,u5,gp
      DOUBLE PRECISION :: R13=1./3.
      DOUBLE PRECISION :: R23=2./3.

      ub=u0+u2-r13*u3+r23*u5+(u1+u4)*gp+u3*gp*gp  !from bottom       
      RETURN 
      END FUNCTION

!********************************************************************
!*******************************************************************

      FUNCTION xpuxr(xdx,u1,u3,u4,gp)
      IMPLICIT NONE
      DOUBLE PRECISION :: xpuxr
      DOUBLE PRECISION :: xdx,u1,u3,u4,gp
      xpuxr=2.0/xdx*(u1-2.0*u3+u4*gp)    
      RETURN 
      END FUNCTION

      FUNCTION xpuxl(xdx,u1,u3,u4,gp)
      IMPLICIT NONE
      DOUBLE PRECISION :: xpuxl
      DOUBLE PRECISION :: xdx,u1,u3,u4,gp
      xpuxl=2.0/xdx*(u1+2.0*u3+u4*gp)    
      RETURN 
      END FUNCTION

      FUNCTION xpuxu(xdx,u1,u3,u4,gp)
      IMPLICIT NONE
      DOUBLE PRECISION :: xpuxu
      DOUBLE PRECISION :: xdx,u1,u3,u4,gp
      xpuxu=2.0/xdx*(u1-u4+2*u3*gp)    
      RETURN 
      END FUNCTION

      FUNCTION xpuxb(xdx,u1,u3,u4,gp)
      IMPLICIT NONE
      DOUBLE PRECISION :: xpuxb
      DOUBLE PRECISION :: xdx,u1,u3,u4,gp
      xpuxb=2.0/xdx*(u1+u4+2*u3*gp)    
      RETURN 
      END FUNCTION

      FUNCTION xpuyu(xdy,u2,u4,u5,gp)
      IMPLICIT NONE
      DOUBLE PRECISION :: xpuyu
      DOUBLE PRECISION :: xdy,u2,u4,u5,gp
      xpuyu=2.0/xdy*(u2-2.0*u5+u4*gp)   
      RETURN 
      END FUNCTION

      FUNCTION xpuyb(xdy,u2,u4,u5,gp)
      IMPLICIT NONE
      DOUBLE PRECISION :: xpuyb
      DOUBLE PRECISION :: xdy,u2,u4,u5,gp
      xpuyb=2.0/xdy*(u2+2.0*u5+u4*gp)   
      RETURN 
      END FUNCTION

      FUNCTION xpuyr(xdy,u2,u4,u5,gp)
      IMPLICIT NONE
      DOUBLE PRECISION :: xpuyr
      DOUBLE PRECISION :: xdy,u2,u4,u5,gp
      xpuyr=2.0/xdy*(u2-u4+2.0*u5*gp)  
      RETURN 
      END FUNCTION

      FUNCTION xpuyl(xdy,u2,u4,u5,gp)
      IMPLICIT NONE
      DOUBLE PRECISION :: xpuyl
      DOUBLE PRECISION :: xdy,u2,u4,u5,gp
      xpuyl=2.0/xdy*(u2+u4+2.0*u5*gp)  
      RETURN 
      END FUNCTION

!************************************************************
!************************************************************
!******** DERIVATIVES OF DoFs inside the element

      FUNCTION xpux(xdx,u1,u3,u4,gpn,gpm)
      IMPLICIT NONE
      DOUBLE PRECISION :: xpux
      DOUBLE PRECISION :: xdx,u1,u3,u4,gpn,gpm
      xpux=2.0/xdx*(u1+2.0*u3*gpn+u4*gpm)
      RETURN 
      END FUNCTION

      FUNCTION xpuy(xdy,u2,u4,u5,gpn,gpm)
      IMPLICIT NONE
      DOUBLE PRECISION :: xpuy
      DOUBLE PRECISION :: xdy,u2,u4,u5,gpn,gpm
      xpuy=2.0/xdy*(u2+u4*gpn+2.0*u5*gpm)
      RETURN 
      END FUNCTION
