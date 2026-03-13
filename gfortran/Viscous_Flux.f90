
!****************************************************************
!*******VISCOUS FLUX FUNCTIONS IN X-DIRECTIONS
!****************************************************************
    
FUNCTION f1v(xpixx,xpixy,xpiyy,xu,xv,xqx,xqy)
IMPLICIT NONE
DOUBLE PRECISION :: f1v
DOUBLE PRECISION :: xpixx,xpixy,xpiyy,xu,xv,xqx,xqy
f1v=0.0        
RETURN 
END FUNCTION
      	
FUNCTION f2v(xpixx,xpixy,xpiyy,xu,xv,xqx,xqy,bx,bulk)     ! SATYA   
IMPLICIT NONE
DOUBLE PRECISION :: f2v
DOUBLE PRECISION :: xpixx,xpixy,xpiyy,xu,xv,xqx,xqy,bx,bulk
f2v=xpixx + bulk*bx  
RETURN 
END FUNCTION

FUNCTION f3v(xpixx,xpixy,xpiyy,xu,xv,xqx,xqy)
IMPLICIT NONE
DOUBLE PRECISION :: f3v
DOUBLE PRECISION :: xpixx,xpixy,xpiyy,xu,xv,xqx,xqy
f3v=xpixy        
RETURN 
END FUNCTION

FUNCTION f4v(xpixx,xpixy,xpiyy,xu,xv,xqx,xqy,bx,bulk)
IMPLICIT NONE
DOUBLE PRECISION :: f4v
DOUBLE PRECISION :: xpixx,xpixy,xpiyy,xu,xv,xqx,xqy,bx,bulk
 f4v=(xpixx + bulk*bx)*xu+xpixy*xv+xqx     
 RETURN 
 END FUNCTION
      

!****************************************************************
!*******VISCOUS FLUX FUNCTIONS IN Y-DIRECTIONS
!****************************************************************

FUNCTION g1v(xpixx,xpixy,xpiyy,xu,xv,xqx,xqy)
IMPLICIT NONE
DOUBLE PRECISION :: g1v
DOUBLE PRECISION :: xpixx,xpixy,xpiyy,xu,xv,xqx,xqy
g1v=0.0      
RETURN 
END FUNCTION

FUNCTION g2v(xpixx,xpixy,xpiyy,xu,xv,xqx,xqy)
IMPLICIT NONE
DOUBLE PRECISION :: g2v
DOUBLE PRECISION :: xpixx,xpixy,xpiyy,xu,xv,xqx,xqy
g2v=xpixy     
RETURN 
END FUNCTION

FUNCTION g3v(xpixx,xpixy,xpiyy,xu,xv,xqx,xqy,bx,bulk)
IMPLICIT NONE
DOUBLE PRECISION :: g3v
DOUBLE PRECISION :: xpixx,xpixy,xpiyy,xu,xv,xqx,xqy,bx,bulk
g3v=xpiyy + bulk*bx      
RETURN 
END FUNCTION

FUNCTION g4v(xpixx,xpixy,xpiyy,xu,xv,xqx,xqy,bx,bulk)
IMPLICIT NONE
DOUBLE PRECISION :: g4v
DOUBLE PRECISION :: xpixx,xpixy,xpiyy,xu,xv,xqx,xqy,bx,bulk
g4v=xpixy*xu+(xpiyy+bulk*bx)*xv+xqy    
RETURN 
END FUNCTION