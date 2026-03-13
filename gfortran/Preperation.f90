
       SUBROUTINE PREPARATION
      
      USE GRID_INFO
      USE FILE_INFO
      USE TIME_INFO
      USE COEFFICIENTS
      USE MACROSCOPIC_INFO

!**************************************************!
!****               ARRAY ALLOCATION           ****!
!**************************************************!
      ALLOCATE ( x(0:nx+1) )
      ALLOCATE ( y(0:ny+1) )

      ALLOCATE ( u0 (0:nx+1,0:ny+1,4) )
      ALLOCATE ( ux (0:nx+1,0:ny+1,4) )
      ALLOCATE ( uy (0:nx+1,0:ny+1,4) )
      ALLOCATE ( uxx(0:nx+1,0:ny+1,4) )
      ALLOCATE ( uxy(0:nx+1,0:ny+1,4) )
      ALLOCATE ( uyy(0:nx+1,0:ny+1,4) )

! OLD STEP DATA
      ALLOCATE ( u0_OLD(0:nx+1,0:ny+1,4) )
      ALLOCATE ( ux_OLD(0:nx+1,0:ny+1,4) )
      ALLOCATE ( uy_OLD(0:nx+1,0:ny+1,4) )
      ALLOCATE ( uxx_OLD(0:nx+1,0:ny+1,4) )
      ALLOCATE ( uxy_OLD(0:nx+1,0:ny+1,4) )
      ALLOCATE ( uyy_OLD(0:nx+1,0:ny+1,4) )

! RIGHT HAND SIDE VALUE
      ALLOCATE (  right0(0:nx+1, 0:ny+1, 4) )
      ALLOCATE (  rightx(0:nx+1, 0:ny+1, 4) )
      ALLOCATE (  righty(0:nx+1, 0:ny+1, 4) )
      ALLOCATE (  rightxx(0:nx+1, 0:ny+1,4) )
      ALLOCATE (  rightxy(0:nx+1, 0:ny+1,4) )
      ALLOCATE (  rightyy(0:nx+1, 0:ny+1,4) )

! FLUX VALUES ON EDGE
      ALLOCATE( fv_vedge(0:nx+1, 0:ny+1, -1:1, 4) )
      ALLOCATE( fv_hedge(0:nx+1, 0:ny+1, -1:1, 4) )

! HIGH ORDER MOMENT VARIABLES
	  ALLOCATE( pixx(0:nx+1,0:ny+1,-1:1,-1:1) )
      ALLOCATE( pixy(0:nx+1,0:ny+1,-1:1,-1:1) )
      ALLOCATE( piyy(0:nx+1,0:ny+1,-1:1,-1:1) )

      ALLOCATE( QX(0:nx+1,0:ny+1,-1:1,-1:1) )
      ALLOCATE( QY(0:nx+1,0:ny+1,-1:1,-1:1) )


      ALLOCATE( bx(0:nx+1,0:ny+1,-1:1,-1:1) )
   
!**************************************************!

!**************************************************!     
!****        SECOND ORDER RUNG-KUTTA SCHEME    ****!
!**************************************************!
      IF(RGK_STEP.EQ.3)THEN
        RGK_COEFF (1)=0.1918
        RGK_COEFF (2)=0.4929
        RGK_COEFF (3)=1.0000
      ELSEIF(RGK_STEP.EQ.4)THEN
        RGK_COEFF (1)=0.1084
        RGK_COEFF (2)=0.2602
        RGK_COEFF (3)=0.5052
        RGK_COEFF (4)=1.0000
      ELSEIF(RGK_STEP.EQ.5)THEN
        RGK_COEFF (1)=0.0695
        RGK_COEFF (2)=0.1602
        RGK_COEFF (3)=0.2898
        RGK_COEFF (4)=0.5060
        RGK_COEFF (5)=1.0000
      ELSE
!*** FIRST ORDER SIMPLEST SCHEME
        RGK_STEP=1
        RGK_COEFF (1)=1.00
      END IF
!**************************************************!

      END SUBROUTINE