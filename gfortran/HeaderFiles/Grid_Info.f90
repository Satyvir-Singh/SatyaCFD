
!=============================================================================================!
      MODULE GRID_INFO
      IMPLICIT NONE

      INTEGER :: NX,NY

      DOUBLE PRECISION,ALLOCATABLE:: x(:),y(:)
      DOUBLE PRECISION:: DX,DY
      DOUBLE PRECISION:: FH,FW

      DOUBLE PRECISION:: Xl,Xr
      DOUBLE PRECISION:: Yb,Yu
      
      
      CHARACTER(LEN=40):: LEFT_BOUNDARY
      CHARACTER(LEN=40):: RIGHT_BOUNDARY
      CHARACTER(LEN=40):: BOTTOM_BOUNDARY
      CHARACTER(LEN=40):: UPPER_BOUNDARY

      END MODULE