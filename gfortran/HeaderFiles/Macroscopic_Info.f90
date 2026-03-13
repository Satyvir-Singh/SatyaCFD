
!=============================================================================================!
MODULE MACROSCOPIC_INFO
IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE,DIMENSION (:,:,:,:):: pixx,pixy,piyy

DOUBLE PRECISION, ALLOCATABLE,DIMENSION (:,:,:,:):: QX,QY
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:,:,:):: bx

END MODULE