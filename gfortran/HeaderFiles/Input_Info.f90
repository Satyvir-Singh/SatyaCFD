!**********************************************************************
!***********       SHOCK - VORTEX INTERACTION             ************!
!********   2D- DISCONTINUOUS GALERKIN METHOD              ***********!
!*****                   SATYVIR SINGH                         *******!
!**                        2017                                    ***!
!**********************************************************************

 MODULE INPUT_INFO
 IMPLICIT NONE


 INTEGER:: NEW_SIMULATION
 CHARACTER(LEN=40):: GOV_EQ_SWITCH
 CHARACTER(LEN=40):: GAS_TYPE

!**** THE R13, R23 ARE CONSTANT VALUES 
 DOUBLE PRECISION :: R13=1.0/3.0
 DOUBLE PRECISION :: R23=2.0/3.0

 DOUBLE PRECISION :: MAX_ERROR
 INTEGER          :: MAX_ITERATION
 INTEGER          :: NPRINT
 DOUBLE PRECISION :: SHOCK_MACH, SHOCK_POSITION

END MODULE