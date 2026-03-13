!**********************************************************************
!***********       SHOCK - VORTEX INTERACTION             ************!
!********   2D- DISCONTINUOUS GALERKIN METHOD              ***********!
!*****                   SATYVIR SINGH                         *******!
!**                        2017                                    ***!
!**********************************************************************

 MODULE TIME_INFO
 
 IMPLICIT NONE

 INTEGER:: ITERATION
 INTEGER:: RGK_STEP
 DOUBLE PRECISION:: RGK_COEFF(5)

 DOUBLE PRECISION:: DT
 DOUBLE PRECISION:: CURRENT_TIME

 DOUBLE PRECISION:: CFL_NUMBER
 DOUBLE PRECISION:: C_TIME_STEP

 DOUBLE PRECISION :: FINAL_FLOW_TIME
 DOUBLE PRECISION :: FLOW_TIME_STAR, REF_TIME

 END MODULE