!**********************************************************************
!***********       SHOCK - VORTEX INTERACTION             ************!
!********   2D- DISCONTINUOUS GALERKIN METHOD              ***********!
!*****                   SATYVIR SINGH                         *******!
!**                        2017                                    ***!
!**********************************************************************

MODULE DIMENSIONLESS_INFO
IMPLICIT NONE
      
DOUBLE PRECISION:: Reynolds       
DOUBLE PRECISION:: Eckert       
DOUBLE PRECISION:: Prandtl      
DOUBLE PRECISION:: xndelta, N_delta


DOUBLE PRECISION:: MACH_REF
DOUBLE PRECISION:: REF_LENGTH

DOUBLE PRECISION:: REF_DENSITY, REF_VELOCITY

DOUBLE PRECISION:: REF_SHEAR_STRESS, REF_HEAT_FLUX
DOUBLE PRECISION:: pref,rhoref


END MODULE