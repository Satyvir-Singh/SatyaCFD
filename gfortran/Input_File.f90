!**********************************************************************
!***********       SHOCK - VORTEX INTERACTION             ************!
!********   2D- DISCONTINUOUS GALERKIN METHOD              ***********!
!*****                   SATYVIR SINGH                         *******!
!**                        2017                                    ***!
!**********************************************************************

SUBROUTINE INPUT_FILE

USE INPUT_INFO
USE TIME_INFO
USE FILE_INFO
USE GRID_INFO
USE GAS_INFO
USE FLOWFIELD_INFO
USE DIMENSIONLESS_INFO

IMPLICIT NONE
INTEGER :: I
CHARACTER(LEN=40) :: DUMMY

!**************************************************!
!**** DEFINE THE FILE NAME FOR POST PROCESSING ****!
!**************************************************!
FILE_DIRECTORY1 ='RESULT/'
FILE_DIRECTORY2 ='RESTART/'
DIRECTORY_2 ='INPUT/'
FILENAME3="RESIDUAL"
FILENAME2="DG2D"
FILETYPE=".plt"
FILETYPE2=".RES"
SUFFIX2 = ".txt"

INPUT_FILENAME = "INPUT_FILE"
!================================================================
! READ INPUT FILE INFORMATION
!================================================================
!****** OPEN GLOBAL INPUT FILE
OPEN(1,FILE=TRIM(DIRECTORY_2)//TRIM(INPUT_FILENAME)//TRIM(SUFFIX2) )

DO I = 1, 6
READ(1,*) 
END DO

!******* NUMBER OF NODES IN COMPUTATIONAL DOMAIN 
READ(1,*) DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, nx
READ(1,*) DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, ny

!****** COMPUTATIONAL DOMAIN INFROMATION
DO I = 1,2
READ(1,*) 
END DO
READ(1,*) DUMMY, DUMMY, DUMMY, xl
READ(1,*) DUMMY, DUMMY, DUMMY, xr
READ(1,*) DUMMY, DUMMY, DUMMY, yb
READ(1,*) DUMMY, DUMMY, DUMMY, yu

!****** BOUNDARY CONDITIONS
DO I = 1,2
READ(1,*) 
END DO
READ(1,*) DUMMY, DUMMY, DUMMY, LEFT_BOUNDARY
READ(1,*) DUMMY, DUMMY, DUMMY, RIGHT_BOUNDARY
READ(1,*) DUMMY, DUMMY, DUMMY, BOTTOM_BOUNDARY
READ(1,*) DUMMY, DUMMY, DUMMY, UPPER_BOUNDARY

!****** STABILITY INFROMATION
DO I = 1,2
READ(1,*) 
END DO
READ(1,*) DUMMY, DUMMY, DUMMY, CFL_NUMBER 
READ(1,*) DUMMY, DUMMY, DUMMY, FINAL_FLOW_TIME
READ(1,*) DUMMY, DUMMY, DUMMY, NPRINT
   
!****** PHYSIACAL QUANTITIES INFROMATION 
DO I = 1,2
READ(1,*) 
END DO  
READ(1,*) DUMMY, DUMMY, DUMMY, GAS_TYPE
READ(1,*) DUMMY, DUMMY, DUMMY, GOV_EQ_SWITCH
READ(1,*) DUMMY, DUMMY, DUMMY, DUMMY,SHOCK_MACH
READ(1,*) DUMMY, DUMMY, DUMMY, SHOCK_POSITION
!READ(1,*) DUMMY, DUMMY, DUMMY, NEW_SIMULATION
                                            
CLOSE(1)

!***********************************************************

RGK_STEP=1

!***** The Maximum number of iteration 
MAX_ITERATION=100000

!***** The Maximum desired accuracy for terminating the simulation 
MAX_ERROR = 1.e-8

NEW_SIMULATION=0            ! NEW_SIMULATION =0 for new simulation,   NEW_SIMULATION =1 FOR OLD SIMULATION

!***** FLOWFIELD INFORMATION
Ref_Pressure=101325.0          !FARFIELD PRESSURE
Ref_Temperature=298.0           !FARFIELD TEMPERATURE
    

END SUBROUTINE