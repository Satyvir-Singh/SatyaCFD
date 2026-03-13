!**********************************************************************
!***********       SHOCK - VORTEX INTERACTION             ************!
!********   2D- DISCONTINUOUS GALERKIN METHOD              ***********!
!*****                   SATYVIR SINGH                         *******!
!**                        2017                                    ***!
!**********************************************************************

SUBROUTINE GAS_INFORMATION
USE GAS_INFO
USE TIME_INFO
USE CONSTANTS
USE GRID_INFO
USE INPUT_INFO
USE FLOWFIELD_INFO
USE DIMENSIONLESS_INFO
USE COEFFICIENTS
IMPLICIT NONE


IF(GAS_TYPE .EQ. "ARGON") THEN
gasr=207.0      
gama=1.667          
Viscosity=2.117e-5     
Thermal_Conductivity=0.016   
Molecular_Diameter=4.17e-10        
Cp=520.0           
Bulk_Ratio=0.0            
VISCOSITY_INDEX =0.81       
     
ELSEIF(GAS_TYPE .EQ. "N2") THEN
gasr=297.0       
gama=1.4         
Viscosity=1.656e-5    
Thermal_Conductivity=0.024  
Molecular_Diameter=4.17e-10        
Cp=1040.0          
Bulk_Ratio=0.8 
VISCOSITY_INDEX =0.74
    
ELSEIF(GAS_TYPE .EQ. "CH4") THEN    
gasr=518.0            
gama=1.313          
Viscosity=1.024e-5     
Thermal_Conductivity=0.035   
Molecular_Diameter=4.83e-10       
Cp=2260.0           
Bulk_Ratio=1.33         
VISCOSITY_INDEX =0.84        
    
ELSEIF(GAS_TYPE .EQ. "H2") THEN          
gasr=4126.0         
gama=1.405         
Viscosity=0.845e-5    
Thermal_Conductivity=0.182   
Molecular_Diameter=2.92e-10        
Cp=14310.0           
Bulk_Ratio=35.0           ! 35           
VISCOSITY_INDEX =0.67     !0.67         
    
ELSEIF(GAS_TYPE .EQ. "CO2") THEN
gasr=189.0          
gama=1.289          
Viscosity=1.38e-5     
Thermal_Conductivity=0.0146    
Molecular_Diameter=5.62e-10        
Cp=2260.0          
Bulk_Ratio=1000.0            !1000.0            
VISCOSITY_INDEX =0.0   !0.93          

END IF

IF(Bulk_Ratio .EQ. 0.0)THEN
GAMA_BULK =0.0
ELSE
GAMA_BULK = (5.0 - 3.0*GAMA)/Bulk_Ratio
END IF

!========================================================
!========================================================

END SUBROUTINE                       

