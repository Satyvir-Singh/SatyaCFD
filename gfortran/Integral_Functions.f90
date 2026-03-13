  
!****************************************************************
!*******integral functions
!****************************************************************  

      FUNCTION rint_edge(fv_1,fv0,fv1)
      IMPLICIT NONE
      DOUBLE PRECISION :: rint_edge
      DOUBLE PRECISION :: fv_1, fv0, fv1
       rint_edge=5./9.*fv_1+8./9.*fv0+5./9.*fv1 
      RETURN 
      END FUNCTION

      FUNCTION rint_area(f_1_1,f_11,f1_1,f11, f0_1,f01,f_10,f10, f00)
      IMPLICIT NONE
      DOUBLE PRECISION :: rint_area
      DOUBLE PRECISION :: f_1_1,f_11,f1_1,f11, f0_1,f01,f_10,f10, f00
       rint_area= 25./81.*(f_1_1+f_11+f1_1+f11)+ 40./81.*(f0_1+f01+f_10+f10)+ 64./81.*f00 
      RETURN 
      END FUNCTION