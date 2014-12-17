c ** LMTREORT -> maximum number of reortoganallizations in lanczos aproach
c ** MXCST -> maximum potential constants
c ** MXDCT - > maximum dicretization points 
c ** MXDM -> maximum dimension (
c ** MXAUX - This varibable is determined by LMTREORT and MXDCT
c ** ([167772,400]; [223696,300])
      INTEGER       LMTREORT, MXCST, MXDCT, MXDM, MXAUX
      PARAMETER     (
     &     LMTREORT = 200, 
     &     MXCST = 14, 
     &     MXDCT = 85000,
     &     MXDM = 3, 
     &     MXAUX = LMTREORT*MXDCT 
     &     ) 

     
 
