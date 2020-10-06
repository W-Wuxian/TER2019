
set xlabel "DEFORMATION"
      set ylabel "S_YY (PA)"
   m="./SUM_ELAST"
 set terminal x11 0
 set grid
 set nokey
 set title "S_YY=f(DEFORMATION) "
plot m using 1:2 

