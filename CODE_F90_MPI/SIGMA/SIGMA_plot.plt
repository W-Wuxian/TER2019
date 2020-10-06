set xlabel "DISTANCE PAR RAPPORT AU BORDS DU TROU (m)"
set ylabel "S_YY (PA)"
m="./DELLTRACE_SIGMA_YY.txt"
set terminal x11 0
set grid
set nokey
set xtics 0.001
set ytics 100000000
set title "EVO DE S_YY "
plot m using 1:2 
