 set xlabel "#iteration"
      set ylabel "||r||"
         set logscale y 2
   m="./DELLres.txt"
 set terminal x11 0
 set grid
 set nokey
 set title "convergence du GC pre-cond"
plot m using 1:2 with lines
