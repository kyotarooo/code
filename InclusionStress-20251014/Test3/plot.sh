plot "log" u (($1-0.0000005)/0.0000001):(-$2/(80000000000*0.001)) w lp title "Srr"
replot "log" u (($1-0.0000005)/0.0000001):(-$3/(80000000000*0.001)) w lp title "Spp"
replot "log" u (($1-0.0000005)/0.0000001):(-$4/(80000000000*0.001)) w lp title "Szz"
replot "log" u (($1-0.0000005)/0.0000001):(-$7/(80000000000*0.001)) w lp title "Srz"
set xrange [0:4]
set yrange [-1:3]
set xtics 1
set grid
rep
pause -1
