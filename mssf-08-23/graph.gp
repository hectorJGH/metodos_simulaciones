set term png
set out "graph.png"

tau(x) = sqrt(25-x*x)

plot "data.txt" w lp, tau(x) w l