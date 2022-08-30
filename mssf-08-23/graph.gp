set term png
set out "graph.png"

F(x) = sqrt(25-x*x)

plot "data.txt" w lp, F(x) w l