ex_euler.x: Pendulo_Euler.cpp
	g++ -g -fsanitize=address -fsanitize=undefined -fsanitize=leak $< -o $@

ex_verlet.x: Pendulo_omegaerlet.cpp
	g++ -g -fsanitize=address -fsanitize=undefined -fsanitize=leak $< -o $@

data_euler.txt: ex_euler.x
	./$< > $@

data_verlet.txt: ex_verlet.x
	./$< > $@

graph_euler: graph.gp data_euler.txt
	mv data_euler.txt data.txt;\
	gnuplot $<

graph_verlet: graph.gp data_verlet.txt
	mv data_verlet.txt data.txt;\
	gnuplot $<

plot: graph.png
	eog $< &

clean:
	rm -f *.x *.txt *.out