all: clean plot

plot: graph.png
	eog $< &

data.txt: a.out
	./$< > $@

graph.png: graph.gp data.txt
	gnuplot $<

clean:
	rm -f *.x *data.txt *.png