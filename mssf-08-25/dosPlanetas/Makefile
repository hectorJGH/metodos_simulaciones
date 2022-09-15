all: clean plot

plot: graph.png
	eog $< &

graph.png: graph.gp data.txt
	gnuplot $<

data.txt: a.out
	./$< > $@

clean:
	rm -f *.x *data.txt *.png

#En construccion
compile:
	@read -p "Enter cpp file:" filec;\
	filecpp=$$filec;\
	g++ -g -fsanitize=address -fsanitize=undefined -fsanitize=leak $(filec)