planeta.x: Planeta.cpp
	g++ -g -fsanitize=address -fsanitize=undefined -fsanitize=leak $< -o $@

planeta.txt: Planeta.x
	$< > $@

orbita.png

#Vector Librry https://github.com/cdelv/Vector3D