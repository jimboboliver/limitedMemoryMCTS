BOOST=F:\boost\include\boost-1_66\

default: main.cpp
	g++ -std=c++0x -o mcts main.cpp -lm -I$(BOOST) -O2

dots:
	python dots.py

debug: max.cpp
	g++ -std=c++0x -o debug max.cpp -lm -g -I$(BOOST) -O2

clean:
	del mcts.exe
	del debug.exe