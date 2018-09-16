BOOST=F:\boost\include\boost-1_66\

default: main.cpp
	g++ -std=c++0x -Wall -o mcts main.cpp -lm -I$(BOOST)

dots:
	python dots.py

debug: main.cpp
	g++ -std=c++0x -Wall -o debug main.cpp -lm -g -I$(BOOST)

clean:
	del mcts.exe
	del debug.exe