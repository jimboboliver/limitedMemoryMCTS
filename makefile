BOOST=F:\boost\include\boost-1_66\

default: mcts.cpp
	g++ -std=c++0x -o mcts mcts.cpp -lm -I$(BOOST) -O2

main: main.cpp mcts.hpp mcts.cpp
	g++ -std=c++0x -o mcts main.cpp mcts.hpp mcts.cpp -lm -I$(BOOST) -O2

dots:
	python dots.py

# debug: main.cpp
# 	g++ -std=c++0x -o debug main.cpp -lm -g -I$(BOOST) -O2

debug: max.cpp
	g++ -std=c++0x -o debug max.cpp -lm -g -I$(BOOST) -O2

clean:
	del mcts.exe
	del debug.exe