BOOST=F:\boost\include\boost-1_66\

default: main.cpp
	g++ -std=c++0x -o mcts main.cpp -lm -I$(BOOST)

dot:
	dot graph0.dot -Tjpg -o graph0.jpg

debug: main.cpp
	g++ -std=c++0x -o debug main.cpp -lm -g -I$(BOOST)

clean:
	del mcts.exe
	del debug.exe