BOOST=F:\boost\include\boost-1_66

default: main.cpp
	g++ -std=c++0x -Wall -o mcts main.cpp -lm -I$(BOOST) -O2

dots:
	python dots.py

debug: main.cpp
	g++ -std=c++0x -Wall -o debug main.cpp -lm -g -I$(BOOST) -O2

profile: main.cpp
	g++ -std=c++0x -Wall -o profile main.cpp -lm -pg -I$(BOOST) -O2

all: profile debug default

clean:
	rm mcts debug