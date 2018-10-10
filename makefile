BOOST=F:\boost\include\boost-1_66
SYMPHONY=SYMPHONY-5.6/SYMPHONY/include
OSI=SYMPHONY-5.6/include/coin

default: main.cpp
	export LD_LIBRARY_PATH=/home/james/limitedMemoryMCTS/SYMPHONY-5.6/lib
	g++ -std=c++0x -Wall -o mcts main.cpp -lm -I$(BOOST) -I$(SYMPHONY) -I$(OSI) -L/home/james/limitedMemoryMCTS/SYMPHONY-5.6/lib -lOsiSym -lSym 

dots:
	python dots.py

debug: main.cpp
	g++ -std=c++0x -Wall -o debug main.cpp -lm -g -I$(BOOST) -I$(SYMPHONY) -I$(OSI)

profile: main.cpp
	g++ -std=c++0x -Wall -o profile main.cpp -lm -pg -I$(BOOST) -O2

all: profile debug default

clean:
	rm mcts debug