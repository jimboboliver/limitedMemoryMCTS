BOOST=F:\boost\include\boost-1_66
LP_SOLVE_INCLUDE=lp_solve
LP_SOLVE_LIB=lp_solve

default: main.cpp
	g++ -std=c++0x -Wall -o mcts main.cpp -lm -I$(BOOST) -I$(LP_SOLVE_INCLUDE) -L$(LP_SOLVE_LIB) -llpsolve55 -O2

dots:
	python dots.py

debug: main.cpp
	g++ -std=c++0x -Wall -o debug main.cpp -lm -g -I$(BOOST) -I$(LP_SOLVE_INCLUDE) -L$(LP_SOLVE_LIB) -llpsolve55

profile: main.cpp
	g++ -std=c++0x -Wall -o profile main.cpp -lm -pg -I$(BOOST) -I$(LP_SOLVE_INCLUDE) -L$(LP_SOLVE_LIB) -llpsolve55 -O2

test: test.cpp
	g++ -std=c++0x -Wall -o test test.cpp -I$(BOOST) -I$(LP_SOLVE_INCLUDE) -L$(LP_SOLVE_LIB) -llpsolve55

all: profile debug default

clean:
	rm mcts debug
