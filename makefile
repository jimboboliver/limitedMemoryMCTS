BOOST=F:\boost\include\boost-1_66
LP_SOLVE_INCLUDE=lp_solve
LP_SOLVE_LIB=lp_solve

COIN_HAS_PKGCONFIG = TRUE

# Linker flags
ifeq ($(COIN_HAS_PKGCONFIG), TRUE)
  LIBS = `PKG_CONFIG_PATH=/home/james/limitedMemoryMCTS/Bonmin-1.8.6/lib64/pkgconfig:/home/james/limitedMemoryMCTS/Bonmin-1.8.6/lib/pkgconfig:/home/james/limitedMemoryMCTS/Bonmin-1.8.6/share/pkgconfig: pkg-config --libs bonmin`
else
  ifeq ($(COIN_CXX_IS_CL), TRUE)
    LIBS = -link -libpath:`$(CYGPATH_W) /home/james/limitedMemoryMCTS/Bonmin-1.8.6/lib` libbonmin.lib 
  else
    LIBS = -L/home/james/limitedMemoryMCTS/Bonmin-1.8.6/lib -lbonmin 
  endif
endif

default: main.cpp
	g++ -std=c++17 -o mcts main.cpp -lm -I$(BOOST) -I$(LP_SOLVE_INCLUDE) -L$(LP_SOLVE_LIB) -llpsolve55 -IBonmin-1.8.6/include/coin $(LIBS) -O3

dots:
	python dots.py

debug: main.cpp
	g++ -std=c++17 -o debug main.cpp -lm -g -I$(BOOST) -I$(LP_SOLVE_INCLUDE) -L$(LP_SOLVE_LIB) -llpsolve55 -IBonmin-1.8.6/include/coin $(LIBS)

profile: main.cpp
	g++ -std=c++17 -o profile main.cpp -lm -pg -I$(BOOST) -I$(LP_SOLVE_INCLUDE) -L$(LP_SOLVE_LIB) -llpsolve55 -IBonmin-1.8.6/include/coin $(LIBS) -O3

test: test.cpp
	g++ -std=c++17 -Wall -o test test.cpp -I$(BOOST) -I$(LP_SOLVE_INCLUDE) -L$(LP_SOLVE_LIB) -llpsolve55

all: profile debug default

clean:
	rm mcts debug
