BOOST=F:\boost\include\boost-1_66

COIN_HAS_PKGCONFIG = TRUE

LIMITEDMEMORYMCTS_SRC = FullOptimiseTMINLP.hpp FullOptimiseTMINLP.cpp LimitedMemoryMCTS.hpp
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

default: main.cpp $(LIMITEDMEMORYMCTS_SRC)
	g++ -std=c++17 -o mcts $(LIMITEDMEMORYMCTS_SRC) main.cpp -lm -I$(BOOST) -IBonmin-1.8.6/include/coin $(LIBS) -O3 -DITERATIONS=300 -DMAX_STATES=10

visual: main.cpp
	g++ -std=c++17 -o visualmcts main.cpp $(LIMITEDMEMORYMCTS_SRC) -lm -I$(BOOST) -IBonmin-1.8.6/include/coin $(LIBS) -O3 -DITERATIONS=30 -DDISPLAY_MODE -DMAX_STATES=10

dots:
	python dots.py

debug: main.cpp
	g++ -std=c++17 -o debug main.cpp $(LIMITEDMEMORYMCTS_SRC) -lm -g -I$(BOOST) -IBonmin-1.8.6/include/coin $(LIBS) -DITERATIONS=100 -DMAX_STATES=10

profile: main.cpp
	g++ -std=c++17 -o profile main.cpp $(LIMITEDMEMORYMCTS_SRC) -lm -pg -I$(BOOST) -IBonmin-1.8.6/include/coin $(LIBS) -O3 -DITERATIONS=100 -DMAX_STATES=10

test: test.cpp
	g++ -std=c++17 -Wall -o test test.cpp -I$(BOOST)
all: profile debug default

clean:
	rm mcts debug
