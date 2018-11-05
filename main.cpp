#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <ctime>
#include <vector>
#include <float.h>

#include <iomanip>
#include <fstream>

#include "LimitedMemoryMCTS.hpp"

/* Game stuff */

typedef struct {
    unsigned x:2; // 0 = nothing, 1 = X, 2 = O
} Spot;

typedef struct {
    Spot spots[9];
} Board;

class VertexProperties {
    public:
        VertexProperties(Board* state_stored, Spot last_player);
        VertexProperties();

        bool has_state;
        Spot player;

    public: // Needed for MCTS interface
        float wins;
        int visits;
        Board* state;
        unsigned int possible_children;
        int terminal;

        void determine_terminal();
        VertexProperties generate_child(int playNum);
        void regenerate_child(int playNum, VertexProperties* child);
        int playout();
        bool equals(VertexProperties otherVP);

    private:
        void num_possible_moves();
};

VertexProperties::VertexProperties() {
    state = 0;
    player = Spot{2};
    has_state = false;
    wins = 0;
    visits = 0;
    possible_children = 0;
    terminal = 0;
}

VertexProperties::VertexProperties(Board* state_stored, Spot last_player) {
    state = state_stored;
    player = last_player;
    has_state = true;
    wins = 0;
    visits = 0;
    possible_children = 0;
    num_possible_moves();
    determine_terminal();
}

void VertexProperties::num_possible_moves() {
    for (int i = 0; i < 9; i++) {
        if (!state->spots[i].x) { // If the spot is empty it's a possible move
            possible_children++;
        }
    }
}

VertexProperties VertexProperties::generate_child(int playNum) {
    int plays[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    int available_plays[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    for (int j = 0; j < 9; j++) {
        if (state->spots[j].x) {
            plays[j] = 1;
        }
    }

    int num_plays = 0;
    for (int i = 0; i < 9; i++) {
        if (!plays[i]) {
            available_plays[num_plays++] = i;
        }
    }
    int chosen_play = available_plays[playNum];

    Board* new_state = new Board;
    *new_state = *state;

    if (player.x == 1) {
        new_state->spots[chosen_play].x = 2;
        return VertexProperties(new_state, Spot{2});
    } else {
        new_state->spots[chosen_play].x = 1;
        return VertexProperties(new_state, Spot{1});
    }
}

void VertexProperties::regenerate_child(int playNum, VertexProperties* child) {
    int plays[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    int available_plays[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    for (int j = 0; j < 9; j++) {
        if (state->spots[j].x) {
            plays[j] = 1;
        }
    }

    int num_plays = 0;
    for (int i = 0; i < 9; i++) {
        if (!plays[i]) {
            available_plays[num_plays++] = i;
        }
    }
    int chosen_play = available_plays[playNum];

    Board* new_state = new Board;
    *new_state = *state;

    if (player.x == 1) {
        new_state->spots[chosen_play].x = 2;
    } else {
        new_state->spots[chosen_play].x = 1;
    }
    child->state = new_state;
    child->has_state = true;
}

bool VertexProperties::equals(VertexProperties otherVP) {
    return state->spots[0].x == otherVP.state->spots[0].x &&
           state->spots[1].x == otherVP.state->spots[1].x &&
           state->spots[2].x == otherVP.state->spots[2].x &&
           state->spots[3].x == otherVP.state->spots[3].x &&
           state->spots[4].x == otherVP.state->spots[4].x &&
           state->spots[5].x == otherVP.state->spots[5].x &&
           state->spots[6].x == otherVP.state->spots[6].x &&
           state->spots[7].x == otherVP.state->spots[7].x &&
           state->spots[8].x == otherVP.state->spots[8].x;
}

Board simulate(Board board, Spot player) {
    Board options[9];
    int numOptions = 0;

    for (int i = 0; i < 9; i++) {
        if (!board.spots[i].x) {
            options[numOptions] = board;
            options[numOptions++].spots[i].x = player.x;
        }
    }

    return options[rand() % numOptions];
}

int no_moves_left(Board board) {
    return board.spots[0].x && board.spots[1].x && board.spots[2].x && board.spots[3].x && board.spots[4].x && 
            board.spots[5].x && board.spots[6].x && board.spots[7].x && board.spots[8].x;
}

/** Return whether tic-tac-toe state is gameover and if so, who won. 0 = not gameover, 1 = X, 2 = O, 3 = draw */
int terminal_board(Board state) {
    if (state.spots[0].x && state.spots[0].x == state.spots[1].x && state.spots[1].x == state.spots[2].x) {
        return state.spots[0].x;
    } else if (state.spots[3].x && state.spots[3].x == state.spots[4].x && state.spots[4].x == state.spots[5].x) {
        return state.spots[3].x;
    } else if (state.spots[6].x && state.spots[6].x == state.spots[7].x && state.spots[7].x == state.spots[8].x) {
        return state.spots[6].x;
    } else if (state.spots[0].x && state.spots[0].x == state.spots[3].x && state.spots[3].x == state.spots[6].x) {
        return state.spots[0].x;
    } else if (state.spots[1].x && state.spots[1].x == state.spots[4].x && state.spots[4].x == state.spots[7].x) {
        return state.spots[1].x;
    } else if (state.spots[2].x && state.spots[2].x == state.spots[5].x && state.spots[5].x == state.spots[8].x) {
        return state.spots[2].x;
    } else if (state.spots[0].x && state.spots[0].x == state.spots[4].x && state.spots[4].x == state.spots[8].x) {
        return state.spots[0].x;
    } else if (state.spots[2].x && state.spots[2].x == state.spots[4].x && state.spots[4].x == state.spots[6].x) {
        return state.spots[2].x;
    } else if (no_moves_left(state)) {
        return 3;
    }
    return 0;
}

int VertexProperties::playout() {
    int term;
    Board currentBoard = *state;
    Spot currentPlayer = player;
    currentPlayer.x ^= 3;
    
    while (!(term = terminal_board(currentBoard))) {
        currentBoard = simulate(currentBoard, currentPlayer);
        currentPlayer.x ^= 3;
    }

    return term;
}

/** Return whether tic-tac-toe state is gameover and if so, who won. 0 = not gameover, 1 = AI, 2 = player, 3 = draw */
void VertexProperties::determine_terminal() {
    terminal = terminal_board(*state);
}

char get_spot(Spot spot) {
    if (spot.x == 1) {
        return 'X';
    } else if (spot.x == 2) {
        return 'O';
    } else {
        return '*';
    }
}

std::string print_board(Board board) {
    std::ostringstream ss;
    ss << get_spot(board.spots[0]) << get_spot(board.spots[1]) << get_spot(board.spots[2]) << '\n' << 
          get_spot(board.spots[3]) << get_spot(board.spots[4]) << get_spot(board.spots[5]) << '\n' << 
          get_spot(board.spots[6]) << get_spot(board.spots[7]) << get_spot(board.spots[8]) << '\n';
    return ss.str();
}

Spot who_played(Board state) {
    int num_played = 0;
    for (int i = 0; i < 9; i++) {
        if (state.spots[i].x) {
            num_played++;
        }
    }

    if (num_played % 2 == 1) {
        return Spot{1};
    } else if (num_played == 0) {
        return Spot{0};
    }
    return Spot{2};
}

std::string get_colour(int who) {
    if (who == 1) {
        return "green";
    } else if (who == 2) {
        return "blue";
    }
    return "black";
}

typedef LimitedMemoryMCTS<VertexProperties> MCTS;

MCTS::vertex_t lastAdded;

template <typename Map>
struct my_node_writer {
    my_node_writer(Map& g_) : g (g_) {};
    template <class Vertex>
    void operator()(std::ostream& out, Vertex v) {
        MCTS::vertex_t vertex = v;
        out << " [label=\"" << g[vertex].visits << "\"]" << std::endl;
        Board board = *(g[vertex].state);
        if (g[vertex].has_state) {
            out << " [label=\"" << print_board(board) << "\"]" << std::endl;
            if (vertex == lastAdded) {
                out << " [color=\"" << "purple" << "\"]" << std::endl;
            } else {
                if (terminal_board(board)) {
                    out << " [color=\"" << "black" << "\"]" << std::endl;
                } else {
                    out << " [color=\"" << get_colour(who_played(board).x) << "\"]" << std::endl;
                }
            }
        } else {
            if (vertex == lastAdded) {
                out << " [color=\"" << "purple" << "\"]" << std::endl;
            } else {
                out << " [color=\"" << "red" << "\"]" << std::endl;
            }
        }
    };
    Map g;
};

template <typename Map>
my_node_writer<Map> node_writer(Map& map) { 
    return my_node_writer<Map>(map); 
}

template <typename Map>
struct my_edge_writer {
    my_edge_writer(Map& g_) : g (g_) {};
    template <class Edge>
    void operator()(std::ostream& out, Edge e) {
        // MCTS::vertex_t vertex = boost::source(e, g);
        MCTS::vertex_t pointing_at = boost::target(e, g);
        // out << " [label=\""<< g[pointing_at].wins / g[pointing_at].visits + 0.5 * sqrt(log(g[vertex].visits) / g[pointing_at].visits) << '\n' << "\"]" << std::endl;
        out << " [label=\""<< g[pointing_at].wins << '\n' << g[pointing_at].visits << "\"]" << std::endl;
    };
    Map g;
};

template <typename Map>
my_edge_writer<Map> edge_writer(Map& map) { 
    return my_edge_writer<Map>(map); 
}

void write_dot(MCTS::Graph* graph, int write_iteration) {
    // Make ID map
    std::map<MCTS::vertex_t, size_t> ids;

    for (auto u : boost::make_iterator_range(boost::vertices(*graph))) {
        ids[u] = ids.size();
    }

    boost::default_writer w;
    // represent graph in DOT format
    std::ofstream myfile;
    std::ostringstream ss;
    ss << "graph/graph" << write_iteration << ".dot";
    myfile.open(ss.str().c_str());

    boost::write_graphviz(myfile, *graph, node_writer(*graph), edge_writer(*graph), w, boost::make_assoc_property_map(ids));
    myfile.close();
}

int main(int argc, char** argv) {
    std::srand(time(NULL));

    MCTS mcts = MCTS(VertexProperties(new Board{0}, Spot{2}));

    MCTS::vertex_t vertex;

    int ITERATIONS = std::stoi(argv[1]);

    while (1) {
        int term;

         // MCTS iterations
        for (int it = 0; it < ITERATIONS; it++) {
            #ifndef DISPLAY_MODE
            std::cout << it << '\n';
            #endif

            vertex = mcts.get_root();

            // Select
            vertex = mcts.select(vertex);

            // If not terminal
            if (!(term = mcts.terminal(vertex))) {
                // Expand
                if (mcts.has_unborn(vertex)) {
                    if (!mcts.has_state(vertex)) {
                        mcts.regenerate(vertex);
                    }
                    vertex = mcts.add_child(vertex);
                    lastAdded = vertex;
                }

                // Simulate
                term = mcts.simulate(vertex);
            }

            // Backpropagate
            mcts.backpropagate(vertex, term);

            #ifdef DISPLAY_MODE
            int nums[2];
            mcts.get_num_regenerated(nums);
            std::cout << mcts.get_num_states() << ' ' << nums[0] << ' ' << nums[1] << '\n';
            #endif
            // if (it % 50 == 0) {
            // mcts.optimise_states(false);
            // } else {
                mcts.optimise_states(true); // mini optimise
            // }
            #ifdef DISPLAY_MODE
            write_dot(mcts.get_graph(), it);
            mcts.reset_num_regenerated();
            while(std::getchar() != '\n');
            #endif
        }
        
        vertex = mcts.make_best_play();
        term = mcts.terminal(vertex);
        Board* rootBoard = new Board;
        *rootBoard = *(mcts.get_vertex_properties(vertex).state);
        mcts.root_changeover(vertex);

        std::cout << print_board(*rootBoard);
        std::cout.flush();

        term = mcts.terminal(vertex);

        if (term == 1) {
            std::cout << "\nGame over!\nAI wins\n";
            delete rootBoard;
            break;
        } else if (term == 3) {
            std::cout << "\nGame over!\nDraw\n";
            delete rootBoard;
            break;
        }

        std::cout << "\nEnter play: ";

        int humanPlay;
        int ch;
        while (1) {
            humanPlay = std::cin.get() - 48;
            while ((ch = std::cin.get()) != '\n' && ch != EOF);
            std::cout << '\n';
            if (humanPlay > 8 || humanPlay < 0 || rootBoard->spots[humanPlay].x) {
                std::cout << "Cannot make play. Enter play: ";
            } else {
                break;
            }
        }

        rootBoard->spots[humanPlay].x = 2;

        std::cout << print_board(*rootBoard) << '\n';
        std::cout.flush();

        term = terminal_board(*rootBoard);

        vertex = mcts.make_play(vertex, VertexProperties(rootBoard, Spot{2}));
        mcts.root_changeover(vertex);

        if (term == 2) {
            std::cout << "\nGame over!\nPlayer wins\n";
            delete rootBoard;
            break;
        } else if (term == 3) {
            std::cout << "\nGame over!\nDraw\n";
            delete rootBoard;
            break;
        }
    }

    int nums[2];
    mcts.get_num_regenerated(nums);
    std::cout << "Total states regenerated by this implementation: " << nums[0] << '\n' << "Total states regenerated by reference implementation: " <<  nums[1] << '\n';

    return 0;
}