#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <ctime>
#include <vector>
#include "mcts.hpp"

#define ITERATIONS 1000

typedef struct {
    unsigned x:2; // 0 = nothing, 1 = X, 2 = O
} Spot;

typedef struct {
    Spot spots[9];
} Board;

struct VertexProperties {
    float wins;
    int visits;
};

Board rootBoard;

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
            get_spot(board.spots[3]) << get_spot(board.spots[4]) << get_spot(board.spots[5]) << '\n' << get_spot(board.spots[6])
            << get_spot(board.spots[7]) << get_spot(board.spots[8]) << '\n';
    return ss.str();
}

int num_possible_moves(Board board) {
    int num_moves = 0;
    for (int i = 0; i < 9; i++) {
        if (!board.spots[i].x) { // If the spot is empty it's a possible move
            num_moves++;
        }
    }
    return num_moves;
}

bool has_unborn(Board board, LimitedMemoryMCTS<VertexProperties>::vertex_t vertex, LimitedMemoryMCTS<VertexProperties>::Graph* graph) {
    return num_possible_moves(board) != boost::out_degree(vertex, (*graph));
}

int no_moves_left(Board* board) {
    return board->spots[0].x && board->spots[1].x && board->spots[2].x && board->spots[3].x && board->spots[4].x && 
            board->spots[5].x && board->spots[6].x && board->spots[7].x && board->spots[8].x;
}

/** Return whether tic-tac-toe state is gameover and if so, who won. 0 = not gameover, 1 = X, 2 = O, 3 = draw */
int terminal(Board* board) {
    if (board->spots[0].x && board->spots[0].x == board->spots[1].x && board->spots[1].x == board->spots[2].x) {
        return board->spots[0].x;
    } else if (board->spots[3].x && board->spots[3].x == board->spots[4].x && board->spots[4].x == board->spots[5].x) {
        return board->spots[3].x;
    } else if (board->spots[6].x && board->spots[6].x == board->spots[7].x && board->spots[7].x == board->spots[8].x) {
        return board->spots[6].x;
    } else if (board->spots[0].x && board->spots[0].x == board->spots[3].x && board->spots[3].x == board->spots[6].x) {
        return board->spots[0].x;
    } else if (board->spots[1].x && board->spots[1].x == board->spots[4].x && board->spots[4].x == board->spots[7].x) {
        return board->spots[1].x;
    } else if (board->spots[2].x && board->spots[2].x == board->spots[5].x && board->spots[5].x == board->spots[8].x) {
        return board->spots[2].x;
    } else if (board->spots[0].x && board->spots[0].x == board->spots[4].x && board->spots[4].x == board->spots[8].x) {
        return board->spots[0].x;
    } else if (board->spots[2].x && board->spots[2].x == board->spots[4].x && board->spots[4].x == board->spots[6].x) {
        return board->spots[2].x;
    } else if (no_moves_left(board)) {
        return 3;
    }
    return 0;
}

bool same_board(Board board1, Board board2) {
    for (int i = 0; i < 9; i++) {
        if (board1.spots[i].x != board2.spots[i].x) {
            return false;
        }
    }
    return true;
}

LimitedMemoryMCTS<VertexProperties>::vertex_t get_parent(LimitedMemoryMCTS<VertexProperties>::vertex_t vertex, LimitedMemoryMCTS<VertexProperties>::Graph* graph) {
    LimitedMemoryMCTS<VertexProperties>::in_edge_iterator in_begin, in_end;
    
    for (boost::tie(in_begin, in_end) = boost::in_edges(vertex, *graph); in_begin != in_end; ++in_begin) {
        return boost::source(*in_begin, *graph);
    }
}

/** Add child and return it. */
LimitedMemoryMCTS<VertexProperties>::vertex_t add_child(LimitedMemoryMCTS<VertexProperties>::vertex_t vertex, LimitedMemoryMCTS<VertexProperties>::Graph* graph) {
    LimitedMemoryMCTS<VertexProperties>::vertex_t child = add_vertex(VertexProperties{.wins = 0, .visits = 0}, (*graph)); // Add node and return the vertex descriptor

    add_edge(vertex, child, (*graph));

    return child;
}

Board simulate(Board* board, Spot player) {
    Board options[9];
    int numOptions = 0;

    for (int i = 0; i < 9; i++) {
        if (!board->spots[i].x) {
            options[numOptions] = *board;
            options[numOptions++].spots[i].x = player.x;
        }
    }

    return options[rand() % numOptions];
}

int get_play(Board parentBoard, LimitedMemoryMCTS<VertexProperties>::vertex_t vertex, LimitedMemoryMCTS<VertexProperties>::Graph* graph) {
    int options[9];
    int numOptions = 0;
    for (int i = 0; i < 9; i++) { // Find all the possible children of the parent
        if (!parentBoard.spots[i].x) {
            options[numOptions++] = i;
        }
    }

    // printf("%d %d\n", boost::out_degree(get_parent(vertex, graph), *graph), boost::out_degree(vertex, *graph));

    int childIndex = 0;
    LimitedMemoryMCTS<VertexProperties>::out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(get_parent(vertex, graph), *graph); ei != ei_end; ++ei) {
        if (vertex == boost::target(*ei, *graph)) { // If this is the right child, set the childIndex
            break;
        }
        childIndex++;
    }

    return options[childIndex]; // Return the child's action
}

Board make_play(Board board, LimitedMemoryMCTS<VertexProperties>::vertex_t vertex, LimitedMemoryMCTS<VertexProperties>::Graph* graph, Spot player) {
    board.spots[get_play(board, vertex, graph)].x = player.x;
    return board;
}

void backpropagate(LimitedMemoryMCTS<VertexProperties>::vertex_t vertex, LimitedMemoryMCTS<VertexProperties>::Graph* graph, int result, Spot player) {
    std::list<LimitedMemoryMCTS<VertexProperties>::vertex_t> path;
    while (boost::in_degree(vertex, *graph)) { // While we are not at the root
        path.push_front(vertex);
        if (result == player.x) {
            (*graph)[vertex].wins += 1;
        } else if (result == 3) {
            (*graph)[vertex].wins += 0.5;
        }
        (*graph)[vertex].visits++;
        player.x = player.x ^ 3; // Switch player
        vertex = get_parent(vertex, graph);
    }
    // Update the root node
    if (result == player.x) {
        (*graph)[vertex].wins += 1;
    } else if (result == 3) {
        (*graph)[vertex].wins += 0.5;
    }    
    (*graph)[vertex].visits += 1;
}

bool is_leaf(LimitedMemoryMCTS<VertexProperties>::vertex_t vertex, LimitedMemoryMCTS<VertexProperties>::Graph* graph) {
    return !boost::in_degree(vertex, *graph) || !boost::out_degree(vertex, *graph);
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

template <typename Map>
struct my_node_writer {
    my_node_writer(Map& g_) : g (g_) {};
    template <class Vertex>
    void operator()(std::ostream& out, Vertex v) {
        std::list<LimitedMemoryMCTS<VertexProperties>::vertex_t> path;
        int player = 2;
        LimitedMemoryMCTS<VertexProperties>::vertex_t vertex = v;
        Board board = rootBoard;

        while (boost::in_degree(vertex, g)) { // While we are not at the root
            path.push_front(vertex);
            vertex = get_parent(vertex, &g);
        }

        LimitedMemoryMCTS<VertexProperties>::vertex_t parent = vertex;
        for (LimitedMemoryMCTS<VertexProperties>::vertex_t testVertex : path) {
            player = player ^ 3;
            board = make_play(board, testVertex, &g, Spot{(unsigned int) player});
            parent = testVertex;
        }

        out << " [label=\"" << print_board(board) << "\"]" << std::endl;
        if (terminal(&board)) {
            out << " [color=\"" << "black" << "\"]" << std::endl;
            out << " [fontcolor=\"" << "black" << "\"]" << std::endl;
        } else {
            out << " [color=\"" << get_colour(player) << "\"]" << std::endl;
            out << " [fontcolor=\"" << get_colour(player) << "\"]" << std::endl;
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
        LimitedMemoryMCTS<VertexProperties>::vertex_t pointing_at = boost::target(e, g);
        out << " [label=\""<< g[pointing_at].wins / g[pointing_at].visits << "\"]" << std::endl;
    };
    Map g;
};

template <typename Map>
my_edge_writer<Map> edge_writer(Map& map) { 
    return my_edge_writer<Map>(map); 
}

void write_dot(LimitedMemoryMCTS<VertexProperties>::Graph* graph, int write_iteration) {
    // Make ID map
    std::map<LimitedMemoryMCTS<VertexProperties>::vertex_t, size_t> ids;

    for (auto u : boost::make_iterator_range(boost::vertices(*graph))) {
        ids[u] = ids.size();
    }

    boost::default_writer w;
    // represent graph in DOT format
    std::ofstream myfile;
    std::ostringstream ss;
    ss << "graph" << write_iteration << ".dot";
    myfile.open(ss.str().c_str());

    boost::write_graphviz(myfile, *graph, node_writer(*graph), edge_writer(*graph), w, boost::make_assoc_property_map(ids));
    myfile.close();
}

std::list<LimitedMemoryMCTS<VertexProperties>::vertex_t> delete_branch(LimitedMemoryMCTS<VertexProperties>::vertex_t vertex, LimitedMemoryMCTS<VertexProperties>::Graph* graph) { // recursively delete all descendant nodes and itself
    std::list<LimitedMemoryMCTS<VertexProperties>::vertex_t> garbage;
    LimitedMemoryMCTS<VertexProperties>::out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, (*graph)); ei != ei_end; ++ei) {
        LimitedMemoryMCTS<VertexProperties>::vertex_t target = boost::target(*ei, *graph);
        garbage.merge(delete_branch(target, graph));
    }

    garbage.push_back(vertex);

    return garbage;
}

LimitedMemoryMCTS<VertexProperties>::vertex_t root_changeover(LimitedMemoryMCTS<VertexProperties>::vertex_t chosen, LimitedMemoryMCTS<VertexProperties>::Graph* graph) {
    // Make the play the new root by first removing all irrelevant branches
    std::list<LimitedMemoryMCTS<VertexProperties>::vertex_t> garbage;
    LimitedMemoryMCTS<VertexProperties>::out_edge_iterator ei, ei_end;
    LimitedMemoryMCTS<VertexProperties>::vertex_t root = LimitedMemoryMCTS<VertexProperties>::get_root(graph);
    for (boost::tie(ei, ei_end) = boost::out_edges(root, *graph); ei != ei_end; ++ei) {
        LimitedMemoryMCTS<VertexProperties>::vertex_t target = boost::target(*ei, *graph);
        if (target != chosen) {
            garbage.merge(delete_branch(target, graph));
        }
    }

    for (LimitedMemoryMCTS<VertexProperties>::vertex_t to_remove : garbage) {
        boost::clear_vertex(to_remove, *graph);
        boost::remove_vertex(to_remove, *graph);
    }

    boost::clear_vertex(root, *graph); // Remove old root. Remove all edges otherwise undefined behaviour occurs after remove_vertex
    boost::remove_vertex(root, *graph);
}

std::list<LimitedMemoryMCTS<VertexProperties>::vertex_t> get_children(LimitedMemoryMCTS<VertexProperties>::vertex_t vertex, LimitedMemoryMCTS<VertexProperties>::Graph* graph) {
    std::list<LimitedMemoryMCTS<VertexProperties>::vertex_t> children;
    LimitedMemoryMCTS<VertexProperties>::out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, (*graph)); ei != ei_end; ++ei) {
        children.push_back(boost::target(*ei, *graph));
    }

    return children;
}

LimitedMemoryMCTS<VertexProperties>::vertex_t make_best_play(LimitedMemoryMCTS<VertexProperties>::Graph* graph) {
    LimitedMemoryMCTS<VertexProperties>::vertex_t root = LimitedMemoryMCTS<VertexProperties>::get_root(graph);
    float max_value = -1;
    float new_value;
    int best_index = 0;
    int i = 0;
    LimitedMemoryMCTS<VertexProperties>::vertex_t best;
    LimitedMemoryMCTS<VertexProperties>::out_edge_iterator ei, ei_end;

    // Calculate best play
    for (boost::tie(ei, ei_end) = boost::out_edges(root, (*graph)); ei != ei_end; ++ei) {
        LimitedMemoryMCTS<VertexProperties>::vertex_t target = boost::target(*ei, *graph);

        new_value = (*graph)[target].wins / (*graph)[target].visits;

        if (new_value > max_value) {
            max_value = new_value;
            best = target;
            best_index = i;
        }
        i++;
    }

    return best;
}

LimitedMemoryMCTS<VertexProperties>::vertex_t make_human_play(Board parentBoard, LimitedMemoryMCTS<VertexProperties>::vertex_t root, LimitedMemoryMCTS<VertexProperties>::Graph* graph, unsigned int action) {
    int options[9];
    int numOptions = 0;
    for (int i = 0; i < 9; i++) { // Find all the possible children of the parent
        if (!parentBoard.spots[i].x) {
            options[numOptions++] = i;
        }
        if (i == action) {
            break;
        }
    }

    LimitedMemoryMCTS<VertexProperties>::vertex_t target;
    int numCurrentChildren = boost::out_degree(root, *graph);
    if (numOptions <= numCurrentChildren) {
        int childIndex = 0;
        LimitedMemoryMCTS<VertexProperties>::out_edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::out_edges(root, *graph); ei != ei_end; ++ei) {
            if (childIndex++ == numOptions - 1) { // If this is the right child, set the childIndex
                target = boost::target(*ei, *graph);
                break;
            }
        }
    } else {
        for (int i = 0; i < numOptions - numCurrentChildren; i++) {
            VertexProperties childProperties{.wins = 0, .visits = 0};
            target = add_vertex(childProperties, (*graph));

            add_edge(root, target, (*graph));
        }
    }

    return target;
}

int main() {
    std::srand(time(NULL));

    int write_iteration = 0;

    LimitedMemoryMCTS<VertexProperties>::Graph graph;

    rootBoard = Board{0};

    LimitedMemoryMCTS<VertexProperties>::vertex_t vertex = boost::add_vertex(VertexProperties{.wins = 0, .visits = 0}, graph); // Add root node

    while (1) {
        int term;

         // MCTS iterations
        for (int it = 0; it < ITERATIONS; it++) {
            Board currentBoard = rootBoard;

            vertex = LimitedMemoryMCTS<VertexProperties>::get_root(&graph);
            Spot player = {1}; // Start with AI player

            // Select
            while (boost::out_degree(vertex, graph) && !has_unborn(currentBoard, vertex, &graph)) {
                vertex = select_child(vertex, &graph);
                currentBoard = make_play(currentBoard, vertex, &graph, player);
                player.x = player.x ^ 3; // Switch player
            }

            // Expand
            if (!(term = terminal(&currentBoard))) {
                if (has_unborn(currentBoard, vertex, &graph)) {
                    vertex = add_child(vertex, &graph);
                    currentBoard = make_play(currentBoard, vertex, &graph, player);
                    player.x = player.x ^ 3; // Switch player
                }

                // Simulate
                Board simBoard = currentBoard;
                while (!(term = terminal(&simBoard))) {
                    simBoard = simulate(&simBoard, player);
                    player.x = player.x ^ 3; // Switch player
                }
            }

            // Backpropagate
            player.x = player.x ^ 3; // Switch player
            backpropagate(vertex, &graph, term, player);
        }
        
        write_dot(&graph, write_iteration++);

        vertex = make_best_play(&graph);
        rootBoard = make_play(rootBoard, vertex, &graph, Spot{1});
        root_changeover(vertex, &graph);

        std::cout << print_board(rootBoard);
        std::cout.flush();

        term = terminal(&rootBoard);

        if (term == 1) {
            printf("\nGame over! AI wins\n");
            break;
        } else if (term == 3) {
            printf("\nGame over! Draw\n");
            break;
        }

        printf("\nEnter play: ");

        int humanPlay;
        int ch;
        while (1) {
            humanPlay = std::cin.get() - 48;
            while ((ch = std::cin.get()) != '\n' && ch != EOF);
            printf("\n");
            if (rootBoard.spots[humanPlay].x || humanPlay > 8 || humanPlay < 0) {
                printf("Cannot make play. Enter play: ");
            } else {
                break;
            }
        }

        vertex = make_human_play(rootBoard, vertex, &graph, humanPlay);
        rootBoard = make_play(rootBoard, vertex, &graph, Spot{2});
        root_changeover(vertex, &graph);

        std::cout << print_board(rootBoard);
        std::cout.flush();

        printf("\n");

        term = terminal(&rootBoard);

        if (term == 2) {
            printf("\nGame over! Player wins\n");
            break;
        } else if (term == 3) {
            printf("\nGame over! Draw\n");
            break;
        }
    }

    return 0;  
}