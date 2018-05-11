#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <ctime>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#define ITERATIONS 1000

typedef struct {
    unsigned x:2; // 0 = nothing, 1 = X, 2 = O
} Spot;

typedef struct {
    Spot spots[9];
} Board;

struct VertexProperties {
    bool has_state;
    Board board;
    float wins;
    int visits;
};

/* define the graph type
        listS: selects the STL list container to store 
                the OutEdge list
        listS: selects the STL vector container to store 
            the vertices
        directedS: selects directed edges
*/
typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, VertexProperties> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
typedef boost::graph_traits<Graph>::edge_descriptor edge_t;
typedef boost::graph_traits <Graph>::edge_iterator edge_iterator;
typedef boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
typedef boost::graph_traits<Graph>::in_edge_iterator in_edge_iterator;
typedef std::map<vertex_t, size_t> IndexMap;

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

vertex_t get_root(Graph* graph) {
    std::pair<vertex_iterator,vertex_iterator> it = boost::vertices(*graph);
    return *it.first;
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

bool has_unborn(vertex_t vertex, Graph* graph) {
    return num_possible_moves((*graph)[vertex].board) != boost::out_degree(vertex, (*graph));
}

vertex_t select_child(vertex_t vertex, Graph* graph) {
    float max_value = 0;
    float new_value;
    vertex_t best;
    out_edge_iterator ei, ei_end;

    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, (*graph)); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, *graph);

        new_value = (*graph)[target].wins / (*graph)[target].visits + sqrt(2 * log((*graph)[vertex].visits) / (*graph)[target].visits);

        if (new_value > max_value) {
            max_value = new_value;
            best = target;
        }
    }

    return best;
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

bool child_does_not_exist(vertex_t vertex, Graph* graph, Board board) {
    out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, (*graph)); ei != ei_end; ++ei) {
        if (same_board((*graph)[boost::target(*ei, *graph)].board, board)) {
            return false;
        }
    }
    return true;
}

bool should_forget(vertex_t vertex, Graph* graph) {
    return get_root(graph) != vertex && boost::out_degree(vertex, *graph) > 1;
}

/** Add a random available move from this leaf as a child and return it to simulate from. */
vertex_t add_child(vertex_t vertex, Graph* graph, Spot player) {
    Board options[9];
    int numOptions = 0;

    Board new_board;
    Board current = (*graph)[vertex].board;

    for (int i = 0; i < 9; i++) {
        if (current.spots[i].x) {
            continue;
        }
        new_board = (*graph)[vertex].board;
        new_board.spots[i] = player;
        if (!(*graph)[vertex].board.spots[i].x && child_does_not_exist(vertex, graph, new_board)) {
            options[numOptions++] = new_board;
        }
    }

    Board chosen = options[rand() % numOptions];

    vertex_t child = add_vertex(VertexProperties{.has_state = true, .board = chosen, .wins = 0, .visits = 0}, (*graph)); // Add node and return the vertex descriptor

    add_edge(vertex, child, (*graph));

    // if not the root or a leaf and all children are added, mark state as forgotten
    if (should_forget(vertex, graph)) {
        (*graph)[vertex].has_state = false;
    }

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

vertex_t get_parent(vertex_t vertex, Graph* graph) {
    in_edge_iterator in_begin, in_end;
    
    for (boost::tie(in_begin, in_end) = boost::in_edges(vertex, *graph); in_begin != in_end; ++in_begin) {   
        return boost::source(*in_begin, *graph);
    }
}

void backpropagate(vertex_t vertex, Graph* graph, int result, Spot player) {
    while (boost::in_degree(vertex, *graph)) { // While we are not at the root
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

bool is_leaf(vertex_t vertex, Graph* graph) {
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

std::string get_colour(bool has_state, vertex_t vertex, Graph* graph) {
    if (has_state) {
        int who = who_played((*graph)[vertex].board).x;
        if (who == 1) {
            return "green";
        } else if (who == 2) {
            return "blue";
        }
        return "black";
    } else {
        return "red";
    }
}

template <typename Map>
struct my_node_writer {
    my_node_writer(Map& g_) : g (g_) {};
    template <class Vertex>
    void operator()(std::ostream& out, Vertex v) {
        out << " [label=\"" << print_board(g[v].board) << "\"]" << std::endl;
        out << " [color=\"" << get_colour(g[v].has_state, v, &g) << "\"]" << std::endl;
        out << " [fontcolor=\"" << get_colour(g[v].has_state, v, &g) << "\"]" << std::endl;
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
        vertex_t pointing_at = boost::target(e, g);
        out << " [label=\"" << g[pointing_at].wins / g[pointing_at].visits << "\"]" << std::endl;
    };
    Map g;
};

template <typename Map>
my_edge_writer<Map> edge_writer(Map& map) { 
    return my_edge_writer<Map>(map); 
}

void write_dot(Graph* graph, int write_iteration) {
    // Make ID map
    std::map<vertex_t, size_t> ids;

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

std::list<vertex_t> delete_branch(vertex_t vertex, Graph* graph) { // recursively delete all descendant nodes and itself
    std::list<vertex_t> garbage;
    out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, (*graph)); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, *graph);
        garbage.merge(delete_branch(target, graph));
    }

    garbage.push_back(vertex);

    return garbage;
}

vertex_t root_changeover(vertex_t chosen, Graph* graph) {
    // Make the play the new root by first removing all irrelevant branches
    std::list<vertex_t> garbage;
    out_edge_iterator ei, ei_end;
    vertex_t root = get_root(graph);
    for (boost::tie(ei, ei_end) = boost::out_edges(root, *graph); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, *graph);
        if (target != chosen) {
            garbage.merge(delete_branch(target, graph));
        }
    }

    for (vertex_t to_remove : garbage) {
        boost::clear_vertex(to_remove, *graph);
        boost::remove_vertex(to_remove, *graph);
    }

    boost::clear_vertex(root, *graph); // Remove old root. Remove all edges otherwise undefined behaviour occurs after remove_vertex
    boost::remove_vertex(root, *graph);
}

std::list<vertex_t> get_leaves(vertex_t vertex, Graph* graph) {
    std::list<vertex_t> leaves;
    out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, (*graph)); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, *graph);
        if (boost::out_degree(target, *graph) == 0) { // If the target is a leaf, just add it to the vector of leaves.
            leaves.push_back(target);
        } else { // Else search for leaves down the target
            leaves.merge(delete_branch(target, graph));
        }
    }

    return leaves;
}

void regenerate_state(vertex_t root, vertex_t vertex, Graph* graph) {
    std::list<vertex_t> leaves = get_leaves(vertex, graph);
    Board board{0};

    if (leaves.size() > 1) {
        board = (*graph)[leaves.front()].board;
        leaves.pop_front();
        for (vertex_t leaf : leaves) {
            for (int i = 0; i < 9; i++) {
                if (board.spots[i].x != (*graph)[leaf].board.spots[i].x) {
                    board.spots[i].x = 0;
                }
            }
        }
    } else {
        board = (*graph)[root].board;
        for (int i = 0; i < 9; i++) {
            if (board.spots[i].x == 0 && (*graph)[leaves.front()].board.spots[i].x == 1) {
                board.spots[i].x = 1;
                break;
            }
        }   
    }

    if (!same_board(board, (*graph)[vertex].board)) {
        printf("Regenerate failed.\n");
        
        std::cout << print_board((*graph)[vertex].board);
        std::cout.flush();
        std::cout << print_board(board);
        std::cout.flush();
        exit(1);
    }

    (*graph)[vertex].has_state = true;
    // (*graph)[vertex].board = board;
}

vertex_t make_best_play(Graph* graph) {
    vertex_t root = get_root(graph);
    float max_value = -1;
    float new_value;
    vertex_t best;
    out_edge_iterator ei, ei_end;

    // Calculate best play
    for (boost::tie(ei, ei_end) = boost::out_edges(root, (*graph)); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, *graph);

        new_value = (*graph)[target].wins / (*graph)[target].visits;

        if (new_value > max_value) {
            max_value = new_value;
            best = target;
        }
    }

    // Regenerate the state if it has been forgotten
    if (!(*graph)[best].has_state) {
        regenerate_state(root, best, graph);
    }

    root_changeover(best, graph);

    return best;
}

vertex_t make_human_play(vertex_t root, Graph* graph, Board board) {
    vertex_t target;
    out_edge_iterator ei, ei_end;
    bool found = false;
    Board test;
    for (boost::tie(ei, ei_end) = boost::out_edges(root, (*graph)); ei != ei_end; ++ei) {
        target = boost::target(*ei, *graph);
        test = (*graph)[target].board;
        if (same_board((*graph)[target].board, board)) {
            found = true;
            break;
        }
    }

    if (!found) {
        VertexProperties childProperties{.has_state = true, .board = board, .wins = 0, .visits = 0};
        target = add_vertex(childProperties, (*graph));

        add_edge(root, target, (*graph));
    }

    (*graph)[target].has_state = true;

    root_changeover(target, graph);

    return target;
}

int main() {
    std::srand(std::time(NULL));

    int write_iteration = 0;

    Graph graph;

    VertexProperties rootProperties{.has_state = true, .board = Board{0}, .wins = 0, .visits = 0};

    vertex_t vertex = boost::add_vertex(rootProperties, graph); // Add root node

    while (1) {
        int term;

         // MCTS iterations
        for (int it = 0; it < ITERATIONS; it++) {
            // printf("Iteration %d\n", it);

            vertex = get_root(&graph);
            Spot player = {1}; // Start with AI player

            // Select
            while (out_degree(vertex, graph) && !has_unborn(vertex, &graph)) {
                vertex = select_child(vertex, &graph);
                player.x = player.x ^ 3; // Switch player
            }

            // Expand
            if (!(term = terminal(&graph[vertex].board))) {
                if (has_unborn(vertex, &graph)) {
                    vertex = add_child(vertex, &graph, player);
                    player.x = player.x ^ 3; // Switch player
                }

                // Simulate
                Board simBoard = graph[vertex].board;
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

        std::cout << print_board(graph[vertex].board);
        std::cout.flush();

        term = terminal(&graph[vertex].board);

        if (term == 1) {
            printf("Game over! AI wins\n");
            break;
        } else if (term == 3) {
            printf("Game over! Draw\n");
            break;
        }

        printf("Enter play: ");

        int humanPlay;
        int ch;
        Board human_board;
        while (1) {
            humanPlay = std::cin.get() - 48;
            while ((ch = std::cin.get()) != '\n' && ch != EOF);
            printf("\n");
            human_board = graph[vertex].board;
            if (human_board.spots[humanPlay].x || humanPlay > 8 || humanPlay < 0) {
                printf("Cannot make play. Enter play: ");
            } else {
                break;
            }
        }

        human_board.spots[humanPlay].x = 2;
        vertex = make_human_play(vertex, &graph, human_board);

        std::cout << print_board(graph[vertex].board);
        std::cout.flush();

        printf("\n");

        term = terminal(&graph[vertex].board);

        if (term == 2) {
            printf("Game over! Player wins\n");
            break;
        } else if (term == 3) {
            printf("Game over! Draw\n");
            break;
        }
    }

    return 0;  
}