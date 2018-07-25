#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include "game.hpp"

#define ITERATIONS 1000

struct VertexProperties {
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

vertex_t get_root(Graph* graph) {
    std::pair<vertex_iterator,vertex_iterator> it = boost::vertices(*graph);
    return *it.first;
}

bool has_unborn(vertex_t vertex, Graph* graph) {
    return num_possible_moves() != boost::out_degree(vertex, (*graph));
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

vertex_t get_parent(vertex_t vertex, Graph* graph) {
    in_edge_iterator in_begin, in_end;
    
    for (boost::tie(in_begin, in_end) = boost::in_edges(vertex, *graph); in_begin != in_end; ++in_begin) {
        return boost::source(*in_begin, *graph);
    }
}

/** Add child and return it. */
vertex_t add_child(vertex_t vertex, Graph* graph) {
    vertex_t child = add_vertex(VertexProperties{.wins = 0, .visits = 0}, (*graph)); // Add node and return the vertex descriptor

    add_edge(vertex, child, (*graph));

    return child;
}

void backpropagate(vertex_t vertex, Graph* graph, int result, Spot player) {
    std::list<vertex_t> path;
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

bool is_leaf(vertex_t vertex, Graph* graph) {
    return !boost::in_degree(vertex, *graph) || !boost::out_degree(vertex, *graph);
}

Spot who_played(State state) {
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
        std::list<vertex_t> path;
        int player = 2;
        vertex_t vertex = v;
        State board = rootState;

        while (boost::in_degree(vertex, g)) { // While we are not at the root
            path.push_front(vertex);
            vertex = get_parent(vertex, &g);
        }

        vertex_t parent = vertex;
        for (vertex_t testVertex : path) {
            player = player ^ 3;
            board = make_play(board, testVertex, &g, Spot{(unsigned int) player});
            parent = testVertex;
        }

        out << " [label=\"" << print_environment(board) << "\"]" << std::endl;
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
        vertex_t pointing_at = boost::target(e, g);
        out << " [label=\""<< g[pointing_at].wins / g[pointing_at].visits << "\"]" << std::endl;
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

std::list<vertex_t> get_children(vertex_t vertex, Graph* graph) {
    std::list<vertex_t> children;
    out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, (*graph)); ei != ei_end; ++ei) {
        children.push_back(boost::target(*ei, *graph));
    }

    return children;
}

vertex_t make_best_play(Graph* graph) {
    vertex_t root = get_root(graph);
    float max_value = -1;
    float new_value;
    int best_index = 0;
    int i = 0;
    vertex_t best;
    out_edge_iterator ei, ei_end;

    // Calculate best play
    for (boost::tie(ei, ei_end) = boost::out_edges(root, (*graph)); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, *graph);

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

vertex_t make_human_play(State parentState, vertex_t root, Graph* graph, unsigned int action) {
    int options[9];
    int numOptions = 0;
    for (int i = 0; i < 9; i++) { // Find all the possible children of the parent
        if (!parentState.spots[i].x) {
            options[numOptions++] = i;
        }
        if (i == action) {
            break;
        }
    }

    vertex_t target;
    int numCurrentChildren = boost::out_degree(root, *graph);
    if (numOptions <= numCurrentChildren) {
        int childIndex = 0;
        out_edge_iterator ei, ei_end;
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

    Graph graph;

    rootState = State{0};

    vertex_t vertex = boost::add_vertex(VertexProperties{.wins = 0, .visits = 0}, graph); // Add root node

    while (1) {
        int term;

         // MCTS iterations
        for (int it = 0; it < ITERATIONS; it++) {
            State currentState = rootState;

            vertex = get_root(&graph);
            Spot player = {1}; // Start with AI player

            // Select
            while (boost::out_degree(vertex, graph) && !has_unborn(currentState, vertex, &graph)) {
                vertex = select_child(vertex, &graph);
                currentState = make_play(currentState, vertex, &graph, player);
                player.x = player.x ^ 3; // Switch player
            }

            // Expand
            if (!(term = terminal(&currentState))) {
                if (has_unborn(currentState, vertex, &graph)) {
                    vertex = add_child(vertex, &graph);
                    currentState = make_play(currentState, vertex, &graph, player);
                    player.x = player.x ^ 3; // Switch player
                }

                // Simulate
                State simState = currentState;
                while (!(term = terminal(&simState))) {
                    simState = simulate(&simState, player);
                    player.x = player.x ^ 3; // Switch player
                }
            }

            // Backpropagate
            player.x = player.x ^ 3; // Switch player
            backpropagate(vertex, &graph, term, player);
        }
        
        write_dot(&graph, write_iteration++);

        vertex = make_best_play(&graph);
        rootState = make_play(rootState, vertex, &graph, Spot{1});
        root_changeover(vertex, &graph);

        std::cout << print_environment(rootState);
        std::cout.flush();

        term = terminal(&rootState);

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
            if (rootState.spots[humanPlay].x || humanPlay > 8 || humanPlay < 0) {
                printf("Cannot make play. Enter play: ");
            } else {
                break;
            }
        }

        vertex = make_human_play(rootState, vertex, &graph, humanPlay);
        rootState = make_play(rootState, vertex, &graph, Spot{2});
        root_changeover(vertex, &graph);

        std::cout << print_environment(rootState);
        std::cout.flush();

        printf("\n");

        term = terminal(&rootState);

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