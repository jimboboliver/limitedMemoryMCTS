#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <ctime>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#define ITERATIONS 10000
#define MAX_STATES 5

template<typename vertex_properties>
class LimitedMemoryMCTS {
    public:
        /* define the graph type
            listS: selects the STL list container to store the OutEdge list
            listS: selects the STL list container to store the vertices
            bidirectionalS: selects bidirectional edges

            template allows for vertex_properties to be written for specific use
        */
        typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, vertex_properties> Graph;
        typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename boost::graph_traits<Graph>::edge_descriptor edge_t;
        typedef typename boost::graph_traits<Graph>::edge_iterator edge_iterator;
        typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
        typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
        typedef typename boost::graph_traits<Graph>::in_edge_iterator in_edge_iterator;
        typedef std::map<vertex_t, size_t> IndexMap;

        LimitedMemoryMCTS(vertex_properties rootVP);

        vertex_t add_vertex(vertex_properties vp);

        vertex_t get_root();
        
        vertex_t select_child(vertex_t vertex);

        vertex_t select(vertex_t vertex);

        vertex_t get_parent(vertex_t vertex);

        std::list<std::pair<float, vertex_t>> leaf_probabilities(vertex_t vertex, float currentProbability);

        void remove_worst_state();

        int get_child_num(vertex_t parent, vertex_t child);

        void regenerate(vertex_t vertex);

        /* Add child and return it. */
        vertex_t add_child(vertex_t vertex);

        bool is_leaf(vertex_t vertex);

        std::list<vertex_t> delete_branch(vertex_t vertex);

        void root_changeover(vertex_t chosen);

        std::list<vertex_properties> get_child_properties(vertex_t vertex);

        std::list<vertex_t> get_children(vertex_t vertex);

        vertex_t make_best_play();

        int terminal(vertex_t vertex);
        
        bool has_unborn(vertex_t vertex);

        int simulate(vertex_t vertex);

        void backpropagate(vertex_t vertex, int result);
    
        Graph* get_graph();

        vertex_properties get_vertex_properties(vertex_t vertex);

        vertex_t make_play(vertex_t root, vertex_properties vp);

        void make_new_graph();

        int get_num_states();

        void forget_state(vertex_t vertex);

        void optimise_states();

        int expected_regeneration(std::map<vertex_t, bool> states);

    private:
        Graph graph = Graph();
        int num_states = 0;
        vertex_t root;
};

/* Game stuff */

typedef struct {
    unsigned x:2; // 0 = nothing, 1 = X, 2 = O
} Spot;

typedef struct {
    Spot spots[9];
} Board;

class VertexProperties {
    public:
        VertexProperties(Board state_stored, Spot last_player);
        VertexProperties();

        bool has_state;
        Board state;
        Spot player;

    public: // Needed for MCTS interface
        float wins;
        int visits;

        unsigned int num_possible_moves();
        int terminal();
        VertexProperties generate_child(int playNum);
        void regenerate_child(int playNum, VertexProperties* child);
        int playout();
        bool equals(VertexProperties otherVP);
};

VertexProperties::VertexProperties() {
    state = Board{0};
    player = Spot{2};
    has_state = true;
    wins = 0;
    visits = 0;
}

VertexProperties::VertexProperties(Board state_stored, Spot last_player) {
    state = state_stored;
    player = last_player;
    has_state = true;
    wins = 0;
    visits = 0;
}

unsigned int VertexProperties::num_possible_moves() {
    unsigned int num_moves = 0;
    for (int i = 0; i < 9; i++) {
        if (!state.spots[i].x) { // If the spot is empty it's a possible move
            num_moves++;
        }
    }
    return num_moves;
}

VertexProperties VertexProperties::generate_child(int playNum) {
    int plays[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    int available_plays[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    for (int j = 0; j < 9; j++) {
        if (state.spots[j].x) {
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

    Board new_state = state;

    if (player.x == 1) {
        new_state.spots[chosen_play].x = 2;
        return VertexProperties(new_state, Spot{2});
    } else {
        new_state.spots[chosen_play].x = 1;
        return VertexProperties(new_state, Spot{1});
    }
}

void VertexProperties::regenerate_child(int playNum, VertexProperties* child) {
    int plays[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    int available_plays[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    for (int j = 0; j < 9; j++) {
        if (state.spots[j].x) {
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

    Board new_state = state;

    if (player.x == 1) {
        new_state.spots[chosen_play].x = 2;
    } else {
        new_state.spots[chosen_play].x = 1;
    }
    child->state = new_state;
    child->has_state = true;
}

bool VertexProperties::equals(VertexProperties otherVP) {
    return state.spots[0].x == otherVP.state.spots[0].x &&
           state.spots[1].x == otherVP.state.spots[1].x &&
           state.spots[2].x == otherVP.state.spots[2].x &&
           state.spots[3].x == otherVP.state.spots[3].x &&
           state.spots[4].x == otherVP.state.spots[4].x &&
           state.spots[5].x == otherVP.state.spots[5].x &&
           state.spots[6].x == otherVP.state.spots[6].x &&
           state.spots[7].x == otherVP.state.spots[7].x &&
           state.spots[8].x == otherVP.state.spots[8].x;
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
    Board currentBoard = state;
    Spot currentPlayer = player;
    currentPlayer.x ^= 3;
    
    while (!(term = terminal_board(currentBoard))) {
        currentBoard = simulate(currentBoard, currentPlayer);
        currentPlayer.x ^= 3;
    }

    return term;
}

/** Return whether tic-tac-toe state is gameover and if so, who won. 0 = not gameover, 1 = AI, 2 = player, 3 = draw */
int VertexProperties::terminal() {
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

typedef LimitedMemoryMCTS<VertexProperties> MCTS;

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

template <typename Map>
struct my_node_writer {
    my_node_writer(Map& g_) : g (g_) {};
    template <class Vertex>
    void operator()(std::ostream& out, Vertex v) {
        MCTS::vertex_t vertex = v;
        out << " [label=\"" << g[vertex].visits << "\"]" << std::endl;
        if (g[vertex].has_state) {
            Board board = g[vertex].state;

            out << " [label=\"" << print_board(board) << "\"]" << std::endl;
            if (terminal_board(board)) {
                out << " [color=\"" << "black" << "\"]" << std::endl;
            } else {
                out << " [color=\"" << get_colour(who_played(board).x) << "\"]" << std::endl;
            }
        } else {
            Board board = g[vertex].state;

            out << " [label=\"" << print_board(board) << "\"]" << std::endl;
            out << " [color=\"" << "red" << "\"]" << std::endl;
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
        MCTS::vertex_t pointing_at = boost::target(e, g);
        // out << " [label=\""<< g[pointing_at].wins / g[pointing_at].visits << "\"]" << std::endl;
        out << " [label=\""<< g[pointing_at].wins << '\n' << g[pointing_at].visits << '\n' << "\"]" << std::endl;
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
    ss << "graph" << write_iteration << ".dot";
    myfile.open(ss.str().c_str());

    boost::write_graphviz(myfile, *graph, node_writer(*graph), edge_writer(*graph), w, boost::make_assoc_property_map(ids));
    myfile.close();
}

template<typename vertex_properties>
LimitedMemoryMCTS<vertex_properties>::LimitedMemoryMCTS(vertex_properties vp) {
    root = add_vertex(vp);
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::add_vertex(vertex_properties vp) {
    if (vp.has_state) {
        num_states++;
    }
    return boost::add_vertex(vp, graph);
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::get_root() {
    return root;
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::select_child(vertex_t vertex) {
    float max_value = 0;
    float new_value;
    vertex_t best;
    out_edge_iterator ei, ei_end;

    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, graph); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, graph);     

        new_value = graph[target].wins / graph[target].visits + 0.5 * sqrt(log(graph[vertex].visits) / graph[target].visits);

        if (new_value > max_value) {
            max_value = new_value;
            best = target;
        }
    }

    return best;
}

template<typename vertex_properties>
bool LimitedMemoryMCTS<vertex_properties>::has_unborn(vertex_t vertex) {
    return graph[vertex].num_possible_moves() != boost::out_degree(vertex, graph);
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::select(vertex_t vertex) {
    while (boost::out_degree(vertex, graph) && !has_unborn(vertex)) {
        vertex = select_child(vertex);
    }

    return vertex;
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::get_parent(vertex_t vertex) {
    in_edge_iterator in_begin, in_end;
    
    for (boost::tie(in_begin, in_end) = boost::in_edges(vertex, graph); in_begin != in_end; ++in_begin) {
        return boost::source(*in_begin, graph);
    }

    printf("Shouldn't get here 2\n");

    exit(1); // Shouldn't get here
}

template<typename vertex_properties>
std::list<std::pair<float, typename LimitedMemoryMCTS<vertex_properties>::vertex_t>> LimitedMemoryMCTS<vertex_properties>::leaf_probabilities(vertex_t vertex, float currentProbability) {
    std::list<std::pair<float, vertex_t>> probabilities;

    out_edge_iterator ei, ei_end;
    std::list<std::pair<float, vertex_t>> UCTs;
    float UCTsum = 0;
    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, graph); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, graph);
        float uct = graph[target].wins / graph[target].visits + sqrt(log(graph[vertex].visits) / graph[target].visits);
        if (std::isnan(uct)) {
            uct = 0;
        }
        UCTsum += uct;
        UCTs.push_back(std::make_pair(uct, target));
    }

    for (typename std::list<std::pair<float, vertex_t>>::iterator it = UCTs.begin(); it != UCTs.end(); it++) {
        float new_probability = currentProbability * ((*it).first / UCTsum);
        if (std::isnan(new_probability)) {
            new_probability = 0;
        }
        if (is_leaf((*it).second) && graph[(*it).second].has_state) { // Add probability to list of leaf probabilities
            probabilities.push_back(std::make_pair(new_probability, (*it).second));
        } else { // Keep walking tree
            probabilities.merge(leaf_probabilities((*it).second, new_probability));
        }
    }

    return probabilities;
}

template<typename vertex_properties>
void LimitedMemoryMCTS<vertex_properties>::forget_state(vertex_t vertex) {
    graph[vertex].has_state = false;
    num_states--;
}

template<typename vertex_properties>
void LimitedMemoryMCTS<vertex_properties>::remove_worst_state() {
    std::list<std::pair<float, vertex_t>> probabilities = leaf_probabilities(get_root(), 1);
    float min_value = 9999999;
    vertex_t min_vertex = 0;
    
    
    for (typename std::list<std::pair<float, vertex_t>>::iterator it = probabilities.begin(); it != probabilities.end(); it++) {
        if ((*it).first < min_value) {
            min_value = (*it).first;
            min_vertex = (*it).second;
        }
    }
    if (min_vertex != 0) {
        forget_state(min_vertex);
    }
}

template<typename vertex_properties>
int LimitedMemoryMCTS<vertex_properties>::get_child_num(vertex_t parent, vertex_t child) {
    out_edge_iterator ei, ei_end;
    int i = 0;
    for (boost::tie(ei, ei_end) = boost::out_edges(parent, graph); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, graph);
        if (target == child) {
            return i;
        }
        i++;
    }

    printf("Shouldn't get here\n");

    exit(1); // Shouldn't get here
}

/* Requires 2 spare states */
template<typename vertex_properties>
void LimitedMemoryMCTS<vertex_properties>::regenerate(vertex_t vertex) {
    std::list<vertex_t> path = std::list<vertex_t>();
    vertex_t parent = get_parent(vertex);
    
    path.push_front(vertex);
    path.push_front(parent);
    while (!graph[parent].has_state) {
        parent = get_parent(parent);
        path.push_front(parent);
    }
    typename std::list<vertex_t>::iterator test = path.end();
    --test;
    vertex_t child;
    for (typename std::list<vertex_t>::iterator it = path.begin(); it != test; it++) {
        parent = *it;
        it++;
        child = *it;
        it--;
        graph[parent].regenerate_child(get_child_num(parent, child), &graph[child]);
        num_states++;
        if (parent != root) { // Forget the parent again if it's not the root
            forget_state(parent);
        }
    }
}

/* Returns number of nodes that will be generated next iteration due to a given set of saved states */
template<typename vertex_properties>
int LimitedMemoryMCTS<vertex_properties>::expected_regeneration(std::map<vertex_t, bool> states) {
    vertex_t vertex = root;
    while (boost::out_degree(vertex, graph) && !has_unborn(vertex)) { // Follow path to find next select
        vertex = select_child(vertex);
    }

    int num_to_regenerate = 0;
    while (vertex != root && !states[vertex]) {
        vertex = get_parent(vertex);
        num_to_regenerate++;
    }

    return num_to_regenerate;
}
 
template<typename vertex_properties>
void LimitedMemoryMCTS<vertex_properties>::optimise_states() {
    std::map<vertex_t, bool> states = std::map<vertex_t, bool>();
    vertex_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::vertices(graph); ei != ei_end; ++ei) {
        vertex_t target = *ei;
        if (target != root) {
            states[target] = false;
        }
    }
    // std::cout << expected_regeneration(states) << '\n';
}

/* Requires 2 spare states */
/* Add child and return it. */
template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::add_child(vertex_t parent) {
    if (!graph[parent].has_state) { // If parent has no state, regenerate it
        regenerate(parent);
    }

    vertex_properties vp = graph[parent].generate_child(boost::out_degree(parent, graph));
    vertex_t child = boost::add_vertex(vp, graph); // Add node and return the vertex descriptor
    boost::add_edge(parent, child, graph);

    num_states++;

    return child;
}

template<typename vertex_properties>
bool LimitedMemoryMCTS<vertex_properties>::is_leaf(vertex_t vertex) {
    return !boost::in_degree(vertex, graph) || !boost::out_degree(vertex, graph);
}

template<typename vertex_properties>
std::list<typename LimitedMemoryMCTS<vertex_properties>::vertex_t> LimitedMemoryMCTS<vertex_properties>::delete_branch(vertex_t vertex) { // recursively delete all descendant nodes and itself
    std::list<vertex_t> probabilities;
    out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, graph); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, graph);
        probabilities.merge(delete_branch(target));
    }

    probabilities.push_back(vertex);

    return probabilities;
}

template<typename vertex_properties>
void LimitedMemoryMCTS<vertex_properties>::root_changeover(vertex_t chosen) {
    // Make the play the new root by first removing all irrelevant branches
    std::list<vertex_t> probabilities;
    out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(root, graph); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, graph);
        if (target != chosen) {
            probabilities.merge(delete_branch(target));
        }
    }

    for (vertex_t to_remove : probabilities) {
        boost::clear_vertex(to_remove, graph);
        boost::remove_vertex(to_remove, graph);
    }

    boost::clear_vertex(root, graph); // Remove old root. Remove all edges otherwise undefined behaviour occurs after remove_vertex
    boost::remove_vertex(root, graph);

    root = chosen;
}

template<typename vertex_properties>
std::list<vertex_properties> LimitedMemoryMCTS<vertex_properties>::get_child_properties(vertex_t vertex) {
    std::list<vertex_properties> children;

    for (auto vd : boost::make_iterator_range(adjacent_vertices(vertex, graph))) {
        children.push_back(graph[vd]);
    }
    return children;
}

template<typename vertex_properties>
std::list<typename LimitedMemoryMCTS<vertex_properties>::vertex_t> LimitedMemoryMCTS<vertex_properties>::get_children(vertex_t vertex) {
    std::list<vertex_t> children;
    out_edge_iterator ei, ei_end;

    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, graph); ei != ei_end; ++ei) {
        children.push_back(boost::target(*ei, graph));
    }
    return children;
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::make_best_play() {
    vertex_t root = get_root();
    float max_value = -1;
    float new_value;
    int i = 0;
    vertex_t best = 0;
    out_edge_iterator ei, ei_end;

    // Calculate best play
    for (boost::tie(ei, ei_end) = boost::out_edges(root, graph); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, graph);

        new_value = graph[target].wins / graph[target].visits;

        if (new_value > max_value) {
            max_value = new_value;
            best = target;
        }
        i++;
    }

    if (!graph[best].has_state) {
        regenerate(best);
    }

    return best;
}

template<typename vertex_properties>
int LimitedMemoryMCTS<vertex_properties>::terminal(vertex_t vertex) {
    return graph[vertex].terminal();
}

template<typename vertex_properties>
int LimitedMemoryMCTS<vertex_properties>::simulate(vertex_t vertex) {
    return graph[vertex].playout();
}

template<typename vertex_properties>
void LimitedMemoryMCTS<vertex_properties>::backpropagate(vertex_t vertex, int result) {
    std::list<vertex_t> path;
    while (boost::in_degree(vertex, graph)) { // While we are not at the root
        path.push_front(vertex);
        vertex = get_parent(vertex);
    }
    int player = 1; // Start with AI player
    for (typename std::list<vertex_t>::iterator it = path.begin(); it != path.end(); ++it) {
        if (result == player) {
            graph[*it].wins += 1;
        } else if (result == 3) {
            graph[*it].wins += 0.5;
        }
        graph[*it].visits += 1;
        player ^= 3; // Switch player
    }
    // Add statistics to root
    if (result == 2) {
        graph[vertex].wins += 1;
    } else if (result == 3) {
        graph[vertex].wins += 0.5;
    }
    graph[vertex].visits += 1;
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::Graph* LimitedMemoryMCTS<vertex_properties>::get_graph() {
    return &graph;
}

template<typename vertex_properties>
vertex_properties LimitedMemoryMCTS<vertex_properties>::get_vertex_properties(vertex_t vertex) {
    return graph[vertex];
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::make_play(vertex_t root, vertex_properties vp) {
    std::list<vertex_t> children = get_children(root);

    for (typename std::list<vertex_t>::iterator it = children.begin(); it != children.end(); ++it) {
        if (!graph[*it].has_state) {
            regenerate(*it);
        }
        if (graph[*it].equals(vp)) {
            return *it;
        }
    }

    // If we got to this point, we need to add the new node
    vertex_t child = boost::add_vertex(vp, graph); // Add node and return the vertex descriptor

    boost::add_edge(root, child, graph);

    return child;
}

template<typename vertex_properties>
int LimitedMemoryMCTS<vertex_properties>::get_num_states() {
    return num_states;
}

int main() {
    int write_iteration = 0;
    
    std::srand(time(NULL));

    MCTS mcts = MCTS(VertexProperties(Board{0}, Spot{2}));

    MCTS::vertex_t vertex = mcts.get_root();

    while (1) {
        int term;

         // MCTS iterations
        for (int it = 0; it < ITERATIONS; it++) {
            vertex = mcts.get_root();

            // Select
            vertex = mcts.select(vertex);

            // If not terminal
            if (!(term = mcts.terminal(vertex))) {
                // Expand
                if (mcts.has_unborn(vertex)) {
                    vertex = mcts.add_child(vertex);
                }

                // Simulate
                term = mcts.simulate(vertex);
            }

            // Backpropagate
            mcts.backpropagate(vertex, term);

            mcts.optimise_states();
        }
        
        write_dot(mcts.get_graph(), write_iteration++);

        vertex = mcts.make_best_play();
        Board rootBoard = mcts.get_vertex_properties(vertex).state;
        mcts.root_changeover(vertex);

        std::cout << print_board(rootBoard);
        std::cout.flush();

        term = mcts.terminal(vertex);

        if (term == 1) {
            std::cout << "\nGame over! AI wins\n";
            break;
        } else if (term == 3) {
            std::cout << "\nGame over! Draw\n";
            break;
        }

        std::cout << "\nEnter play: ";

        int humanPlay;
        int ch;
        while (1) {
            humanPlay = std::cin.get() - 48;
            while ((ch = std::cin.get()) != '\n' && ch != EOF);
            std::cout << '\n';
            if (humanPlay > 8 || humanPlay < 0 || rootBoard.spots[humanPlay].x) {
                std::cout << "Cannot make play. Enter play: ";
            } else {
                break;
            }
        }

        rootBoard.spots[humanPlay].x = 2;

        vertex = mcts.make_play(vertex, VertexProperties(rootBoard, Spot{2}));
        mcts.root_changeover(vertex);

        std::cout << print_board(rootBoard) << '\n';
        std::cout.flush();

        term = terminal_board(rootBoard);

        if (term == 2) {
            std::cout << "\nGame over! Player wins\n";
            break;
        } else if (term == 3) {
            std::cout << "\nGame over! Draw\n";
            break;
        }
    }

    return 0;  
}