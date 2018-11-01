#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <ctime>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <float.h>

#include <iomanip>
#include <fstream>

#include "CoinPragma.hpp"
#include "CoinTime.hpp"
#include "CoinError.hpp"

#include "BonOsiTMINLPInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "BonCbc.hpp"
#include "BonBonminSetup.hpp"

#include "BonOACutGenerator2.hpp"
#include "BonEcpCuts.hpp"
#include "BonOaNlpOptim.hpp"

#include "BonAmplInterface.hpp"

#include "BonTMINLP.hpp"

#define MAX_STATES 10

// #define DISPLAY_MODE

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

        ~LimitedMemoryMCTS();

        vertex_t add_vertex(vertex_properties vp);

        vertex_t get_root();
        
        vertex_t select_child(vertex_t vertex);

        vertex_t select(vertex_t vertex);

        vertex_t get_parent(vertex_t vertex);

        std::map<vertex_t, std::pair<double, double>> prop_UCBs(vertex_t vertex, double currentUCB, int currentDepth);

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

        void optimise_states(bool mini);

        bool has_state(vertex_t vertex);

        double UCB(double wins, int visits, int parent_visits);

        void reset_num_regenerated();

        void get_num_regenerated(int* nums);

    private:
        Graph graph = Graph();
        int num_states = 0;
        vertex_t root;
        int num_regenerated = 0;
        int ref_num_regenerated = 0;
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

template<typename vertex_properties>
LimitedMemoryMCTS<vertex_properties>::LimitedMemoryMCTS(vertex_properties vp) {
    root = add_vertex(vp);
}

/* Destructor, deletes all remaining states */
template<typename vertex_properties>
LimitedMemoryMCTS<vertex_properties>::~LimitedMemoryMCTS() {
    for (auto u : boost::make_iterator_range(boost::vertices(graph))) {
        if (graph[u].has_state) {
            delete graph[u].state;
        }
    }
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
double LimitedMemoryMCTS<vertex_properties>::UCB(double wins, int visits, int parent_visits) {
    return wins / visits + 0.5 * sqrt(log(parent_visits) / visits);
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::select_child(vertex_t vertex) {
    double max_value = 0;
    double new_value;
    vertex_t best = 0;
    out_edge_iterator ei, ei_end;

    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, graph); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, graph);     

        new_value = UCB(graph[target].wins, graph[target].visits, graph[vertex].visits);

        if (new_value > max_value) {
            max_value = new_value;
            best = target;
        }
    }

    return best;
}

template<typename vertex_properties>
bool LimitedMemoryMCTS<vertex_properties>::has_unborn(vertex_t vertex) {
    return graph[vertex].possible_children != boost::out_degree(vertex, graph);
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::select(vertex_t vertex) {
    while (boost::out_degree(vertex, graph) && !has_unborn(vertex)) {
        vertex = select_child(vertex);
        ref_num_regenerated++;
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
void LimitedMemoryMCTS<vertex_properties>::forget_state(vertex_t vertex) {
    if (graph[vertex].has_state) {
        delete graph[vertex].state;
        graph[vertex].has_state = false;
        num_states--;
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
    std::vector<vertex_t> path = std::vector<vertex_t>();
    vertex_t parent = get_parent(vertex);
    
    path.push_back(vertex);
    path.push_back(parent);
    while (!graph[parent].has_state) {
        parent = get_parent(parent);
        path.push_back(parent);
    }
    typename std::vector<vertex_t>::iterator test = path.begin();
    ++test;
    vertex_t child;
    for (typename std::vector<vertex_t>::iterator it = path.end(); it-- != test;) {
        parent = *it;
        it--;
        child = *it;
        it++;
        graph[parent].regenerate_child(get_child_num(parent, child), &graph[child]);
        num_regenerated++;
        num_states++;
        if (parent != path[path.size() - 1]) { // Forget the parent again if it's not the node which we are regenerating from
            forget_state(parent);
        }
    }
}

/* Proportional UCB = UCB as a percentage * depth */
template<typename vertex_properties>
std::map<typename LimitedMemoryMCTS<vertex_properties>::vertex_t, std::pair<double, double>> LimitedMemoryMCTS<vertex_properties>::prop_UCBs(vertex_t vertex, double currentUCB, int currentDepth) {
    std::map<vertex_t, std::pair<double, double>> proportional_UCBs;

    out_edge_iterator ei, ei_end;
    std::list<std::pair<double, vertex_t>> UCTs;
    double UCTsum = 0;
    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, graph); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, graph);
        double uct = UCB(graph[target].wins, graph[target].visits, graph[vertex].visits);
        if (std::isnan(uct) || uct == 0) {
            uct = 0.00001; // Really small number more than 0
        }
        UCTsum += uct;
        UCTs.push_back(std::make_pair(uct, target));
    }

    for (typename std::list<std::pair<double, vertex_t>>::iterator it = UCTs.begin(); it != UCTs.end(); it++) {
        double new_UCB = currentUCB * ((*it).first / UCTsum);
        // double new_UCB = (*it).first;
        if (std::isnan(new_UCB)) {
            new_UCB = 0;
        }
        proportional_UCBs[(*it).second] = std::make_pair(new_UCB, new_UCB * currentDepth);
        if (!is_leaf((*it).second)) { // Keep walking tree
            std::map<vertex_t, std::pair<double, double>> other_UCBs = prop_UCBs((*it).second, new_UCB, currentDepth + 1);
            proportional_UCBs.insert(other_UCBs.begin(), other_UCBs.end());
        }
    }

    return proportional_UCBs;
}

using namespace  Ipopt;
using namespace Bonmin;

class MyTMINLP : public TMINLP {
    public:
        /// Default constructor.
        MyTMINLP(std::vector<int> nodesofinterest, std::vector<double> propucb, std::vector<double> propucbdepth, std::vector<std::vector<int>> p, double* r, std::vector<Index> hr, std::vector<Index> hc, int nh, std::vector<int> sp, bool* nonlin);
        
        /// virtual destructor.
        virtual ~MyTMINLP(){}

            // /** Copy constructor.*/
        // MyTMINLP(const MyTMINLP &other):printSol_(other.printSol_){}
        /** Assignment operator. no data = nothing to assign*/
        //MyTMINLP& operator=(const MyTMINLP&) {}

        
        /** \name Overloaded functions specific to a TMINLP.*/
        //@{
        /** Pass the type of the variables (INTEGER, BINARY, CONTINUOUS) to the optimizer.
             \param n size of var_types (has to be equal to the number of variables in the problem)
        \param var_types types of the variables (has to be filled by function).
        */
        virtual bool get_variables_types(Index n, VariableType* var_types);
        
        /** Pass info about linear and nonlinear variables.*/
        virtual bool get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types);

        /** Pass the type of the constraints (LINEAR, NON_LINEAR) to the optimizer.
         \param m size of const_types (has to be equal to the number of constraints in the problem)
        \param const_types types of the constraints (has to be filled by function).
        */
        virtual bool get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types);
        //@}  
            
        /** \name Overloaded functions defining a TNLP.
             * This group of function implement the various elements needed to define and solve a TNLP.
             * They are the same as those in a standard Ipopt NLP problem*/
        //@{
        /** Method to pass the main dimensions of the problem to Ipopt.
                \param n number of variables in problem.
                \param m number of constraints.
                \param nnz_jac_g number of nonzeroes in Jacobian of constraints system.
                \param nnz_h_lag number of nonzeroes in Hessian of the Lagrangean.
                \param index_style indicate wether arrays are numbered from 0 (C-style) or
                from 1 (Fortran).
                \return true in case of success.*/
        virtual bool get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,
                                    Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style);
        
        /** Method to pass the bounds on variables and constraints to Ipopt. 
             \param n size of x_l and x_u (has to be equal to the number of variables in the problem)
            \param x_l lower bounds on variables (function should fill it).
            \param x_u upper bounds on the variables (function should fill it).
            \param m size of g_l and g_u (has to be equal to the number of constraints in the problem).
            \param g_l lower bounds of the constraints (function should fill it).
            \param g_u upper bounds of the constraints (function should fill it).
        \return true in case of success.*/
        virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                    Index m, Number* g_l, Number* g_u);
        
        /** Method to to pass the starting point for optimization to Ipopt.
            \param init_x do we initialize primals?
            \param x pass starting primal points (function should fill it if init_x is 1).
            \param m size of lambda (has to be equal to the number of constraints in the problem).
            \param init_lambda do we initialize duals of constraints? 
            \param lambda lower bounds of the constraints (function should fill it).
            \return true in case of success.*/
        virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                        bool init_z, Number* z_L, Number* z_U,
                                        Index m, bool init_lambda,
                                        Number* lambda);
        
        /** Method which compute the value of the objective function at point x.
            \param n size of array x (has to be the number of variables in the problem).
            \param x point where to evaluate.
            \param new_x Is this the first time we evaluate functions at this point? 
            (in the present context we don't care).
            \param obj_value value of objective in x (has to be computed by the function).
            \return true in case of success.*/
        virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

        /** Method which compute the gradient of the objective at a point x.
            \param n size of array x (has to be the number of variables in the problem).
            \param x point where to evaluate.
            \param new_x Is this the first time we evaluate functions at this point? 
            (in the present context we don't care).
            \param grad_f gradient of objective taken in x (function has to fill it).
            \return true in case of success.*/
        virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

        /** Method which compute the value of the functions defining the constraints at a point
            x.
            \param n size of array x (has to be the number of variables in the problem).
            \param x point where to evaluate.
            \param new_x Is this the first time we evaluate functions at this point? 
            (in the present context we don't care).
            \param m size of array g (has to be equal to the number of constraints in the problem)
            \param grad_f values of the constraints (function has to fill it).
            \return true in case of success.*/
        virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

        /** Method to compute the Jacobian of the functions defining the constraints.
            If the parameter values==NULL fill the arrays iCol and jRow which store the position of
            the non-zero element of the Jacobian.
            If the paramenter values!=NULL fill values with the non-zero elements of the Jacobian.
            \param n size of array x (has to be the number of variables in the problem).
            \param x point where to evaluate.
            \param new_x Is this the first time we evaluate functions at this point? 
            (in the present context we don't care).
            \param m size of array g (has to be equal to the number of constraints in the problem)
            \param grad_f values of the constraints (function has to fill it).
            \return true in case of success.*/
        virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                                Index m, Index nele_jac, Index* iRow, Index *jCol,
                                Number* values);
        
        virtual bool eval_h(Index n, const Number* x, bool new_x,
                            Number obj_factor, Index m, const Number* lambda,
                            bool new_lambda, Index nele_hess, Index* iRow,
                            Index* jCol, Number* values);

        
        /** Method called by Ipopt at the end of optimization.*/  
        virtual void finalize_solution(TMINLP::SolverReturn status,
                                        Index n, const Number* x, Number obj_value);
        
        //@}

        virtual const SosInfo * sosConstraints() const{return NULL;}
        virtual const BranchingInfo* branchingInfo() const{return NULL;}

    private:
        int num_vertices;
        std::vector<int> nodes_of_interest;
        std::vector<double> prop_UCB;
        std::vector<double> prop_UCB_Depth;
        std::vector<std::vector<int>> paths;
        double* result;
        std::vector<Index> hessRow;
        std::vector<Index> hessCol;
        int num_hess;
        std::vector<int> start_points;
        bool* is_non_linear;
};

/* Bonmin stuff */
MyTMINLP::MyTMINLP(std::vector<int> nodesofinterest, std::vector<double> propucb, std::vector<double> propucbdepth, std::vector<std::vector<int>> p, double* r, std::vector<Index> hr, std::vector<Index> hc, int nh, std::vector<int> sp, bool* nonlin) {
    num_vertices = propucb.size();
    nodes_of_interest = nodesofinterest;
    paths = p;
    prop_UCB = propucb;
    prop_UCB_Depth = propucbdepth;
    result = r;
    hessRow = hr;
    hessCol = hc;
    num_hess = nh;
    start_points = sp;
    is_non_linear = nonlin;
}

bool MyTMINLP::get_variables_types(Index n, VariableType* var_types) {
    for (int i = 0; i < num_vertices; i++) {
        var_types[i] = BINARY;
    }

    return true;
}

bool MyTMINLP::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types) {
    for (int i = 0; i < num_vertices; i++) {
        if (is_non_linear[i]) {
            var_types[i] = Ipopt::TNLP::NON_LINEAR;
        } else {
            var_types[i] = Ipopt::TNLP::LINEAR;
        }
    }

    return true;
}

bool MyTMINLP::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types) {
    const_types[0] = Ipopt::TNLP::LINEAR;
    return true;
}
bool MyTMINLP::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style) {
    n = num_vertices; // number of variables
    m = 1; // number of constraints
    nnz_jac_g = num_vertices; // number of non zeroes in Jacobian
    nnz_h_lag = num_hess; // number of non zeroes in Hessian of Lagrangean (symmetric, only one triangle)

    index_style = TNLP::C_STYLE;
    return true;
}

bool MyTMINLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
    for (int i = 0; i < num_vertices; i++) {
        x_l[i] = 1.;
        x_u[i] = 1.;  
    }

    for (int index : nodes_of_interest) { // activate decision variables that are of interest
        x_l[index] = 0;
    }

    g_l[0] = -DBL_MAX;
    g_u[0] = MAX_STATES - 4;

    return true;
}

bool MyTMINLP::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
    // start at the current configuration
    for (int i = 0; i < num_vertices; i++) {
        x[i] = start_points[i];
    }

    return true;
}

bool MyTMINLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
    obj_value = 0;
    for (auto path : paths) { // non-linear path iteration part
        double states = 0;
        double prev_prod = 1;
        for (int i : path) {
            prev_prod *= x[i];
            states += prev_prod;
        }
        obj_value -= prop_UCB[path[0]] * (path.size() - states); // multiply proportional UCB of leaf by number of nodes needed to regenerate in path
    }

    for (int i = 0; i < num_vertices; i++) { // linear whole tree part
        obj_value -= prop_UCB_Depth[i] * (1 - x[i]); // multiply (proportional UCB * depth) by state
    }

    return true;
}

bool MyTMINLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
    for (int j = 0; j < num_vertices; j++) {
        grad_f[j] = 0;
        for (auto path : paths) { // non-linear path iteration part
            double states = 0;
            double prev_prod = 1;
            bool in_path = false;
            for (int i : path) {
                if (i == j) { // term contains the variable we are differentiating by
                    in_path = true;
                } else {
                    prev_prod *= x[i];
                    if (prev_prod == 0) {
                        break;
                    }
                }

                if (in_path) { // subsequent terms all will have the variable so add them to the sum
                    states += prev_prod;
                }
            }
            grad_f[j] += prop_UCB[path[0]] * states; // add this component of the gradient (non-linear component)
        }
        grad_f[j] += prop_UCB_Depth[j]; // linear component
    }

    return true;
}

bool MyTMINLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
    g[0] = num_vertices;

    for (int i = 0; i < num_vertices; i++) {
        g[0] -= x[i];
    }
    
    return true;
}

bool MyTMINLP::eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nnz_jac, Index* iRow, Index *jCol, Number* values) {
    if (values == NULL) {

        for (int i = 0; i < num_vertices; i++) {
            iRow[i] = 0;
            jCol[i] = i;
        }
    } else {

        for (int i = 0; i < num_vertices; i++) {
            values[i] = -1;
        }
    }

    return true;
}

bool MyTMINLP::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values) {
    Index idx;
    if (values == NULL) {
        for (idx = 0; idx < hessRow.size(); idx++) {
            iRow[idx] = hessRow[idx];
            jCol[idx] = hessCol[idx];
        }
    } else {
        for (idx = 0; idx < hessRow.size(); idx++) {
            values[idx] = 0;
            int row = hessRow[idx];
            int col = hessCol[idx];

            for (auto path : paths) { // non-linear path iteration part
                bool found_row = false; // First differentiation
                bool found_col = false; // Second differentiation
                double states = 0;
                double prev_prod = 1;
                for (int i : path) {
                    if (i == row) { // term contains the variable we are differentiating by
                        found_row = true;
                    } else if (i == col) {
                        found_col = true;
                    } else {
                        prev_prod *= x[i];
                        if (prev_prod == 0) {
                            break;
                        }
                    }

                    if (found_row && found_col) { // subsequent terms all will have the variables so add them to the sum
                        states += prev_prod;
                    }
                }
                if (found_row && found_col) { // There is at least one term that will differentiate twice into this matrix index
                    values[idx] += prop_UCB[path[0]] * states;
                }
            }
        }
    }

    return true;
}

void MyTMINLP::finalize_solution(TMINLP::SolverReturn status, Index n, const Number* x, Number obj_value) {
    for (int i = 0; i < num_vertices; i++) {
        result[i] = x[i];
    }
}

template<typename vertex_properties>
void LimitedMemoryMCTS<vertex_properties>::optimise_states(bool mini) {
    std::map<vertex_t, std::pair<double, double>> prop_UCB_Map = prop_UCBs(root, 1, 1);
    std::vector<double> prop_UCB;
    std::vector<double> prop_UCB_Depth;
    std::map<vertex_t, int> vertices_to_index = std::map<vertex_t, int>();
    std::map<int, vertex_t> index_to_vertex = std::map<int, vertex_t>();
    std::vector<vertex_t> leaves = std::vector<vertex_t>();
    int j = 0;
    std::vector<int> nodes_of_interest = std::vector<int>();
    std::vector<int> start_points = std::vector<int>();

    for (auto u : boost::make_iterator_range(boost::vertices(graph))) {
        if (u != root) {
            if (!(mini && !graph[u].has_state)) {
                nodes_of_interest.push_back(j);
            }
            vertices_to_index[u] = j;
            index_to_vertex[j] = u;
            std::pair<double, double> ucbs = prop_UCB_Map[u];
            prop_UCB.push_back(ucbs.first);
            prop_UCB_Depth.push_back(ucbs.second);
            start_points.push_back(!graph[u].has_state);
            j++;
        }
        if (is_leaf(u)) {
            leaves.push_back(u);
        }
    }
    std::vector<std::vector<int>> paths = std::vector<std::vector<int>>();
    for (auto it = leaves.begin(); it != leaves.end(); it++) {
        std::vector<int> path = std::vector<int>();
        vertex_t vertex = *it;
        while (vertex != root) { // create path to root
            path.push_back(vertices_to_index[vertex]);
            vertex = get_parent(vertex);
        }
        
        paths.push_back(path);
    }

    double* results = new double[vertices_to_index.size()];

    int num_hess = 0;
    std::vector<Index> hessRow;
    std::vector<Index> hessCol;

    bool* is_non_linear = new bool[vertices_to_index.size()]; // vector of whether each variable is in at least one non-linear term

    for (int i = 0; i < vertices_to_index.size(); i++) {
        is_non_linear[i] = false;
    }

    // Find non-zero elements in Hessian of Lagrangian
    for (int row = 1; row < prop_UCB.size(); row++) {
        for (int col = 0; col < row + 1; col++) { // Only look at lower triangle (no non-zeroes will be in diagonal for this problem)
            for (auto path : paths) { // non-linear path iteration part
                bool found_row = false; // First differentiation
                bool found_col = false; // Second differentiation
                for (int i : path) {
                    if (path.size() > 1) {
                        is_non_linear[i] = true;
                    }
                    if (i == row) { // term contains the variable we are differentiating by
                        found_row = true;
                    } else if (i == col) {
                        found_col = true;
                    }
                }
                if (found_row && found_col) { // There is at least one term that will differentiate twice into this matrix index
                    num_hess++;
                    hessRow.push_back(row);
                    hessCol.push_back(col);
                    break; // We can continue onto the next matrix index now
                }
            }
        }
    }

    using namespace Ipopt;
    using namespace Bonmin;
    SmartPtr<MyTMINLP> tminlp = new MyTMINLP(nodes_of_interest, prop_UCB, prop_UCB_Depth, paths, results, hessRow, hessCol, num_hess, start_points, is_non_linear);

    BonminSetup bonmin;
    bonmin.initializeOptionsAndJournalist();
    bonmin.options()->SetIntegerValue("print_level", 0); // don't print debug from IP Opt
    bonmin.options()->SetIntegerValue("bonmin.bb_log_level", 0); // don't print debug from Bonmin
    bonmin.options()->SetIntegerValue("bonmin.fp_log_level", 0);
    bonmin.options()->SetIntegerValue("bonmin.lp_log_level", 0);
    bonmin.options()->SetIntegerValue("bonmin.milp_log_level", 0);
    bonmin.options()->SetIntegerValue("bonmin.nlp_log_at_root", 0);
    bonmin.options()->SetIntegerValue("bonmin.nlp_log_level", 0);
    bonmin.options()->SetIntegerValue("bonmin.oa_cuts_log_level", 0);
    bonmin.options()->SetIntegerValue("bonmin.oa_log_level", 0);
    bonmin.options()->SetStringValue("sb","yes");

    //Now initialize from tminlp
    bonmin.initialize(GetRawPtr(tminlp));

    //Set up done, now let's branch and bound
    try {
        Bab bb;
        bb(bonmin);//process parameter file using Ipopt and do branch and bound using Cbc
    }
    catch(TNLPSolver::UnsolvedError *E) {
        //There has been a failure to solve a problem with Ipopt.
        std::cerr<<"Ipopt has failed to solve a problem"<<std::endl;
    }
    catch(OsiTMINLPInterface::SimpleError &E) {
        std::cerr<<E.className()<<"::"<<E.methodName()
            <<std::endl
            <<E.message()<<std::endl;
    }
    catch(CoinError &E) {
        std::cerr<<E.className()<<"::"<<E.methodName()
            <<std::endl
            <<E.message()<<std::endl;
    }

    std::map<vertex_t, int> states = std::map<vertex_t, int>();
    j = 0;
    for (auto u : boost::make_iterator_range(boost::vertices(graph))) {
        if (u != root) {
            states[u] = !results[j++]; // Inverted as model uses 1 = no state
        }
    }
    states[root] = 1;

    // For each path follow it down and regenerate the states required
    std::vector<vertex_t> path = std::vector<vertex_t>();
    for (auto it = leaves.begin(); it != leaves.end(); it++) {
        vertex_t vertex = *it;
        while (vertex != root) { // create path to root
            path.push_back(vertex);
            vertex = get_parent(vertex);
        }
        path.push_back(root);

        vertex_t child;
        // Follow path and regenerate states
        for (int i = path.size() - 1; i > 0; i--) {
            vertex_t parent = path[i];
            child = path[i - 1];

            if (i == 1 && !states[child]) { // ensure last child's state is forgotten if needed
                forget_state(child);
            } else if (!graph[child].has_state) { // Regenerate child state
                graph[parent].regenerate_child(get_child_num(parent, child), &graph[child]);
                num_states++;
            }

            if (!states[parent]) { // parent should not keep state
                forget_state(parent);
            }
        }
        path.clear();
    }

    delete[] results;
    delete[] is_non_linear;
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

    lastAdded = child;

    return child;
}

template<typename vertex_properties>
bool LimitedMemoryMCTS<vertex_properties>::is_leaf(vertex_t vertex) {
    return !boost::out_degree(vertex, graph);
}

template<typename vertex_properties>
std::list<typename LimitedMemoryMCTS<vertex_properties>::vertex_t> LimitedMemoryMCTS<vertex_properties>::delete_branch(vertex_t vertex) { // recursively delete all descendant nodes and itself
    std::list<vertex_t> garbage;
    out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, graph); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, graph);
        garbage.merge(delete_branch(target));
    }

    garbage.push_back(vertex);

    return garbage;
}

template<typename vertex_properties>
void LimitedMemoryMCTS<vertex_properties>::root_changeover(vertex_t chosen) {
    // Make the play the new root by first removing all irrelevant branches
    std::list<vertex_t> garbage;
    out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(root, graph); ei != ei_end; ++ei) {
        vertex_t target = boost::target(*ei, graph);
        if (target != chosen) {
            garbage.merge(delete_branch(target));
        }
    }

    for (vertex_t to_remove : garbage) {
        boost::clear_vertex(to_remove, graph);
        if (graph[to_remove].has_state) {
            delete graph[to_remove].state;
        }
        boost::remove_vertex(to_remove, graph);
    }

    boost::clear_vertex(root, graph); // Remove old root. Remove all edges otherwise undefined behaviour occurs after remove_vertex
    delete graph[root].state;
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
    double max_value = -1;
    double new_value;
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

    int term = terminal(best);

    if (!graph[best].has_state) {
        regenerate(best);
    }

    return best;
}

template<typename vertex_properties>
int LimitedMemoryMCTS<vertex_properties>::terminal(vertex_t vertex) {
    return graph[vertex].terminal;
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
            delete vp.state; // Dellocate the unused vp state
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

template<typename vertex_properties>
bool LimitedMemoryMCTS<vertex_properties>::has_state(vertex_t vertex) {
    return graph[vertex].has_state;
}

template<typename vertex_properties>
void LimitedMemoryMCTS<vertex_properties>::reset_num_regenerated() {
    num_regenerated = 0;
    ref_num_regenerated = 0;
}

template<typename vertex_properties>
void LimitedMemoryMCTS<vertex_properties>::get_num_regenerated(int* nums) {
    nums[0] = num_regenerated;
    nums[1] = ref_num_regenerated;
}

int main() {
    // int write_iteration = 0; Each function gets its own color in the visualization graphs: e.g. fun
    
    std::srand(time(NULL));

    MCTS mcts = MCTS(VertexProperties(new Board{0}, Spot{2}));

    MCTS::vertex_t vertex;

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
                mcts.optimise_states(false);
            // } else {
                // mcts.optimise_states(true); // mini optimise
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