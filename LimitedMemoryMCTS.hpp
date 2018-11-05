#ifndef LMMCTS_INC
#define LMMCTS_INC

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include "OptimiseTMINLP.hpp"

#endif

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