#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

template<typename vertex_properties>
class LimitedMemoryMCTS {
    public:
        /* define the graph type
            listS: selects the STL list container to store 
                    the OutEdge list
            listS: selects the STL vector container to store 
                the vertices
            directedS: selects directed edges

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

        static vertex_t add_vertex(vertex_properties vp);

        static vertex_t get_root();
        
        static vertex_t select_child(vertex_t vertex);

        static vertex_t select(vertex_t vertex);

        static vertex_t get_parent(vertex_t vertex);

        /* Add child and return it. */
        static vertex_t add_child(vertex_t vertex, vertex_properties vp);

        bool is_leaf(vertex_t vertex);

        static std::list<vertex_t> delete_branch(vertex_t vertex);

        static vertex_t root_changeover(vertex_t chosen);

        std::list<vertex_t> get_children(vertex_t vertex);

        static vertex_t make_best_play();
    
    private:
        Graph graph;
};

#include "mcts.cpp" // Done like this to allow for separation of implementation from template declaration