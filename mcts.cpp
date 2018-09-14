template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::add_vertex(vertex_properties vp) {
    boost::add_vertex(vp, LimitedMemoryMCTS<vertex_properties>::graph);
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::get_root() {
    std::pair<vertex_iterator, vertex_iterator> it = boost::vertices(LimitedMemoryMCTS<vertex_properties>::graph);
    return *it.first;
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::select_child(LimitedMemoryMCTS<vertex_properties>::vertex_t vertex) {
    float max_value = 0;
    float new_value;
    LimitedMemoryMCTS<vertex_properties>::vertex_t best;
    LimitedMemoryMCTS<vertex_properties>::out_edge_iterator ei, ei_end;

    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, (LimitedMemoryMCTS<vertex_properties>::graph)); ei != ei_end; ++ei) {
        LimitedMemoryMCTS<vertex_properties>::vertex_t target = boost::target(*ei, LimitedMemoryMCTS<vertex_properties>::graph);

        new_value = (LimitedMemoryMCTS<vertex_properties>::graph)[target].wins / (LimitedMemoryMCTS<vertex_properties>::graph)[target].visits + sqrt(2 * log((LimitedMemoryMCTS<vertex_properties>::graph)[vertex].visits) / (LimitedMemoryMCTS<vertex_properties>::graph)[target].visits);

        if (new_value > max_value) {
            max_value = new_value;
            best = target;
        }
    }

    return best;
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::select(LimitedMemoryMCTS<vertex_properties>::vertex_t vertex) {
    while (boost::out_degree(vertex, LimitedMemoryMCTS<vertex_properties>::graph) && LimitedMemoryMCTS<vertex_properties>::graph[vertex].num_possible_moves() == boost::out_degree(vertex, LimitedMemoryMCTS<vertex_properties>::graph)) {
        vertex = LimitedMemoryMCTS<VertexProperties>::select_child(vertex);
    }
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::get_parent(LimitedMemoryMCTS<vertex_properties>::vertex_t vertex) {
    LimitedMemoryMCTS<vertex_properties>::in_edge_iterator in_begin, in_end;
    
    for (boost::tie(in_begin, in_end) = boost::in_edges(vertex, LimitedMemoryMCTS<vertex_properties>::graph); in_begin != in_end; ++in_begin) {
        return boost::source(*in_begin, LimitedMemoryMCTS<vertex_properties>::graph);
    }
}

/** Add child and return it. */
template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::add_child(LimitedMemoryMCTS<vertex_properties>::vertex_t vertex, vertex_properties vp) {
    LimitedMemoryMCTS<vertex_properties>::vertex_t child = add_vertex(vp, (LimitedMemoryMCTS<vertex_properties>::graph)); // Add node and return the vertex descriptor

    add_edge(vertex, child, LimitedMemoryMCTS<vertex_properties>::graph);

    return child;
}

template<typename vertex_properties>
bool LimitedMemoryMCTS<vertex_properties>::is_leaf(LimitedMemoryMCTS<vertex_properties>::vertex_t vertex) {
    return !boost::in_degree(vertex, LimitedMemoryMCTS<vertex_properties>::graph) || !boost::out_degree(vertex, LimitedMemoryMCTS<vertex_properties>::graph);
}

template<typename vertex_properties>
std::list<typename LimitedMemoryMCTS<vertex_properties>::vertex_t> LimitedMemoryMCTS<vertex_properties>::delete_branch(LimitedMemoryMCTS<vertex_properties>::vertex_t vertex) { // recursively delete all descendant nodes and itself
    std::list<LimitedMemoryMCTS<vertex_properties>::vertex_t> garbage;
    LimitedMemoryMCTS<vertex_properties>::out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, (LimitedMemoryMCTS<vertex_properties>::graph)); ei != ei_end; ++ei) {
        LimitedMemoryMCTS<vertex_properties>::vertex_t target = boost::target(*ei, LimitedMemoryMCTS<vertex_properties>::graph);
        garbage.merge(LimitedMemoryMCTS<vertex_properties>::delete_branch(target));
    }

    garbage.push_back(vertex);

    return garbage;
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::root_changeover(LimitedMemoryMCTS<vertex_properties>::vertex_t chosen) {
    // Make the play the new root by first removing all irrelevant branches
    std::list<LimitedMemoryMCTS<vertex_properties>::vertex_t> garbage;
    LimitedMemoryMCTS<vertex_properties>::out_edge_iterator ei, ei_end;
    LimitedMemoryMCTS<vertex_properties>::vertex_t root = LimitedMemoryMCTS<vertex_properties>::get_root();
    for (boost::tie(ei, ei_end) = boost::out_edges(root, LimitedMemoryMCTS<vertex_properties>::graph); ei != ei_end; ++ei) {
        LimitedMemoryMCTS<vertex_properties>::vertex_t target = boost::target(*ei, LimitedMemoryMCTS<vertex_properties>::graph);
        if (target != chosen) {
            garbage.merge(LimitedMemoryMCTS<vertex_properties>::delete_branch(target));
        }
    }

    for (LimitedMemoryMCTS<vertex_properties>::vertex_t to_remove : garbage) {
        boost::clear_vertex(to_remove, LimitedMemoryMCTS<vertex_properties>::graph);
        boost::remove_vertex(to_remove, LimitedMemoryMCTS<vertex_properties>::graph);
    }

    boost::clear_vertex(root, LimitedMemoryMCTS<vertex_properties>::graph); // Remove old root. Remove all edges otherwise undefined behaviour occurs after remove_vertex
    boost::remove_vertex(root, LimitedMemoryMCTS<vertex_properties>::graph);
}

template<typename vertex_properties>
std::list<typename LimitedMemoryMCTS<vertex_properties>::vertex_t> LimitedMemoryMCTS<vertex_properties>::get_children(LimitedMemoryMCTS<vertex_properties>::vertex_t vertex) {
    std::list<LimitedMemoryMCTS<vertex_properties>::vertex_t> children;
    LimitedMemoryMCTS<vertex_properties>::out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(vertex, LimitedMemoryMCTS<vertex_properties>::graph); ei != ei_end; ++ei) {
        children.push_back(boost::target(*ei, LimitedMemoryMCTS<vertex_properties>::graph));
    }

    return children;
}

template<typename vertex_properties>
typename LimitedMemoryMCTS<vertex_properties>::vertex_t LimitedMemoryMCTS<vertex_properties>::make_best_play() {
    LimitedMemoryMCTS<vertex_properties>::vertex_t root = LimitedMemoryMCTS<vertex_properties>::get_root();
    float max_value = -1;
    float new_value;
    int best_index = 0;
    int i = 0;
    LimitedMemoryMCTS<vertex_properties>::vertex_t best;
    LimitedMemoryMCTS<vertex_properties>::out_edge_iterator ei, ei_end;

    // Calculate best play
    for (boost::tie(ei, ei_end) = boost::out_edges(root, LimitedMemoryMCTS<vertex_properties>::graph); ei != ei_end; ++ei) {
        LimitedMemoryMCTS<vertex_properties>::vertex_t target = boost::target(*ei, LimitedMemoryMCTS<vertex_properties>::graph);

        new_value = LimitedMemoryMCTS<vertex_properties>::graph[target].wins / LimitedMemoryMCTS<vertex_properties>::graph[target].visits;

        if (new_value > max_value) {
            max_value = new_value;
            best = target;
            best_index = i;
        }
        i++;
    }

    return best;
}