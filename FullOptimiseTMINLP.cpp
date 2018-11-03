#include "FullOptimiseTMINLP.hpp"

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