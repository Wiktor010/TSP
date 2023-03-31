#include "TSP.hpp"

#include <algorithm>
#include <stack>
#include <optional>

std::ostream& operator<<(std::ostream& os, const CostMatrix& cm) {
    for (std::size_t r = 0; r < cm.size(); ++r) {
        for (std::size_t c = 0; c < cm.size(); ++c) {
            const auto& elem = cm[r][c];
            os << (is_inf(elem) ? "INF" : std::to_string(elem)) << " ";
        }
        os << "\n";
    }
    os << std::endl;

    return os;
}

/* PART 1 */

/**
 * Create path from unsorted path and last 2x2 cost matrix.
 * @return The vector of consecutive vertex.
 */
path_t StageState::get_path() {
    cost_t lower_bound = reduce_cost_matrix();
    update_lower_bound(lower_bound);

    for (std::size_t i = 0; i < matrix_.size(); i++) {
        for (std::size_t j = 0; j < matrix_[i].size(); j++) {
            if (is_inf(matrix_[i][j])) {
                continue;
            }
            if (matrix_[i][j] == 0) {
                vertex_t vector_to_patch(i, j);
                append_to_path(vector_to_patch);
            }
        }
    }

    unsorted_path_t sorted_path;
    sorted_path.push_back(unsorted_path_[0]);
    for (std::size_t i = 0; i < unsorted_path_.size() - 1; i++) {
        for (std::size_t j = 0; j < unsorted_path_.size(); j++) {
            if (sorted_path[i].col == unsorted_path_[j].row) {
                sorted_path.push_back(unsorted_path_[j]);
            }
        }
    }

    path_t sorted_to_print_path;
    for(std::size_t i = 0; i < sorted_path.size(); i++){
        sorted_to_print_path.push_back(sorted_path[i].col + 1);
    }
    return sorted_to_print_path;
}

/**
 * Get minimum values from each row and returns them.
 * @return Vector of minimum values in row.
 */
std::vector<cost_t> CostMatrix::get_min_values_in_rows() const {
    std::vector<cost_t> min_values;
    int p;
    cost_t min;
    for (std::size_t i = 0; i < matrix_.size(); i++) {
        p = 0;
        min = 0;
        for (std::size_t j = 0; j < matrix_.size(); j++) {
            if (is_inf(matrix_[i][j])){
                continue;
            }
            else if (p != 0 and not is_inf(matrix_[i][j])) {
                if (min > matrix_[i][j]) {
                    min = matrix_[i][j];
                }
                else{
                    continue;
                }
            }
            else if (p == 0 and not is_inf(matrix_[i][j])) {
                min = matrix_[i][j];
                p += 1;
            }
        }
        min_values.push_back(min);
    }
    return min_values;
}


/**
 * Reduce rows so that in each row at least one zero value is present.
 * @return Sum of values reduced in rows.
 */
cost_t CostMatrix::reduce_rows() {
    std::vector<cost_t> min_values = get_min_values_in_rows();
    for(std::size_t i = 0; i < matrix_.size();i++){
        for(std::size_t j = 0; j < matrix_.size(); j++){
            if(!is_inf(matrix_[i][j])){
                matrix_[i][j] = matrix_[i][j] - min_values[i];
            }
        }
    }
    return std::accumulate(std::begin(min_values), std::end(min_values), 0);

}
/**
 * Get minimum values from each column and returns them.
 * @return Vector of minimum values in columns.
 */
std::vector<cost_t> CostMatrix::get_min_values_in_cols() const {
    std::vector<cost_t> min_values;
    int p;
    cost_t min;
    for (std::size_t i = 0; i < matrix_.size(); i++) {
        p = 0;
        min = 0;
        for (auto& j : matrix_) {
            if (is_inf(j[i])){
                continue;
            }
            else if (p != 0 and not is_inf(j[i])) {
                if (min > j[i]) {
                    min = j[i];
                }
                else{
                    continue;
                }
            }
            else if (p == 0 and not is_inf(j[i])) {
                min = j[i];
                p += 1;
            }
        }
        min_values.push_back(min);
    }
    return min_values;
}


/**
 * Reduces rows so that in each column at least one zero value is present.
 * @return Sum of values reduced in columns.
 */
cost_t CostMatrix::reduce_cols() {
    std::vector<cost_t>min_values = get_min_values_in_cols();
    for(auto&i : matrix_){
        for (std::size_t j = 0; j < i.size(); ++j){
            if(!is_inf(i[j])){
                i[j] = i[j] - min_values[j];
            }
        }
    }
    return std::accumulate(std::begin(min_values),std::end(min_values),0);

}

/**
 * Get the cost of not visiting the vertex_t (@see: get_new_vertex())
 * @param row
 * @param col
 * @return The sum of minimal values in row and col, excluding the intersection value.
 */
cost_t CostMatrix::get_vertex_cost(std::size_t row, std::size_t col) const {
    CostMatrix matrix(matrix_);
    cost_t vertex_cost;
    matrix.matrix_[row][col] = INF;
    vertex_cost = matrix.get_min_values_in_rows()[row] + matrix.get_min_values_in_cols()[col];
    return vertex_cost;

}

/* PART 2 */

/**
 * Choose next vertex to visit:
 * - Look for vertex_t (pair row and column) with value 0 in the current cost matrix.
 * - Get the vertex_t cost (calls get_vertex_cost()).
 * - Choose the vertex_t with maximum cost and returns it.
 * @param cm
 * @return The coordinates of the next vertex.
 */
/* Wybierz następny wierzchołek do odwiedzenia:
* - Poszukaj wierzchołka_t (para wiersz i kolumna) o wartości 0 w aktualnej macierzy kosztów.
* - Pobierz koszt vertex_t (wywołuje get_vertex_cost()).
* - Wybierz wierzchołek_t z maksymalnym kosztem i zwróć go.
*/
NewVertex StageState::choose_new_vertex() {
    std::vector<NewVertex> next_vertex = {};

    for(std::size_t r = 0; r < matrix_.size(); r++) {
        for(std::size_t c = 0; c < matrix_.size(); c++) {
            if (matrix_.get_matrix()[r][c] == 0) {
                vertex_t current_vertex(r,c);
                NewVertex new_vertex(current_vertex, matrix_.get_vertex_cost(r, c));
                next_vertex.push_back(new_vertex);

            }

        }
    }
    NewVertex* vertex_ptr = &next_vertex[0];
    for(auto &i : next_vertex) {
        if(i.cost > (*vertex_ptr).cost) {
            vertex_ptr = &i;
        }
    }
    return *vertex_ptr;
}

/**
 * Update the cost matrix with the new vertex.
 * @param new_vertex
 */
void StageState::update_cost_matrix(vertex_t new_vertex) {

    std::size_t row_id = new_vertex.row;
    std::size_t col_id = new_vertex.col;
    for(std::size_t i = 0; i < matrix_.size(); i++) {
        for(std::size_t j = 0; j < matrix_[i].size(); j++) {
            if(i == row_id or j == col_id){
                matrix_[i][j] = INF;
            }
        }
    }
    matrix_[col_id][row_id] = INF;

    unsorted_path_t value;
    for(std::size_t i = 0; i < unsorted_path_.size(); i++) {
        value.push_back(unsorted_path_[i]);
    }
    std::size_t  size_value = value.size();
    std::size_t  size_value_2 = size_value;
    for(std::size_t z = 0; z < size_value_2; z++) {
        size_value = value.size();
        for(std::size_t i = 0; i < size_value; i++) {
            for(std::size_t j = 0; j < size_value; j++) {
                if(value[i].col == value[j].row and i != j) {
                    value.push_back({value[i].row, value[j].col});
                    matrix_[value[j].col][value[i].row] = INF;
                }
            }
        }
    }

}

/**
 * Reduce the cost matrix.
 * @return The sum of reduced values.
 */
cost_t StageState::reduce_cost_matrix() {
    cost_t lower_bound = matrix_.reduce_rows() + matrix_.reduce_cols();
    return lower_bound;
}

/**
 * Given the optimal path, return the optimal cost.
 * @param optimal_path
 * @param m
 * @return Cost of the path.
 */
cost_t get_optimal_cost(const path_t& optimal_path, const cost_matrix_t& m) {
    cost_t cost = 0;

    for (std::size_t idx = 1; idx < optimal_path.size(); ++idx) {
        cost += m[optimal_path[idx - 1] - 1][optimal_path[idx] - 1];
    }

    // Add the cost of returning from the last city to the initial one.
    cost += m[optimal_path[optimal_path.size() - 1] - 1][optimal_path[0] - 1];

    return cost;
}

/**
 * Create the right branch matrix with the chosen vertex forbidden and the new lower bound.
 * @param m
 * @param v
 * @param lb
 * @return New branch.
 */
StageState create_right_branch_matrix(cost_matrix_t m, vertex_t v, cost_t lb) {
    CostMatrix cm(m);
    cm[v.row][v.col] = INF;
    return StageState(cm, {}, lb);
}

/**
 * Retain only optimal ones (from all possible ones).
 * @param solutions
 * @return Vector of optimal solutions.
 */
tsp_solutions_t filter_solutions(tsp_solutions_t solutions) {
    cost_t optimal_cost = INF;
    for (const auto& s : solutions) {
        optimal_cost = (s.lower_bound < optimal_cost) ? s.lower_bound : optimal_cost;
    }

    tsp_solutions_t optimal_solutions;
    std::copy_if(solutions.begin(), solutions.end(),
                 std::back_inserter(optimal_solutions),
                 [&optimal_cost](const tsp_solution_t& s) { return s.lower_bound == optimal_cost; }
    );

    return optimal_solutions;
}

/**
 * Solve the TSP.
 * @param cm The cost matrix.
 * @return A list of optimal solutions.
 */
tsp_solutions_t solve_tsp(const cost_matrix_t& cm) {

    StageState left_branch(cm);

    // The branch & bound tree.
    std::stack<StageState> tree_lifo;

    // The number of levels determines the number of steps before obtaining
    // a 2x2 matrix.
    std::size_t n_levels = cm.size() - 2;

    tree_lifo.push(left_branch);   // Use the first cost matrix as the root.

    cost_t best_lb = INF;
    tsp_solutions_t solutions;

    while (!tree_lifo.empty()) {

        left_branch = tree_lifo.top();
        tree_lifo.pop();

        while (left_branch.get_level() != n_levels && left_branch.get_lower_bound() <= best_lb) {
            // Repeat until a 2x2 matrix is obtained or the lower bound is too high...

            if (left_branch.get_level() == 0) {
                left_branch.reset_lower_bound();
            }

            // 1. Reduce the matrix in rows and columns.
            cost_t new_cost = left_branch.reduce_cost_matrix(); // @TODO (KROK 1)

            // 2. Update the lower bound and check the break condition.
            left_branch.update_lower_bound(new_cost);
            if (left_branch.get_lower_bound() > best_lb) {
                break;
            }

            // 3. Get new vertex and the cost of not choosing it.
            NewVertex new_vertex = left_branch.choose_new_vertex(); // @TODO (KROK 2)

            // 4. @TODO Update the path - use append_to_path method.
            left_branch.append_to_path(new_vertex.coordinates);
            // 5. @TODO (KROK 3) Update the cost matrix of the left branch.
            left_branch.update_cost_matrix(new_vertex.coordinates);
            // 6. Update the right branch and push it to the LIFO.
            cost_t new_lower_bound = left_branch.get_lower_bound() + new_vertex.cost;
            tree_lifo.push(create_right_branch_matrix(cm, new_vertex.coordinates,
                                                      new_lower_bound));
        }

        if (left_branch.get_lower_bound() <= best_lb) {
            // If the new solution is at least as good as the previous one,
            // save its lower bound and its path.
            best_lb = left_branch.get_lower_bound();
            path_t new_path = left_branch.get_path();
            solutions.push_back({get_optimal_cost(new_path, cm), new_path});
        }
    }

    return filter_solutions(solutions); // Filter solutions to find only optimal ones.
}
