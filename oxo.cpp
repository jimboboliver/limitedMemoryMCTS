#include "game.hpp"

typedef struct {
    unsigned x:2; // 0 = nothing, 1 = X, 2 = O
} Spot;

class Board {
    public:
        Spot spots[9] = {0};
        Spot nextPlayer = Spot{1};

        void switch_player() {
            nextPlayer.x = nextPlayer.x ^ 2;
        }
};

Board currentBoard = {0};

char get_spot(Spot spot) {
    if (spot.x == 1) {
        return 'X';
    } else if (spot.x == 2) {
        return 'O';
    } else {
        return '*';
    }
}

std::string print_environment(Board board) {
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

int no_moves_left(Board* board) {
    return board->spots[0].x && board->spots[1].x && board->spots[2].x && board->spots[3].x && board->spots[4].x && 
            board->spots[5].x && board->spots[6].x && board->spots[7].x && board->spots[8].x;
}

/** Return whether tic-tac-toe state is gameover and if so, who won. 0 = not gameover, 1 = AI, 2 = Human, 3 = draw */
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

/**
 * Runs simulation from current board and returns terminal(final_board_state)
 */
int simulate() {
    Board options[9];
    int numOptions = 0;
    int term;

    Board board = currentBoard;
    while (!(term = terminal(&board))) {
        for (int i = 0; i < 9; i++) {
            if (!board.spots[i].x) {
                options[numOptions] = board;
                options[numOptions++].spots[i].x = board.nextPlayer.x;
            }
        }

        board = options[rand() % 9];
        board.switch_player();
    }

    return term;
}

int get_play() {
    int options[9];
    int numOptions = 0;
    for (int i = 0; i < 9; i++) { // Find all the possible children of the parent
        if (!parentBoard.spots[i].x) {
            options[numOptions++] = i;
        }
    }

    // printf("%d %d\n", boost::out_degree(get_parent(vertex, graph), *graph), boost::out_degree(vertex, *graph));

    int childIndex = 0;
    out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(get_parent(vertex, graph), *graph); ei != ei_end; ++ei) {
        if (vertex == boost::target(*ei, *graph)) { // If this is the right child, set the childIndex
            break;
        }
        childIndex++;
    }

    return options[childIndex]; // Return the child's action
}

void make_play(int index) {
    currentBoard.spots[get_play(index)].x = player.x;
}