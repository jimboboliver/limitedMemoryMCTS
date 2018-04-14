#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ITERATIONS 1000

typedef struct {
    unsigned x:2; // 0 = nothing, 1 = X, 2 = O
} Spot;

typedef struct {
    Spot spots[9];
} Board;

typedef struct Node {
    float wins; // draw counts for 0.5 therefore float
    int visits;
    Board board;
    int numChild; // Number of child nodes
    struct Node* children;
} Node;

/** Use UCB to select optimal child */
Node* select_child(Node* node) {
    float max_value = 0;
    float new_value;
    Node* best;

    for (int i = 0; i < node->numChild; i++) {
        new_value = node->children[i].wins / node->children[i].visits + sqrt(2 * log(node->visits) / node->children[i].visits);

        if (new_value > max_value) {
            max_value = new_value;
            best = &node->children[i];
        }
    }

    return best;
}

int no_moves_left(Board* board) {
    return board->spots[0].x && board->spots[1].x && board->spots[2].x && board->spots[3].x && board->spots[4].x && board->spots[5].x && board->spots[6].x && board->spots[7].x && board->spots[8].x;
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

int has_unborn(Board* board) {
    for (int i = 0; i < 9; i++) {
        if (!board->spots[i].x) {
            return 1;
        }
    }
    return 0;
}

/** Add a random available move from this leaf as a child and return it to simulate from. */
Node* add_child(Node* node, Spot player) {
    Board options[9];
    int numOptions = 0;

    for (int i = 0; i < 9; i++) {
        if (!node->board.spots[i].x) {
            options[numOptions] = node->board;
            options[numOptions++].spots[i].x = player.x;
        }
    }

    int index = node->numChild;
    // Choose a random move
    if (!node->numChild++) {
        node->children = (Node*) calloc(1, sizeof(Node));
    } else {
        node->children = (Node*) realloc(node->children, sizeof(Node) * node->numChild);
        node->children[index].wins = 0; // Zero everything else
        node->children[index].visits = 0;
        node->children[index].numChild = 0;
        node->children[index].children = 0; // TODO: remove redundant set to null pointer
    }
    node->children[index].board = options[rand() % numOptions]; // Set the board

    return &node->children[index];
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

void deallocate_branch(Node* node) {
    if (node->numChild) {
        for (int i = 0; i < node->numChild; i++) {
            deallocate_branch(&node->children[i]);
        }
    }
    free(node);
}

void print_board(Node* node) {
    printf("%d%d%d\n%d%d%d\n%d%d%d\n", node->board.spots[0].x, node->board.spots[1].x, node->board.spots[2].x, node->board.spots[3].x, node->board.spots[4].x, node->board.spots[5].x, node->board.spots[6].x, node->board.spots[7].x, node->board.spots[8].x);
}

Node* root_changeover(Node* root, int chosen) {
    for (int i = 0; i < root->numChild; i++) { // Deallocate the irrelevant children
        if (i != chosen) {
            deallocate_branch(&root->children[i]);
        }
    }
    Node* oldRoot = root;
    root = &root->children[chosen];
    free(oldRoot); // Deallocate the old root
    return root;
}

void print_tree(Node* root) {
    char* buffer;
}

int main() {

    srand(time(NULL));

    Node* root = (Node*) calloc(1, sizeof(Node));

    while (1) { // Until game ends

        int i;

        for (int it = 0; it < ITERATIONS; it++) { // Until move must be made
            printf("Iteration %d\n", it);
            // Select
            int term;
            Node** visited = (Node**) malloc(sizeof(Node*));
            int numVisited = 1;
            Node* node = root;
            visited[0] = node;
            Spot player = {1}; // Start with AI player
            while (node->numChild && !has_unborn(&node->board)) { // While the current node is not either a leaf or has unadded children
                node = select_child(node);
                visited = (Node**) realloc(visited, sizeof(Node*) * ++numVisited);
                visited[numVisited - 1] = node;
                player.x = player.x ^ 3; // Switch player
            }

            // Expand
            if (!terminal(&node->board)) {
                node = add_child(node, player);
                visited = (Node**) realloc(visited, sizeof(Node*) * ++numVisited);
                visited[numVisited - 1] = node;
                player.x = player.x ^ 3; // Switch player

                // Simulate
                Board simBoard = node->board;
                while (!(term = terminal(&simBoard))) {
                    simBoard = simulate(&simBoard, player);
                    player.x = player.x ^ 3; // Switch player
                }
            }

            // Backpropagate
            if (term == 1) {
                for (i = 0; i < numVisited; i++) {
                    visited[i]->wins += 1;
                    visited[i]->visits += 1;
                }
            } else if (term == 3) {
                for (i = 0; i < numVisited; i++) {
                    visited[i]->wins += 0.5;
                    visited[i]->visits += 1;
                }
            } else {
                for (i = 0; i < numVisited; i++) {
                    visited[i]->visits += 1;
                }
            }
            free(visited);
        }

        int chosen = 0;

        for (i = 1; i < root->numChild; i++) {
            if (root->children[chosen].wins / root->children[chosen].visits < root->children[i].wins / root->children[i].visits) {
                chosen = i;
            }
        }
        root = root_changeover(root, chosen);
        print_board(root);

        if (terminal(&root->board)) {
            printf("Game over\n");
            break;
        }

        printf("Enter play: ");

        int humanPlay = getc(stdin) - 48;

        printf("\n");

        for (i = 0; i < root->numChild; i++) {
            if (root->children[i].board.spots[humanPlay].x == 'O') {
                chosen = i;
                break;
            }
        }
        root = root_changeover(root, chosen);

        print_board(root);

        if (terminal(&root->board)) {
            printf("Game over\n");
            break;
        }
    }
    
    deallocate_branch(root);

    return 0;
}