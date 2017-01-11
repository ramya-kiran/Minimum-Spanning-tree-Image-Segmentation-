import numpy as np


class Board:
    # initializing the board with a stack to keep track of the remaining moves
    # and a list to keep the board state.
    def __init__(self):
        self.stack = []
        for i in range(7):
            self.stack.append([])
            for j in range(6):
                self.stack[i].append([j, i])
        self.lst = np.zeros((6, 7))
        self.turn = 1
        self.value = None
        self.list_moves = []
        self.count = 0


# Generates all the allowed possible moves using the stack.
    def generate_moves(self):
        place = []
        for i in range(7):
            if self.stack[i]:
                place.append(self.stack[i][-1][1])
        return place


# makes a move which is given to this function, and changes the turn of
# the player.
    def make_move(self, move):
        # implement allowed moved concept
        if self.stack:
            row, col = self.stack[move].pop()
            self.lst[row, col] = self.turn
            self.list_moves.append([row, col])
            self.turn = -self.turn
        else:
            pass


# unmake a move that has been made, and changes the board states to the
# previous move state.
    def unmake_last_move(self):
        if self.list_moves:
            r, c = self.list_moves[-1]
            self.lst[r, c] = 0
            self.stack[c].extend([[r, c]])
            self.turn = -self.turn
            self.list_moves.pop()
        else:
            pass


# This function checks if the last move won. It generates the horizontal, vertical, and diagonal elements and scans it
# to find 4 plays together. 
    def last_move_won(self):
        if not self.list_moves:
            return False
        else:
            r, c = self.list_moves[-1]
            play = self.lst[r, c]
            # horizontal +ve direction
            hx1 = self.scan(play, r, c, 0, 1)
            # horizontal -ve direction
            hx2 = self.scan(play, r, c, 0, -1)
            # vertical +ve direction
            vy1 = self.scan(play, r, c, -1, 0)
            # vertical -ve direction
            vy2 = self.scan(play, r, c, 1, 0)
            # diagonal in 45 degrees
            dd1 = self.scan(play, r, c, -1, 1)
            # diagonal in 135 degrees
            dd2 = self.scan(play, r, c, -1, -1)
            # diagonal in 225 degrees
            dd3 = self.scan(play, r, c, 1, -1)
            # diagonal in 315 degrees
            dd4 = self.scan(play, r, c, 1, 1)
            horizontal = hx1 + hx2
            vertical = vy1 + vy2
            diagonal1 = dd1 + dd3
            diagonal2 = dd2 + dd4
            if play == 1:
                if (sum(horizontal) >= 5 * play) or (sum(vertical) >= 5 * play) or (sum(diagonal1) >= 5 * play) \
                        or (sum(diagonal2) >= 5 * play):
                    return True
                else:
                    return False
            elif play == -1:
                if (sum(horizontal) <= 5 * play) or (sum(vertical) <= 5 * play) or (sum(diagonal1) <= 5 * play) \
                        or (sum(diagonal2) <= 5 * play):
                    return True
                else:
                    return False
            else:
                print("No valid player")

# This function is to just display the board state
    def __str__(self):
        # return self.lst, self.stack
        return str(self.lst) + ":" + str(self.stack)

# scanning is used by last_move_won function.
    def scan(self, play, x, y, dx, dy):
        checking = []
        while (0 <= x < 6) and (0 <= y < 7) and self.lst[x, y] == play:
            checking.append(self.lst[x, y])
            x = x + dx
            y = y + dy
        return checking
