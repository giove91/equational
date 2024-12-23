import argparse
import ast

from gurobipy import Model, GRB


# def int_to_bit_tuple(n):
#     # Ensure the input is between 0 and 15
#     if not (0 <= n <= 15):
#         raise ValueError("Input must be an integer between 0 and 15 inclusive.")
#     return tuple(reversed(((n >> 3) & 1, (n >> 2) & 1, (n >> 1) & 1, n & 1)))


def compute_inverse_product(x, y):
    """
    Compute x^{-1}y using the multiplication table of a group.

    :param x: First element (0-based indexing).
    :param y: Second element (0-based indexing).
    :return: Index of the result (0-based indexing).
    """
    return x ^ y


def construct_L_y(diamond_table, y):
    return diamond_table[y]


def compose_permutations(perm1, perm2):
    """
    Compose two permutations.

    :param perm1: List[int], the first permutation
    :param perm2: List[int], the second permutation
    :return: List[int], the composed permutation
    """
    n = len(perm1)
    return [perm1[perm2[i]] for i in range(n)]


def inverse_permutation(perm):
    """
    Compute the inverse of a given permutation with a validity check.

    :param perm: List[int], the permutation to invert (a list of size n where each number in {0, ..., n-1} appears exactly once)
    :return: List[int], the inverse permutation
    :raises ValueError: If the input is not a valid permutation
    """
    n = len(perm)

    # Check if the input is a valid permutation
    if sorted(perm) != list(range(n)):
        raise ValueError("The input is not a valid permutation.")

    # Compute the inverse
    inv = [0] * n
    for i in range(n):
        inv[perm[i]] = i
    return inv


def check_pairwise_commutativity(operators):
    """
    Check if all operators pairwise commute.

    :param operators: List[List[int]], a list of permutations (operators)
    :return: bool, True if all operators commute, False otherwise
    """
    n = len(operators)
    for i in range(n):
        for j in range(i + 1, n):
            # Compose operators in both orders
            L1_L2 = compose_permutations(operators[i], operators[j])
            L2_L1 = compose_permutations(operators[j], operators[i])
            # Check if they are equal
            if L1_L2 != L2_L1:
                print(f"Operators at index {i} and {j} do not commute.")
                return False
    return True


def check_linear(diamond_table):
    """
    Linearity test by Terence Tao:
        One consequence of linearity is that the operators LyLy′−1, being translations,
        commute with each other and with the RxRx′−1 (assuming right invertibility in the latter case).
        So this may be one test to verify if a given magma is linear or not.
    """
    n = len(diamond_table)
    operators = []
    for y1 in range(n):
        for y2 in range(n):
            operators.append(
                compose_permutations(
                    construct_L_y(diamond_table, y1),
                    inverse_permutation(construct_L_y(diamond_table, y2)),
                )
            )

    operators = list(set(tuple(op) for op in operators))

    print(len(operators))
    for operator in operators:
        print(operator)

    return check_pairwise_commutativity(operators)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--derangement', action='store_true')
    parser.add_argument('-s', '--num_solutions', type=int, default=None)

    args = parser.parse_args()

    # Set the dimension size
    n = 16

    # Create a new model
    m = Model("Equation677_translation_invariant")

    # Create variables w[i,j] for i,j in {0, ..., n-1}
    # w[i,j] = 1 means that f(i) = j
    w = m.addVars(n, n, vtype=GRB.BINARY, name="w")

    # Add constraints: For every i, sum of w[i,j] over j is 1
    m.addConstrs(
        (sum(w[i,j] for j in range(n)) == 1 for i in range(n)),
        name="row_col_sum"
    )

    # Add constraints: For every j, sum of w[i,j] over i is 1
    m.addConstrs(
        (sum(w[i,j] for i in range(n)) == 1 for j in range(n)),
        name="row_col_sum"
    )

    if args.derangement:
        # Derangement
        # m.addConstrs(
        #     w[i,i] == 0 for i in range(n)
        # )

        # Impose f(0) = 1 (still no solutions, even if weaker)
        m.addConstr(w[0, 1] == 1)

    # Add the constraint for equation 677:
    # For every x, y, a, b, c in {0, ..., n-1},
    # v[y,c,x] >= v[y,x,a] + v[a,y,b] + v[x,b,c] - 2
    # where v[i,j,k] = w[j-i, k-i], assuming a group structure on {0, ..., n-1}
    # (we try Z/nZ to start...)

    def get_v(ii, jj, kk):
        """
        Return the variable that tells us if ii diamond jj = kk in the magma, i.e.,
        kk = ii f(ii^{-1} jj), i.e.,
        ii^{-1} kk = f(ii^{-1} jj)
        """
        # return w[(jj-ii) % n, (kk-ii) % n]
        return w[
            compute_inverse_product(ii, jj),
            compute_inverse_product(ii, kk)
        ]

    for x in range(n):
        for y in range(n):
            for a in range(n):
                for b in range(n):
                    for c in range(n):
                        m.addConstr(
                            get_v(y,c,x) >= get_v(y,x,a) + get_v(a,y,b) + get_v(x,b,c) - 2,
                            name=f"eq_677_{x}_{y}_{a}_{b}_{c}"
                        )

    m.setObjective(sum(w[i, i] for i in range(n)), GRB.MINIMIZE)



    # Update the model to ensure variables and constraints are integrated
    m.update()

    if args.num_solutions is not None:
        m.setParam('PoolSearchMode', 2)
        m.setParam('PoolSolutions', args.num_solutions)

    # (Optional) Optimize the model
    m.optimize()

    print(f'Number of solutions: {m.SolCount}')

    # Print the solution (if feasible)
    if m.status == GRB.OPTIMAL:
        diamond_table = []

        for i in range(n):
            row = []
            for j in range(n):
                for k in range(n):
                    if get_v(i,j,k).X > 0.5:
                        # print(f'{i} * {j} = {k}')
                        # print(f"x[{i},{j},{k}] = {x[i,j,k].X}")
                        row.append(k)

            print(row)
            diamond_table.append(row)

        # print(diamond_table)

        for i in range(n):
            for j in range(n):
                if w[i, j].X > 0.5:
                    print(f'f({i}) = {j}')

                # print(f'w[{i}, {j}] = {w[i, j].X}')

        print(f'Linear check: {check_linear(diamond_table)}')

        # # Go through all solutions
        # for i in range(m.SolCount):
        #     m.setParam(GRB.Param.SolutionNumber, i)
        #     print(f'Solution #{i}')


if __name__ == '__main__':
    main()