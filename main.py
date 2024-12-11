import argparse
from gurobipy import Model, GRB


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('n', type=int)
    parser.add_argument('--y', type=int, choices=[0, 1], default=0,
                        help='Element y for which we attempt 1 * y = 2 * y = z')
    parser.add_argument('--z', type=int, choices=[0, 1, 3], default=0,
                        help='Element z for which we attempt 1 * y = 2 * y = z')

    args = parser.parse_args()

    # Set the dimension size
    n = args.n

    # Create a new model
    m = Model("Equation677")

    # Create variables x[i,j,k] for i,j,k in {0, ..., n-1}
    # x[i,j,k] = 1 means that i * j = k
    v = m.addVars(n, n, n, vtype=GRB.BINARY, name="v")

    # Add constraints: For every i, j, sum of x[i,j,k] over k is 1
    m.addConstrs(
        (sum(v[i,j,k] for k in range(n)) == 1 for i in range(n) for j in range(n)),
        name="row_col_sum"
    )

    # Add constraints: left multiplication is invertible
    m.addConstrs(
        sum(v[i,j,k] for j in range(n)) == 1
        for i in range(n) for k in range(n)
    )

    # Add the constraint for equation 677:
    # For every x, y, a, b, c in {0, ..., n-1},
    # v[y,c,x] >= v[y,x,a] + v[a,y,b] + v[x,b,c] - 2
    for x in range(n):
        for y in range(n):
            for a in range(n):
                for b in range(n):
                    for c in range(n):
                        m.addConstr(
                            v[y,c,x] >= v[y,x,a] + v[a,y,b] + v[x,b,c] - 2,
                            name=f"eq_677_{x}_{y}_{a}_{b}_{c}"
                        )

    # (Optional) Set an objective, for example, maximize the sum of all variables
    # m.setObjective(sum(v[i,j,k] for i in range(n) for j in range(n) for k in range(n)), GRB.MAXIMIZE)

    # Right multiplication is not injective:
    # There exist y and z1 != z2 such that z1 * y = z2 * y
    # Wlog z1 = 1, z2 = 2, and y in {0, 1}

    # m.addConstrs(
    #     v[1, y_special, k] == v[2, y_special, k]
    #     for k in range(n)
    # )

    m.setObjective(v[1, args.y, args.z] + v[2, args.y, args.z], GRB.MAXIMIZE)



    # Update the model to ensure variables and constraints are integrated
    m.update()

    # (Optional) Optimize the model
    m.optimize()

    # Print the solution (if feasible)
    if m.status == GRB.OPTIMAL:
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if v[i,j,k].X > 0.5:
                        print(f'{i} * {j} = {k}')
                        # print(f"x[{i},{j},{k}] = {x[i,j,k].X}")


if __name__ == '__main__':
    main()