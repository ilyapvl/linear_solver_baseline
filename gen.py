#!/usr/bin/env python3
from pathlib import Path


def generate(n: int = 20,
             c: float = 10.0,
             b_val: float = 1.0,
             d_val: float = 1.0,
             a_val: float = 0.1,
             e_val: float = 0.1,
             outfile: str = "input.txt") -> None:
    

    A = [[0.0] * n for _ in range(n)]
    rhs = [0.0] * n

    for i in range(1, n + 1):           # i — 1-based, как в условии
        row = i - 1                     # индекс строки в массиве
        A[row][row] = c                 # c_i на главной диагонали

        if i >= 2:
            A[row][row - 1] = b_val     # b_i (поддиагональ)
        if i <= n - 1:
            A[row][row + 1] = d_val     # d_i (наддиагональ)
        if i >= 3:
            A[row][row - 2] = a_val     # a_i (вторая поддиагональ)
        if i <= n - 2:
            A[row][row + 2] = e_val     # e_i (вторая наддиагональ)

        rhs[row] = float(i)             # f_i = i

    path = Path(outfile)
    with path.open("w") as f:
            f.write(f"{n}\n")
            for row in A:
                f.write(" ".join(f"{v:g}" for v in row) + "\n")
            f.write(" ".join(f"{v:g}" for v in rhs) + "\n")

    print(f"ok")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("-n", "--size", type=int, default=20)
    p.add_argument("-o", "--output", default="input.txt")
    args = p.parse_args()

    generate(n=args.size, outfile=args.output)
