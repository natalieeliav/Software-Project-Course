import numpy as np
import pandas as pd
import sys
import mysymnmf
import random
np.random.seed(0)

def print_matrix(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[i]) - 1):
            print("%.4f" % matrix[i][j], end=",")
        print("%.4f" % matrix[i][-1])


def init_H(n, k, m):
    upper_limit = 2 * np.sqrt(m / k)
    matrix = np.random.uniform(0, upper_limit , size=(n, k))
    return matrix.tolist()


def matrix_avg(matrix) :
    total_sum = np.sum(matrix)
    total_elements = np.size(matrix)
    return total_sum / total_elements
    
def init_points(file_name) :
    try:
        df = pd.read_csv(file_name, header=None)
    except FileNotFoundError:
        print("Invalid input!")
        return

    points = df.values
    n = len(points)
    d = len(points[0])
    points = points.tolist()
    return points, n, d

def final_H(points, n, d, k) :
        norm = mysymnmf.norm(points, n, d, k)
        m = matrix_avg(np.array(norm))
        h = init_H(n, k, m)
        result = mysymnmf.symnmf(points, h, n, d, k)
        return result


def main():
    file_name = sys.argv[3]
    k = int(sys.argv[1])
    points, n, d = init_points(file_name)

    result = None
    if sys.argv[2] == "symnmf":
        print_matrix(final_H(points, n, d, k))
    elif sys.argv[2] == "sym":
        result = mysymnmf.sym(points, n, d, k)
    elif sys.argv[2] == "ddg":
        result = mysymnmf.ddg(points, n, d, k)
    elif sys.argv[2] == "norm":
        result = mysymnmf.norm(points, n, d, k)
        print_matrix(result)
    return

if __name__ == "__main__":
    main()
