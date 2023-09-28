# Software Project Course - symNMF algorithm

The project revolves around the implementation of a clustering algorithm known as symmetric Non-negative Matrix Factorization (symNMF).
The primary goal is to apply this algorithm to various datasets and compare its performance with Kmeans clustering.
SymNMF is utilized to transform a set of data points into a similarity matrix, which is subsequently used to cluster the data into a specified number of clusters (k).

contains the following files:
1. symnmf.py: Python interface of the code.
2. symnmf.h: C header file.
3. symnmf.c: C interface of the code.
4. symnmfmodule.c: Python C API wrapper.
5. analysis.py: Analyze of the algorithm.
6. setup.py: The setup file.
7. Makefile: script to build the C interface.

The project was done as part of the course "software project" in Tel Aviv university.

# How It Works

SymNMF begins by constructing a similarity matrix (A) based on the Euclidean distances between data points.
The degree matrix (D) is calculated, representing the degrees of each vertex in the graph.
Then, a normalized similarity matrix (W) is derived.
The core of the algorithm lies in optimizing the matrix H, which is initialized randomly and updated iteratively.
The convergence of H is achieved when a maximum iteration number is reached or a specific threshold is met.
This process effectively identifies clusters within the data.

# Usage

If you want to utilize this project to perform clustering or similarity analysis on your data, follow these steps:

1. Parameter Setup: Specify the number of clusters (k) you want to identify and choose the goal of your analysis:
    If you want to perform the full symNMF clustering algorithm, set the goal as 'symnmf'.
    If you need to calculate the similarity matrix, set the goal as 'sym'.
    To obtain the Diagonal Degree Matrix, use 'ddg' as your goal.
    If you're interested in the normalized similarity matrix, set the goal to 'norm'.

2. Prepare Input Data: Organize your data points into an input file with a '.txt' format.
   This file will serve as the source data for your analysis.

3. Run the Program: Execute the 'symnmf.py' Python program, providing the required arguments, including the chosen goal and the name of your input data file.

4. Retrieve Results: The program will generate output based on your specified goal.
   You will receive matrices, such as the similarity matrix, diagonal degree matrix, or clustering matrix H, depending on your chosen goal.
   These matrices will be presented in a comma-separated format.

5. For further analysis and comparison of clustering results, consider using 'analysis.py.'
   It allows you to assess the performance of SymNMF versus Kmeans on your dataset using the silhouette score.



