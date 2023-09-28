#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define beta 0.5
#define epsilon 0.0001
#define max_iter 300

double distance(double* x, double* y, int d);
double similarity(double* x, double* y, int d);
double** get_similarity_matrix(double** points, int n, int d);
double** get_diag_deg_matrix(double** A, int n);
double** negative_root_matrix(double **D, int n);
double** matrixMultiplication(int m, int n, int l, double** A, double** B);
double** get_W(double** D, double** A, int n);
void update_H(double** H, int n, int k, double** W);
int check_convergence(double** A, double** B,int n,int k);
double** transpose(double** A, int n, int k);
double** create_matrix(int n, int d);
void free_matrix(double** matrix, int n);
void print_matrix(double** matrix, int n, int k);
void copy_matrix_inplace(double** A, double** B, int n, int k);
void get_dimensions(FILE* file, int* dim);
void get_datapoints(FILE* file, double** datapoints, int n, int d);


/* returns the distance between two vectors*/
double distance(double* x, double* y, int d){
    double sum;
    int i;
    sum = 0;
    for (i = 0 ; i < d ; i++){
        sum += ((x[i] - y[i]) * (x[i] - y[i]));
    }
    return sum;
}

/* calculats similarity */
double similarity(double* x, double* y, int d){
    return exp( (-0.5) * (distance(x,y,d)) );
}

/* Returns the similarity matrix */
double** get_similarity_matrix(double** points, int n, int d){
    double** A;
    int i,j;
    A = create_matrix(n,n);
    for (i = 0 ; i < n ; i++){
        for (j = 0 ; j < n ; j++){
            if (i == j){
                A[i][j] = 0;
            }
            else{
                A[i][j] = similarity(points[i],points[j],d);
            }
        }
    }
    return A;
}

/* Calculates nodes degree matrix */
double** get_diag_deg_matrix(double** A, int n){
    double** D;
    int i, j;
    D = create_matrix(n,n);
    for (i = 0 ; i < n ; i++){
        for (j = 0 ; j < n ; j++){
            D[i][i] += A[i][j];
        }
    }
    return D;
}

/* returns a negative squared root for given matrix */
double** negative_root_matrix(double **D, int n){
    double** root;
    int i;
    root = create_matrix(n, n);
    for (i = 0; i < n; i++){
        root[i][i] = 1 / sqrt(D[i][i]);
    }
    return root;
}

/* Matrix Multiplication,ret value = A(m*n) * B(n*l) */
double** matrixMultiplication(int m, int n, int l, double** A, double** B) {
    double** C;
    int i, j, k;
    C = create_matrix(m, l);
    for (i = 0; i < m; i++) {
        for (j = 0; j < l; j++) {
            for (k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j]; 
            }
        }
    }
    return C;
}


/* return normalized sim matrix W  */
double** get_W(double** D, double** A, int n){
    double **D_root, **W, **tmp;
    tmp = create_matrix(n, n);
    D_root = negative_root_matrix(D, n);
    tmp = matrixMultiplication(n, n, n, D_root, A);
    W = matrixMultiplication(n, n, n, tmp, D_root);
    free_matrix(tmp, n);
    free_matrix(D_root, n);
    return W;
}


int check_convergence(double** A, double** B,int n,int k){
    double sum, curr;
    int i, j;
    sum = 0.0;
    curr = 0;
    for (i = 0; i < n; i++){
        for (j = 0; j < k; j++){
            curr = pow((A[i][j]-B[i][j]),2);
            sum += curr;
        }
    }


    if (sum < 0.0001){
        return 1;
    }
    return 0;
}



/* Transpose a Matrix inplace*/
double** transpose(double** A, int n, int k) {
    int i,j;
    double** transposed;
    transposed = create_matrix(k, n);
    for(i = 0; i < n; i++) {
        for(j = 0; j < k; j++) {
            transposed[j][i] = A[i][j];
        }
    }
    return transposed;
}

/* Creates a n X n size matrix with zeroes */
double** create_matrix(int n, int d) {
    double** mat;
    int i;
    mat = (double**)calloc(n, sizeof(double*));
    if (mat == NULL){
        printf("An Error Has Occurred");
        return NULL;
    }
    for (i = 0; i < n; i++){
        mat[i] = (double*)calloc(d, sizeof(double));
        if (mat[i] == NULL){
            printf("An Error Has Occurred");
            return NULL;
        }
    }
    return mat;
}
/* Free space for matrix of length n */
void free_matrix(double** matrix, int n) {
    int i;
    for (i=0; i<n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

/* Prints matrix of length n*/

void print_matrix(double** matrix, int n, int k) {
    int i,j;
    for (i = 0; i < n; ++i) {
        printf("%.4f", matrix[i][0]);
        for (j = 1; j < k; j++) {
            printf(",%.4f", matrix[i][j]);
        }
        printf("\n");
    }

}

/* gets mats A and B, and copys B to A */
void copy_matrix_inplace(double** A, double** B, int n, int k) {
    int i,j;
    for(i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            A[i][j] = B[i][j];
        }
    }
}

/* gets file with matrix, returns array of dimensions */
void get_dimensions(FILE* file, int* dim) {
    double a;
    char c;
    int n = 0, d = 0;
    while (fscanf(file ,"%lf",&a) == 1){
        c = fgetc(file);
        if (c == ',') {
            d++;
        } else {
            n++;
        }
    }
    d = d/n + 1;
    dim[0] = n;
    dim[1] = d;
}

/* reads matrix in file and writes the matrix into datapoints */
void get_datapoints(FILE* file, double** datapoints, int n, int d) {
    int i,j;
    int check;
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            check = fscanf(file,"%lf",&datapoints[i][j]);
            fgetc(file);
            if (check != 1){
                return ;
            }
        }
    }
}


void update_H(double** H, int n, int k, double** W){
    double **wh, **hhh, **hh;
    int i, j;
    double tmp1, tmp2, curr;
    tmp1 = 0;
    tmp2 = 0;
    curr = 0;
    wh = create_matrix(n, k);
    hh = create_matrix(n, n);
    hhh = create_matrix(n, k);
    wh = matrixMultiplication(n, n, k, W, H);
    hh = matrixMultiplication(n, k, n, H, transpose(H, n, k));
    hhh = matrixMultiplication(n, n, k, hh, H);
    for (i = 0; i < n; i++){
        for (j = 0; j < k; j ++){
            tmp1 = wh[i][j] / hhh[i][j];
            tmp2 = beta + tmp1*beta;
            curr = H[i][j];
            H[i][j] = curr * tmp2;
        }
    }
    free_matrix(wh, n);
    free_matrix(hhh, n);
    free_matrix(hh, n);

}

int main(int argc, char* argv[]) {
    FILE* file;
    double** points;
    int n, d, args_amt;
    char* goal;
    double** sym_mat;
    double** ddg;
    double** W;
    int* dim = calloc(sizeof(int), 2);
    args_amt = argc;
    if (args_amt != 3){
        return 0;
    }
    if (dim == NULL){
        printf("An Error Has Occurred1");
        return 0;
    }

    goal = argv[1];
    if (strcmp(goal, "sym") && strcmp(goal, "ddg") && strcmp(goal, "norm")) {
        printf("Invalid Input!");
        return 0;
    }

    file = fopen(argv[2], "r");
    if (file == NULL){
        printf("An Error Has Occurred");
        return 0;
    }
    get_dimensions(file, dim);
    n = dim[0];
    d = dim[1];

    points = create_matrix(n, d);
    fseek(file, 0 , 0);
    get_datapoints(file, points, n, d);
    fclose(file);
    sym_mat = get_similarity_matrix(points, n, d);



    if (!strcmp(goal, "sym")) {
        print_matrix(sym_mat, n, n);
        free_matrix(sym_mat, n);
    } else if (!strcmp(goal, "ddg")) {
        ddg = get_diag_deg_matrix(sym_mat, n);
        print_matrix(ddg, n, n);
        free_matrix(sym_mat, n);
        free_matrix(ddg, n);
    } else if (!strcmp(goal, "norm")) {
        ddg = get_diag_deg_matrix(sym_mat, n);
        W = get_W(ddg, sym_mat, n);
        print_matrix(W, n, n);
        free_matrix(sym_mat, n);
        free_matrix(ddg, n);
        free_matrix(W, n);
    } 
    free(dim);
    free_matrix(points,n);
    return 1;
}
