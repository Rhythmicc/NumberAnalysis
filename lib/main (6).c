#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void printmat(float *A, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            printf("%f ", A[i * n + j]);
        printf("\n");
    }
}

void printmat_trans(float *A, int n)
{
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
            printf("%f ", A[i * n + j]);
        printf("\n");
    }
}

void printvec(float *x, int n)
{
    for (int i = 0; i < n; i++)
        printf("%f\n", x[i]);
}

void matvec(float *A, float *x, float *y, int n)
{
    for (int i = 0; i < n; i++)
    {
        y[i] = 0;
        for (int j = 0; j < n; j++)
            y[i] += A[i * n + j] * x[j];
    }
}

float dotproduct(float *vec1, float *vec2, int n)
{
    float result = 0;
    for (int i = 0; i < n; i++)
        result += vec1[i] * vec2[i];
    return result;
}

void gramschmidt(float *V, float *X, int m, int n)
{
    //X[0 * m + 0] = 3;
    //X[1 * m + 1] = -1;
    //X[1 * m + 0] = -1;
    //X[1 * m + 1] = 2;
    // m vectors of length n
    float *sum = (float *)malloc(sizeof(float) * n);

    for (int i = 0; i < m; i++)
    {
        float *xi = &X[i * n];
        float *vi = &V[i * n];
        memset(sum, 0, sizeof(float) * n);

        for (int j = 0; j < i; j++)
        {
            float *vj = &V[j * n];
            float dp1 = dotproduct(xi, vj, n);
            float dp2 = dotproduct(vj, vj, n);

            for (int k = 0; k < n; k++)
                sum[k] += (dp1 / dp2) * vj[k];
        }

        for (int k = 0; k < n; k++)
            vi[k] = xi[k] - sum[k];
    }

    free(sum);
}

// A is m x n, AT is n x m
void transpose(float *AT, float *A, int m, int n)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            AT[j * m + i] = A[i * n + j];
}

float vec2norm(float *x, int n)
{
    float sum = 0;
    for (int i = 0; i < n; i++)
        sum += x[i] * x[i];
    return sqrt(sum);
}

void matmat(float *C, float *A, float *B, int m, int k, int n)
{
    memset(C, 0, sizeof(float) * m * n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            for (int kk = 0; kk < k; kk++)
                C[i * n + j] += A[i * k + kk] * B[kk * n + j];
}

void matmat_transA(float *C, float *AT, float *B, int m, int k, int n)
{
    memset(C, 0, sizeof(float) * m * n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            for (int kk = 0; kk < k; kk++)
                C[i * n + j] += AT[kk * m + i] * B[kk * n + j];
}

void matmat_transB(float *C, float *A, float *BT, int m, int k, int n)
{
    memset(C, 0, sizeof(float) * m * n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            for (int kk = 0; kk < k; kk++)
                C[i * n + j] += A[i * k + kk] * BT[j * k + kk];
}

void qr(float *Q, float *R, float *A, int m, int n)
{
    gramschmidt(Q, A, n, n);

    for (int i = 0; i < m; i++)
    {
        float *qi = &Q[i * m];
        float norm = vec2norm(qi, n);

        for (int j = 0; j < n; j++)
            qi[j] = qi[j] / norm;
    }

    float *QT = (float *)malloc(sizeof(float) * m * n);
    transpose(QT, Q, m, n);

    matmat_transA(R, QT, A, n, m, n);

    free(QT);
}

void qreigensolver(float *eigenvalue, float *eigenvector, float *A, int n, int maxiter, float threshold)
{
    float *Ai = (float *)malloc(sizeof(float) * n * n);
    float *Ai_new = (float *)malloc(sizeof(float) * n * n);
    float *Q = (float *)malloc(sizeof(float) * n * n);
    float *R = (float *)malloc(sizeof(float) * n * n);
    float *error = (float *)malloc(sizeof(float) * n);
    
    float *Qaccum = (float *)malloc(sizeof(float) * n * n);
    memset(Qaccum, 0, sizeof(float) * n * n);
    for (int i = 0; i < n; i++)
        Qaccum[i * n + i] = 1;

    memcpy(Ai, A, sizeof(float) * n * n);

    int iter = 0;
    float errornorm = 0;
    do
    {
        //printf("\niter = %i\n", iter);
        memset(Q, 0, sizeof(float) * n * n);
        memset(R, 0, sizeof(float) * n * n);

        qr(Q, R, Ai, n, n);
        matmat_transB(eigenvector, Qaccum, Q, n, n, n);
        //matmat(eigenvector, Q, Qaccum, n, n, n);
        memcpy(Qaccum, eigenvector, sizeof(float) * n * n);
        
        matmat_transB(Ai_new, R, Q, n, n, n);
        for (int i = 0; i < n; i++)
            error[i] = Ai_new[i * n + i] - Ai[i * n + i];
        errornorm = vec2norm(error, n);
        //printf("errornorm = %f\n", errornorm);

        memcpy(Ai, Ai_new, sizeof(float) * n * n);
        iter++;
    }
    while (iter < maxiter && errornorm > threshold);

    for (int i = 0; i < n; i++)
        eigenvalue[i] = Ai_new[i * n + i];
    
    //printf("Ai_new = \n");
    //printmat(Ai_new, n);

    free(Ai);
    free(Ai_new);
    free(Q);
    free(R);
    free(error);
}

void tridiagonalization(float *A, float *T, float *V, int n)
{
    memset(T, 0, sizeof(float) * n * n);

    float *va = (float *)malloc(sizeof(float) * n);
    float *vb = (float *)malloc(sizeof(float) * n);
    float *wp = (float *)malloc(sizeof(float) * n);
    float *w = (float *)malloc(sizeof(float) * n);
    memset(va, 0, sizeof(float) * n);
    va[0] = 1.0;

    matvec(A, va, wp, n);
    float alpha = dotproduct(wp, va, n);
    T[0 * n + 0] = alpha;
    for (int i = 0; i < n; i++)
        w[i] = wp[i] - alpha * va[i];
    for (int i = 0; i < n; i++)
        V[i * n + 0] = va[i];

    for (int j = 1; j < n; j++)
    {
        float beta = vec2norm(w, n);
        for (int i = 0; i < n; i++)
            vb[i] = w[i] / beta;
        matvec(A, vb, wp, n);
        float alpha = dotproduct(wp, vb, n);
        T[j * n + j] = alpha;
        for (int i = 0; i < n; i++)
            w[i] = wp[i] - alpha * vb[i] - beta * va[i];

        for (int i = 0; i < n; i++)
            V[i * n + j] = vb[i];
        memcpy(vb, va, sizeof(float) * n);
        T[(j - 1) * n + j] = beta;
        T[j * n + (j - 1)] = beta;
    }

    printf("\nT = \n");
    printmat(T, n);

    printf("\nV = \n");
    printmat(V, n);

    free(va);
    free(vb);
    free(wp);
    free(w);
}

int main(int argc, char **argv)
{
    // method: power
    char *method = argv[1];
    printf("\n");

    char *filename = argv[2];
    printf ("filename = %s\n", filename);

    int n;

    FILE *file;
    file = fopen(filename, "r+");
    fscanf(file, "%i", &n);
    printf("The order of the matrix is %i\n", n);
    
    float *A = (float *)malloc(sizeof(float) * n * n);


    // read A
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            fscanf(file, "%f", &A[i * n + j]);

    fclose(file);

    // print A
    printf("\nA = \n");
    printmat(A, n);

    if (strcmp(method, "qrsvd") == 0)
    {
        printf ("\nQR method for finding eigenvalues.\n");
        float *eigenvalues = (float *)malloc(sizeof(float) * n);
        float *eigenvector = (float *)malloc(sizeof(float) * n * n);

        qreigensolver(eigenvalues, eigenvector, A, n, 20, 0.000001);

        printf("\neigenvalues = \n");
        printvec(eigenvalues, n);
        
        printf("\neigenvectors = \n");
        printmat(eigenvector, n);

        free(eigenvalues);
        free(eigenvector);
    }
    else if (strcmp(method, "lanczos") == 0)
    {
        printf ("\nLanczos method for SVD.\n");

        float *T = (float *)malloc(sizeof(float) * n * n);
        float *V = (float *)malloc(sizeof(float) * n * n);
        float *VT = (float *)malloc(sizeof(float) * n * n);

        tridiagonalization(A, T, V, n);

        float *eigenvalues = (float *)malloc(sizeof(float) * n);
        float *eigenvector = (float *)malloc(sizeof(float) * n * n);

        qreigensolver(eigenvalues, eigenvector, T, n, 20, 0.000001);

        printf("\neigenvalues = \n");
        printvec(eigenvalues, n);

        printf("\neigenvectors = \n");
        printmat(eigenvector, n);

        float *F = (float *)malloc(sizeof(float) * n * n);
        float *FT = (float *)malloc(sizeof(float) * n * n);
        float *Ap = (float *)malloc(sizeof(float) * n * n);
        float *golden = (float *)malloc(sizeof(float) * n * n);
        float *EV = (float *)malloc(sizeof(float) * n * n);
        memset(EV, 0, sizeof(float) * n * n);
        for (int i = 0; i < n-1; i++)
            EV[i * n + i] = eigenvalues[i];

        matmat(F, V, eigenvector, n, n, n);
        transpose(FT, F, n, n);
        matmat(Ap, F, EV, n, n, n);
        matmat(golden, Ap, FT, n, n, n);

        printf("\nAp = \n");
        printmat(golden, n);

        free(eigenvalues);
        free(eigenvector);

        free(T);
        free(V);
        free(VT);
    }

    free(A);


}
