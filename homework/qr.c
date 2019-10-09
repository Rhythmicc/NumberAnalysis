#pragma GCC optimize("O3")
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define pos(i,j) (i)*n+(j)
#define max(a,b) (a)>(b)?(a):(b)
#define min(a,b) (a)<(b)?(a):(b)
#define swap(a,b) tmp = *(a), *(a)=*(b),*(b)=tmp
#define MALLOC(p,type,n) type*p = (type*)malloc(sizeof(type)*n)


void matvec(const float *A, const float *x, float *y, int n) {
    for (int i = 0; i < n; ++i) {
        y[i] = 0;
        for (int j = 0; j < n; ++j)y[i] += A[pos(i, j)] * x[j];
    }
}

float dotproduct(const float *a, const float *b, int n) {
    float res = 0;
    for (int i = 0; i < n; ++i)res += a[i] * b[i];
    return res;
}

float vec2norm(const float *x, int n) {
    float sum = 0;
    for (int i = 0; i < n; ++i)sum += x[i] * x[i];
    return (float)sqrt((double)sum);
}

void power(float*eigenvalues, float *A, int n, int maxiter, float threshold) {
    size_t sz = sizeof(float) * n;
    float*eigenvector = (float*)malloc(sz);
    float*eigenvector_new = (float*)malloc(sz);
    for(int i=0;i<n;++i)eigenvector[i] = 1;
    float evalue = eigenvector[0];
    float evalue_new;
    int iter = 0;
    float error = 0;
    do{
        matvec(A, eigenvector, eigenvector_new, n);
        evalue_new = eigenvector_new[0];
        for(int i=0;i<n;++i)eigenvector_new[i] /= evalue_new;
        memcpy(eigenvector, eigenvector_new, sz);
        error = fabsf((evalue_new - evalue) / evalue_new);
        evalue = evalue_new;
        ++iter;
    }while(iter < maxiter && error > threshold);
    eigenvalues[0] = evalue;
    free(eigenvector);
}

void gramschmidt(float*V,float*X,int m,int n){
    MALLOC(sum,float,n);
    for(int i=0;i<m;++i){
        float*vi = V + pos(i,0);
        float*xi = X + pos(i,0);
        memset(sum,0, sizeof(float)*n);
        for(int j=0;j<i;++j){
            float*vj = V + pos(j,0);
            float dp1 = dotproduct(xi,vj,n);
            float dp2 = dotproduct(vj,vj,n);
            for(int k=0;k<n;++k){
                sum[k] += (dp1/dp2)*vj[k];
            }
        }
        for(int k=0;k<n;++k)vi[k] = xi[k] - sum[k];
    }
    free(sum);
}

void T(float*AT,float*A,int m,int n){
    for(int i=0;i<m;++i)for(int j=0;j<n;++j)AT[j*m+i] = A[i*n+j];
}

void matmat_transA(float*C,float*AT,float*B,int m,int n,int k){
    memset(C,0, sizeof(float)*n*m);
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int l=0;l<k;++l)C[pos(i,j)] += AT[l*m+i] * B[l*n+j];
}

void matmat_transB(float*C,float*A,float*BT,int m,int n,int k) {
    memset(C,0, sizeof(float)*n*m);
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int l=0;l<k;++l)C[pos(i,j)] += A[i*k+l] * BT[j*k+l];
}

void qr(float*Q,float*R,float*A,int m,int n) {
    gramschmidt(Q, A, m, n);
    for (int i = 0; i < n; ++i) {
        float *qi = Q + pos(i, 0);
        float norm = vec2norm(qi, n);
        for (int j = 0; j < n; ++j)qi[j] /= norm;
    }
    MALLOC(QT, float, m * n);
    T(QT, Q, m, n);
    matmat_transA(R, QT, A, m, n, n);
    free(QT);
}

void qreigensolver(float*eigenvalues, float *A, int n, int maxiter, float threshold) {
    MALLOC(Ai, float, n * n);
    MALLOC(Ai_new, float, n * n);
    MALLOC(Q, float, n * n);
    MALLOC(R, float, n * n);
    MALLOC(error, float, n);
    memcpy(Ai, A, sizeof(float) * n * n);
    int iter = 0;
    float errornorm = 0;
    do {
        qr(Q, R, Ai, n, n);
        matmat_transB(Ai_new, R, Q, n, n, n);
        for (int i = 0; i < n; ++i)error[i] = Ai_new[i * n + i] - Ai[i * n + i];
        errornorm = vec2norm(error, n);
        memcpy(Ai, Ai_new, sizeof(float) * n * n);
        ++iter;
    } while (iter < maxiter && errornorm > threshold);
    for (int i = 0; i < n; ++i)eigenvalues[i] = Ai_new[pos(i, i)];
    free(Ai), free(Ai_new), free(Q), free(R), free(error);
}

void (*func[])(float*eigenvalues, float *A, int n, int maxiter, float threshold) = {power,qreigensolver};

int main(int argc, char **argv) {
    int mode = argv[1][0] - '1';
    int n;
    scanf("%d",&n);
    float*A = (float*)malloc(sizeof(float)*n*n);
    float *eigenvalues = (float*)malloc(sizeof(float)*n);
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)scanf("%f",A+pos(i,j));
    func[mode](eigenvalues,A,n,1000,1e-5f);
    for(int i=0;i<n;++i)printf("%f%c",eigenvalues[i],i==n-1?'\n':' ');ß
    return 0;
}
