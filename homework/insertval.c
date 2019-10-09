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

float qpow(float dest, int k) {
    float res = 1;
    while (k) {
        if (k & 1)res *= dest;
        dest *= dest;
        k >>= 1;
    }
    return res;
}

void lu(float*A,float*x,float*b,int n){
    size_t sz = sizeof(float) * 1ll * n * n;
    float *L = (float *)malloc(sz), *U = (float *)malloc(sz), *y = (float *)malloc(sz / n);
    for (int i = 0; i < n; ++i)L[pos(i, i)] = 1;
    for (int i = 0; i < n; ++i) {
        U[pos(0, i)] = A[pos(0, i)];
        L[pos(i, 0)] = A[pos(i, 0)] / U[0];
    }
    for (int k = 1; k < n; ++k) {
        for (int i = k; i < n; ++i) {
            float sum = 0;
            for (int j = 0; j < k; ++j)sum += L[pos(k, j)] * U[pos(j, i)];
            U[pos(k, i)] = A[pos(k, i)] - sum;
        }
        for (int i = k; i < n; ++i) {
            float sum = 0;
            for (int j = 0; j < k; ++j)sum += L[pos(i, j)] * U[pos(j, k)];
            if (U[pos(k, k)] != 0)L[pos(i, k)] = (A[pos(i, k)] - sum) / U[pos(k, k)];
        }
    }
    for (int i = 0; i < n; ++i) {
        float sum = 0;
        for (int j = 0; j < i; ++j)sum += L[pos(i, j)] * y[j];
        y[i] = b[i] - sum;
    }
    for (int i = n - 1; i >= 0; --i) {
        float sum = 0;
        for (int j = i + 1; j < n; ++j)sum += U[pos(i, j)] * x[j];
        if (U[pos(i, i)] != 0)
            x[i] = (y[i] - sum) / U[pos(i, i)];
    }
    free(L), free(U), free(y);
}

float polynomial(float*in, float*fin, float x, int n) {
    size_t sz = sizeof(float) * n;
    float *matA = (float *) malloc(sz * n);
    float *vecx = (float *) malloc(sz);
    float *vecb = (float *) malloc(sz);
    memcpy(vecb, fin, sz);
    for (int i = 0; i < n; ++i)for (int j = 0; j < n; ++j)matA[pos(i, j)] = qpow(in[i], j);
    lu(matA, vecx, vecb, n);
    float res = 0;
    for (int i = 0; i < n; ++i)res += vecx[i] * qpow(x, i);
    free(matA), free(vecx), free(vecb);
    return res;
}

float lagrange(float*in, float*fin, float x, int n){
    float res=0;
    for(int i=0;i<n;++i){
        float prod = 1;
        for(int j=0;j<n;++j)if(i!=j)prod *= (x-in[j]) / (in[i]-in[j]);
        res += fin[i] * prod;
    }
    return res;
}

float newton(float*in, float*fin, float x, int n) {
    float res = 0;
    float *diff = (float *) malloc(sizeof(float) * n * n);
    memcpy(diff, fin, sizeof(float) * n);
    for (int i = 1; i < n; ++i)
        for (int j = 0; j < n - i; ++j)
            diff[pos(i, j)] =
                    (diff[pos((i - 1), j)] - diff[pos((i - 1), (j + 1))]) /
                    (in[j] - in[j + i]);
    for(int i=0;i<n;++i){
        float prod = 1;
        for(int j=0;j<i;++j)prod *= x-in[j];
        res += prod*diff[pos(i,0)];
    }
    free(diff);
    return res;
}

float cubic_spline(float*in, float*fin, float x, int n) {
    /**
     * 算法思路:
     * 1. 先利用二分查找确定x的位置，如果在in中找到了x，立即返回
     * 2. 利用带状矩阵求解线性方程组，对于解集a,b,c:
     *      a = (b[i+1]-b[i]) / (in[i+1]-in[i]) / 3
     *      b = 方程组的解
     *      c = (fin[i + 1] - fin[i]) / (in[i + 1] - in[i])
               -  (2.0f * b[i] + b[i + 1]) * (in[i + 1] - in[i]) / 3
     * 3. 解为:
     *      h = x - in[l];
     *      res = ((a[l]*h + b[l])*h + c[l])*h + fin[l]
     */
    int l = 0, r = n;
    while (l < r) {
        int mid = (l + r) / 2;
        if (in[mid] < x) l = mid + 1;
        else if (in[mid] > x) r = mid - 1;
        else return fin[mid];
    }
    l = min(l,n-1);
    l = in[l]<x?l:l-1;
    l = max(l,0);
    size_t sz = sizeof(float) * n;
    float *upper = (float *) malloc(sz * 2);
    float *lower = (float *) malloc(sz * 2);
    float *rhs = (float *) malloc(sz);
    float *b = (float *) malloc(sz), *vecy = (float *) malloc(sz);
    float *a = (float*)malloc(sz), *c=(float*)malloc(sz);
    memset(upper,0,sz),memset(lower,0,sz);
    memset(b, 0, sz),memset(vecy, 0, sz);
    for (int i = 1; i < n - 1; ++i) {
        lower[pos(1, i)] = (in[i] - in[i - 1]) / 3.0f;
        upper[i] = (in[i + 1] - in[i - 1]) * 2 / 3;
        upper[pos(1, i)] = (in[i + 1] - in[i]) / 3;
        rhs[i] = (fin[i + 1] - fin[i]) / (in[i + 1] - in[i]) - (fin[i] - fin[i - 1]) / (in[i] - in[i - 1]);
    }
    upper[0] = 2;
    upper[n] = 0;
    rhs[0] = 0;
    upper[n - 1] = 2;
    lower[pos(1, n - 1)] = 0;
    rhs[n - 1] = 0;
    int i_mx, j_mx, j_mn;
    for (int i = 0; i < n; ++i) {
        lower[i] = 1.0f / upper[i];
        j_mn = max(0, i - 1);
        j_mx = min(n - 1, i + 1);
        for (int j = j_mn; j <= j_mx; ++j) {
            int k = j - i;
            if (k >= 0)upper[pos(k, i)] *= lower[i];
            else lower[pos(-k, i)] *= lower[i];
        }
        upper[i] = 1;
    }
    for (int k = 0; k < n; ++k) {
        i_mx = min(n - 1, k + 1);
        float tx = 0;
        for (int i = k + 1; i <= i_mx; ++i) {
            if (upper[k] != 0.0) {
                tx =- lower[pos(i - k, i)] / upper[k];
                lower[pos(i - k, i)] = -tx;
                j_mx = min(n - 1, k + 1);
                for (int j = k + 1; j <= j_mx; ++j) {
                    int tk = j - i;
                    if (tk >= 0)upper[pos(tk, i)] += tx * upper[pos(j - k, k)];
                    else lower[pos(-tk, i)] += tx * upper[pos(j - k, k)];
                }
            }
        }
    }
    int j_start, j_stop;
    float sum;
    for (int i = 0; i < n; ++i) {
        sum = 0;
        j_start = max(0, i - 1);
        for (int j = j_start; j < i; ++j)sum += lower[pos(i - j, i)] * vecy[j];
        vecy[i] = (rhs[i] * lower[i]) - sum;
    }
    for (int i = n - 1; i >= 0; --i) {
        sum = 0;
        j_stop = min(n - 1, i + 1);
        for (int j = i + 1; j <= j_stop; ++j)sum += upper[pos(j - i, i)] * b[j];
        b[i] = (vecy[i] - sum) / upper[i];
    }
    for(int i=0;i<n-1;++i){
        a[i] = (b[i+1]-b[i]) / (in[i+1]-in[i]) / 3;
        c[i] = (fin[i + 1] - fin[i]) / (in[i + 1] - in[i])
               -  (2.0f * b[i] + b[i + 1]) * (in[i + 1] - in[i]) / 3;
    }
    float h=in[n-1]-in[n-2];
    a[n-1] = 0;
    c[n-1] = 3.0f * a[n - 2] * h * h + 2.0f * b[n - 2] * h + c[n - 2];
    h = x - in[l];
    float res;
    if(x<in[0])res = (b[0]*h + c[0])*h + fin[0];
    else if(x>in[n-1])res = (b[n-1]*h+c[n-1]*h+fin[n-1]);
    else res = ((a[l]*h + b[l])*h + c[l])*h + fin[l];
    printf("a:");
    for(int i=0;i<n;++i)printf("%f%s",a[i],i==n-1?"\n":", ");
    printf("b:");
    for(int i=0;i<n;++i)printf("%f%s",b[i],i==n-1?"\n":", ");
    printf("c:");
    for(int i=0;i<n;++i)printf("%f%s",c[i],i==n-1?"\n":", ");
    free(upper), free(lower), free(rhs), free(vecy), free(a), free(b), free(c);
    return res;
}

float(*func[])(float*in, float*fin, float x, int n) = {polynomial, lagrange, newton, cubic_spline};

int main(int argc, char **argv) {
    int mode = argv[1][0] - '1';
    int n;
    scanf("%d", &n);
    float *in = (float *) malloc(sizeof(float) * n);
    float *fin = (float *) malloc(sizeof(float) * n);
    float x;
    for (int i = 0; i < n; ++i)scanf("%f", in + i);
    for (int i = 0; i < n; ++i)scanf("%f", fin + i);
    scanf("%f", &x);
    printf("%f\n",func[mode](in,fin,x,n));
    free(in), free(fin);
    return 0;
}
