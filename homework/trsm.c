#pragma GCC optimize("O3")
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define pos(i,j) (i)*n+(j)


const double eps = 1e-5;

void swap(float *a, float *b) {
    float tmp = *a;
    *a = *b;
    *b = tmp;
}

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

float difference(const float *x, const float *xx, int n) {
    float sum = 0;
    for (int i = 0; i < n; ++i)sum += (x[i] - xx[i]) * (x[i] - xx[i]);
    return (float)sqrt((double)sum);
}

float vec2norm(const float *x, int n) {
    float sum = 0;
    for (int i = 0; i < n; ++i)sum += x[i] * x[i];
    return (float)sqrt((double)sum);
}

int check_solution(const float*A,const float*b,const float*x,int n,int print){
    float*residual = malloc(sizeof(float)*n);
    for(int i=0;i<n;++i){
        float sum = -b[i];
        for(int j=0;j<n;++j)sum += A[pos(i,j)] * x[j];
        residual[i] = sum;
    }
    float res = vec2norm(residual,n);
    if(print)printf("        residual: %f\n", res);
    return res<1e-4?1:0;
}

void gauss_2017011344(const float *A, float *x, const float *b, int n) {
    size_t sz = sizeof(float) * n * n;
    float *At = malloc(sz);
    memcpy(At, A, sz);
    memcpy(x, b, sz / n);
    for (int k = 0, col = 0; k < n && col < n; ++k, ++col) {
        int max_r = k;
        for (int i = k + 1; i < n; ++i)
            if (fabsf(At[pos(i, col)]) > fabsf(At[pos(max_r, col)]))max_r = i;
        if (fabsf(At[pos(max_r, col)]) < eps)return; /// 无解
        if (k != max_r) {
            for (int j = col; j < n; ++j)swap(At + pos(k, j), At + pos(max_r, j));
            swap(x + k, x + max_r);
        }
        //if(At[pos(k,col,n)]!=0)
        x[k] /= At[pos(k, col)];
        for (int j = col + 1; j < n; ++j)if(At[pos(k,col)]!=0)At[pos(k, j)] /= At[pos(k, col)];
        At[pos(k, col)] = 1;
        for (int i = 0; i < n; ++i)
            if (i != k) {
                x[i] -= x[k] * At[pos(i, k)];
                for (int j = col + 1; j < n; ++j)At[pos(i, j)] -= At[pos(k, j)] * At[pos(i, col)];
                At[pos(i, col)] = 0;
            }
    }
    free(At);
}

void Lu_doolittle_2017011344(const float *A, float *x, const float *b, int n) {
    size_t sz = sizeof(float) * 1ll * n * n;
    float *L = malloc(sz), *U = malloc(sz), *y = malloc(sz / n);
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

void Lu_crout_2017011344(const float *A, float *x, const float *b, int n) {
    size_t sz = sizeof(float) * n * n;
    float *L = malloc(sz), *U = malloc(sz), *y = malloc(sz / n);
    memset(L, 0, sz);
    memset(U, 0, sz);
    for (int i = n; i >= 1; --i) {
        for (int j = n - i; j < n; ++j) {
            float sum = 0;
            for (int k = 0; k < n - i; ++k)sum += L[pos(j, k)] * U[pos(k, n - i)];
            L[pos(j, n - i)] = A[pos(j, n - i)] - sum;
        }
        for (int j = n - i + 1; j < n; ++j) {
            float sum = 0;
            for (int k = 0; k < n - i; ++k)sum += L[pos(n - i, k)] * U[pos(k, j)];
            if (L[pos(n - i, n - i)] != 0)
                U[pos(n - i, j)] = (A[pos(n - i, j)] - sum) / L[pos(n - i, n - i)];
        }
    }
    for (int i = 0; i < n; ++i)U[pos(i, i)] = 1;
    y[0] = b[0] / L[0];
    for (int i = 1; i < n; ++i) {
        float sum = 0;
        for (int j = 0; j < i; ++j)sum += L[pos(i, j)] * y[j];
        if (L[pos(i, i)] != 0)y[i] = (b[i] - sum) / L[pos(i, i)];
    }
    x[n - 1] = y[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        float sum = 0;
        for (int j = i + 1; j < n; ++j)sum += U[pos(i, j)] * x[j];
        x[i] = y[i] - sum;
    }
    free(L), free(U), free(y);
}

void cholesky_2017011344(const float *A, float *x, const float *b, int n) {
    size_t sz = sizeof(float) * n * n;
    float *L = malloc(sz);
    memcpy(x, b, sz / n);
    for (int k = 0; k < n; ++k) {
        float sum = 0;
        for (int i = 0; i < k; ++i)sum += L[pos(k, i)] * L[pos(k, i)];
        sum = A[pos(k, k)] - sum;
        L[pos(k, k)] = (float) sqrt(sum > 0 ? sum * 1.0 : 0);
        for (int i = k + 1; i < n; ++i) {
            sum = 0;
            for (int j = 0; j < k; ++j)sum += L[pos(i, j)] * L[pos(k, j)];
            if (L[pos(k, k)] != 0)
                L[pos(i, k)] = (A[pos(i, k)] - sum) / L[pos(k, k)];
        }
        for (int j = 0; j < k; ++j)L[pos(j, k)] = 0;
    }
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < k; ++i)x[k] -= x[i] * L[pos(k, i)];
        if (L[pos(k, k)] != 0)
            x[k] /= L[pos(k, k)];
    }
    free(L);
}

void jacobi_2017011344(const float *A, float *x, const float *b, int n, int *iter, int maxiter, float threshold) {
    size_t sz = sizeof(float) * n;
    float *xx = malloc(sz);
    float *At = malloc(sz * n);
    memcpy(At, A, sz * n);
    float tmp[n];
    for (int i = 0; i < n; ++i) {
        if (fabsf(At[pos(i, i)]) < threshold) {
            int flag = 0;
            for (int j = 0; j < n; ++j)
                if (fabsf(At[pos(i, j)]) > threshold) {
                    memcpy(tmp, At + pos(i, 0), sz);
                    memcpy(At + pos(i, 0), At + pos(j, 0), sz);
                    memcpy(At + pos(j, 0), tmp, sz);
                    flag = 1;
                    break;
                }
            if (!flag) {
                puts("    * Matrix is not legal.");
                free(xx), free(At);
                return;
            }
            i = 0;
        }
    }
    memset(x, 0, sz);
    for (*iter = 0; *iter < maxiter; ++*iter) {
        for (int i = 0; i < n; ++i) {
            xx[i] = b[i];
            for (int j = 0; j < n; ++j)if (i != j)xx[i] -= At[pos(i, j)] * x[j];
            xx[i] /= At[pos(i, i)];
        }
        if (difference(x, xx, n) < threshold) break;
        memcpy(x, xx, sz);
    }
    free(xx), free(At);
}

void gs_2017011344(const float *A, float *x, const float *b, int n, int *iter, int maxiter, float threshold) {
    size_t sz = sizeof(float) * n;
    float *xx = malloc(sz);
    for (*iter = 0; *iter < maxiter; ++*iter) {
        for (int i = 0; i < n; ++i) {
            xx[i] = b[i];
            for (int j = 0; j < i; ++j)xx[i] -= A[pos(i, j)] * xx[j];
            for (int j = i + 1; j < n; ++j)xx[i] -= A[pos(i, j)] * x[j];
            if(isnormal(A[pos(i, i)]))xx[i] /= A[pos(i, i)];
        }
        if (difference(x, xx, n) < threshold) {
            free(xx);
            return;
        }
        memcpy(x, xx, sz);
    }
    free(xx);
}

void sor_2017011344(const float *A, float *x, const float *b, int n, int *iter, int maxiter, float threshold) {
    size_t sz = sizeof(float) * n;
    float w = 1.0f;
    memset(x,0,sz);
    for (*iter = 0; *iter < maxiter; ++*iter) {
        float norm = 0;
        for (int i = 0; i < n; ++i) {
            float xx = x[i], sum = 0;
            for (int j = 0; j < n; ++j)if (i != j) sum += A[pos(i, j)] * x[j];
            float tmp = (1 - w) * xx + (b[i] - sum) / A[pos(i, i)];
            if(isnormal(tmp))x[i] = tmp;
            norm = fabsf(xx-x[i]) > norm?fabsf(xx-x[i]):norm;
        }
        if (norm < threshold)break;
    }
}

void cg_2017011344(const float *A, float *x, const float *b, int n, int *iter, int maxiter, float threshold) {
    size_t sz = sizeof(float) * n;
    memset(x, 0, sz);
    float *residual = malloc(sz), *y = malloc(sz), *p = malloc(sz), *q = malloc(sz);
    *iter = 0;
    float rho = 0, rho_1 = 0;
    matvec(A, x, y, n);
    for (int i = 0; i < n; ++i)residual[i] = b[i] - y[i];
    do {
        rho = dotproduct(residual, residual, n);
        if (*iter) {
            float beta = rho / rho_1;
            for (int i = 0; i < n; ++i)p[i] = residual[i] + beta * p[i];
        } else
            memcpy(p, residual, sz);
        matvec(A, p, q, n);
        float alpha = rho / dotproduct(p, q, n);
        if(isnormal(alpha)) {
            for (int i = 0; i < n; ++i)x[i] += alpha * p[i];
            for (int i = 0; i < n; ++i)residual[i] -= alpha * q[i];
        }
        rho_1 = rho;
        if (vec2norm(residual,n) < threshold)break;
    } while (++(*iter) < maxiter);
    free(residual), free(y), free(p), free(q);
}

void (*func_direct[4])(const float *, float *, const float *, int) ={
        gauss_2017011344, Lu_doolittle_2017011344,
        Lu_crout_2017011344, cholesky_2017011344
};

void (*func_iter[4])(const float *A, float *x, const float *b, int n, int *iter, int maxiter, float threshold) = {
        jacobi_2017011344, gs_2017011344,
        sor_2017011344, cg_2017011344
};

char *names[8] = {"gauss", "Lu_doolittle", "Lu_crout", "cholesky", "jacobi", "GS", "sor", "conjugate"};

void print_solution(const float*A,const float*b,float *x, int n) {
    printf("    * \033[1;31mSolution:\033[0m\n");
    int flag = 0;
    for (int i = 0; i < n; ++i)
        if (isnormal(x[i])) {
            flag = 1;
            break;
        }
    if (!flag) {
        puts("        No solution!");
        return;
    }
    puts(check_solution(A,b,x,n,1)?"        match condition":"        not match condition");
    printf("        ");
    if (n <= 10)for (int i = 0; i < n; ++i)printf("%.2f%c", x[i], i == n - 1 ? '\n' : '\t');
    else puts("\033[1;31msolution is too long!\033[0m");
}

int main(int argc, char **argv) {
    FILE *fp = fopen("func_time", "a");
    int n, t;
    struct timeval s,e;
    int mod = argc > 1 ? argv[1][0] - '1' : 0;
    if (mod > 7)return -1;
    printf("\033[1;31mUse function: %s\033[0m\n", names[mod]);
    //scanf("%d", &t);
    t = 1;
    while (t--) {
        scanf("%d", &n);
        size_t sz = sizeof(float) * n;
        float *A = malloc(sz * n);
        for (int i = 0; i < n; ++i)for (int j = 0; j < n; ++j)scanf("%f", A + pos(i, j));
        float *b = malloc(sz), *x = malloc(sz);
        for (int i = 0; i < n; ++i)scanf("%f", b + i);
        int iter = 0;
        gettimeofday(&s,NULL);
        if (mod < 4)func_direct[mod](A, x, b, n);
        else func_iter[mod - 4](A, x, b, n, &iter, 1000, (float) 1e-4);
        gettimeofday(&e,NULL);
        double second = (double)(e.tv_sec - s.tv_sec) + (double)(e.tv_usec - s.tv_usec) / 1e6;
        printf("    * Calculate Matrix A with width:\033[1;31m%4d\033[0m and use time: \033[1;31m%.2lfs\033[0m.",
               n, second);
        if(mod<4)puts("");
        else printf("Iterations:\033[1;31m%d\033[0m\n",iter);
        print_solution(A, b, x, n);
        fprintf(fp, "%.2lf\n", second);
        fflush(fp);
        free(A), free(x), free(b);
    }
    fclose(fp);
    return 0;
}
