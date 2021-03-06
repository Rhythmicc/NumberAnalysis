#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "csr_sparse.h"

void matvec(float *A, float *x, float *y, int m, int n) {
    for (int i = 0; i < m; ++i) {
        y[i] = 0;
        for (int j = 0; j < n; ++j)
            y[i] += A[i * n + j] * x[j];
    }
}

float dotproduct(float *vec1, float *vec2, int n) {
    float result = 0;
    for (int i = 0; i < n; ++i)
        result += vec1[i] * vec2[i];
    return result;
}

float vec2norm(float *x, int n) {
    float sum = 0;
    for (int i = 0; i < n; ++i)
        sum += x[i] * x[i];
    return sqrt(sum);
}

void cg(float *A, float *x, float *b, float *residual, float *y, float *p, float *q, int n, int maxiter,
        float threshold) {
    size_t sz = sizeof(float) * n;
    float vecb = vec2norm(b, n) * threshold;
    memset(x, 0, sz);
    int iter = 0;
    float rho = 0, rho_1 = 0;
    matvec(A, x, y, n, n);
    for (int i = 0; i < n; ++i)residual[i] = b[i] - y[i];
    do {
        rho = dotproduct(residual, residual, n);
        if (iter) {
            float beta = rho / rho_1;
            for (int i = 0; i < n; ++i)p[i] = residual[i] + beta * p[i];
        } else
            memcpy(p, residual, sz);
        matvec(A, p, q, n, n);
        float alpha = rho / dotproduct(p, q, n);
        if (isnormal(alpha)) {
            for (int i = 0; i < n; ++i)x[i] += alpha * p[i];
            for (int i = 0; i < n; ++i)residual[i] -= alpha * q[i];
        }
        rho_1 = rho;
        if (vec2norm(residual, n) < vecb)break;
    } while (++iter < maxiter);
}

void update_line(float*sm,float*yp, float val, int f,int k) {
    for (int l = k; l < f; ++l) *(sm++) += val * *(yp++);
}

const ul block = 160;

void update(csr_matrix *csr, float *X, float *Y, int f, float lamda) {
    MALLOC(residual_all, float, f * block);
    MALLOC(y_all, float, f * block);
    MALLOC(p_all, float, f * block);
    MALLOC(q_all, float, f * block);
    MALLOC(smat_all, float, f * f * block);
    MALLOC(svec_all, float, f * block);
    int *ia = csr->ia, *ja = csr->ja, m = csr->m;
    float *v = csr->val;

#pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < m; ++u) {
        int tid = omp_get_thread_num(), st=ia[u], end=ia[u+1];
        size_t sz = tid * f;
        float *smat = smat_all + sz * f;
        float *svec = svec_all + sz;
        float *residual = residual_all + sz;
        float *y = y_all + sz;
        float *p = p_all + sz;
        float *q = q_all + sz;
        float *xu = X + u * f;
        memset(smat, 0, sizeof(float) * f * f);
        memset(svec, 0, sizeof(float) * f);
        for (int j = st; j < end; ++j) { /// nnz * f * f * 0.5
            float *yn = Y + ja[j] * f, val = *yn, *sm = smat;
            for (int k = 0; k < f; ++k) {
                update_line(sm + k, yn, val, f, k);
                sm += f;
                ++yn;
            }
        }
        for (int i = 0; i < f; ++i) {
            smat[i * f + i] += lamda;
            for (int j = i + 1; j < f; ++j) smat[j * f + i] = smat[i * f + j];
        }
        for (int j = st; j < end; ++j) { /// nnz * f
            float *yn = Y + ja[j] * f, val = v[j];
            for (int i = 0; i < f; ++i)svec[i] += val * *(yn++);
        }
        cg(smat, xu, svec, residual, y, p, q, f, 100, 1e-4);
    }
    free(smat_all), free(svec_all), free(residual_all), free(y_all), free(p_all), free(q_all);
}

void als_recsys(csr_matrix *csr, float *X, float *Y, int f, float lamda) {
    // create YT and Rp
    int iter = 0, nnz = csr->nnz, *ia = csr->ia, *ja = csr->ja, m = csr->m;
    float error, *v = csr->val;
    csr_matrix *csc = new_csr(csr->n, csr->m, csr->nnz);
    csr_transpose(csr, csc);

    double error_old = 0.0;
    double error_new;
    struct timeval t1, t2;
    double time_updatex = 0;
    double time_updatey = 0;
    double time_validate = 0;
    do {
        gettimeofday(&t1, NULL);
        update(csr, X, Y, f, lamda);
        gettimeofday(&t2, NULL);
        time_updatex += (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
        gettimeofday(&t1, NULL);
        update(csc, Y, X, f, lamda);
        gettimeofday(&t2, NULL);
        time_updatey += (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
        gettimeofday(&t1, NULL);
        error_new = 0;
#pragma omp parallel for reduction(+:error_new) schedule(dynamic)
        for (int i = 0; i < m; ++i) {
            double sum, cal=0;
            for (int j = ia[i]; j < ia[i + 1]; ++j) {
                sum = 0;
                for (int k = 0; k < f; ++k)sum += X[i * f + k] * Y[ja[j] * f + k];
                cal += (sum - v[j]) * (sum - v[j]);
            }
            error_new += cal;
        }
        error_new = sqrt(error_new / nnz);
        gettimeofday(&t2, NULL);
        time_validate += (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
        error = fabs(error_new - error_old) / error_new;
        error_old = error_new;
        printf("iter = %i, error = %f, error_new = %f\n", iter, error, error_new);
        ++iter;
    } while (iter < 1000 && error > 0.000092);

    printf("\nUpdate X\t%5.2f ms\n", time_updatex);
    printf("Update Y\t%5.2f ms\n", time_updatey);
    printf("Validate\t%5.2f ms\n", time_validate);
    del_csr(csc);
}

int main(int argc, char **argv) {
    // parameters
    omp_set_num_threads(block);
    int f = atoi(argv[2]);
    char *filename = argv[1];
    printf("filename = %s\n", filename);
    csr_matrix *csr = load_mtx(filename);
    printf("The order of the rating matrix R is %d by %d, #nonzeros = %d\n",
           csr->m, csr->n, csr->nnz);

    printf("The latent feature is %i \n", f);
    float *X = (float *) malloc(sizeof(float) * csr->m * f);
    memset(X, 0, sizeof(float) * csr->m * f);
    float *Y = (float *) malloc(sizeof(float) * csr->n * f);
    for (int i = 0; i < csr->n * f; ++i) Y[i] = 1;
    float lamda = 0.1;
    struct timeval t1,t2;
    gettimeofday(&t1,NULL);
    als_recsys(csr, X, Y, f, lamda);
    gettimeofday(&t2,NULL);
    printf("Total\t%5.2f ms\n", (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0);
    free(X), free(Y);
    del_csr(csr);
}
