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
    return (float) sqrt((double) sum);
}

void transpose(float *AT, float *A, int m, int n) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            AT[j * m + i] = A[i * n + j];
}

void tridiagonalization(float *A, float *T, float *V, int n) {
    memset(T, 0, sizeof(float) * n * n);
    float *va = (float *) malloc(sizeof(float) * n);
    float *vb = (float *) malloc(sizeof(float) * n);
    float *wp = (float *) malloc(sizeof(float) * n);
    float *w = (float *) malloc(sizeof(float) * n);
    memset(va, 0, sizeof(float) * n);
    va[0] = 1.0;
    matvec(A, va, wp, n);
    float alpha = dotproduct(wp, va, n);
    T[0 * n + 0] = alpha;
    for (int i = 0; i < n; i++)
        w[i] = wp[i] - alpha * va[i];
    for (int i = 0; i < n; i++)
        V[i * n + 0] = va[i];
    for (int j = 1; j < n; j++) {
        float beta = vec2norm(w, n);
        for (int i = 0; i < n; i++) vb[i] = w[i] / beta;
        matvec(A, vb, wp, n);
        alpha = dotproduct(wp, vb, n);
        T[j * n + j] = alpha;
        for (int i = 0; i < n; i++) w[i] = wp[i] - alpha * vb[i] - beta * va[i];
        for (int i = 0; i < n; i++) V[i * n + j] = vb[i];
        memcpy(vb, va, sizeof(float) * n);
        T[(j - 1) * n + j] = beta;
        T[j * n + (j - 1)] = beta;
    }
    free(va), free(vb), free(wp), free(w);
}

void gramschmidt(float *V, float *X, int m, int n) {
    MALLOC(sum, float, n);
    for (int i = 0; i < m; ++i) {
        float *vi = V + pos(i, 0);
        float *xi = X + pos(i, 0);
        memset(sum, 0, sizeof(float) * n);
        for (int j = 0; j < i; ++j) {
            float *vj = V + pos(j, 0);
            float dp1 = dotproduct(xi, vj, n);
            float dp2 = dotproduct(vj, vj, n);
            for (int k = 0; k < n; ++k) {
                sum[k] += (dp1 / dp2) * vj[k];
            }
        }
        for (int k = 0; k < n; ++k)vi[k] = xi[k] - sum[k];
    }
    free(sum);
}

void matmat(float *C, float *A, float *B, int m, int k, int n) {
    memset(C, 0, sizeof(float) * m * n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            for (int kk = 0; kk < k; kk++)
                C[i * n + j] += A[i * k + kk] * B[kk * n + j];
}

void matmat_transA(float *C, float *AT, float *B, int m, int k, int n) {
    memset(C, 0, sizeof(float) * m * n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            for (int kk = 0; kk < k; kk++)
                C[i * n + j] += AT[kk * m + i] * B[kk * n + j];
}

void matmat_transB(float *C, float *A, float *BT, int m, int k, int n) {
    memset(C, 0, sizeof(float) * m * n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            for (int kk = 0; kk < k; kk++)
                C[i * n + j] += A[i * k + kk] * BT[j * k + kk];
}

void qr(float *Q, float *R, float *A, int m, int n) {
    gramschmidt(Q, A, m, n);
    for (int i = 0; i < n; ++i) {
        float *qi = Q + pos(i, 0);
        float norm = vec2norm(qi, n);
        for (int j = 0; j < n; ++j)qi[j] /= norm;
    }
    MALLOC(QT, float, m * n);
    transpose(QT, Q, m, n);
    matmat_transA(R, QT, A, m, n, n);
    free(QT);
}

void qreigensolver(float *eigenvalue, float *eigenvector, float *A, int n, int maxiter, float threshold) {
    MALLOC(Ai, float, n * n);
    MALLOC(Ai_new, float, n * n);
    MALLOC(Q, float, n * n);
    MALLOC(R, float, n * n);
    MALLOC(error, float, n);
    MALLOC(Qaccum, float, n * n);
    memset(Qaccum, 0, sizeof(float) * n * n);
    for (int i = 0; i < n; i++) Qaccum[pos(i, i)] = 1;
    memcpy(Ai, A, sizeof(float) * n * n);
    int iter = 0;
    float errornorm = 0;
    do {
        memset(Q, 0, sizeof(float) * n * n);
        memset(R, 0, sizeof(float) * n * n);
        qr(Q, R, Ai, n, n);
        matmat_transB(eigenvector, Qaccum, Q, n, n, n);
        memcpy(Qaccum, eigenvector, sizeof(float) * n * n);
        matmat_transB(Ai_new, R, Q, n, n, n);
        for (int i = 0; i < n; i++) error[i] = Ai_new[pos(i, i)] - Ai[pos(i, i)];
        errornorm = vec2norm(error, n);
        memcpy(Ai, Ai_new, sizeof(float) * n * n);
        iter++;
    } while (iter < maxiter && errornorm > threshold);
    for (int i = 0; i < n; i++) eigenvalue[i] = Ai_new[i * n + i];
    free(Ai), free(Ai_new), free(Q), free(R), free(error);
}

void lanczos(float *A, float *U, float *sigma, float *V, int n) {
    MALLOC(T, float, n * n);
    MALLOC(VV, float, n * n);
    tridiagonalization(A, T, VV, n);
    MALLOC(eigenvalues, float, n); /// sigma
    MALLOC(eigenvector, float, n * n); /// -U
    qreigensolver(eigenvalues, eigenvector, T, n, 20, 1e-5f);
    MALLOC(F, float, n * n);
    MALLOC(FT, float, n * n); /// -V
    matmat(F, VV, eigenvector, n, n, n);
    transpose(FT, F, n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            eigenvector[pos(i, j)] = -eigenvector[pos(i, j)];
            FT[pos(i, j)] = -FT[pos(i, j)];
        }
    memcpy(U, eigenvector, sizeof(float) * n * n);
    memcpy(sigma, eigenvalues, sizeof(float) * n);
    memcpy(V, FT, sizeof(float) * n * n);
    free(eigenvalues), free(eigenvector), free(T);
    free(VV), free(F), free(FT);
}

void read_RGBA_to_float(unsigned *raw_img, int n) {
    unsigned val = 0;
    unsigned pixel;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < 4; ++k) {
                scanf("%u", &pixel);
                val = (val << 8) | pixel;
            }
            raw_img[pos(i, j)] = val;
        }
}

void float_to_RGBA(unsigned *img, float *ret, int n) {
    int chan = n * n;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            unsigned pixel = img[pos(i, j)];
            for (int k = 3; k >= 0; --k) {
                ret[k * chan + pos(i, j)] = pixel & 0xff;
                pixel >>= 8;
            }
        }
}

void write_float_to_RGBA(unsigned *img, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            unsigned pixel = img[pos(i, j)];
            unsigned vector[4];
            for (int k = 3; k >= 0; --k) {
                vector[k] = pixel & 0xff;
                pixel >>= 8;
            }
            for (int k = 0; k < 4; ++k)printf("%d%c", vector[k], j == n - 1 && k == 3 ? '\n' : ' ');
        }
}

void lanczos_2017011344(unsigned *raw_img, unsigned *ret_img, int n, float rate) {
    MALLOC(RGBA, float, 4 * n * n);
    float_to_RGBA(raw_img, RGBA, n);
    MALLOC(chan, float, n * n);
    MALLOC(result, float, n * n);
    MALLOC(ans, float, n * n);
    MALLOC(U, float, n * n);
    MALLOC(sigma, float, n);
    MALLOC(V, float, n * n);
    for (int k = 0; k < 4; ++k) {
        memcpy(chan, RGBA + k * n * n, sizeof(float) * n * n);
        if (k == 3)memcpy(ans, chan, sizeof(float) * n * n);
        else {
            matmat_transB(result, chan, chan, n, n, n);
            lanczos(chan, U, sigma, V, n);
            MALLOC(tmp_sigma, float, n);
            memcpy(tmp_sigma, sigma, sizeof(float) * n);
            lanczos(result, V, sigma, U, n);
            //memcpy(sigma, tmp_sigma, sizeof(float)*n);
            float sum = 0, temp = 0;
            int n_sigma = 0;
            for (int i = 0; i < n; ++i)sum += sigma[i];
            while (temp / sum < rate)temp += sigma[n_sigma++];
            MALLOC(tmp, float, n_sigma * n_sigma);
            memset(tmp, 0, sizeof(float) * n_sigma * n_sigma);
            for (int i = 0; i < n_sigma; ++i)tmp[i * n_sigma + i] = sigma[i];
            matmat(result, U, tmp, n, n_sigma, n_sigma);
            matmat(ans, result, V, n, n_sigma, n);
            float mx = -1, mn = MAXFLOAT;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    mx = max(mx, ans[pos(i, j)]), mn = min(mn, ans[pos(i, j)]);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    ans[pos(i, j)] = (ans[pos(i, j)] - mn) / (mx - mn) * 255;
            free(tmp_sigma), free(tmp);
        }
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                ret_img[pos(i, j)] = (ret_img[pos(i, j)] << 8) + (int) ans[pos(i, j)];
    }
    free(RGBA), free(chan), free(U), free(sigma), free(V), free(result), free(ans);
}

int main(int argc, char **argv) {
    int n, m;
    scanf("%d%d", &n, &m);
    if (n != m)return -1;
    MALLOC(img, unsigned, n * n);
    MALLOC(ret, unsigned, n * n);
    read_RGBA_to_float(img, n);
    lanczos_2017011344(img, ret, n, 0.8);
    printf("%d %d\n", n, m);
    write_float_to_RGBA(ret, n);
    free(img), free(ret);
    return 0;
}
