#pragma GCC optimize("O3")
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define pos(i,j) (i)*n+(j)
#define range(i,a,b) for(int i=a;i<b;++i)
#define max(a,b) (a)>(b)?(a):(b)
#define min(a,b) (a)<(b)?(a):(b)
#define swap(a,b) tmp = *(a), *(a)=*(b),*(b)=tmp
#define MALLOC(p,type,n) type*p = (type*)malloc(sizeof(type)*(n))
#define MODE_SPGEMM

typedef struct {
    int m, n, nnz, nrow, __REAL_NNZ__, __REAL_NROW__;
    int *row, *col;
    float *val;
}csr_sparse;

void csr_realloc(csr_sparse*a) { /// explore csr
    if (a->nnz == a->__REAL_NNZ__) {
        a->__REAL_NNZ__ += 10;
        a->col = realloc(a->col, sizeof(int) * a->__REAL_NNZ__);
        a->val = realloc(a->val, sizeof(float) * a->__REAL_NNZ__);
    }
    if (a->nrow == a->__REAL_NROW__) {
        a->__REAL_NROW__ += 10;
        a->row = realloc(a->row, sizeof(int) * a->__REAL_NROW__);
    }
}

void csr_malloc(csr_sparse*C) { /// init csr
    C->__REAL_NROW__ = C->nrow + 1;
    C->__REAL_NNZ__ = C->nnz;
    C->row = (int *) malloc(sizeof(int) * C->__REAL_NROW__);
    C->col = (int *) malloc(sizeof(int) * C->nnz);
    C->val = (float *) malloc(sizeof(float) * C->nnz);
}

void csr_add_val(csr_sparse*C, float val, int col, int*nxt_line) { /// call this function iteratively
    csr_realloc(C);
    if (*nxt_line) {
        ++C->nrow;
        C->row[C->nrow] = C->row[C->nrow - 1];
    }
    ++C->row[C->nrow];
    C->col[C->nnz] = col;
    C->val[C->nnz++] = val;
    *nxt_line = 0; /// lazy variable
}

int csr_get_element_rowIndx(csr_sparse*a, int indx) { /// lower_bound
    int l = 0, r = a->nrow;
    while (l < r) {
        unsigned m = (l + r) >> 1;
        if (a->row[m] >= indx)r = m;
        else if (a->row[m] < indx) l = m + 1;
        else return max(0, m - 1);
    }
    return max(0, l - 1);
}

csr_sparse*spgemm_2017011344(csr_sparse*A, csr_sparse*B) {
    /// The implementation of this function is what I came up with on the fly and maybe there's room for optimization
    if (A->n != B->m)return NULL;
    MALLOC(C, csr_sparse, 1);
    C->m = A->m;
    C->n = B->n;
    C->nrow = 0;
    C->nnz = 0;
    csr_malloc(C);
    for (int i = 1; i <= A->nrow; ++i) {
        int B_col = 0;
        int nxt_line = 1;
        while (B_col < B->n) {
            int B_row = 0, A_indx = A->row[i - 1], B_indx = 0;
            float sum = 0;
            while (A_indx < A->row[i]) { /// mul A_col and B_row
                while (A_indx < B_row && A_indx < A->row[i]) ++A_indx;
                if (A_indx >= A->row[i])break;
                while (B->col[B_indx] != B_col && B_indx < B->nnz) ++B_indx;
                if (B_indx >= B->nnz)break;
                B_row = csr_get_element_rowIndx(B, B_indx + 1); /// (lower_bound)
                if (i - 1 == B_row)sum += A->val[A_indx] * B->val[B_indx];
                ++B_indx;
            }
            if (sum != 0)csr_add_val(C, sum, B_col, &nxt_line); /// add new value
            ++B_col;
        }
    }
    return C;
}

void print_csr(csr_sparse*a) {
    printf("row: ");
    range(i, 0, a->nrow + 1)printf("%d%c", a->row[i], i == a->nrow ? '\n' : ' ');
    printf("col: ");
    range(i, 0, a->nnz)printf("%d%c", a->col[i], i == a->nnz - 1 ? '\n' : ' ');
    printf("val: ");
    range(i, 0, a->nnz)printf("%f%c", a->val[i], i == a->nnz - 1 ? '\n' : ' ');
}

void free_csr(csr_sparse*a) {
    free(a->val), free(a->row), free(a->col);
    free(a);
}

int main(int argc, char **argv) {
    MALLOC(A, csr_sparse, 1);
    scanf("%d%d%d%d", &A->m, &A->n, &A->nrow, &A->nnz);
    csr_malloc(A);
    range(i, 0, A->nrow + 1)scanf("%d", A->row + i);
    range(i, 0, A->nnz)scanf("%d", A->col + i);
    range(i, 0, A->nnz)scanf("%f", A->val + i);
    csr_sparse *C = spgemm_2017011344(A, A);
    print_csr(C);
    //range(i,0,nrow)range(j,rowptr[i],rowptr[i+1])y[i] += x[colidx[j]] * val[j];
    free_csr(A);
    free_csr(C);
    return 0;
}
