#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define pos(i,j) ((i)*n+(j))
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
#define swap(a,b) tmp = *(a), *(a)=*(b),*(b)=tmp
#define RMALLOC(type,n) (type*)malloc(sizeof(type)*(n))
#define MALLOC(p,type,n) type*p = RMALLOC(type, n)
#ifndef VALUE_TYPE
#define VALUE_TYPE float
#endif

typedef struct {
    int id;
    VALUE_TYPE val;
}pair;

int pair_cmp(const void*a, const void*b) {
    return ((pair *) a)->id < ((pair *) b)->id ? -1 : 1;
}

typedef struct {
    pair*src;
    unsigned len, _REAL_LEN;
}vector;
#define vector_get_element(v, indx) (v->src[indx])

void clean_vector(vector*v){
    qsort(v->src, v->len, sizeof(pair), pair_cmp);
}

void del_vector(vector*v){
    free(v->src);
    free(v);
}

vector*new_vector(){
    MALLOC(ret, vector, 1);
    ret->len = 0;
    ret->_REAL_LEN = 1;
    ret->src = RMALLOC(pair, 1);
    return ret;
}

void veccpy(vector*dest, vector*src){
    dest->_REAL_LEN = src->_REAL_LEN;
    dest->len = src->len;
    dest->src = (pair*)realloc(dest->src, sizeof(pair) * dest->_REAL_LEN);
    memcpy(dest->src, src->src, sizeof(pair) * dest->len);
}

void vector_push_back(vector*v,pair val){
    if(v->len==v->_REAL_LEN){
        v->_REAL_LEN += 10;
        v->src = (pair*)realloc(v->src, sizeof(pair)*v->_REAL_LEN);
    }
    v->src[v->len].id = val.id;
    v->src[v->len++].val = val.val;
}

typedef struct {
    int m, n;
    vector**src;
}csr_matrix;

void clean_csr_matrix(csr_matrix*csr){
    for(int i=0; i < csr->m; ++i)clean_vector(csr->src[i]);
}

csr_matrix*new_csr_matrix(int m, int n){
    MALLOC(res, csr_matrix, 1);
    res->m = m;
    res->n = n;
    res->src = RMALLOC(vector*, m);
    for(int i=0;i<m;++i)res->src[i] = new_vector();
    return res;
}

void del_csr(csr_matrix*csr){
    for(int i=0; i < csr->m; ++i)del_vector(csr->src[i]);
    free(csr->src);
}

csr_matrix*load_mtx(char*filename) {
    FILE *fp = fopen(filename, "r");
    char line[105];
    while (fgets(line, 100, fp), line[0] == '%');
    unsigned n, m, nnz, row, col;
    VALUE_TYPE val;
    sscanf(line, "%d%d%d", &m, &n, &nnz);
    csr_matrix *ret = new_csr_matrix(m, n);
    for (int i = 0; i < nnz; ++i) {
        fscanf(fp, "%u%u%f", &row, &col, &val);
        vector_push_back(ret->src[row-1], (pair){col-1, val});
    }
    fclose(fp);
    clean_csr_matrix(ret);
    return ret;
}

void csr_transpose(csr_matrix*t, csr_matrix*src) {
    vector *vec = new_vector();
    for (int i = 0; i < src->m; ++i)
        for (int j = 0; j < src->src[i]->len; ++j) {
            pair tmp = vector_get_element(src->src[i], j);
            vector_push_back(vec, (pair) {(i * src->n) + tmp.id, tmp.val});
        }
    clean_vector(vec);
    for(int i=0;i<vec->len;++i){
        pair*tmp = vec->src+i;
        int col = tmp->id / src->n;
        int row = tmp->id % src->n;
        vector_push_back(t->src[row], (pair){col, tmp->val});
    }
    free(vec);
}

void csrvec(csr_matrix*a, vector*v, vector*res) {
    if (vector_get_element(v, v->len - 1).id >= a->n)return;
    for (int i = 0; i < a->m; ++i) {
        int pa = 0, pv = 0;
        VALUE_TYPE sum = 0;
        while (pa < a->src[i]->len && pv < v->len){
            while (vector_get_element(a->src[i], pa).id < vector_get_element(v, pv).id &&
                   pa < a->src[i]->len)
                ++pa;
            if (pa == a->src[i]->len)break;
            while (vector_get_element(v, pv).id < vector_get_element(a->src[i], pa).id &&
                   pv < v->len)
                ++pv;
            if (pv == v->len)break;
            sum += vector_get_element(a->src[i], pa).val * vector_get_element(v, pv).val;
            ++pa, ++pv;
        }
        if(fabsf(sum)>=1e-5)vector_push_back(res, (pair){i, sum});
    }
}

void csr_to_matrix(csr_matrix*a, VALUE_TYPE*dest){
    for(int i=0;i<a->m;++i){
        for(int j=0;j<a->src[i]->len;++j){
            pair tmp = vector_get_element(a->src[i], j);
            dest[i*a->n + tmp.id] = tmp.val;
        }
    }
}

void csrcsr_to_matrix(csr_matrix*a, csr_matrix*b, unsigned*indxs, VALUE_TYPE*res) {
    if (a->n != b->m) return;
    int flag = indxs==NULL?0:1;
    unsigned cscm = b->n, cscn = b->m;
    csr_matrix *csc = new_csr_matrix(cscm, cscn);
    csr_transpose(csc, b);
    memset(res, 0, sizeof(VALUE_TYPE) * (a->m) * (b->n));
    for (int i = 0; i < a->m; ++i) { /// a rows
        for (int j = 0; j < csc->m; ++j) { /// b cols
            int pa = 0, pb = 0;
            while (pa < a->src[i]->len && pb < csc->src[j]->len) {
                while (vector_get_element(a->src[i], pa).id < vector_get_element(csc->src[j], pb).id &&
                       pa < a->src[i]->len)
                    ++pa;
                if (pa == a->src[i]->len)break;
                while (vector_get_element(csc->src[j], pb).id < vector_get_element(a->src[i], pa).id &&
                       pb < csc->src[j]->len)
                    ++pb;
                if (pb == csc->src[j]->len)break;
                res[i * cscm + j] += vector_get_element(a->src[i], pa).val * vector_get_element(csc->src[j], pb).val;
                ++pa, ++pb;
            }
            if(flag && fabsf(res[i*cscm+j])>=1e-5)indxs[++indxs[0]] = i*cscm+j;
        }
    }
    del_csr(csc);
}

void matrix_to_csr(VALUE_TYPE*src, unsigned *idxs, unsigned n, csr_matrix*c) {
    for(int i=1;i<=idxs[0];++i){
        int row = idxs[i] / n;
        int col = idxs[i] % n;
        vector_push_back(c->src[row], (pair){col, src[idxs[i]]});
    }
}

void csrcsr(csr_matrix*c, csr_matrix*a, csr_matrix*b){
    MALLOC(idxs, unsigned, 40000);
    MALLOC(res, VALUE_TYPE, (a->m*b->n));
    csrcsr_to_matrix(a, b, idxs, res);
    matrix_to_csr(res, idxs, b->n, c);
    free(res);
    free(idxs);
}

void printcsr(csr_matrix*csr){
    for(int i=0;i<min(csr->m, 5); ++i){
        for(int j=0; j < csr->src[i]->len; ++j){
            pair tmp = vector_get_element(csr->src[i], j);
            printf("%d %d %f\n", i, tmp.id, tmp.val);
        }
    }
}

void printmat(VALUE_TYPE*mat, int m, int n){
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            printf("%4.2f ", mat[i * n + j]);
        printf("\n");
    }
}

int main(int argc, char **argv) {
    csr_matrix*mtx = load_mtx("case.mtx");
    printcsr(mtx);
    puts("------------");
    csr_matrix*mtxT = new_csr_matrix(mtx->n, mtx->m);
    csr_transpose(mtxT, mtx);
    printcsr(mtxT);
    puts("------------");
    csr_matrix*mat = new_csr_matrix(mtx->m, mtx->m);
    csrcsr(mat, mtx, mtxT);
    printcsr(mat);
    del_csr(mtx);
    del_csr(mtxT);
    del_csr(mat);
    return 0;
}
