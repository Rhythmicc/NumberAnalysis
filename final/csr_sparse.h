#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
#define RMALLOC(type,n) (type*)malloc(sizeof(type)*(n))
#define MALLOC(p,type,n) type*p = RMALLOC(type, n)
#ifndef VALUE_TYPE
#define VALUE_TYPE float
#endif
typedef unsigned long ul;
typedef unsigned long long ull;

typedef struct {
    ul id;
    VALUE_TYPE val;
} pair;

int pair_cmp(const void *a, const void *b) {
    return ((pair *) a)->id < ((pair *) b)->id ? -1 : 1;
}

typedef struct {
    pair *src;
    ul len, _REAL_LEN;
} vector;
#define vector_get_element(v, indx) (v->src[indx])
#define any_of_vector(v, i) for(pair*i=v->src; i<v->src+v->len;++i)
#define end_of_vector(v) (v->src+v->len)

void clean_vector(vector *v) {
    qsort(v->src, v->len, sizeof(pair), pair_cmp);
}

void del_vector(vector *v) {
    free(v->src);
    free(v);
}

vector *new_vector() {
    MALLOC(ret, vector, 1);
    ret->len = 0;
    ret->_REAL_LEN = 1;
    ret->src = RMALLOC(pair, 1);
    return ret;
}

void veccpy(vector *dest, vector *src) {
    dest->_REAL_LEN = src->_REAL_LEN;
    dest->len = src->len;
    dest->src = (pair *) realloc(dest->src, sizeof(pair) * dest->_REAL_LEN);
    memcpy(dest->src, src->src, sizeof(pair) * dest->len);
}

void vector_push_back(vector *v, pair val) {
    if (v->len == v->_REAL_LEN) {
        v->_REAL_LEN += 1000;
        v->src = (pair *) realloc(v->src, sizeof(pair) * v->_REAL_LEN);
    }
    v->src[v->len].id = val.id;
    v->src[v->len++].val = val.val;
}

typedef struct {
    ul m, n, nnz;
    vector **src;
} csr_matrix;
#define anyrow_of_csr(csr, i) for(vector**i=csr->src;i<csr->src+csr->m;++i)
#define indx_of_iter(csr, i) (i-csr->src)

void clean_csr_matrix(csr_matrix *csr) {
    for (int i = 0; i < csr->m; ++i)clean_vector(csr->src[i]);
}

csr_matrix *new_csr_matrix(size_t m, size_t n) {
    MALLOC(res, csr_matrix, 1);
    res->m = m;
    res->n = n;
    res->src = RMALLOC(vector*, m);
    anyrow_of_csr(res, i) (*i) = new_vector();
    return res;
}

void del_csr(csr_matrix *csr) {
    anyrow_of_csr(csr, i)del_vector((*i));
    free(csr->src);
}

csr_matrix *load_mtx(char *filename) {
    FILE *fp = fopen(filename, "r");
    char line[105];
    while (fgets(line, 100, fp), line[0] == '%');
    size_t n, m, nnz, row, col;
    VALUE_TYPE val;
    sscanf(line, "%lu%lu%lu", &m, &n, &nnz);
    csr_matrix *ret = new_csr_matrix(m, n);
    ret->nnz = nnz;
    for (int i = 0; i < nnz; ++i) {
        fscanf(fp, "%lu%lu%f", &row, &col, &val);
        vector_push_back(ret->src[row - 1], (pair) {col - 1, val});
    }
    fclose(fp);
    clean_csr_matrix(ret);
    return ret;
}

void csr_to_matrix(csr_matrix *a, VALUE_TYPE *dest) {
    anyrow_of_csr(a, i)any_of_vector((*i), j)dest[(i - a->src) * a->n + j->id] = j->val;
}

void csr_to_vec(csr_matrix*src, vector*vec) {
    anyrow_of_csr(src, i)any_of_vector((*i), j)vector_push_back(vec, (pair){(i-src->src)*src->n+j->id,j->val});
}

void vec_to_csr(vector *src, unsigned n, csr_matrix *c) {
    c->nnz = src->len;
    any_of_vector(src, i)vector_push_back(c->src[i->id / n], (pair) {i->id % n, i->val});
}

void vecT_to_csr(vector *src, unsigned n, csr_matrix *c) {
    c->nnz = src->len;
    any_of_vector(src, i)vector_push_back(c->src[i->id % n], (pair) {i->id / n, i->val});
}

void csr_transpose(csr_matrix *t, csr_matrix *src) {
    vector *vec = new_vector();
    for (int i = 0; i < src->m; ++i)
        any_of_vector(src->src[i], j)vector_push_back(vec, (pair) {i * src->n + j->id, j->val});
    clean_vector(vec);
    vecT_to_csr(vec, src->n, t);
    free(vec);
}

void csrvec(csr_matrix *a, vector *v, vector *res) {
    if (vector_get_element(v, v->len - 1).id >= a->n)return;
    int indx = 0;
    anyrow_of_csr(a, i) {
        pair *pa = (*i)->src, *pv = v->src;
        pair *enda = end_of_vector((*i)), *endv = end_of_vector(v);
        VALUE_TYPE sum = 0;
        while (pa < enda && pv < endv) {
            while (pa->id < pv->id && pa < enda)++pa;
            if (pa == enda)break;
            while (pv->id < pa->id && pv < endv)++pv;
            if (pv == endv)break;
            sum += pa->val * pv->val;
            ++pa, ++pv;
        }
        if (fabsf(sum) >= 1e-5)vector_push_back(res, (pair) {indx, sum});
        ++indx;
    }
}

void csrcsr_to_vec(csr_matrix *a, csr_matrix *b, vector *res) {
    if (a->n != b->m) return;
    unsigned cscm = b->n, cscn = b->m, indx = 0;
    csr_matrix *csc = new_csr_matrix(cscm, cscn);
    csr_transpose(csc, b);
    anyrow_of_csr(a, i) { /// a rows
        anyrow_of_csr(csc, j) { /// b cols
            pair *pa = (*i)->src, *pb = (*j)->src;
            pair *enda = end_of_vector((*i)), *endb = end_of_vector((*j));
            VALUE_TYPE sum = 0;
            while (pa < enda && pb < endb) {
                while (pa->id < pb->id && pa < enda)++pa;
                if (pa == enda)break;
                while (pb->id < pa->id && pb < endb)++pb;
                if (pb == endb)break;
                sum += pa->val * pb->val;
                ++pa, ++pb;
            }
            if (fabsf(sum) >= 1e-5f)vector_push_back(res, (pair) {indx, sum});
            ++indx;
        }
    }
    del_csr(csc);
}

void csrcsr(csr_matrix *c, csr_matrix *a, csr_matrix *b) {
    vector *res = new_vector();
    csrcsr_to_vec(a, b, res);
    vec_to_csr(res, b->n, c);
    del_vector(res);
}

void matmat_to_vec(vector *res, float *A, float *B, int m, int k, int n, ul*one_base_cal_index) {
    VALUE_TYPE sum;
    res->len = 0;
    if (!one_base_cal_index) {
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                sum = 0;
                for (int kk = 0; kk < k; kk++)//C[i * n + j] += A[i * k + kk] * B[kk * n + j];
                    sum += A[i * k + kk] * B[kk * n + j];
                if (fabsf(sum) >= 1e-5f)vector_push_back(res, (pair) {i * n + j, sum});
            }
    } else {
        ul len = one_base_cal_index[0];
        ul *p = one_base_cal_index + 1;
        while (len--) {
            sum = 0;
            for (int kk = 0; kk < k; kk++)sum += A[(*p / n) * k + kk] * B[kk * n + (*p % n)];
            vector_push_back(res, (pair) {*p, sum});
            ++p;
        }
    }
}

void matmat_transB_to_vec(vector *res, float *A, float *BT, int m, int k, int n, ul*one_base_cal_index) {
    VALUE_TYPE sum;
    res->len = 0;
    if (!one_base_cal_index) {
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                sum = 0;
                for (int kk = 0; kk < k; kk++)//C[i * n + j] += A[i * k + kk] * B[kk * n + j];
                    sum += A[i * k + kk] * BT[j * k + kk];
                if (fabsf(sum) >= 1e-5f)vector_push_back(res, (pair) {i * n + j, sum});
            }
    } else {
        ul len = one_base_cal_index[0];
        ul *p = one_base_cal_index + 1;
        while (len--) {
            sum = 0;
            for (int kk = 0; kk < k; kk++)sum += A[(*p / n) * k + kk] * BT[(*p % n) * k + kk];
            vector_push_back(res, (pair) {*p, sum});
            ++p;
        }
    }
}
