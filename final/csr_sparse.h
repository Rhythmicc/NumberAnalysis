#include <sys/mman.h>
#include <fcntl.h>
#include <ctype.h>
#include <pthread.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define min(a,b) (a)<(b)?(a):(b)
#define RMALLOC(type,n) (type*)malloc(sizeof(type)*(n))
#define MALLOC(p,type,n) type*p = RMALLOC(type, n)
#ifndef VALUE_TYPE
#define VALUE_TYPE float
#endif
typedef long ul;

typedef struct {
    int m, n, nnz;
    int *ia, *ja;
    VALUE_TYPE *val;
}csr_matrix;

void csr_transpose(const csr_matrix*Mtx,csr_matrix *ret){
    ret->m = Mtx->n;
    ret->n = Mtx->m;
    ret->ia = calloc(sizeof(int),(Mtx->n+1));
    ret->ja = malloc(sizeof(int)*(Mtx->nnz));
    ret->val = malloc(sizeof(VALUE_TYPE)*(Mtx->nnz));
    int *cnt = (int*)calloc(sizeof(int),(Mtx->n+1));
    int nnz = Mtx->nnz;
    int row = 0,col;
    for(int i = 0 ; i < nnz ; ++i)ret->ia[Mtx->ja[i]+1]++;
    for(int i = 1 ; i <= Mtx->n ; ++i) ret->ia[i] += ret->ia[i - 1];
    for(int i = 0 ; i < nnz ; ++i){
        col = Mtx->ja[i];
        ret->ja[ret->ia[col]+cnt[col]] = row;
        ret->val[ret->ia[col]+cnt[col]] = Mtx->val[i];
        ++cnt[col];
        if(i+1==Mtx->ia[row+1])++row;
    }
    free(cnt);
}

void del_csr(csr_matrix *P){
    free(P->val);
    free(P->ja);
    free(P->ia);
}

csr_matrix* new_csr(int m,int n,int nnz){
    MALLOC(res,csr_matrix,1);
    res->m = m;
    res->n = n;
    res->nnz = nnz;
    res->ia = RMALLOC(int,m+1);
    res->ja = RMALLOC(int, nnz);
    res->val = RMALLOC(VALUE_TYPE, nnz);
    return res;
}

float readDouble(char **buffer) {
    while (!isdigit(**buffer)) ++*buffer;
    int flag = 1;
    if (*(*buffer - 1) == '-') flag = -flag;
    int ret = 0, dep = 0, dec = 1;
    while (isdigit(**buffer)) {
        ret = (ret << 1) + (ret << 3) + **buffer - '0';
        ++*buffer;
    }
    if (**buffer == '.') {
        ++*buffer;
        while (isdigit(**buffer)) {
            dep = (dep << 1) + (dep << 3) + **buffer - '0';
            ++*buffer;
            dec = (dec << 1) + (dec << 3);
        }
    }
    return flag * (ret + 1.0 * dep / dec);
}

ul readInt(char **buffer) {
    while (!isdigit(**buffer)) ++*buffer;
    ul flag = 1;
    if (*(*buffer - 1) == '-') flag = -flag;
    ul ret = 0;
    while (isdigit(**buffer)) {
        ret = (ret << 1) + (ret << 3) + **buffer - '0';
        ++*buffer;
    }
    return ret * flag;
}

csr_matrix *load_mtx(char *filename) {
#define MAX_BUFFER_SIZE 1024ll*1024*1024*2
    int m, n, nnz,idxi, idxj, fd = 0;
    if ((fd = open(filename, O_RDONLY)) < 0) printf("error\n");
    char *fp = mmap(NULL, MAX_BUFFER_SIZE, PROT_READ, MAP_SHARED, fd, 0);
    if (fp == NULL) {
        printf("error in mmap\n");
        return NULL;
    }
    unsigned totSz = sizeof(char) * (strlen(fp) + 1);
    char *bufferdIn = (char *) malloc(totSz);
    if (bufferdIn == NULL) {
        printf("Cannot alloc memory\n");
        return NULL;
    }
    memcpy(bufferdIn, fp, totSz);
    char *iter = bufferdIn;
    char line[105];
    do {
        sscanf(iter, "%[^\n]", line);
        iter += strlen(line) + 1;
    } while (line[0] == '%');
    sscanf(line, "%d %d %d", &m, &n, &nnz);
    int *cnt_row = (int *) calloc(sizeof(int), (m + 1));
    int *id_row = (int *) malloc(nnz * sizeof(int));
    int *id_col = (int *) malloc(nnz * sizeof(int));
    double fval;
    VALUE_TYPE *val = (VALUE_TYPE *) malloc(nnz * sizeof(VALUE_TYPE));
    for (int i = 0; i < nnz; i++) {
        idxi = readInt(&iter)-1;
        idxj = readInt(&iter)-1;
        fval = readDouble(&iter);
        ++cnt_row[idxi + 1];
        id_row[i] = idxi;
        id_col[i] = idxj;
        val[i] = fval;
    }
    for(int i = 0; i < m ; ++i)cnt_row[i + 1]+=cnt_row[i];
    csr_matrix*res = new_csr(m, n, nnz);
    nnz = cnt_row[m];
    memcpy(res->ia, cnt_row, (m + 1) * sizeof(int));
    memset(cnt_row, 0, (m + 1) * sizeof(int));

    for (int i = 0; i < nnz; ++i) {
        int offset = res->ia[id_row[i]] + cnt_row[id_row[i]];
        res->ja[offset] = id_col[i];
        res->val[offset] = val[i];
        cnt_row[id_row[i]]++;
    }
    // free tmp space
    free(id_col),free(val),free(id_row),free(cnt_row),free(bufferdIn);
    return res;
}
