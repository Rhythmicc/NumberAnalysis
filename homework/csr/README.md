## CSR GEMM

 - data structure: csr_sparse
 
 | variable | means |
 |:---|:---|
 |m,n| size of matrix |
 |nrow|compressd row number|
 |nnz|non-zero element number|
 |\_\_REAL_NROW__|real space of row|
 |\_\_REAL_NNZ__|real space of col and val|
 |row,col,val|basic data with csr|
 
 - function based on csr_sparse:
 
 |name|means|
 |:---|:---|
 |csr_malloc|init csr|
 |csr_realloc|expand csr for more space|
 |csr_add_val|call this function iteratively to insert val for csr|
 |spgemm_2017011344|function to calculate csr * csr|
 |print_csr|show csr|
 |free_csr|free csr|

## probable problem:
  
  - The implementation of method `spgemm_2017011344` is what I came up with on the fly and maybe there's room for optimization 