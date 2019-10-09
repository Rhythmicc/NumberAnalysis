#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
void vecadd_2017011344(float*c, float*a,float*b,int len){
	for(int i=0;i<len;++i)c[i] = a[i]+b[i];
}

void dotprod_2017011344(float*result,float*a,float*b,int len){
	*result = 0;
	for(int i=0;i<len;++i)*result += a[i]*b[i];
}

void matvec_2017011344(float*y,float*A,float*x,int m,int n){
	memset(y,0,sizeof(float)*m);
	for(int i=0;i<m;++i)for(int j=0;j<n;++j)y[i] += A[i*m+j]*x[j];
}

void matmat_2017011344(float*C,float*A,float*B,int m,int k,int n){
	memset(C,0,sizeof(float)*(m*n));
	for(int i=0;i<m;++i)for(int j=0;j<n;++j)for(int p=0;p<k;++p)C[i*m+j] += A[i*m+p]*B[p+j*n];
}

float rand_float(){
	return rand()%1000+rand()%100*1.0/100;
}

int main(int argc, char *argv[]) {
	srand(time(NULL));
	float*a,*b,*c;
	puts("+------------------------------------------------------------------------------------------------------------------------------------+");
	puts("|                                                   test for func: vecadd_2017011344                                                 |");
	puts("+------------------------------------------------------------------------------------------------------------------------------------+");
	int len = 10;
	a = (float*) malloc(sizeof(float)*len);
	b = (float*) malloc(sizeof(float)*len);
	c = (float*) malloc(sizeof(float)*len);
	printf("vector A with length %d:\n",len);
	for(int i=0;i<len;++i){
		a[i] = rand_float();
		printf("%6.2f%c",a[i],i==len-1?'\n':'\t');
	}
	puts("--------------------------------------------------------------------------------------------------------------------------------------");
	printf("vector B with length %d:\n",len);
	for(int i=0;i<len;++i){
		b[i] = rand_float();
		printf("%6.2f%c",b[i],i==len-1?'\n':'\t');
	}
	puts("--------------------------------------------------------------------------------------------------------------------------------------");
	vecadd_2017011344(c, a, b, len);
	printf("vector C with length %d:\n",len);
	for(int i=0;i<len;++i)printf("%6.2f%c",c[i],i==len-1?'\n':'\t');
	puts("+------------------------------------------------------------------------------------------------------------------------------------+");
	puts("|                                                   test for func: dotprod_2017011344                                                |");
	puts("+------------------------------------------------------------------------------------------------------------------------------------+");
	printf("vector A with length %d:\n",len);
	for(int i=0;i<len;++i){
		a[i] = rand_float();
		printf("%6.2f%c",a[i],i==len-1?'\n':'\t');
	}
	puts("--------------------------------------------------------------------------------------------------------------------------------------");
	printf("vector B with length %d:\n",len);
	for(int i=0;i<len;++i){
		b[i] = rand_float();
		printf("%6.2f%c",b[i],i==len-1?'\n':'\t');
	}
	puts("--------------------------------------------------------------------------------------------------------------------------------------");
	c = (float*)realloc(c, sizeof(float));
	dotprod_2017011344(c, a, b, len);
	printf("dotprod_2017011344 result:%.2f\n",*c);
	puts("+------------------------------------------------------------------------------------------------------------------------------------+");
	puts("|                                                   test for func: matvec_2017011344                                                 |");
	puts("+------------------------------------------------------------------------------------------------------------------------------------+");
	int m=10,n=10;
	len = m*n;
	a = (float*)realloc(a, sizeof(float)*len);
	b = (float*)realloc(b, sizeof(float)*n);
	c = (float*)realloc(c, sizeof(float)*m);
	printf("Matrix A with row %d and col %d:\n",m,n);
	for(int i=0;i<len;++i){
		a[i] = rand_float();
		printf("%6.2f%c",a[i],i%n==n-1?'\n':'\t');
	}
	puts("--------------------------------------------------------------------------------------------------------------------------------------");
	printf("vector x with length %d:\n",n);
	for(int i=0;i<n;++i){
		b[i] = rand_float();
		printf("%6.2f%c",b[i],i==n-1?'\n':'\t');
	}
	puts("--------------------------------------------------------------------------------------------------------------------------------------");
	matvec_2017011344(c, a, b, m, n);
	printf("vector y with length %d:\n",m);
	for(int i=0;i<m;++i)printf("%6.2f%c",c[i],i==m-1?'\n':' ');
	puts("+------------------------------------------------------------------------------------------------------------------------------------+");
	puts("|                                                   test for func: matmat_2017011344                                                 |");
	puts("+------------------------------------------------------------------------------------------------------------------------------------+");
	b = (float*)realloc(b, sizeof(float)*len);
	c = (float*)realloc(c, sizeof(float)*len);
	printf("Matrix A with row %d and col %d:\n",m,n);
	for(int i=0;i<len;++i){
		a[i] = rand_float();
		printf("%6.2f%c",a[i],i%n==n-1?'\n':'\t');
	}
	puts("--------------------------------------------------------------------------------------------------------------------------------------");
	printf("Matrix B with row %d and col %d:\n",m,n);
	for(int i=0;i<len;++i){
		b[i] = rand_float();
		printf("%6.2f%c",b[i],i%n==n-1?'\n':'\t');
	}
	puts("--------------------------------------------------------------------------------------------------------------------------------------");
	matmat_2017011344(c, a, b, m, n, n);
	printf("Matrix C with row %d and col %d:\n",m,n);
	for(int i=0;i<len;++i)printf("%6.2f%c",c[i],i%n==n-1?'\n':'\t');
	puts("--------------------------------------------------------------------------------------------------------------------------------------");
	return 0;
}