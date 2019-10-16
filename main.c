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

float func(float x){
    return sqrt(1.0 - x*x);
}

float midpoint_2017011344(float(*func)(float), float a, float b, int nseg) {
    float integral = 0;
    float interval = (b - a) / (float) nseg;
    for (int i = 0; i < nseg; ++i) {
        float l = a + i * interval;
        float r = a + (i + 1) * interval;
        float m = (l + r) / 2.0;
        float y = func(m);
        integral += y * interval;
    }
    return integral;
}

float trapezoid_2017011344(float(*func)(float), float a, float b, int nseg) {
    float integral = 0;
    float interval = (b - a) / (float) nseg;
    for (int i = 0; i < nseg; ++i) {
        float l = a + i * interval;
        float r = a + (i + 1) * interval;
        float ly = func(l), ry = func(r);
        integral += (ly + ry) * interval / 2.0;
    }
    return integral;
}

float simpson13_2017011344(float(*func)(float), float a, float b, int nseg) {
    float integral = 0;
    float interval = (b - a) / (float) nseg;
    for (int i = 0; i <= nseg; ++i) {
        float y = func(a + i * interval);
        if (i == 0 || i == nseg)integral += y;
        else if (i & 1)integral += 4.0 * y;
        else integral += 2.0 * y;
    }
    integral = integral * interval / 3.0;
    return integral;
}

float simpson38_2017011344(float(*func)(float), float a, float b, int nseg) {
    float integral = 0;
    float interval = (b - a) / (float) nseg;
    for (int i = 0; i <= nseg; ++i) {
        float y = func(a + i * interval);
        if (i == 0 || i == nseg)integral += y;
        else if (i % 3 == 1)integral += 3.0 * y;
        else if (i % 3 == 2)integral += 3.0 * y;
        else integral += 2.0 * y;
    }
    integral *= interval / 8.0 * 3;
    return integral;
}

float romberg_2017011344(float(*func)(float), float a, float b, int iter) {
    iter /= 10;
    MALLOC(result, float, iter * iter);
    memset(result, 0, sizeof(float) * iter * iter);
    float interval = (b - a) / 2.0;
    int nseg = 2;
    for (int i = 0; i < iter; ++i) {
        float intergral = 0;
        for (int j = 0; j <= nseg; ++j) {
            float y = func(a + j * interval);
            if (j == 0 || j == nseg)intergral += y;
            else intergral += 2 * y;
        }
        result[i * iter] = intergral * interval / 2;
        interval /= 2.0;
        nseg <<= 1;
    }
    for (int j = 1; j < iter; ++j)
        for (int i = 0; i < iter - j; ++i)
            result[i * iter + j] = (4 * result[(i + 1) * iter + j - 1] - result[(i * iter + j - 1)]) / 3;
    free(result);
    return result[iter-1];
}

float gauss_2017011344(float(*func)(float), float a, float b, int ignore) {
    float w1 = 5.0 / 9, w2 = 8.0 / 9, w3 = 5.0 / 9;
    float z1 = -0.7746, z2 = 0, z3 = 0.7746;
    float x1 = (b - a) * z1 / 2.0 + (b + a) / 2.0;
    float x2 = (b - a) * z2 / 2.0 + (b + a) / 2.0;
    float x3 = (b - a) * z3 / 2.0 + (b + a) / 2.0;
    float y1 = func(x1), y2 = func(x2), y3 = func(x3);
    float intergral = ((b - a) / 2) * (w1 * y1 + w2 * y2 + w3 * y3);
    return intergral;
}

float montecario_2017011344(float(*func)(float), float a, float b, int nseg) {
    float integral = 0, n = nseg;
    while (nseg--) {
        float x = a + (b - a) * (nseg+1)/n;
        integral += func(x);
    }
    return integral * (b - a) / n;
}

float (*funcs[])(float(*)(float), float, float, int) = {
        midpoint_2017011344, trapezoid_2017011344,
        simpson13_2017011344, simpson38_2017011344,
        montecario_2017011344,gauss_2017011344,
        romberg_2017011344
};
char*func_names[]={
        "midpoint","trapezoid","simpson13","simpson38","montecario","gauss","romberg"
};

int main(int argc, char **argv) {
    int nseg = atoi(argv[1]);
    const float PI = acosf(-1);
    printf("| Function\t| Result\t| Error  \t|\n");
    for (int mode = 0; mode < 6; ++mode) {
        float res = 2 * funcs[mode](func, -1, 1, nseg);
        float error = fabsf(PI - res);
        printf("| \033[1;31m%-10s\033[0m\t| \033[1;31m%f\033[0m\t| \033[1;31m%f\033[0m\t|\n", func_names[mode],
               res, error);
    }
    return 0;
}
