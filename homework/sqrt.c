#pragma GCC optimize("O3")
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pos(i,j) i*n+j

float fa(float x) {
    return x * x - 4;
}

float fb(float x) {
    return x + 1 - x * x / 4;
}

float fcx(float x) {
    return x * x - 4;
}

float fcdx(float x) {
    return 2 * x;
}

float bisection_2017011344(float(*func)(float), float left, float right, int maxiter, float threshold) {
    int iter = 0;
    float m = (left + right) / 2;
    while (iter < maxiter && fabsf(func(m)) > threshold) {
        ++iter;
        if (func(left) * func(m) < 0)right = m;
        else if (func(right) * func(m) < 0)left = m;
        m = (left + right) / 2;
    }
    return m;
}

float fixedpoint_2017011344(float(*func)(float), float initguess, int maxiter, float threshold) {
    float pre = 0, x = initguess;
    int iter = 0;
    while (iter++ < maxiter) {
        pre = x;
        x = func(x);
        if (fabsf(x - pre) < threshold)break;
    }
    return x;
}

float newtonraphson_2017011344(float(*func)(float), float(*funcderivative)(float), float initguess, int maxiter, float threshold) {
    float pre = 0, x = initguess;
    int iter = 0;
    while (iter++ < maxiter) {
        pre = x;
        x = x - func(x) / funcderivative(x);
        if (fabsf(x - pre) < threshold)break;
    }
    return x;
}

float secant_2017011344(float(*func)(float), float initguess0, float initguess1, int maxiter, float threshold){
    float pre = initguess0, x = initguess1;
    int iter = 0;
    while(iter++ < maxiter){
        float tmp = x;
        x = x - (x-pre) / (func(x)-func(pre)) * func(x);
        pre = tmp;
        if (fabsf(x - pre) < threshold)break;
    }
    return x;
}

int main(int argc, char **argv) {
    int mode = argv[1][0] - '1';
    switch (mode) {
        case 0:
            printf("%f\n", bisection_2017011344(fa, 0, 15, 100, 1e-4f));
            break;
        case 1:
            printf("%f\n", fixedpoint_2017011344(fb, 4, 100, 1e-4f));
            break;
        case 2:
            printf("%f\n", newtonraphson_2017011344(fcx, fcdx, 4, 100, 1e-4f));
            break;
        case 3:
            printf("%f\n", secant_2017011344(fcx,0, 4, 100, 1e-4f));
            break;
        default:
            puts("Not a legal mode!");
            break;
    }
    return 0;
}
