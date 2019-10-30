#pragma GCC optimize("O3")
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

float func(float x, float y) {
    return y + exp(2 * x);
}

float eular(float(*func)(float, float), float h, float x0, float y0, float xn) {
    float x_old = x0;
    float y_old = y0;
    float x_new, y_new;
    do {
        x_new = x_old + h;
        y_new = y_old + h * func(x_old, y_old);
        x_old = x_new;
        y_old = y_new;
    } while (x_new < xn);
    return y_new;
}

float midpoint(float(*func)(float, float), float h, float x0, float y0, float xn) {
    float x_old = x0;
    float y_old = y0;
    float x_new, y_new;
    do {
        x_new = x_old + h;
        float k1 = func(x_old, y_old);
        float k2 = func(x_old + 0.5 * h, y_old + 0.5 * k1 * h);
        y_new = y_old + k2 * h;
        x_old = x_new;
        y_old = y_new;
    } while (x_new < xn);
    return y_new;
}

float rungekutta2nd(float(*func)(float, float), float h, float x0, float y0, float xn) {
    float x_old = x0;
    float y_old = y0;
    float x_new, y_new;
    do {
        x_new = x_old + h;
        float k1 = func(x_old, y_old);
        float k2 = func(x_old + h, y_old + k1 * h);
        y_new = y_old + (k1 + k2) / 2 * h;
        x_old = x_new;
        y_old = y_new;
    } while (x_new < xn);
    return y_new;
}

float rungekutta4th(float(*func)(float, float), float h, float x0, float y0, float xn) {
    float x_old = x0;
    float y_old = y0;
    float x_new, y_new;
    do {
        x_new = x_old + h;
        float k1 = func(x_old, y_old);
        float k2 = func(x_old + 0.5 * h, y_old + 0.5 * k1 * h);
        float k3 = func(x_old + 0.5 * h, y_old + 0.5 * k2 * h);
        float k4 = func(x_old + h, y_old + k3 * h);
        y_new = y_old + (k1 + 2 * k2 + 2 * k3 + k4) / 6 * h;
        x_old = x_new;
        y_old = y_new;
    } while (x_new < xn);
    return y_new;
}

float (*f[])(float(*func)(float, float), float h, float x0, float y0, float xn) ={eular, midpoint, rungekutta2nd,
                                                                                  rungekutta4th};

char *name[] = {"eular", "midpoint", "rungekutta2nd", "rungekutta4th"};

int main(int argc, char **argv) {
    int id = atoi(argv[1]) - 1;
    float h = 0.1, x0 = 0, y0 = 1, xn = 0.3;
    printf("%s %f\n", name[id], f[id](func, h, x0, y0, xn));
    return 0;
}
