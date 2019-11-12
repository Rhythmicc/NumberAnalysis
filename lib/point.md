## 1.直接法解法(LU分解, cholesky分解)

- LU分解

  - 4*4矩阵:

  - $$
    A = 
    \begin{bmatrix}
    4 & 2 & 2 & 1 \\
    2 & 18 & 3 & 2 \\
    5 & 3 & 8 & 1 \\
    6 & 2 & 2 & 5 \\
    \end{bmatrix}
    
    B = \begin{bmatrix}
    9 \\
    27 \\
    12 \\
    8 \\
    \end{bmatrix}
    $$

  - Ax=b => LUx = b => Ly=b, Ux = y

  - $$
    L = \begin{bmatrix}
    1 \\
    0.5 & 1 \\
    1.25 & 0.03 & 1 \\
    1.5 & -0.06 & -0.16 & 1\\
    \end{bmatrix}
    
    U = \begin{bmatrix}
    4 & 2 & 2 & 1\\
    0 & 17 & 2 & 1.5\\
    0 & 0 & 5.44 & -0.27 \\
    0 & 0 & 0 & 3.54\\
    \end{bmatrix}
    $$

  - $$
    y = \begin{bmatrix}
    9 \\
    22.5 \\
    0.07 \\
    -4.16 \\
    \end{bmatrix}
    
    x = \begin{bmatrix}
    1.85 \\
    1.43 \\
    -0.05 \\
    -0.18 \\
    \end{bmatrix}
    $$

- Cholesky分解
$$
  A = 
\begin{bmatrix}
  4 & 2 & 2 & 1 \\
2 & 18 & 3 & 2 \\
  2 & 3 & 8 & 1 \\
1 & 2 & 1 & 5 \\
  \end{bmatrix}
$$
  
- $A = LL^T, Ly = b, L^Tx = y$
  
  - $l_{11} = \sqrt{a11}=2, l_{2...4|1} = a_{2...4|1} / l_{11} = [1,1,0.5]$
  
  - $l_{22} = \sqrt{a_{22}-l_{21}^2} = 4.12, l_{32} = \frac{a_{32}-l_{31}*l_{21}}{l_{22}} = 0.49, l_{42} = \frac{a_{42}-l_{41}*l_{21}}{l_{22}} = 4.36$
  
  - $l_{33}=\sqrt{a_{33}-l_{31}^2-l_{32}^2} = 2.6, l_{43} = \frac{a_{43}-l_{41}*l_{31}-l_{42}*l_{32}}{l_{33}} = 0.12$
  
  - $l_{44} = \sqrt{a_{44}-l_{41}^2-l_{42}^2-l_{43}^2} = 2.15$
  
  - $$
    L = \begin{bmatrix}
    2 & 0 & 0&0\\
    1&4.12&0&0\\
    1&0.49&2.6&0\\
    0.5&0.36&0.12&2.15\\
    \end{bmatrix}
    
    x = \begin{bmatrix}
    1.13\\
    1.18\\
    0.64\\
    0.72\\
    \end{bmatrix}
    $$


## 2.迭代法解法(Jacobi, G-S)

- Jacobi:
  - 转换迭代公式, XJB往下迭代就行了.
- G-S:
  - 同jacobi, 迭代时用最新值

## 3.非线性方程求根(二分,固定点法,牛顿法)

- 二分, XJB算就行了
- 固定点, 把方程变换, 迭代求值($x^2 = 4 => x = x + 1 - \frac{x^2}{4}$
- 牛顿法(切线法)

## 4.插值(多项式插值, Lagrange, 牛顿插值)

- polynomial: 

  $A[i][j] = in[i]^j, b[i] = fin[i], res = sum(x[i] * x^i)$
  
- Lagrange:

  $x = [0, 2, 5]$

  $y = [0,3,1]$

  Prod1 = (x-2)*(x-5) / (0-2)(0-5)

- Newton:

  $A[i][j] = (A[i-1][j] - A[i-1][j+1]) / (in[j] - in[j+1])$

  $prodi = \prod_0^{i-1}{x-in[j]}, res += prodi * A[i][0]$

## 5.特征值与特征向量(power法)

$$
A = \begin{bmatrix}
4&2&2&1\\
2&18&3&2\\
2&3&8&1\\
1&2&1&5\\
\end{bmatrix}
$$

设: $x0 = [1,1,1,1]^T$, $x1 = Ax0, val = x1[0], x1/=x1[0]$... -> 收敛


## 6.数值积分(中点, 梯形)

- midpoint
- 梯形
- simpson13, 首尾 + 4 * 奇数点 + 2 * 偶数点

## 7.常微分方程初值问题(欧拉公式, 龙格库兹)

