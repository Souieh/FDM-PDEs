# Computation of convection of vector field using finite difference method

This folder contains the source code of finite difference method for solving the convection vector field of the this two problems  

## Dependices 
    latex 
    matplotlib        3.6.2
    numpy             1.23.4

## first Problem

```math
\left\{
    \begin{array}{ll}
        -div(A)\nabla \zeta_k  = div(A e_k)  &, k=1,2  \;\;,\;\;in\; Y_s \\
        A\nabla \zeta_k \cdot v = A  e_k \cdot v &, k=1,2 \\
        \zeta  : Y-periodic &
    \end{array}
    \right.
```


## second Problem

```math
\left\{
\begin{array}{l}
    -div(A)\nabla \gamma  = 0   \;\;,\;\;in\; Y_s \\
    A\nabla \gamma_k \cdot v = -\alpha(y)  \\
    \gamma  : Y-periodic 
\end{array}
\right.
```

## 

```math
Y=]0,1[^2 \text{ , } T=[1/4,3/4]^2\text{ , } Y_s = Y-T \text{ , } \Sigma=\partial Ys \text{ , } 
```
```math
\alpha = 1 ,A 
\begin{bmatrix}
1 & 0 \\
1& 1
\end{bmatrix}
```
