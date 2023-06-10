# System-Solvers

This repository contains code for solving linear systems using various algorithms such as Gaussian Elimination and pivoting methods in matrices. The provided algorithms are designed to call on each other to solve problems related to linear systems.

## Algorithm 1: Finding the Smallest Integer

The first function is designed to find the smallest integer $p$ such that $p$ is the smallest $i$ and the absolute value of $M[p,i]$ over $s[p]$ is the maximum over $j$ in $\{i, i+1, \ldots, n-1\}$ of $\left|\frac{M[j,i]}{s[j]}\right|$.

## Algorithm 2: Scaled Partial Pivoting

The second function, scaled partial pivoting, is used to convert the augmented matrix $[A \, b]$ into an array representing the augmented matrix $M$. It then performs a forward elimination with scaled partial pivoting until all entries below the main diagonal in the first $c$ columns are 0.

## Algorithm 3: Solving Linear Systems

The third function, spp_solve, is the final step in finding the solution $x$ to $Ax = b$. It computes the solution using forward elimination with partial pivoting followed by backward substitution. The output is an `numpy.ndarray` called `x` with shape $(n,1)$.

## Algorithm 4: PLU Decomposition

The fourth function, PLU, returns arrays representing a permutation matrix $P$, a lower triangular matrix $L$, and an upper triangular matrix $U$ such that $A = PLU$. The parameters are $A$, which is the matrix $A$, and an integer $n$, which is the dimension of matrix $A$ $(n \times n)$.

## Algorithm 5: Jacobi Method

The fifth function, Jacobi, returns an array of approximations to the solution of $Ax=b$ obtained using the Jacobi method. The parameters are $A$, an array representing the matrix $A$ of shape $(n,n)$, $b$, an array representing the vector $b$ of shape $(n,1)$, $n$, an integer that is at least 2, $x_0$, an array representing the initial approximation $x_0$ of shape $(n,1)$, and $N$, the number of iterations to be performed. This returns an `x_approx` array, whose column 0 is $x_0$ and, for $i=1,2,\ldots,N$, whose column $i$ is the approximation to the solution of $Ax=b$ obtained after performing $i$ iterations of the Jacobi method starting from $x_0$. The shape of `x_approx` is $(n,N+1)$.

## Algorithm 6: Convergence Plotting

The sixth function, Jacobi_plot, plots the convergence of $\|x - x^{(k)}\|_{\infty}$ and $\|x-x^{(k)}\|_2$. It uses `matplotlib.pyplot` to create a plot illustrating the convergence behavior of the Jacobi method.

These algorithms provide a set of tools for solving linear systems efficiently and accurately.

