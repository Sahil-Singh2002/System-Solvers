# System-Solvers
Course work 2, has combination of Algorithms such as Gaussian Elimination, and pivoting Method in Matrics.

In this code, we are given a set a task to built algorithms which call on each other to solve problem related to linear-systems.

The First function we built is designed to help find the smallest integer p such that p is the smallest i and the absolute of M[p,i] over s[p] is the maximum over j in {i,i+1,...,n-1} of |M[j,i]|/s[j] 

The second function scaled partial pivoting is used to return an array representing the augmented matrix M at the start from the augmented matrix [A b] and then procced to perform a forwward elimintation with scaled partial pivoting until all of the entries below the main diagonal in th first c columns are 0.

The third function spp_solve is the final step for finding the solution x to Ax = b computed using forward elimination with partial pivoting followed by backward substitution. I till return x which is a numpy.ndarray of shape (n,1) array.  

The fourth function called PLU returns arrays representing a permutation matrix P, a lower triangular matrix L and an upper triangular matrix U such that A=PLU. With the parameters being A which is our matrix A and integer n which is the for the matrix A which dimension (n,n).

The fifth function Jocobi returns an array of approximations to the solution of Ax=b obtained using the Jacobi method. The Parameters A array representing the matrix A of shape (n,n), b array representing the vector b of shape (n,1), n integer that is at least 2, x0 array representing the initial approximation x0 of shape (n,1) and N is the number of iterations to be performed. This returns x_approx array whose column 0 is x0 and, for i=1,2,...,N, whose column i is the approximation to the solution of Ax=b obtained after performing i iterations of the Jacobi method starting from x0 of shape (n,N+1).

The sixth function Jocobi plot plots the convergence of norm(x - x^k)_inf and norm(x-x^k)_2.

    $||x-x^{(k)}||_\infty$
    $||x-x^{(k)}||_2$
