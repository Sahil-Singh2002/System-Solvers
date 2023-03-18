"""
MATH2019 CW2 systemsolvers module

@author: Sahil Singh
"""

import numpy as np
import matplotlib.pyplot as plt
import warmup_solution as ws

def no_pivoting(A,b,n,c):
    
    """
    Returns an array representing the augmented matrix M arrived at by
    starting from the augmented matrix [A b] and performing forward
    elimination without row interchanges until all of the entries below the
    main diagonal in the first c columns are 0.
    
    Parameters
    ----------
    A : numpy.ndarray of shape (n,n)
        array representing the square matrix A.
    b : numpy.ndarray of shape (n,1)
        array representing the column vector b.
    n : integer
        integer that is at least 2.
    c : integer
        positive integer that is at most n-1.
    
    Returns
    -------
    M : numpy.ndarray of shape (n,n+1)
        array representing the augmented matrix M.
    """
    
    # Create the initial augmented matrix
    M=np.hstack((A,b))
    
    # Continue here:...
    for i in range(c):
        for j in range(i+1,n):
            m=M[j,i]/M[i,i]
            M[j,i]=0
            M[j,i+1:n+1]=M[j,i+1:n+1]-m*M[i,i+1:n+1]
    
    return M

def backward_substitution(M,n):
    
    """
    Returns an array representing the solution x of Ux=v computed using
    backward substitution where U is an upper triangular matrix.
    
    Parameters
    ----------
    M : numpy.ndarray of shape (n,n+1)
        Array representing the augmented matrix [U v].
    n : integer
        integer that is at least 2.
        
    Returns
    -------
    x : numpy.ndarray of shape (n,1)
        array representing the solution x.
    """
    
    
    x=np.zeros([n,1])
    
    x[n-1]=M[n-1,n]/M[n-1,n-1]
    
    for i in range(n-2,-1,-1):
        x[i]=(M[i,n]-np.dot(M[i,i+1:n],x[i+1:n]))/M[i,i]
    
    return x

def no_pivoting_solve(A,b,n):
    
    """
    Returns an array representing the solution x to Ax=b computed using
    forward elimination with no pivoting followed by backward substitution.
    
    Parameters
    ----------
    A : numpy.ndarray of shape (n,n)
        array representing the square matrix A.
    b : numpy.ndarray of shape (n,1)
        array representing the column vector b.
    n : integer
        integer that is at least 2.
    
    Returns
    -------
    x : numpy.ndarray of shape (n,1)
        array representing the solution x.
    """
    
    M=no_pivoting(A,b,n,n-1)
    x=backward_substitution(M,n)
    
    return x

def find_max(M,s,n,i):
    """
    
    Returns an integer representing the row number for which
    when we have the augmented matrix [A b] as M we take the coloumn i
    where it take the row of M and s takes the ratio of the absolute and 
    and the biggest ratio will be the row which will be returned as p.    
    
    Parameters
    ----------
    M : numpy.ndarray of shape (n,n+1)
        array representing the augmented matrix M.
    s : numpy.ndarray of shape (n,)
        array representing the column vector s where all it's elements are positive.
    n : integer
        integer that is at least 2.
    i : integer
        positive integer that is at most n-2.
    
    Returns
    -------
    p : integer
        positive integer that is between 0 and n-1
    """
    if np.shape(M) != (n,n+1):
        raise ValueError("The dimensions of the matrix should be n x (n+1) matrix ")
    if n < 2 or type(n) != int:
        raise ValueError("The input digit must be a positive Integer at least 2")
    if i > n-2 or type(i)!=int:
        raise ValueError("The input for i needs to be positive Integer at most n-2")
    
    m,p = 0,0
    p_list = []
    for j in range(n):
        # Iterate and store
        p_list.append(abs(M[j,i]/s[j]))
        if p_list[j] > m:
            m = p_list[j]
            p = j
    return p 

def scaled_partial_pivoting(A,b,n,c):
    """
    Returns an array representing the augmented matrix M arrived at by
    starting from the augmented matrix [A b] and performing forward
    elimination with scaled partial pivoting until all of the entries below the
    main diagonal in the first c columns are 0.
    
    Parameters
    ----------
    A : numpy.ndarray of shape (n,n)
        array representing the square matrix A.
    b : numpy.ndarray of shape (n,1)
        array representing the column vector b.
    n : integer
        integer that is at least 2.
    c : integer
        positive integer that is at most n-1.
    
    Returns
    -------
    M : numpy.ndarray of shape (n,n+1)
        array representing the augmented matrix M.
    """
    if np.shape(A) != (n,n):
        raise ValueError("The dimensions of the matrix should be n x n matrix ")
    if np.shape(b) != (n,1):
        raise ValueError("The dimensions of the vector b should be n x 1 matrix ")
    if n < 2 or type(n) != int:
        raise ValueError("The input digit must be a positive Integer at least 2")
    if c > n-1 or type(n) != int:
        raise ValueError("The input digit must be a positive Integer at most n-1")
    
    s=np.amax(np.abs(A),1)
    M=np.hstack((A,b))
 
    for i in range (c):
        # function to find smallest value and itds row position to do the swap
        p = find_max(M,s,n,i)
        s[[i,p]]=s[[p,i]]
        M[[i,p],:]=M[[p,i],:]
        
        for j in range(i+1,n):
            # Iterate
            m=M[j,i]/M[i,i]
            M[j,i]=0
            M[j,i+1:n+1]=M[j,i+1:n+1]-m*M[i,i+1:n+1]
    return M

def spp_solve(A,b,n):
    
    """
    Returns an array representing the solution x to Ax=b computed using
    forward elimination with scaled partial pivoting followed by backward substitution.
    
    Parameters
    ----------
    A : numpy.ndarray of shape (n,n)
        array representing the square matrix A.
    b : numpy.ndarray of shape (n,1)
        array representing the column vector b.
    n : integer
        integer that is at least 2.
    
    Returns
    -------
    x : numpy.ndarray of shape (n,1)
        array representing the solution x.
    """
    if np.shape(A) != (n,n):
        raise ValueError("The dimensions of the matrix should be n x n matrix ")
    if np.shape(b) != (n,1):
        raise ValueError("The dimensions of the vector b should be n x 1 matrix ")
    if n < 2 or type(n) != int:
        raise ValueError("The input digit must be a positive Integer at least 2")
    M = scaled_partial_pivoting(A,b,n,n-1)
    x = backward_substitution(M,n)
    return x

def PLU(A,n):
    """
    Returns arrays P L U representing .
    
    Parameters
    ----------
    A : numpy.ndarray of shape (n,n)
        array representing the square matrix A.
    n : integer
        integer that is at least 2.
    
    Returns
    -------
    P : numpy.ndarray of shape (n,n)
        array representing permutation matrix P.
    L : numpy.ndarray of shape (n,n)
        array representing lower triangular matrix L.
    U : numpy.ndarray of shape (n,n)
        array representing the upper triangular matrix U.
        
        
    """
    if n < 2 or type(n) != int:
        raise ValueError("The input digit must be a positive Integer at least 2")
    if np.shape(A) != (n,n):
        raise ValueError("The dimensions of the matrix should be n x n matrix ")
    # Set P=I
    P=np.identity(n)
    # Set L to be a zero matrix
    L=np.zeros([n,n])
    # Set U=A
    U=A.copy()
    # Start loop
    for i in range(n):
        # condition just helps us stop the for loop once we get a 
        # value which is > 0. we treat it as more of a boolion
        condition = True 
        for row_n in range(i,n) :
            # Stopping criterion
            if condition is True:
                if abs(U[row_n][i])>10**(-15):
                    condition = False
                    s = row_n
        # row swap criterion
        if s != i:
            P[[i,s],:]=P[[s,i],:]
            L[[i,s],:]=L[[s,i],:]
            U[[i,s],:]=U[[s,i],:] 
        for j in range(i+1,n):
            L[j,i] = U[j,i]/U[i,i]
            U[j,i] = 0
            for k in range(i+1,n):
                U[j,k] = U[j,k] - L[j,i]*U[i,k] 
    P = np.transpose(P)
    L = L + np.identity(n)   
    return P, L, U

def Jacobi(A,b,n,x0,N):
    
    # Continue here:...    
    x_approx = np.zeros([n,N+1])
    
    for k in np.arange(0,n):
        x_approx[k,0]= x0[k]
    for i in range(0,N):   
        for j in range(0,n):
            x_approx[j,i+1] = (b[j] - np.dot( a[j,0:j], x_approx[0:j,i] ) - np.dot(A[j,j+1:n],x_approx[j+1:n,i])) /A[j,j]
    return x_approx

def Jacobi_plot(A,b,n,x0,N):

    # Create array of k values
    k_array = np.arange(N+1)
    # Prepare figure
    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.set_xlabel("$k$")
    ax.grid(True)
    
    # Continue here:...
    x_approx = Jacobi(A,b,n,x0,N)
    x_array = no_pivoting_solve(A,b,n)
    #create an array with zeros inside of it so its just empty when we start
    list_i, list_2 = np.zeros(N+1), np.zeros(N+1)
      
    for k in range(N+1):
        h = x[0:n,0] -x_approx[0:n,i]
        error_2[i] = np.linalg.norm(h,2)
        error_i[i] = np.linalg.norm(h,np.inf)
        
    
    
    # Plot
    ax.plot(k_array,error_i,"s",label="$||x-x^{(k)}||_\infty$",linestyle ='--')
    ax.plot(k_array,error_2,"o",label="$||x-x^{(k)}||_2$",linestyle ='-.')
    # Add legend
    ax.legend()
    return fig, ax
