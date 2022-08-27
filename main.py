import numpy as np
import math
from scipy.linalg import solve

# Checking if the entered matrix is minimum of 10 by 10 matrix.
def isLessthan10x10(arr):
    if(len(arr) < 10):
        return True;
    return False

def gaussElimination(a,c,x ): # a is coefficient matrix, c is right-hand side vector, x is solution.

   #pivotting part
    a=np.array((a),dtype=float)
    c=np.array((c),dtype=float)
    n = len(c)
    for i in range(0,n-1):     # Looping through the columns
        
        if np.abs(a[i,i])==0:
            for k in range(i+1,n):
                if np.abs(a[k,i])>np.abs(a[i,i]):
                    a[[i,k]]=a[[k,i]]             # Swapping rows.
                    c[[i,k]]=c[[k,i]]
                    break
                    
        for j in range(i+1,n):     # Looping through diagonals.
            m = a[j,i]/a[i,i]
            a[j,:] = a[j,:] - m*a[i,:]
            c[j] = c[j] - m*c[i]
    

    #back-subs part
          
    x[n-1] = c[n-1]/a[n-1,n-1]    # Solve for the last entry
    for i in range(n-2,-1,-1):      # Looping from last value to initial value.
        summed = 0
        for j in range(i+1,n):        # For X's, sum and move to right hand side.
            summed = summed + a[i,j]*x[j]
        x[i] = (c[i] - summed)/a[i,i]
    return x


def isDominant(a): # Checks whether the matrix is diagonally dominant or not for the Gauss-Siedel Method.
    abs = np.abs(a)
    d= np.diag(abs)
    wo_d = np.sum(abs, axis= 1) - d
    if np.all(d > wo_d):
        return True
    else:
        return False


def gaussSeidel(a, x ,c): # a is coefficient matrix, c is right-hand side vector, x is solution.
    a=np.array((a),dtype=float)
    c=np.array((c),dtype=float) 
    n = len(a)

    if isDominant(a) == False:
        print("This matrix is not diagonally dominant thus Gauss-Siedel Method cannot be applied for the value you have entered.")
   
    else:
        for j in range(0, n):        #Looping n times.
            d = c[j]                 #Temporarily defining an array for C array.
            for i in range(0, n):    #Calculating ith value of array
                if(j != i):
                    d-=a[j][i] * x[i]
            x[j] = d / a[j][j]       #Complete the first iteration process
        return x                     #Returning the X value


#Command part using the functions defined above:

a = eval (input("Please enter the coefficient array (2d array form): ")) #taking inputs and creating related arrays                         
c = eval (input("Please enter the right-hand side vector (1d array form ): "))    
n = len(a)                
x =  np.zeros(n, float)    # Creating initial x which has the same size with coefficient array and full of zeros in order to initialize the iteration.
iteration_number = 100 # The value in which the Gauss-Siedel method is iterated.

if isLessthan10x10(a) or isLessthan10x10(c):
    print("Please enter a matrix at least 10x10. Re-run the script.") #check the 10x10 condition of both A and C arrays.
    exit()
else:
    for i in range(0, iteration_number):            
        x = gaussSeidel(a, x, c)        
    print("Gauss Seidel result: ", x)
    x =  np.zeros(n, float)  # Reinitialize x for other matrix operation.
    x =  gaussElimination(a, c, x)
    print("Gauss Elimination method with pivoting result: ", x)
    x = np.zeros(n,float) # Reinitialize x for other matrix operation.
    x = np.linalg.solve(a, c) #Numpy Solver
    print ("Solution with Numpy solver to check the correctness of the algoritm coded by me: \n", x)