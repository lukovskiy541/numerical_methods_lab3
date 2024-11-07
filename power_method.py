import math


A = [[5.18, 1.12, 0.95, 1.32, 0.83],
     [1.12, 6.28, 2.12, 0.57, 0.91],
     [0.95, 2.12, 6.13, 1.29, 1.57],
     [1.32, 0.57, 1.29, 5.67, 1.25],
     [0.83, 0.91, 1.57, 1.25, 5.21]
]

b = [1]*len(A)

def check_eigenvector(A, eigenvalue, eigenvector):
    Av = matrix_vector_multiply(A, eigenvector)  
    scaled_vector = [eigenvalue * x for x in eigenvector]  
    
    print("A * v = ", Av)
    print("lambda * v = ", scaled_vector)
    
  
    difference = [Av[i] - scaled_vector[i] for i in range(len(Av))]
   
    
    return difference


def dot_product(v1, v2):
    return sum(x * y for x, y in zip(v1, v2))

def vector_norm(v):
    return math.sqrt(sum(x * x for x in v))

def matrix_vector_multiply(matrix, vector):
    return [dot_product(row, vector) for row in matrix]

def power_method(A, b, epsilon=1e-4, max_iterations=1000):
    lambda_k = 0
    b_k = b[:]
    
    for i in range(max_iterations):
        print(i)
        b_k1 = matrix_vector_multiply(A, b_k)
  
        m = 0
     
        lambda_k1 = b_k1[m] / b_k[m]
        print('lambda_k1: ', lambda_k1)
        norm_b_k1 = vector_norm(b_k1)
        
        b_k1 = [x / norm_b_k1 for x in b_k1]
 
        if abs(lambda_k1 - lambda_k) <= epsilon:
            return lambda_k1, b_k1
        
        lambda_k = lambda_k1
        
        b_k = b_k1[:]
    
    return lambda_k1, b_k1


def matrix_infinity_norm(A):
    return max(sum(abs(x) for x in row) for row in A)

def build_matrix_B(A, norm_A_inf):
    size = len(A)
    B = [[(norm_A_inf if i == j else 0) - A[i][j] for j in range(size)] for i in range(size)]
    return B

def check_determinants(A):
    n = len(A)
    for i in range(1, n + 1):
        submatrix = [row[:i] for row in A[:i]]
        det = determinant(submatrix)
        print(f"Det({i}) = {det:.2f}")
        if det <= 0:
            return False
    return True

def determinant(A):
    n = len(A)
    if n == 1:
        return A[0][0]
    if n == 2:
        return A[0][0] * A[1][1] - A[0][1] * A[1][0]
    det = 0
    for j in range(n):
        submatrix = [row[:j] + row[j+1:] for row in A[1:]]
        sign = (-1) ** j
        det += sign * A[0][j] * determinant(submatrix)
    return det




print("\nПеревірка визначників головних мінорів:")
if check_determinants(A):
    print("Всі головні мінори додатні")
else:
    print("Не всі головні мінори додатні")



assessment_above = matrix_infinity_norm(A)


B = build_matrix_B(A,assessment_above)

print('B: ', B)


print('Assessment_above: ',assessment_above)
    
lambda_max_B, vector_max_B = power_method(B, b)

print('Lambda max B: ', lambda_max_B)

lambda_max_A, vector_max_A = power_method(A,b)

print('Lambda max A: ', lambda_max_A)
print('Eigenvector: ', vector_max_A  )
print('\n')

print("Check for eigenvector lambda max A")
check_eigenvector(A, lambda_max_A, vector_max_A)
print('\n')
lambda_min_A = assessment_above - lambda_max_B

print('Lambda min A: ', lambda_min_A, )
print('Eigenvector: ', vector_max_B )
print('\n')
print("Check for eigenvector lambda min A")
check_eigenvector(A, lambda_min_A, vector_max_B)



