import math

def identity_matrix(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]

def matrix_multiply(A, B):
    n = len(A)
    result = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                result[i][j] += A[i][k] * B[k][j]
    return result

def transpose_matrix(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A))]

def matrix_vector_multiply(matrix, vector):
    return [sum(matrix[i][j] * vector[j] for j in range(len(vector))) for i in range(len(matrix))]

def max_off_diagonal(A):
    n = len(A)
    max_val = 0
    p, q = 0, 0
    for i in range(n):
        for j in range(i+1, n):
            if abs(A[i][j]) > abs(max_val):
                max_val = A[i][j]
                p, q = i, j
    return p, q, max_val

def jacobi_method(A, epsilon=1e-4):
    n = len(A)
    U = identity_matrix(n)
    iterations = 0
    
    while True:
        p, q, max_val = max_off_diagonal(A)
        
        if abs(max_val) < epsilon:
            break
        
        if A[p][p] == A[q][q]:
            phi = math.pi / 4
        else:
            phi = 0.5 * math.atan(2 * A[p][q] / (A[p][p] - A[q][q]))
        
        U_pq = identity_matrix(n)
        U_pq[p][p] = math.cos(phi)
        U_pq[q][q] = math.cos(phi)
        U_pq[p][q] = -math.sin(phi)
        U_pq[q][p] = math.sin(phi)
        
        U_pq_T = transpose_matrix(U_pq)
        A = matrix_multiply(matrix_multiply(U_pq_T, A), U_pq)
        
        U = matrix_multiply(U, U_pq)
        
        iterations += 1

        # Виведення на кожній ітерації
        print(f"\nIteration {iterations}:")
        print(f"Max off-diagonal element: A[{p}][{q}] = {max_val}")
        print(f"Rotation angle phi: {phi}")
        print("Updated matrix A:")
        for row in A:
            print(["{}".format(x) for x in row])
        print("Current eigenvalues :", [A[i][i] for i in range(n)])
        
        if iterations > 10000:
            break

    eigenvalues = [A[i][i] for i in range(n)]
    eigenvectors = [[U[i][j] for i in range(n)] for j in range(n)]
    
    return eigenvalues, eigenvectors

A = [[5.18, 1.12, 0.95, 1.32, 0.83],
     [1.12, 6.28, 2.12, 0.57, 0.91],
     [0.95, 2.12, 6.13, 1.29, 1.57],
     [1.32, 0.57, 1.29, 5.67, 1.25],
     [0.83, 0.91, 1.57, 1.25, 5.21]
]

eigenvalues, eigenvectors = jacobi_method(A)

print("\nEigenvalues and eigenvectors:")
for i in range(len(eigenvalues)):
    transposed_vector = eigenvectors[i][::-1]
    print(f"Eigenvalue {eigenvalues[i]:.4f} - Eigenvector: {transposed_vector}")
