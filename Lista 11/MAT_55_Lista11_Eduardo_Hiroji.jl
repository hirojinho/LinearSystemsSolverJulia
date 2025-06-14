# =====================================================================
# *********************************************************************
#                          MAT-55 2025 - Lista 11
# *********************************************************************
# =====================================================================
# Dupla:
#
#
# =====================================================================
#                   Método de Gram-Schmidt Clássico
# =====================================================================
# Dados de entrada:
# a     matriz, mxn com posto completo
# tol	escalar, tolerância 
#
# Saída:
# Q,R   fatoração QR reduzida da matriz a

function clgs(coefficient_matrix::Array{Float64,2},tol::Float64)
    ##################
    #Digite seu código aqui
    rank = size(coefficient_matrix,2)
    m = size(coefficient_matrix,1)
    Q = zeros(m,rank)
    R = zeros(rank,rank)

    for j in 1:rank
        v = coefficient_matrix[:,j]
        for i in 1:j-1
            R[i,j] = Q[:,i]'*coefficient_matrix[:,j]
            v = v - R[i,j]*Q[:,i]
        end
        R[j,j] = norm(v)
        if R[j,j] < tol
            error("A matriz não é de posto completo")
        end
        Q[:,j] = v/R[j,j]
    end
    ##################
    return Q, R
end

# =====================================================================
#                   Método de Gram-Schmidt Modificado 
# =====================================================================
# Dados de entrada:
# a     matriz, mxn com posto completo
# tol	escalar, tolerância 
#
# Saída:
# Q,R   fatoração QR reduzida da matriz a

function mgs(coefficient_matrix::Array{Float64,2},tol::Float64)
    ##################
    #Digite seu código aqui
    rank = size(coefficient_matrix,2)   
    m = size(coefficient_matrix,1)
    Q = zeros(m,rank)
    R = zeros(rank,rank)

    for j in 1:rank
        v = coefficient_matrix[:,j]
        for i in 1:j-1
            R[i,j] = Q[:,j]'*v
            v = v - R[i,j]*Q[:,j]
        end
        R[j,j] = norm(v)
        if R[j,j] < tol
            error("A matriz não é de posto completo")
        end
        Q[:,j] = v/R[j,j]
    end
    ##################
    return Q,R
end

function householder(coefficient_matrix::Array{Float64,2})
    m = size(coefficient_matrix,1)
    n = size(coefficient_matrix,2)

    for j = 1:n
        x = coefficient_matrix[j:m,j]
        v = x
        v[1] = v[1] + sign(v[1])*norm(v)
        v = v/norm(v)
        coefficient_matrix[j:m,j:n] = coefficient_matrix[j:m,j:n] - 2*v'*coefficient_matrix[j:m,j:n]*v
    end

    return coefficient_matrix
end

function norm(v::Array{Float64,1})
    return sqrt(sum(v.^2))
end

# =====================================================================
# Exercício 7
# =====================================================================
#Dados do problema

coefficient_matrix = [
    0.0  0.7  1.0;
   -0.7  0.0  0.7;
   -1.0 -0.7  0.0;
   -0.7 -1.0 -0.7;
    0.0 -0.7 -1.0;
    0.7  0.0 -0.7;
    1.0  0.7  0.0;
    0.7  1.0  0.7;
    0.0 -0.7  1.0;
    0.7  0.0 -0.7;
   -1.0  0.7  0.0;
    0.7 -1.0  0.7;
    0.0 -0.7 -1.0;
   -0.7  0.0  0.7;
    1.0 -0.7  0.0;
   -0.7  1.0 -0.7;
]

rhs_vector = [
    0.7;
    0.0;
   -0.7;
   -1.0;
   -0.7;
    0.0;
    0.7;
    1.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
]

# Clássico
Q,R = clgs(coefficient_matrix, 1e-10)

x = R\Q'*rhs_vector

println("\nclgs")
println("For the classical Gram-Schmidt, the solution is: ", x)

# Modificado
Q,R = mgs(coefficient_matrix, 1e-10)

x = R\Q'*rhs_vector

println("\nmgs")
println("For the modified Gram-Schmidt, the solution is: ", x)

# Increasing the approximation
coefficient_matrix = replace(coefficient_matrix, 0.7 => 0.707, -0.7 => -0.707)

# Clássico
Q,R = clgs(coefficient_matrix, 1e-10)

x = R\Q'*rhs_vector

println("\nclgs")
println("For the classical Gram-Schmidt, increasing the approximation, the solution is: ", x)

# Modificado
Q,R = mgs(coefficient_matrix, 1e-10)

x = R\Q'*rhs_vector

println("\nmgs")
println("For the modified Gram-Schmidt, increasing the approximation, the solution is: ", x)

# Householder
coefficient_matrix = householder(coefficient_matrix)

println("\nhouseholder")
println("For the Householder, increasing the approximation, the solution is: ", coefficient_matrix)

# =====================================================================
# Comente os resultados encontrados





