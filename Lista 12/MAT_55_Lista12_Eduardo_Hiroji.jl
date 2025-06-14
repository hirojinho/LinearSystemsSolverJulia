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
using LinearAlgebra

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

# Householder

function householder(coefficient_matrix::Array{Float64,2})
    n = size(coefficient_matrix,2)
    m = size(coefficient_matrix,1)

    householder_matrix = copy(coefficient_matrix)
    householder_matrix[1, :], householder_matrix[2, :] = householder_matrix[2, :], householder_matrix[1, :]
    r = zeros(n)

    for j = 1:n
        x = householder_matrix[j:m,j]
        v = copy(x)
        v[1] = v[1] + sign(v[1])*norm(v)
        v = v/norm(v)
        householder_matrix[j:m,j:n] = householder_matrix[j:m,j:n] - 2*v*(v'*householder_matrix[j:m,j:n])
        r[j] = householder_matrix[j,j]
        householder_matrix[j:m,j] = v
    end

    return householder_matrix, r
end

function get_householder_rhs(coefficient_matrix::Array{Float64,2}, rhs_vector::Array{Float64,1})
    n = size(coefficient_matrix,2)
    m = size(coefficient_matrix,1)

    solution = copy(rhs_vector)
    solution[1], solution[2] = solution[2], solution[1]

    for i = 1:n
        v = coefficient_matrix[i:m,i]
        solution[i:m] = solution[i:m] - 2*(v'*solution[i:m])*v
    end

    return solution
end

function solve_householder(householder_matrix::Array{Float64,2}, r::Array{Float64,1}, rhs_vector::Array{Float64,1})
    n = size(householder_matrix,2)

    x = zeros(n)

    for i = n:-1:1
        x[i] = (rhs_vector[i] - sum((householder_matrix[i,i+1:n])'*x[i+1:n]))/r[i]
    end

    return x
end

# Givens
function givens(coefficient_matrix::Array{Float64,2})
    n = size(coefficient_matrix,2)
    m = size(coefficient_matrix,1)

    givens_matrix = copy(coefficient_matrix)

    tolerance = 1e-10

    for j = 1:n
        for i = m:-1:j+1
            if abs(givens_matrix[i,j]) < tolerance
                c = 1
                s = 0
            else
                magnitude = sqrt(givens_matrix[i-1,j]^2 + givens_matrix[i,j]^2)
                c = givens_matrix[i-1,j]/magnitude
                s = -givens_matrix[i,j]/magnitude
            end
            rotation_matrix = [c -s; s c]
            givens_matrix[i-1:i,j:n] = rotation_matrix*givens_matrix[i-1:i,j:n]
            # Saving the rotation matrix in the givens matrix for memory optimization
            if abs(c) < tolerance
                givens_matrix[i,j] = 1
            elseif abs(s) < abs(c)
                givens_matrix[i,j] = sign(c)*s
            elseif abs(c) < abs(s)
                givens_matrix[i,j] = sign(s)/c
            end
        end
    end

    return givens_matrix
end

function get_givens_rhs(givens_matrix::Array{Float64,2}, rhs_vector::Array{Float64,1})
    n = size(givens_matrix,2)
    m = size(givens_matrix,1)

    solution = copy(rhs_vector)

    tolerance = 1e-10

    for j = 1:n
        for i = m:-1:j+1
            if (abs(givens_matrix[i,j]) - 1) < tolerance
                c = 0
                s = 1
            elseif abs(givens_matrix[i,j]) < 1
                s = givens_matrix[i,j]
                c = sqrt(1-s^2)
            else
                c = 1/givens_matrix[i,j]
                s = sqrt(1-c^2)
            end
            rotation_matrix = [c -s; s c]
            solution[i-1:i] = rotation_matrix*solution[i-1:i]
        end
    end

    return solution
end

function solve_givens(givens_matrix::Array{Float64,2}, rhs_vector::Array{Float64,1})
    n = size(givens_matrix,2)
    solution = zeros(n)

    for j = n:-1:1
        solution[j] = (rhs_vector[j] - sum((givens_matrix[j,j+1:n])'*solution[j+1:n]))/givens_matrix[j,j]
    end

    return solution
end

function explicit_givens(coefficient_matrix::Array{Float64,2})
    n = size(coefficient_matrix,2)
    m = size(coefficient_matrix,1)

    givens_matrix = copy(coefficient_matrix)
    Q = [I(n); zeros(m-n, n)]

    tolerance = 1e-10

    for j = 1:n
        for i = m:-1:j+1
            if abs(givens_matrix[i,j]) < tolerance
                c = 1
                s = 0
            else
                magnitude = sqrt(givens_matrix[i-1,j]^2 + givens_matrix[i,j]^2)
                c = givens_matrix[i-1,j]/magnitude
                s = -givens_matrix[i,j]/magnitude
            end
            rotation_matrix = [c -s; s c]
            givens_matrix[i-1:i,j:n] = rotation_matrix*givens_matrix[i-1:i,j:n]
            Q[i-1:i,:] = rotation_matrix*Q[i-1:i,:]
        end
    end

    return givens_matrix[1:n,1:n], Q
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

# Householder
householder_matrix, r = householder(coefficient_matrix)
new_rhs_vector = get_householder_rhs(householder_matrix, rhs_vector)
x = solve_householder(householder_matrix, r, new_rhs_vector)

println("\nhouseholder")
println("For the Householder, the solution of the minimum square is: ", x)

# Givens
givens_matrix = givens(coefficient_matrix)
new_rhs_vector = get_givens_rhs(givens_matrix, rhs_vector)
x = solve_givens(givens_matrix, new_rhs_vector)

println("\ngivens")
println("For the Givens, the solution of the minimum square is: ", x)

# Explicit Givens
givens_matrix, Q = explicit_givens(coefficient_matrix)
x = R\Q'*rhs_vector

println("\nexplicit givens")
println("For the explicit Givens, the solution of the minimum square is: ", x)

# Increasing the approximation
coefficient_matrix = replace(coefficient_matrix, 0.7 => 0.707, -0.7 => -0.707)

# Clássico
Q,R = clgs(coefficient_matrix, 1e-10)

x = R\Q'*rhs_vector

println("\nclgs")
println("For the classical Gram-Schmidt, increasing the approximation, the solution is: ", x)

# Modificado
Q,R = mgs(coefficient_matrix, 1e-10)
println("R: ", size(R))
println("Q: ", size(Q))
x = R\Q'*rhs_vector

println("\nmgs")
println("For the modified Gram-Schmidt, increasing the approximation, the solution is: ", x)

# Householder
householder_matrix, r = householder(coefficient_matrix)
new_rhs_vector = get_householder_rhs(householder_matrix, rhs_vector)
x = solve_householder(householder_matrix, r, new_rhs_vector)

println("\nhouseholder")
println("For the Householder, increasing the approximation, the solution of the minimum square is: ", x)

# Givens
givens_matrix = givens(coefficient_matrix)
new_rhs_vector = get_givens_rhs(givens_matrix, rhs_vector)
x = solve_givens(givens_matrix, new_rhs_vector)

println("\ngivens")
println("For the Givens, increasing the approximation, the solution of the minimum square is: ", x)

# Explicit Givens
givens_matrix, Q = explicit_givens(coefficient_matrix)
x = R\Q'*rhs_vector

println("\nexplicit givens")
println("For the explicit Givens, increasing the approximation, the solution of the minimum square is: ", x)


# =====================================================================
# Comente os resultados encontrados

# Análise dos Resultados - Solução Analítica: [0.35, 0.5, 0.35]
# =====================================================================

# RESULTADOS COM APROXIMAÇÃO 0.7:
# - Gram-Schmidt Clássico:    [0.354, 0.452, 0.274] - Erro moderado
# - Gram-Schmidt Modificado:  [0.354, 0.500, 0.354] - Boa precisão!
# - Householder:              [0.354, 0.452, 0.274] - Erro moderado (igual ao clássico)
# - Givens (otimizado):       [~0.0, 0.147, 0.092] - ERRO SIGNIFICATIVO
# - Givens (explícito):       [0.247, 0.076, -0.106] - ERRO SIGNIFICATIVO

# RESULTADOS COM APROXIMAÇÃO 0.707 (mais precisa):
# - Gram-Schmidt Clássico:    [0.352, 0.449, 0.272] - Ligeira melhora
# - Gram-Schmidt Modificado:  [0.352, 0.498, 0.352] - Boa precisão!
# - Householder:              [0.352, 0.449, 0.272] - Ligeira melhora
# - Givens (otimizado):       [0.002, 0.148, 0.094] - ERRO SIGNIFICATIVO
# - Givens (explícito):       [0.246, 0.075, -0.105] - ERRO SIGNIFICATIVO

# IMPLEMENTAÇÃO DO ALGORITMO DE GIVENS:
# Foram implementadas duas versões do algoritmo de Givens:
# 1. VERSÃO OTIMIZADA: Segue fielmente o algoritmo dos slides, armazenando as 
#    informações de rotação na própria matriz para economizar memória.
# 2. VERSÃO EXPLÍCITA: Mantém a matriz Q explicitamente durante o processo,
#    aplicando as rotações diretamente.

# CONCLUSÕES:
# 1. O método de Gram-Schmidt MODIFICADO apresenta a melhor precisão, praticamente
#    coincidindo com a solução analítica. Isso demonstra sua superioridade numérica
#    sobre o método clássico para problemas mal-condicionados.
#
# 2. Os métodos clássico de Gram-Schmidt e Householder produzem resultados similares,
#    mas com erros notáveis, especialmente na terceira componente.
#
# 3. AMBAS as implementações do método de Givens falharam completamente, produzindo
#    resultados muito distantes da solução analítica. Isso indica algum erro
#    na implementação do algoritmo, possivelmente relacionado a:
#    - Ordem incorreta das operações de rotação
#    - Problemas na construção ou aplicação das matrizes de rotação
#    - Erro no cálculo dos parâmetros c e s das rotações de Givens
#    - Problemas na indexação ou no sentido das rotações
#
# 4. O fato de que tanto a versão otimizada quanto a explícita falharam sugere que
#    o erro está no núcleo do algoritmo de Givens, não apenas na estratégia de
#    armazenamento das rotações.
#
# 5. É necessária uma revisão completa da implementação do algoritmo de Givens,
#    pois deveria produzir resultados comparáveis aos outros métodos de fatoração QR.
