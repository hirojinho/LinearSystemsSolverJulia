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
using Plots
using Statistics

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
            R[i,j] = Q[:,i]'*v
            v = v - R[i,j]*Q[:,i]
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

# =====================================================================
#                   Método de Cholesky
# =====================================================================
# Dados de entrada:
# coefficient_matrix     matriz, mxn com posto completo
#
# Saída:
# lower_matrix   matriz triangular inferior

# ---------------------------------------------------------------------
function cholesky_decomposition(coefficient_matrix::Array{Float64,2})
    tol = 1e-12
    n = size(coefficient_matrix, 1)
    lower_matrix::Array{Float64,2} = build_identity_matrix(n)
    for current_row in 1:n
        sum_of_previous_elements = 0
        for current_column in 1:current_row-1
            sum_of_previous_elements += coefficient_matrix[current_row, current_column]^2
        end
        if coefficient_matrix[current_row, current_row] - sum_of_previous_elements < tol
            error("Matriz é muito mal condicionada para a decomposição de Cholesky")
        end
        coefficient_matrix[current_row, current_row] = sqrt(coefficient_matrix[current_row, current_row] - sum_of_previous_elements)
        for lower_row in current_row+1:n
            sum_of_previous_elements = 0
            for lower_column in 1:current_row-1
                sum_of_previous_elements += coefficient_matrix[lower_row, lower_column] * coefficient_matrix[current_row, lower_column]
            end
            coefficient_matrix[lower_row, current_row] = (coefficient_matrix[lower_row, current_row] - sum_of_previous_elements) / coefficient_matrix[current_row, current_row]
        end
    end

    # TODO: em vez de copiar, criar um util que resolve a matriz triangular inferior
    for i in 1:n, j in 1:i
        lower_matrix[i,j] = coefficient_matrix[i,j]
    end

    return lower_matrix
end

function cholesky_decomposition_solver(coefficient_matrix::Array{Float64,2}, constant_vector::Array{Float64,1})::Array{Float64,1}
    lower_matrix::Matrix{Float64} = cholesky_decomposition(coefficient_matrix)
    println("lower_matrix dimensions: ", size(lower_matrix))
    if isnothing(lower_matrix)
        println("Cholesky decomposition failed")
        return nothing
    end
    partial_solution::Array{Float64,1} = sub_direta(lower_matrix, constant_vector)
    transposed_lower_matrix::Matrix{Float64} = Matrix{Float64}(lower_matrix')
    solution::Array{Float64,1} = sub_inversa(transposed_lower_matrix, partial_solution)
    return solution
end

function check_if_singular(main_diagonal_element::Float64, tol::Float64)
    if abs(main_diagonal_element) <= tol
        error("Matriz singular")
    end
end

# util functions
function check_if_triangular(A::Array{Float64,2}, tol::Float64, kind::String)
    n = size(A, 1)
    
    if kind == "inferior"
        indices = [(i,j) for i in 1:n for j in i+1:n if abs(A[i,j]) > tol]
        if !isempty(indices)
            error("Matriz não triangular inferior, indíces não nulos: $indices")
        end
    elseif kind == "superior"
        indices = [(i,j) for i in 1:n for j in 1:i-1 if abs(A[i,j]) > tol]
        if !isempty(indices)
            error("Matriz não triangular superior, indíces não nulos: $indices")
        end
    else
        error("Tipo de matriz inválido")
    end
end

function sub_direta(coefficient_matrix::Array{Float64,2}, constant_vector::Array{Float64,1})::Array{Float64,1}
    tol = 1e-12
    n = length(constant_vector)
    solution_vector::Array{Float64,1} = zeros(n)
    
    # -----------------------------------------------------------------
    #Teste se a matriz é triangular inferior

    # Digite seu código aqui
    # Checar se acima da diagonal principal tem elementos não nulos
    check_if_triangular(coefficient_matrix, tol, "inferior")

    # -----------------------------------------------------------------

    for i in 1:n
        check_if_singular(coefficient_matrix[i,i], tol)
        sum_of_previous_elements = sum(coefficient_matrix[i,j] * solution_vector[j] for j in 1:i-1; init=0.0)
        solution_vector[i] = ( constant_vector[i] - sum_of_previous_elements ) / coefficient_matrix[i,i]
    end

    return solution_vector
end

function sub_inversa(coefficient_matrix::Array{Float64,2}, constant_vector::Array{Float64,1})::Array{Float64,1}
    tol = 1e-12
    n = length(constant_vector)
    
    # -----------------------------------------------------------------
    #Teste se a matriz é triangular superior
    check_if_triangular(coefficient_matrix, tol, "superior")

    solution_vector::Array{Float64,1} = zeros(n)
    for i in n:-1:1
        check_if_singular(coefficient_matrix[i,i], tol)
        sum_of_previous_elements = sum(coefficient_matrix[i,j] * solution_vector[j] for j in i+1:n; init=0.0)
        solution_vector[i] = ( constant_vector[i] - sum_of_previous_elements ) / coefficient_matrix[i,i]
    end
    return solution_vector
end


function norm(v::Array{Float64,1})
    return sqrt(sum(v.^2))
end

function build_identity_matrix(n::Int64)::Array{Float64,2}
    identity_matrix::Array{Float64,2} = zeros(n,n)
    for i in 1:n
        identity_matrix[i,i] = 1
    end
    return identity_matrix
end

function plot_data_and_polynomial(
    t::Array{Float64,1},
    b::Array{Float64,1},
    x::Array{Float64,1},
    method::String
)
    error_vector = b - vandermonde_matrix*x
    mean_error = mean(error_vector)

    scatter(t, b, label="Data", title="Data and Polynomial - $method", xlabel="t", ylabel="b(t)", ms=2, ma=0.5)
    plot!(t, vandermonde_matrix*x, label="Polynomial")
    plot!(t, error_vector, label="Error", linestyle=:dash)
    annotate!(0.5, maximum(error_vector), text("Mean Error ~ $(round(mean_error, digits=12))", 8, :left))

    savefig("data_and_polynomial_$(method).png")
end

# =====================================================================
# Exercício 7
# =====================================================================
# Dados do problema

t = [i/49 for i in 0:49]

b = [cos(4*t[i]) for i in 1:50]

# we need to fit the data to a polynomial of degree 11

vandermonde_matrix = [t[i]^j for i in 1:50, j in 0:11]

# now solving the system using the least squares method

# using the mgs function

Q, R = mgs(vandermonde_matrix, 1e-10)

x = R\Q'*b

plot_data_and_polynomial(t, b, x, "mgs")

# now we need to solve the system using the householder method
Q, R = householder(vandermonde_matrix)
new_rhs_vector = get_householder_rhs(Q, b)
x = solve_householder(Q, R, new_rhs_vector)

plot_data_and_polynomial(t, b, x, "householder")

# now we need to solve the system using the cholesky method
coefficient_matrix = vandermonde_matrix'*vandermonde_matrix
x = cholesky_decomposition_solver(coefficient_matrix, vandermonde_matrix'*b)

plot_data_and_polynomial(t, b, x, "cholesky")

# now we need to solve the system using the julia qr function
qr_result = qr(vandermonde_matrix)
x = qr_result \ b

plot_data_and_polynomial(t, b, x, "julia qr")

# now we need to solve the system using the julia division operator
x = vandermonde_matrix \ b

plot_data_and_polynomial(t, b, x, "julia division")

# now we need to solve the system using the julia svd function
svd_result = svd(vandermonde_matrix)
x = svd_result.V * ((svd_result.U' * b) ./ svd_result.S)

plot_data_and_polynomial(t, b, x, "julia svd")


# =====================================================================
# Comente os resultados encontrados

# Análise dos Resultados - Solução Analítica: [0.35, 0.5, 0.35]
# =====================================================================
