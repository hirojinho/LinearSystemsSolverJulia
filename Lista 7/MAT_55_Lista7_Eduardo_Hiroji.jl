# =====================================================================
# *********************************************************************
#                    MAT-55 2024 - Lista 7
# *********************************************************************
# =====================================================================
# Dupla:
#
#
# =====================================================================
#                        Método de Jacobi =====================================================================
# Dados de entrada:
# a     matriz, nxn não singular
# b	vetor, n
# x	vetor, n, aproximação inicial
# tol	escalar, tolerância para o critério de parada
# kmax	escalar, número máximo de iterações permitido
#
# Saída:
# x     aproximação para a solução do sistema Ax=b
function jacobi(coefficients::Array{Float64,2}, rhs::Array{Float64,1}, initial_guess_original::Array{Float64,1}, tolerance::Float64, max_iterations::Int64)::Array{Float64,1}
##################
#Digite seu código aqui
    num_rows = size(coefficients, 1)
    next_guess = copy(initial_guess_original)
    initial_guess = copy(initial_guess_original)
    for iteration in 1:max_iterations
        for row_index in 1:num_rows
            if coefficients[row_index, row_index] < tolerance throw(ArgumentError("Matrix is singular or too close to singular.")) end
            accumulated_sum = 0.0
            for column_index in 1:num_rows
                if column_index != row_index
                    accumulated_sum += coefficients[row_index, column_index] * initial_guess[column_index]
                end
            end
            next_guess[row_index] = (rhs[row_index] - accumulated_sum) / coefficients[row_index, row_index]
        end
        if stop_criterion(initial_guess, next_guess, tolerance)
            println("Jacobi method converged in $iteration iterations")
            return next_guess
        elseif iteration == max_iterations
            throw(ArgumentError("Maximum number of iterations reached."))
        end
        initial_guess = copy(next_guess)
    end
##################
end

# =====================================================================
#                        Método de Gauss-Seidel =====================================================================
# Dados de entrada:
# a     matriz, nxn não singular
# b	vetor, n
# x	vetor, n, aproximação inicial
# tol	escalar, tolerância para o critério de parada
# kmax	escalar, número máximo de iterações permitido
#
# Saída:
# x     aproximação para a solução do sistema Ax=b

function gauss_seidel(coefficients::Array{Float64,2}, rhs::Array{Float64,1}, initial_guess_original::Array{Float64,1}, tolerance::Float64, max_iterations::Int64)
##################
#Digite seu código aqui
    num_rows = size(coefficients, 1)
    next_guess = copy(initial_guess_original)
    initial_guess = copy(initial_guess_original)
    for iteration in 1:max_iterations
        for row_index in 1:num_rows
            if coefficients[row_index, row_index] < tolerance throw(ArgumentError("Matrix is singular or too close to singular.")) end
            accumulated_sum = 0.0
            for column_index in 1:num_rows
                if column_index != row_index
                    guess_to_sum = (column_index < row_index ? next_guess : initial_guess)[column_index]
                    accumulated_sum += coefficients[row_index, column_index] * guess_to_sum
                end
            end
            next_guess[row_index] = (rhs[row_index] - accumulated_sum) / coefficients[row_index, row_index]
        end
        if stop_criterion(initial_guess, next_guess, tolerance)
            println("Gauss-Seidel method converged in $iteration iterations")
            return next_guess
        elseif iteration == max_iterations
            throw(ArgumentError("Maximum number of iterations reached."))
        end
        initial_guess = copy(next_guess)
    end
##################
end

# =====================================================================
#                        Método SOR =====================================================================
# Dados de entrada:
# a     matriz, nxn não singular
# b	vetor, n
# x	vetor, n, aproximação inicial
# w	escalar, parâmetro do método SOR
# tol	escalar, tolerância para o critério de parada
# kmax	escalar, número máximo de iterações permitido
#
# Saída:
# x     aproximação para a solução do sistema Ax=b

function sor(a::Array{Float64,2}, b::Array{Float64,1}, x::Array{Float64,1}, w::Float64, tol::Float64, kmax::Int64)

##################
#Digite seu código aqui


##################

    return x
end

# =====================================================================
# 			Helper functions
# =====================================================================
function norm(vector::Array{Float64,1})::Float64
    return sqrt(sum(vector.^2))
end

function stop_criterion(initial_guess::Array{Float64,1}, next_guess::Array{Float64,1}, tolerance::Float64)::Bool
    return norm(initial_guess - next_guess)/norm(next_guess) < tolerance
end

# =====================================================================
# 			Dados do problema
# =====================================================================
# Digite aqui os dados do exercício 5
function build_block_matrix(m::Int64)::Matrix{Float64}
    n = (m - 1)^2
    block_matrix = zeros(n, n)

    for row_index in 1:n
        block_matrix[row_index, row_index] = 4
        if row_index - 1 > 0
            block_matrix[row_index, row_index - 1] = -1
        end
        if row_index + 1 <= n
            block_matrix[row_index, row_index + 1] = -1
        end
        if row_index - (m - 1) > 0
            block_matrix[row_index, row_index - (m - 1)] = -1
        end
        if row_index + (m - 1) <= n
            block_matrix[row_index, row_index + (m - 1)] = -1
        end
    end

    return block_matrix
end

function build_solution_vector(m::Int64, h::Float64)::Vector{Float64}
    n = (m - 1)^2
    rhs = zeros(n)

    for row_index in 1:n
        subvector_index = (row_index - 1) % (m - 1) + 1
        # first subvector
        if row_index < (m - 1)
            x_i = subvector_index * h
            last_subvector_index_term = (subvector_index == m - 1 ? h : 0)
            # round to 9 decimal places, as the stopping criterion is 10^-8
            rhs[row_index] = round((x_i - 1) * sin(x_i + last_subvector_index_term), digits=9)
        # last subvector
        elseif row_index > (n - (m - 1))
            x_i = subvector_index * h
            last_subvector_index_term = (subvector_index == m - 1 ? h : 0)
            rhs[row_index] = round(x_i*(2 - x_i) + last_subvector_index_term, digits=9)
        # middle subvectors
        elseif row_index % (m - 1) == 0
            y_j = div(row_index, m - 1) * h
            rhs[row_index] = round(y_j, digits=9)
        end
    end

    return rhs
end

function main()
    m = 4
    h = 0.1
    tolerance = 1e-8
    max_iterations = 1000
    n = (m - 1)^2
    block_matrix::Matrix{Float64} = build_block_matrix(m)
    rhs::Vector{Float64} = build_solution_vector(m, h)
    initial_guess = zeros(n)
    jacobi_solution = jacobi(block_matrix, rhs, initial_guess, tolerance, max_iterations)
    gauss_seidel_solution = gauss_seidel(block_matrix, rhs, initial_guess, tolerance, max_iterations)
end

main()
# =====================================================================
# 		         Resultados obtidos
# =====================================================================
#Exiba os resultados obtidos em uma tabela. Você pode usar o comando @printf para imprimir os dados. No terminal julia, digite ? e na sequência @printf para obter ajuda sobre o comando.

using Printf

# =====================================================================
# 		           Comentários
# =====================================================================
# Digite aqui os seus comentários, e discuta sobre o desempenho dos métodos.










