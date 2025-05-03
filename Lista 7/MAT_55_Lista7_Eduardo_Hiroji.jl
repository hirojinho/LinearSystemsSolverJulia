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

function jacobi(a::Array{Float64,2}, b::Array{Float64,1}, x::Array{Float64,1}, tol::Float64, kmax::Int64)

##################
#Digite seu código aqui


##################

    return x
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

function gauss_seidel(a::Array{Float64,2}, b::Array{Float64,1}, x::Array{Float64,1}, tol::Float64, kmax::Int64)

##################
#Digite seu código aqui


##################

    return x
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
    solution_vector = zeros(n)

    for row_index in 1:n
        subvector_index = (row_index - 1) % (m - 1) + 1
        # first subvector
        if row_index < (m - 1)
            x_i = subvector_index * h
            last_subvector_index_term = (subvector_index == m - 1 ? h : 0)
            # round to 9 decimal places, as the stopping criterion is 10^-8
            solution_vector[row_index] = round((x_i - 1) * sin(x_i + last_subvector_index_term), digits=9)
        # last subvector
        elseif row_index > (n - (m - 1))
            x_i = subvector_index * h
            last_subvector_index_term = (subvector_index == m - 1 ? h : 0)
            solution_vector[row_index] = round(x_i*(2 - x_i) + last_subvector_index_term, digits=9)
        # middle subvectors
        elseif row_index % (m - 1) == 0
            y_j = div(row_index, m - 1) * h
            solution_vector[row_index] = round(y_j, digits=9)
        end
    end

    return solution_vector
end

function main()
    m = 4
    block_matrix::Matrix{Float64} = build_block_matrix(m)
    println(block_matrix)
    h = 0.1
    solution_vector::Vector{Float64} = build_solution_vector(m, h)
    println(solution_vector)
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
#Digite aqui os seus comentários, e discuta sobre o desempenho dos métodos.










