# =====================================================================
# *********************************************************************
#                    MAT-55 2024 - Lista 01 - Exercício 6 
# *********************************************************************
# =====================================================================

# =====================================================================
#                    Algoritmo de Substituição Direta
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# A     matriz nxn triangular inferior, não singular
# b     vetor n
# ---------------------------------------------------------------------
# Saída:
# b     se A é não singular, b é a solução do sistema linear Ax = b

# TODO: Verificar se a matriz é singular junto com o loop de cálculo em vez de um loop separado
function check_if_singular(main_diagonal_element::Float64, tol::Float64)
    if abs(main_diagonal_element) <= tol
        error("Matriz singular")
    end
end

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

function sub_direta(coefficient_matrix::Array{Float64,2}, constant_vector::Array{Float64,1})
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

# =====================================================================
#                    Algoritmo de Substituição Inversa
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# A     matriz nxn triangular superior, não singular
# b     vetor n
# ---------------------------------------------------------------------
# Saída:
# b     se A é não singular, b é a solução do sistema linear Ax = b

function sub_inversa(coefficient_matrix::Array{Float64,2}, constant_vector::Array{Float64,1})::Array{Float64,1}
    tol = 1e-12
    n = length(b)
    
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

function get_multiplier(
    pivot_row::Array{Float64,1}, 
    row_to_reduce::Array{Float64,1}
)
    tol = 1e-12
    first_non_zero_index = findfirst(x -> abs(x) > tol, pivot_row)
    if first_non_zero_index === nothing
        error("Pivot row is zero")
    end
    multiplier = row_to_reduce[first_non_zero_index] / pivot_row[first_non_zero_index]
    return multiplier
end

function get_multiplier(
    pivot_row::Array{Float64,1}, 
    row_to_reduce::Array{Float64,1},
    pivot_row_index::Int
)
    tol = 1e-12
    if pivot_row[pivot_row_index] === 0.0
        error("Pivot row is zero")
    end
    multiplier = row_to_reduce[pivot_row_index] / pivot_row[pivot_row_index]
    return multiplier
end

function reduce_row(
    pivot_row::Array{Float64,1},
    row_to_reduce::Array{Float64,1}
)
    multiplier = get_multiplier(pivot_row, row_to_reduce)
    return row_to_reduce - multiplier * pivot_row
end

# =====================================================================
#                    Algoritmo de Eliminação Gaussiana
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# A     matriz nxn triangular superior, não singular
# b     vetor n
# ---------------------------------------------------------------------
# Saída:
# b     se A é não singular, b é a solução do sistema linear Ax = b

function elim_gauss(coefficient_matrix::Array{Float64,2}, constant_vector::Array{Float64,1})
    n = length(constant_vector)
    reduced_matrix::Array{Float64,2} = copy(coefficient_matrix)
    reduced_constant_vector::Array{Float64,1} = copy(constant_vector)

    # Create an extended matrix with the constant vector
    extended_matrix = hcat(reduced_matrix, reduced_constant_vector)

    #Digite seu código aqui
    for i in 1:n-1
        pivot_row = copy(extended_matrix[i,:])
        rows_to_reduce = i+1:n
        reduced_rows = [reduce_row(pivot_row, extended_matrix[j,:]) for j in rows_to_reduce]
        extended_matrix[rows_to_reduce,:] = vcat(reduced_rows'...)
    end

    # Extract the reduced coefficient matrix and constant vector
    reduced_matrix = extended_matrix[:,1:n]
    reduced_constant_vector = extended_matrix[:,n+1]

    return sub_inversa(reduced_matrix, reduced_constant_vector)
end

# =====================================================================
#                    Decomposição LU
# =====================================================================
# ---------------------------------------------------------------------
function build_identity_matrix(n::Int)
    return [i == j ? 1.0 : 0.0 for i in 1:n, j in 1:n]
end

function build_multiplier_matrix(multipliers_vector::Array{Float64,1}, system_size::Int)
    vector_size = length(multipliers_vector)
    column_to_add = system_size - vector_size
    identity_matrix = build_identity_matrix(system_size)
    for row_index in 1:vector_size
        identity_matrix[row_index + column_to_add, column_to_add] = multipliers_vector[row_index]
    end
    return identity_matrix
end

function build_multiplier_matrix(coefficient_matrix::Array{Float64,2}, pivot_row_index::Int)
    size_of_matrix = size(coefficient_matrix, 1)
    rows_to_reduce = pivot_row_index+1:size_of_matrix
    multipliers_vector = map(
        x -> - get_multiplier(coefficient_matrix[pivot_row_index,:], coefficient_matrix[x,:], pivot_row_index),
        rows_to_reduce
        )
    multiplier_matrix = build_multiplier_matrix(multipliers_vector, size_of_matrix)
    return multiplier_matrix
end

function build_permutation_matrix(max_row_index::Int, row_index::Int, size_of_matrix::Int)
    permutation_matrix = build_identity_matrix(size_of_matrix)
    permutation_matrix[row_index, :], permutation_matrix[max_row_index, :] = permutation_matrix[max_row_index, :], permutation_matrix[row_index, :]
    return permutation_matrix
end

function lower_upper_decomposition(
    coefficient_matrix::Array{Float64,2},
    pivoting::Bool = false
)
    size_of_matrix = size(coefficient_matrix, 1)

    lower_matrix::Array{Float64,2} = build_identity_matrix(size_of_matrix)
    upper_matrix::Array{Float64,2} = zeros(size_of_matrix, size_of_matrix)
    permutation_matrix::Array{Float64,2} = build_identity_matrix(size_of_matrix)

    for row_index in 1:size_of_matrix - 1
        if pivoting
            max_row_index = reduce(
                (max_row_index, row_index) -> abs(coefficient_matrix[row_index, row_index]) > abs(coefficient_matrix[max_row_index, row_index])
                    ? row_index
                    : max_row_index,
                row_index:size_of_matrix,
                init=row_index
            )
            permutation_matrix = permutation_matrix * build_permutation_matrix(max_row_index, row_index, size_of_matrix)
            coefficient_matrix[max_row_index, :], coefficient_matrix[row_index, :] = coefficient_matrix[row_index, :], coefficient_matrix[max_row_index, :]
        end
        multiplier_matrix = build_multiplier_matrix(coefficient_matrix, row_index)
        coefficient_matrix = multiplier_matrix * coefficient_matrix
        coefficient_matrix[row_index+1:end, row_index] = multiplier_matrix[row_index+1:end, row_index]
    end

    # Get lower triangular part from coefficient_matrix
    for i in 1:size_of_matrix
        for j in 1:size_of_matrix
            if i > j
                lower_matrix[i,j] = coefficient_matrix[i,j]
            end
        end
    end

    # Get upper triangular part from coefficient_matrix
    for i in 1:size_of_matrix
        for j in 1:size_of_matrix
            if i <= j
                upper_matrix[i,j] = coefficient_matrix[i,j]
            end
        end
    end

    # Lower matrix is the inverse of the product of the multiplier matrixes
    return inv(lower_matrix), upper_matrix, permutation_matrix
end

function lower_upper_decomposition_solver(
    coefficient_matrix::Array{Float64,2},
    constant_vector::Array{Float64,1},
    pivoting::Bool = false
)
    lower_matrix, upper_matrix, permutation_matrix = lower_upper_decomposition(coefficient_matrix, pivoting)
    if pivoting
        constant_vector = permutation_matrix * constant_vector
    end
    partial_solution::Array{Float64,1} = sub_direta(lower_matrix, constant_vector)
    return sub_inversa(upper_matrix, partial_solution)
end

# =====================================================================
# =====================================================================
#			PROGRAMA PRINCIPAL
# =====================================================================
# =====================================================================

# Implemente um programa para resolver o sistema linear Ax = b, com 
# opção para o usuário fornecer a matriz A, o vetor b e escolher o método 
# utilizado, dentre as opções:
# a: Algoritmo de Substituição Direta.
# b: Algoritmo de Substituição Inversa. 
# c: Eliminação Gaussiana.

# Dados do sistema
# Digite aqui os dados do sistema linear
# A = [1 0 0; 2 1 0; 3 4 1]
# b = [1, 2, 3]
# A = [1 4 3; 0 1 2; 0 0 1]
# b = [3, 2, 1]
function solve_linear_system(A::Array{Int,2}, b::Array{Int,1})
    println("Escolha o método:")
    println("a: Algoritmo de Substituição Direta.")
    println("b: Algoritmo de Substituição Inversa.") 
    println("c: Eliminação Gaussiana.")
    println("d: Decomposição LU sem pivoteamento.")
    println("e: Decomposição LU com pivoteamento.")

    metodo = readline()

    A_float = map(Float64, A)
    b_float = map(Float64, b)

    if metodo == "a"
        println("Método de Substituição Direta selecionado.")
        x = sub_direta(A_float, b_float)
        println("Solução: $x")
        return x
    elseif metodo == "b"
        println("Método de Substituição Inversa selecionado.")
        x = sub_inversa(A_float, b_float)
        println("Solução: $x")
        return x
    elseif metodo == "c"
        println("Método de Eliminação Gaussiana selecionado.")
        x = elim_gauss(A_float, b_float)
        println("Solução: $x")
        return x
    elseif metodo == "d"
        println("Método de Decomposição LU sem pivoteamento selecionado.")
        x = lower_upper_decomposition_solver(A_float, b_float)
        println("Solução: $x")
        return x
    elseif metodo == "e"
        println("Método de Decomposição LU com pivoteamento selecionado.")
        x = lower_upper_decomposition_solver(A_float, b_float, true)
        println("Solução: $x")
        return x
    else
        println("Opção inválida.")
        return nothing
    end
end

# Example usage
A = [2 1 -1; -3 -1 2; -2 1 2]
b = [8, -11, -3]
solve_linear_system(A, b)

