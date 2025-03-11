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

function sub_inversa(coefficient_matrix::Array{Float64,2}, constant_vector::Array{Float64,1})
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
    else
        println("Opção inválida.")
        return nothing
    end
end

# Example usage
A = [2 1 -1; -3 -1 2; -2 1 2]
b = [8, -11, -3]
solve_linear_system(A, b)

