using Printf

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
function jacobi(
    coefficients::Array{Float64,2},
    rhs::Array{Float64,1},
    initial_guess_original::Array{Float64,1},
    tolerance::Float64,
    max_iterations::Int64)::Tuple{Int64, Array{Float64,1}}
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
            println("Jacobi method converged in $iteration iterations for m = $(Int64(sqrt(num_rows) + 1))")
            return iteration, next_guess
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

function gauss_seidel(
    coefficients::Array{Float64,2},
    rhs::Array{Float64,1},
    initial_guess_original::Array{Float64,1},
    tolerance::Float64,
    max_iterations::Int64)::Tuple{Int64, Array{Float64,1}}
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
            println("Gauss-Seidel method converged in $iteration iterations for m = $(Int64(sqrt(num_rows) + 1))")
            return iteration, next_guess
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

function sor(
    coefficients::Array{Float64,2},
    rhs::Array{Float64,1},
    initial_guess_original::Array{Float64,1},
    omega::Float64,
    tolerance::Float64,
    max_iterations::Int64)::Tuple{Int64, Array{Float64,1}}
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
            next_guess[row_index] = (1 - omega)*initial_guess[row_index] +
                omega*(rhs[row_index] - accumulated_sum) / coefficients[row_index, row_index]
        end
        if stop_criterion(initial_guess, next_guess, tolerance)
            println("SOR method converged in $iteration iterations for m = $(Int64(sqrt(num_rows) + 1)) and omega = $omega")
            return iteration, next_guess
        elseif iteration == max_iterations
            throw(ArgumentError("Maximum number of iterations reached."))
        end
        initial_guess = copy(next_guess)
    end
##################
end

# =====================================================================
#                        Método dos Gradientes Conjugados
# =====================================================================
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

function conjugate_gradients(
    coefficients::Array{Float64,2},
    rhs::Array{Float64,1},
    initial_guess_original::Array{Float64,1},
    tolerance::Float64,
    max_iterations::Int64)::Tuple{Int64, Array{Float64,1}}
##################
#Digite seu código aqui
    residue = rhs - coefficients * initial_guess_original
    conjugate_direction = copy(residue)
    guess = copy(initial_guess_original)
    for iteration in 1:max_iterations
        last_guess = copy(guess)
        alpha = (residue')*conjugate_direction/(conjugate_direction'*coefficients*conjugate_direction)
        guess += alpha * conjugate_direction
        residue = rhs - coefficients * guess
        beta = - (residue')*coefficients*conjugate_direction/(conjugate_direction'*coefficients*conjugate_direction)
        conjugate_direction = residue + beta * conjugate_direction
        if stop_criterion(last_guess, guess, tolerance)
            println("Conjugate Gradients method converged in $iteration iterations for m = $(Int64(sqrt(size(coefficients, 1)) + 1))")
            return iteration, guess
        elseif iteration == max_iterations
            throw(ArgumentError("Maximum number of iterations reached."))
        end
    end
##################
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
    h = [0.1, 0.05, 0.025]
    tolerance = 1e-8
    max_iterations = Int64(1e8)
    omega = collect(1:0.1:1.9)

    results = []

    for h_i in h
        m = Int64(1/h_i)
        n = (m - 1)^2

        block_matrix::Matrix{Float64} = build_block_matrix(m)
        rhs::Vector{Float64} = build_solution_vector(m, h_i)
        initial_guess = zeros(n)

        jacobi_iters, _ = jacobi(block_matrix, rhs, initial_guess, tolerance, max_iterations)
        push!(results, (method="Jacobi", m=m, omega=NaN, iterations=jacobi_iters))
        gauss_seidel_iters, _ = gauss_seidel(block_matrix, rhs, initial_guess, tolerance, max_iterations)
        push!(results, (method="Gauss-Seidel", m=m, omega=NaN, iterations=gauss_seidel_iters))

        for omega_i in omega
            sor_iters, _ = sor(block_matrix, rhs, initial_guess, omega_i, tolerance, max_iterations)
            push!(results, (method="SOR", m=m, omega=omega_i, iterations=sor_iters))
        end

        conjugate_gradients_iters, _ = conjugate_gradients(block_matrix, rhs, initial_guess, tolerance, max_iterations)
        push!(results, (method="Conjugate Gradients", m=m, omega=NaN, iterations=conjugate_gradients_iters))

        println("\n")
    end


    # =====================================================================
    # 		         Resultados obtidos
    # =====================================================================
    #Exiba os resultados obtidos em uma tabela. Você pode usar o comando @printf para imprimir os dados. No terminal julia, digite ? e na sequência @printf para obter ajuda sobre o comando.
    @printf("%-25s %-8s %-8s %-12s\n", "Method", "m", "Omega", "Iterations")
    @printf("%s\n", "-"^45)
    for r in results
        if isnan(r.omega)
            @printf("%-25s %-8d %-8s %-12d\n", r.method, r.m, "-", r.iterations)
        else
            @printf("%-25s %-8d %-8.1f %-12d\n", r.method, r.m, r.omega, r.iterations)
        end
    end
    @printf("%s\n", "-"^45)
end

main()

# =====================================================================
# 		           Comentários
# =====================================================================
# Digite aqui os seus comentários, e discuta sobre o desempenho dos métodos.
# Os métodos iterativos Jacobi, Gauss-Seidel e SOR foram aplicados para resolver sistemas lineares provenientes da discretização de um problema,
# variando o parâmetro  m (relacionado ao tamanho da malha) e, no caso do SOR, o parâmetro de relaxação ω.
# =====================================================================
# 1. Método de Jacobi:
# 
# O método de Jacobi apresentou o maior número de iterações para convergência em todos os casos analisados. Isso ocorre porque ele utiliza apenas
# os valores da iteração anterior para atualizar todas as incógnitas, o que geralmente resulta em uma convergência mais lenta.
# À medida que o valor de m aumenta (ou seja, o sistema fica maior), o número de iterações cresce significativamente,
# evidenciando a limitação do método para sistemas de grande porte.
# =====================================================================
# 2. Método de Gauss-Seidel:
# 
# O método de Gauss-Seidel foi mais eficiente que o Jacobi, necessitando de menos iterações para atingir a convergência.
# Isso se deve ao fato de que, a cada iteração, o método já utiliza os valores mais atualizados das incógnitas, acelerando o processo de convergência.
# Ainda assim, o número de iterações cresce com o aumento de m, mas de forma menos acentuada em comparação ao Jacobi.
# =====================================================================
# 3. Método SOR (Successive Over-Relaxation):
# 
# O método SOR mostrou-se o mais eficiente entre os três, especialmente para valores ótimos do parâmetro de relaxação ω.
# Observa-se que, para cada valor de m, existe um valor de ω que minimiza o número de iterações necessárias para a convergência.
# Para valores de ω próximos de 1 (caso particular que recupera o Gauss-Seidel), o desempenho é igual ao do Gauss-Seidel.
# À medida que ω aumenta, o número de iterações diminui até atingir um mínimo, e depois volta a aumentar se ω for elevado demais,
# podendo até prejudicar a convergência. Para os casos analisados, valores de ω entre 1.6 e 1.9 proporcionaram uma redução significativa
# no número de iterações, especialmente para sistemas maiores.
# =====================================================================
# Resumo geral:
#
# O Jacobi é o mais simples, mas o mais lento.
# O Gauss-Seidel já apresenta uma melhora considerável.
# O SOR, com escolha adequada de ω, é o mais eficiente, podendo reduzir drasticamente o número de iterações necessárias para convergência,
# principalmente em sistemas de maior dimensão. Esses resultados ilustram a importância de escolher o método iterativo adequado e, no caso do SOR,
# ajustar corretamente o parâmetro de relaxação para obter o melhor desempenho possível.
# =====================================================================
