# =====================================================================
# *********************************************************************
#                          MAT-55 2024 - Lista 10
# *********************************************************************
# =====================================================================
# Dupla: Eduardo Hiroji
#
#
#Para calcular a SVD você pode usar a SVD do pacote LinearAlgebra
#svd(A)
using LinearAlgebra
using CSV, DataFrames, Statistics

#Dados do problema:
#X = [];
problem_matrix = CSV.read("daily-treasury-rates.csv", DataFrame) |> Matrix |> x -> x[:, 2:end]
difference_matrix = diff(problem_matrix, dims=1)

############################
#Matriz na forma de desvio de média:
#A = [];
mean_array = mean(difference_matrix, dims=1)
mean_difference_matrix = difference_matrix .- mean_array

############################
#Matriz de covariância:
#K = A'*A;
#K = A/(m-1);
num_rows = size(mean_difference_matrix)[1]
covariance_matrix = mean_difference_matrix' * mean_difference_matrix / (num_rows - 1)
U, S, V = svd(covariance_matrix)
component_matrix = U * Diagonal(S)

component_contribution = S ./ sum(S)
for (i, contribution) in enumerate(component_contribution)
    println("Componente $i: $(round(contribution * 100, digits=2))%")
end

# =====================================================================
# 		           Comentários
# =====================================================================
# A análise de componentes principais (PCA) aplicada às taxas de juros do tesouro revelou
# resultados bastante significativos sobre a estrutura e comportamento dessas taxas. O primeiro
# componente principal, responsável por 83.37% da variância total dos dados, demonstra uma
# forte correlação entre as taxas de juros em diferentes prazos. Esta dominância do primeiro
# componente sugere a existência de um fator comum que influencia de maneira coordenada todas
# as taxas de juros do tesouro.
#
# Complementando esta análise, o segundo componente principal, que explica 9.72% da variância,
# representa movimentos independentes nas taxas que não são capturados pelo primeiro
# componente. Juntos, estes dois componentes explicam mais de 93% da variância total dos
# dados, indicando que a maior parte da variabilidade das taxas de juros pode ser
# compreendida através destes dois fatores principais.
#
# A alta concentração de variância no primeiro componente (83.37%) é particularmente
# reveladora. Este resultado indica que as taxas de juros do tesouro em diferentes prazos
# tendem a se mover de forma bastante coordenada, respondendo de maneira similar a
# mudanças nas condições econômicas e financeiras. Esta característica sugere que existe
# um fator comum subjacente que influencia todas as taxas.
#
# A possibilidade de reduzir a dimensionalidade dos dados para apenas dois componentes
# sem perda significativa de informação (mantendo mais de 93% da variância) tem
# implicações práticas importantes. Esta redução simplifica significativamente a análise
# das taxas de juros, permitindo uma compreensão mais clara dos principais fatores que
# influenciam seu comportamento. Os movimentos secundários independentes, capturados
# pelo segundo componente, podem representar ajustes específicos em determinados prazos
# ou respostas a fatores particulares do mercado.
#
# Em conclusão, a análise PCA revela uma estrutura de dados altamente correlacionada,
# onde as taxas de juros do tesouro em diferentes prazos são predominantemente
# influenciadas por um fator comum, com movimentos secundários independentes
# complementando esta dinâmica principal. Esta compreensão é valiosa para o
# entendimento da estrutura de taxas de juros e pode auxiliar na tomada de decisões
# relacionadas a investimentos e análise de risco.