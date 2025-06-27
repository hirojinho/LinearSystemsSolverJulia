# Apresentação - Decomposição em Modos Dinâmicos (DMD)

## Sobre a Apresentação

Esta apresentação foi desenvolvida para uma exposição de **10-15 minutos** sobre o trabalho "Decomposição em Modos Dinâmicos: Teoria Computacional e Conexões com Métodos de Álgebra Linear".

## Estrutura da Apresentação

### 1. Introdução e Motivação (3-4 minutos)
- Contextualização do DMD
- Conexões com métodos de MAT-55
- Motivação para o estudo teórico

### 2. Teoria Computacional (4-5 minutos)
- Formulação matricial do problema
- Papel central da decomposição SVD
- Construção do operador DMD reduzido

### 3. Algoritmos DMD (3-4 minutos)
- Algoritmo DMD Clássico
- Algoritmo DMD Exato
- Comparação com PCA

### 4. Fundamentação Teórica (2-3 minutos)
- Conexão com operador de Koopman
- Justificativa para sistemas não-lineares

### 5. Conclusões (1-2 minutos)
- Principais contribuições
- Reflexões finais

## Como Compilar

### Usando o Makefile
```bash
make -f Makefile_apresentacao apresentacao
```

### Usando pdflatex diretamente
```bash
pdflatex apresentacao_dmd.tex
pdflatex apresentacao_dmd.tex  # Segunda compilação
```

### Visualizar (macOS)
```bash
make -f Makefile_apresentacao view
```

## Dicas para a Apresentação

### Tempo por Slide
- **Slide de título**: 30 segundos
- **Slides de conteúdo**: 1-2 minutos cada
- **Slides de algoritmo**: 2-3 minutos
- **Slide de conclusão**: 1-2 minutos

### Pontos Importantes a Enfatizar
1. **Conexão com MAT-55**: Sempre relacionar com os métodos estudados
2. **Aspecto algorítmico**: Foco na implementação, não apenas teoria
3. **Elegância matemática**: Como métodos clássicos se estendem naturalmente
4. **Operador de Koopman**: Fundamentação teórica rigorosa

### Possíveis Perguntas
- **Sobre estabilidade numérica**: Herança das propriedades da SVD
- **Sobre aplicações**: Mencionar que o foco foi teórico/algorítmico
- **Sobre implementação**: Algoritmos utilizam apenas métodos de MAT-55
- **Sobre diferenças DMD/PCA**: Temporal vs espacial

## Recursos Adicionais

### Se precisar de imagens
O slide sobre formulação matricial menciona uma figura `exemplo_temporal.png`. Se necessário, pode ser removido ou substituído por uma explicação verbal.

### Customização do Tema
Para mudar a aparência:
```latex
\usetheme{Madrid}      % Temas: Warsaw, Berlin, Singapore, etc.
\usecolortheme{seahorse}  % Cores: whale, dolphin, rose, etc.
```

### Notas do Apresentador
Para adicionar notas privadas:
```latex
\note{Sua nota aqui}
```

## Checklist Pré-Apresentação

- [ ] Revisar todos os slides
- [ ] Verificar se as equações estão corretas
- [ ] Testar o tempo de apresentação
- [ ] Preparar respostas para possíveis perguntas
- [ ] Verificar se o PDF foi gerado corretamente

## Backup Plan

Caso haja problemas técnicos:
1. Tenha uma versão impressa dos slides principais
2. Prepare uma versão resumida de 7-10 minutos
3. Saiba explicar os conceitos principais sem slides 