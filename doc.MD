# Implementação/Teste da eficiência do método dos Gradientes Conjugados

## Autoria
Nome: Rubens Zandomenighi Laszlo  GRR:20206147

## Objetivo 


### Métodos 
São utilizadas as funções para geração aleatória da matriz, e posteriormente é transformada em diagonal dominante 
Para amenização de erros, também é transformado o Sistema Linear em, ao início dos métodos : 
(A^T) * A * x = (A^T) * b

#### Método com pré condicionador de Jacobi

Utilização do método de pré condicionamento de jacobi para a amenização de problema de estabilidade numérica.
Ax = b → M−1Ax = M−1b
Sendo M-1 = Matriz diagonal tal que M=D.

### Método sem pré condicionamento 


### Calculo de tempo e erros
É feito o cálculo do tempo do resíduo, tempo média iteração método e tempo para aplicação dos pré condicionadores. 
Os métodos para cálculos de erro utilizados são :

calcNormaMaxRel -> (max (|xi - xi-1| / |xi| ),
seja residuo = b - A * x
calcErroNormaEuc -> sqrt((residuo)^2)

## Arquitetura 
Lista de  Arquivos
     
    Diretório Principal
      |
      |----main.c - Programa principal para execução das funções
      |----utils.c - Contém funções úteis para auxílio na execução, como tratamento de entrada
      |----opmatrizes.c - Contém funções efetuar operações com matrizes
      |----sislin.c - Contém funções lidar com o sistema linear
      |----makefile - arquivo para compilação e correlatos
      |----doc - Documentação do projeto

## Compilação 

all: compila e produz um executável chamado cgSolver;
clean: remove todos os arquivos temporários e os arquivos gerados pelo Makefile (*.o, executável, *.bak, etc.).

## Execução 
cgSolver -n <n> -k <k> -p <p> -i <i> -e <ε> -o <arquivo_saida> 

<n>: (n>10) parâmetro obrigatório definindo a dimensão do Sistema Linear.
<k>: (k>1 e k ímpar)  parâmetro obrigatório definindo o número de diagonais da matriz A.
<p>:  parâmetro obrigatório indicando o pré-condicionador a ser utilizado:
p = 0: sem pré-condicionador
p > 0: pré-condicionador de Jacobi
<i>: parâmetro obrigatório definindo o número máximo de iterações a serem executadas.
<ε>: parâmetro opcional definindo o erro aproximado absoluto máximo, considerando a norma max (relativa) em x ( max ( |xi - xi-1| / |xi| ) < ε ), a ser considerado como critério de parada além do parâmetro <i>.
Na ausência desta opção, usar como critério de parada apenas o parâmetro <i>.
<arquivo_saida>: parâmetro obrigatório no qual arquivo_saida é o caminho completo para o arquivo que vai conter a solução.