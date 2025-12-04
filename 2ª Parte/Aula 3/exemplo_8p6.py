from scipy.optimize import NonlinearConstraint, minimize, fsolve
import numpy as np
# Inicialização de variáveis globais
beta=np.zeros(3)
sig=np.zeros(3)
Dv=0
ve1=0
## Funcao para otimizar o foguete pelo metodo semi-analitico do livro
# BFM - Brute Force Method
def otimizacao_3_estagios(beta,sig,Dv,ve1):
    ## 
    ## Definicao do conjunto de busca
    # Numero de pontos da discretizacao em cada direcao
    N=50
    # Valores maximo e minimo de alfa2 e alfa 3 - os parametros de otimizacao
    alfm=0.1
    alfM=0.5
    
    # Vetores contendo os valores de alfa2 e alfa3 do conjunto de busca
    alf2=np.linspace(alfm,alfM,N)
    alf3=np.linspace(alfm,alfM,N)
    ## Inicializacoes
    # Vetor para armazenar os valores calculados da razao de carga util total -
    # o parametro otimizado
    lamT=np.zeros((N,N))
    # Vetor para armazenar cada valor de razao de carga util do primeiro
    # estagio obtido no processo de busca do maximo
    lamb1=np.zeros((N,N))
    ## Iteracoes
    for i in range(N):
      for j in range(N):
          # Monta o vetor alfa - Razoes de carga util normalizadas - passa como
          # parâmetro para a funcao que encontra lambda 1 resolvendo a
          # equacao de foguete
          alf=np.array([1, alf2[i], alf3[j]])
          # Calcula lambda 1 que satisfaz a equacao de foguete - restricao do
          # metodo de otimizacao
          dados = (alf, beta, sig, Dv, ve1)
          lamb1[i,j] = fsolve(resolve_lamb1,0.9,args=dados)
          # Calcula a razao de carga util total - parametro otimizado
          lamT[i,j]=lamb1[i,j]**3*alf2[i]*alf3[j]
          # Verifica se a solucao eh factivel
          if (np.imag(lamb1[i,j])!=0)or(lamb1[i,j]>=1):
              lamT[i,j]=-1 # Exclui do cálculo do máximo
    # Valor maximo de da razao de carga util total
    lamT_max=np.max(lamT)
    indiceLinear=np.argmax(lamT)
    ijMax=np.unravel_index(indiceLinear,lamT.shape) # Obtém linha e coluna
    imax=ijMax[0];jmax=ijMax[1]
    # Razoes de carga util normalizadas do segundo e terceiro estagios que maximizam a razao
    # de carga util total
    alf23_max=np.array([alf2[imax],alf3[jmax]])
    # Razao de carga util do primeiro estagio associada ao valor maximo
    lamb1_max=lamb1[imax,jmax]
    #
    return lamT,alf2,alf3,lamb1,lamT_max,alf23_max,lamb1_max
#
## 
# Funcao objetivo para encontrar lambda 1 (razao de carga util do primeiro estagio)
# que satisfaz  a equacao de foguete para um dado valor de alfa2 (razao de
# carga util normalizada do segundo estagio)
def resolve_lamb1(lamb1,*dados):
    ## Passagem de parametros por variaveis globais
    alf, beta, sig, Dv, ve1=dados
    if lamb1<0: # Força para não dar resultado negativo
        lamb1=1
    ## O objetivo eh satisfazer a equacao de foguete
    y=-Dv/ve1-np.sum(beta*np.log(sig+(1-sig)*lamb1*alf))
    #
    return y
#
## Funcao objetivo da funcao minimize
def objNEstagios(lam):
    # Entrada
    # lam: razoes de carga util dos N estagios do veiculo
    # Saida
    # lamT: razao de carga util total com o sinal negativo (a funcao fmincon
    # faz minimizacao, para maximizar, eh necessario inverter o sinal do indice
    # de desempenho
    ## A razao de carga util total eh simplesmente o produto das razoes de carga util de cada estagio
    lamT=-np.prod(lam)
    #
    return lamT
#
## Funcao de restricoes nao lineares
def eqFoguete(lam):
    # lam: razoes de carga util dos N estagios do veiculo
    # Saida
    # c: restricoes de desigualdade (nao ha neste problema)
    # ceq: restricoes de igualdade. Eh a equacao de foguete,
    # para que o impulso de velocidade desejado seja satisfeito
    ## beta, sig, Dv e ve1 são variaveis globais
    ## Restricao de igualdade nao linear: a equacao de foguete
    y = -Dv-np.sum(ve1*beta*np.log(sig+(1-sig)*lam))
    #
    return y
#

# Exemplo 8.6 do livro
# TEWARI, A. Atmospheric and Space Flight Dynamics: 
# Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
# Otimizacao da distribuicao de massa entre os estagios para um foguete de
# dois estagios
## Entrada de dados, com base no exemplo 8.6
# Impulso de velocidade requerido para orbita baixa
Dv=9.5; # km/s
Dv=Dv*1000; # m/s
# Razoes estruturais - fixas para o projeto
sig=np.array([0.07, 0.05, 0.05])
# Velocidade de exaustao do primeiro estagio
ve1=200*9.81;   # Propelente solido
# Velocidades de exaustao normalizadas do segundo e terceiro estagios
beta2=1.5; # UDMH/NO4
beta3=1.75;    # querosene/LO2
# Carga util que deve ser coloca em orbita baixa. Soh eh usada ao final do
# processo para determinar o tamanho final do foguete, mas nao impacta na
# distribuicao de massa entre os estagios
mL=5000;    # kg
## Resolve o problema de otimizacao
# Monta o vetor de velocidades de exaustao normalizadas
beta=np.array([1, beta2, beta3])
# Chama a funcao de otimizacao
print("aqui")
lamT,alf2,alf3,lamb1,lamT_max,alf23_max,lamb1_max=otimizacao_3_estagios(beta,sig,Dv,ve1)
## Pos processamento
# Massas iniciais do foguete otimo
m01=mL/lamT_max;
m02=m01*lamb1_max;
m03=m02*lamb1_max*alf23_max[0];
# Mostra os resultados na tela
print('******************************');
print('Otimização com a função de força bruta');
print('******************************');
print('Maxima razao de carga util: lambda_T_max = %f' % (lamT_max))
print('Razoes de carga util normalizadas do segundo e terceiro estagios que maximizam lambda_T:');
print('(alfa_2 = %f, alfa_3 = %f)' % (alf23_max[0],alf23_max[1]))
print('Razoes de carga util do foguete otimo:');
print('lambda_1 = %f, lambda_2 = %f, lambda_3 = %f' % (lamb1_max, alf23_max[0]*lamb1_max, alf23_max[1]*lamb1_max))
print('Massa inicial do foguete antes da queima do primeiro estagio: m_0_1 = %f kg' % (m01))
print('Massa total do primeiro estagio: m_e_1 = %f kg' % (m01-m02))
print('Massa total do segundo estagio: m_e_2 = %f kg' % (m02-m03))
print('Massa total do terceiro estagio: m_e_3 = %f kg' % (m03-mL))
## Verificacao da solucao
# Calculo do impulso total do foguete
alf=np.array([1, alf23_max[0], alf23_max[1]])
lamb1=lamb1_max
Dv=-np.sum(ve1*beta*np.log(sig+(1-sig)*lamb1*alf))
lamT=lamb1**3*np.prod(alf)
print('******************************');
print('Verificacao do Resultado');
print('******************************');
print('Impulso total de velocidade do foguete otimo: Dv = %f km/s' % (Dv/1000))
print('Razao de carga util total do foguete otimo: lambda_T = %f' % (lamT))
#
## Otimizacao usando a funcao minimize
# Chute inicial
lam0=np.array([1/3, 1/3, 1/3])
# Limites inferior e superior das razoes de carga util
minmax1=(0.01, 0.999)
minmax2=(0.01, 0.999)
minmax3=(0.01, 0.999)
# Restrições não lineares
noncon=NonlinearConstraint(eqFoguete, 0, 0)
# Chama a funcao de otimizacao
res= minimize(objNEstagios,lam0,bounds=(minmax1,minmax2,minmax3),constraints=noncon)
lam=res.x
lamT=res.fun
# Troca o sinal de lamT (foi negativado para maximizar)
lamT=-lamT;
# Mostra os resultados
# Massas iniciais do foguete otimo
m01=mL/lamT;
m02=m01*lam[0]
m03=m02*lam[1]
# Impulso total do foguete
Dv=-np.sum(ve1*beta*np.log(sig+(1-sig)*lam))
# Mostra os resultados na tela
print('******************************');
print('Otimização com a função minimize');
print('******************************');
print('Maxima razao de carga util: lambda_T_max = %f' % (lamT))
print('Razoes de carga util do foguete otimo:');
print('lambda_1 = %f, lambda_2 = %f, lambda_3 = %f' % (lam[0],lam[1],lam[2]))
print('Massa inicial do foguete antes da queima do primeiro estagio: m_0_1 = %f kg' % (m01))
print('Massa total do primeiro estagio: m_e_1 = %f kg' % (m01-m02))
print('Massa total do segundo estagio: m_e_2 = %f kg' % (m02-m03))
print('Massa total do terceiro estagio: m_e_3 = %f kg' % (m03-mL))
print('******************************');
print('Verificacao do Resultado');
print('******************************');
print('Impulso total de velocidade do foguete otimo: Dv = %f km/s' % (Dv/1000))
print('Razao de carga util total do foguete otimo: lambda_T = %f' % (lamT))