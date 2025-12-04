from scipy.optimize import NonlinearConstraint, minimize, fsolve
import numpy as np
import matplotlib.pyplot as plt
# Inicialização de variáveis globais
beta=np.zeros(2)
sig=np.zeros(2)
Dv=0
ve1=0
## Funcao para otimizar o foguete pelo metodo semi-analitico do livro
# BFM - Brute Force Method
def otimizacao_2_estagios(beta,sig,Dv,ve1):
    ## 
    ## Definicao do conjunto de busca
    # Numero de pontos da discretizacao
    N=1000;
    # Valores maximo e minimo de alfa2 - o parametro de otimizacao
    alfm=0.001;
    alfM=1;
    # Vetor contendo os valores de alfa2 do conjunto de busca
    alf2=np.linspace(alfm,alfM,N)
    ## Inicializacoes
    # Vetor para armazenar os valores calculados da razao de carga util total -
    # o parametro otimizado
    lamT=np.zeros(N)
    # Vetor para armazenar cada valor de razao de carga util do primeiro
    # estagio obtido no processo de busca do maximo
    lamb1=np.zeros(N)
    ## Iteracoes
    for i in range(N):
      # Monta o vetor alfa - Razoes de carga util normalizadas - passa como
      # parametro para a funcao que encontra lambda 1 resolvendo a
      # equacao de foguete
      alf=np.array([1, alf2[i]])
      # Calcula lambda 1 que satisfaz a equacao de foguete - restricao do
      # metodo de otimizacao
      dados = (alf, beta, sig, Dv, ve1)
      lamb1[i] = fsolve(resolve_lamb1,0.9,args=dados)
      # Calcula a razao de carga util total - parametro otimizado
      lamT[i]=lamb1[i]**2*alf2[i]
    #
    # Valor maximo de da razao de carga util total
    lamT_max=np.max(lamT)
    i_max=np.argmax(lamT)
    # Razao de carga util normalizada do segundo estagio que maximiza a razao
    # de carga util total
    alf2_max=alf2[i_max]
    # Razao de carga util do primeiro estagio associada ao valor maximo
    lamb1_max=lamb1[i_max]
    #
    return lamT,alf2,lamb1,lamT_max,alf2_max,lamb1_max
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
# Exemplo 8.5 do livro
# TEWARI, A. Atmospheric and Space Flight Dynamics: 
# Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
# Otimizacao da distribuicao de massa entre os estagios para um foguete de
# dois estagios
## Entrada de dados, com base no exemplo 8.5
# Impulso de velocidade requerido para orbita baixa
Dv=9.5; # km/s
Dv=Dv*1000; # m/s
# Razoes estruturais - fixas para o projeto
sig=np.array([0.07, 0.05])
# Velocidade de exaustao do primeiro estagio
ve1=200*9.81;   # Propelente solido
# Duas possibilidades de beta2 (velocidade de exaustao normalizada do
# segundo estagio
beta2a=1.5; # UDMH/NO4
beta2b=1.75;    # querosene/LO2
# Carga util que deve ser coloca em orbita baixa. Soh eh usada ao final do
# processo para determinar o tamanho final do foguete, mas nao impacta na
# distribuicao de massa entre os estagios
mL=5000;    # kg
## Resolve para o propelente do segundo estagio UDMH/NO4
# Monta o vetor de velocidades de exaustao normalizadas
beta=np.array([1,beta2a])
# Chama a funcao de otimizacao
lamTa,alf2a,lamb1a,lamT_maxa,alf2_maxa,lamb1_maxa=otimizacao_2_estagios(beta,sig,Dv,ve1)
## Resolve para o propelente do segundo estagio querosene/LO2
beta=np.array([1,beta2b])
# Chama a funcao de otimizacao
lamTb,alf2b,lamb1b,lamT_maxb,alf2_maxb,lamb1_maxb=otimizacao_2_estagios(beta,sig,Dv,ve1)
## Graficos
plt.close('all')
plt.figure(1)
plt.plot(alf2a,lamTa, label = "Propelente UDMH/NO4 - beta_2 = 1,5")
plt.plot(alf2b,lamTb, label = "Propelente querosene/LO2 - beta_2 = 1,75")
plt.plot(alf2_maxa,lamT_maxa,marker='*');plt.plot(alf2_maxb,lamT_maxb,marker='*')
plt.legend();plt.show()
plt.xlabel('alpha_2=lambda_2/lambda_1');plt.ylabel('lambda_T')
#
# Massas iniciais do foguete otimo
m01a=mL/lamT_maxa
m02a=m01a*lamb1_maxa
m01b=mL/lamT_maxb
m02b=m01b*lamb1_maxb
# Verificacao da solucao
# Calculo do impulso total do foguete
alfa=np.array([1, alf2_maxa])
lamb1a=lamb1_maxa
beta=np.array([1, beta2a])
Dva=-np.sum(ve1*beta*np.log(sig+(1-sig)*lamb1a*alfa))
alfb=np.array([1, alf2_maxb])
lamb1b=lamb1_maxb
beta=np.array([1, beta2b])
Dvb=-np.sum(ve1*beta*np.log(sig+(1-sig)*lamb1b*alfb))
# Razoes de carga util de cada estagio 
# Mostra os resultados na tela
print('**************************************************************************');
print('Resultados para o propelente UDMH/NO4 no segundo estagio');
print('**************************************************************************');
print('Maxima razao de carga util total: %f ' % (lamT_maxa))
print('Razao de carga util normalizada do segundo estagio (alfa2) que maximiza lambdaT: %f ' % (alf2_maxa))
print('Razoes de carga util do primeiro e segundo estagio do foguete otimo:');
print([lamb1_maxa,alf2_maxa*lamb1_maxa]);
print('Massa inicial do foguete antes da queima do primeiro estagio: %f kg' % (mL/lamT_maxa))
print('Massa total do primeiro estagio: %f kg' % (m01a-m02a))
print('Massa total do segundo estagio: %f kg' % (m02a-mL))
print('Impulso total de velocidade do foguete otimo: %f km/s' % (Dva/1000))
print('*******************************************************************************');
print('Resultados para o propelente querosene/LO2 no segundo estagio');
print('*******************************************************************************');
print('Maxima razao de carga util: %f ' % (lamT_maxb))
print('Razao de carga util normalizada do segundo estagio (alfa2) que maximiza lambdaT: % f' % (alf2_maxb))
print('Razoes de carga util do primeiro e segundo estagio do foguete otimo:');
print([lamb1_maxb,alf2_maxb*lamb1_maxb]);
print('Massa inicial do foguete antes da queima do primeiro estagio: %f kg' % (mL/lamT_maxb))
print('Massa total do primeiro estagio: %f kg' % (m01b-m02b))
print('Massa total do segundo estagio: %f kg' % (m02b-mL))
print('Impulso total de velocidade do foguete otimo: %f km/s' % (Dvb/1000))
#
## Otimizacao usando a funcao minize
# Chute inicial
lam0=np.array([1/3, 1/3])
# Limites inferior e superior das razoes de carga util
minmax1=(0.01, 0.999)
minmax2=(0.01, 0.999)
#
# Propelente UDMH/NO4: Monta o vetor de velocidades de exaustao normalizadas
beta=np.array([1, beta2a])
# Restrições não lineares
noncon=NonlinearConstraint(eqFoguete, 0, 0)
# Chama a funcao de otimizacao
res= minimize(objNEstagios,lam0,bounds=(minmax1,minmax2),constraints=noncon)
lama=res.x
lamTa=res.fun
# Troca o sinal de lamT (foi negativado para maximizar)
lamTa=-lamTa;
#
# Propelente querosene/LO2: Monta o vetor de velocidades de exaustao normalizadas
beta=np.array([1, beta2b])
# Restrições não lineares
noncon=NonlinearConstraint(eqFoguete, 0, 0)
# Chama a funcao de otimizacao
res= minimize(objNEstagios,lam0,bounds=(minmax1,minmax2),constraints=noncon)
lamb=res.x
lamTb=res.fun
# Troca o sinal de lamT (foi negativado para maximizar)
lamTb=-lamTb;
#
# Mostra os resultados
# Massas iniciais do foguete otimo
m01a=mL/lamTa;
m02a=m01a*lama[0]
m01b=mL/lamTb;
m02b=m01b*lamb[0]
# Impulso total do foguete
beta=np.array([1, beta2a])
Dva=-np.sum(ve1*beta*np.log(sig+(1-sig)*lama))
beta=np.array([1, beta2b])
Dvb=-np.sum(ve1*beta*np.log(sig+(1-sig)*lamb))
print('**********************************************************');
print('Resultados determinados com a funcao minimize');
print('**********************************************************');
print('**************************************************************************');
print('Resultados para o propelente UDMH/NO4 no segundo estagio');
print('**************************************************************************');
print('Maxima razao de carga util: %f ' % (lamTa))
print('Razao de carga util normalizada do segundo estagio (alfa2) que maximiza lambdaT: %f' % (lama[1]/lama[0]))
print('Razoes de carga util do foguete otimo:')
print(lama);
print('Massa inicial do foguete antes da queima do primeiro estagio: %f kg' % (m01a))
print('Massa total do primeiro estagio: %f kg' % (m01a-m02a))
print('Massa total do segundo estagio: %f kg' % (m02a-mL))
print('Impulso total de velocidade do foguete otimo: %f km/s' % (Dva/1000))
print('*******************************************************************************');
print('Resultados para o propelente querosene/LO2 no segundo estagio');
print('*******************************************************************************');
print('Maxima razao de carga util: %f ' % (lamTb))
print('Razao de carga util normalizada do segundo estagio (alfa2) que maximiza lambdaT: %f ' % (lamb[1]/lamb[0]))
print('Razoes de carga util do foguete otimo:');
print(lamb);
print('Massa inicial do foguete antes da queima do primeiro estagio: %f kg' % (m01b))
print('Massa total do primeiro estagio: %f kg' % ((m01b-m02b)))
print('Massa total do segundo estagio: %f kg' % (m02b-mL))
print('Impulso total de velocidade do foguete otimo: %f km/s' % (Dvb/1000))
