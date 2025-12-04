import numpy as np
## Funcao para calculo das derivadas de massa de
# carga util com respeito as massas estruturais
def der_mL_dms(mp,ms,mL,ve):
    # Numero de estagios
    N=mp.size
    # Vetores de massa final e inicial por estagio
    mf=np.zeros(N)
    m0=np.zeros(N)
    for i in range(N):
        m0[i]=mL+np.sum(mp[i:N])+np.sum(ms[i:N])
        mf[i]=m0[i]-mp[i]
    # Derivadas da massa de carga util com respeito as massas estruturais
    dmLds=np.zeros(N)
    for k in range(N):
        dmLds[k]=-np.sum(ve[0:k+1]*(1/m0[0:k+1]-1/mf[0:k+1]))/np.sum(ve*(1/m0-1/mf))
    #
    return dmLds
## Funcao para calculo das derivadas de massa de
# carga util com respeito as massas de propelente
def der_mL_dmp(mp,ms,mL,ve):
    # Numero de estagios
    N=mp.size
    # Vetores de massa final e inicial por estagio
    mf=np.zeros(N)
    m0=np.zeros(N)
    for i in range(N):
        m0[i]=mL+np.sum(mp[i:N])+np.sum(ms[i:N])
        mf[i]=m0[i]-mp[i]
    # Derivadas da massa de carga util com respeito as massas propelente
    dmLdp=np.zeros(N)
    for k in range(N):
        dmLdp[k]=-(ve[k]/m0[k]+np.sum(ve[0:k]*(1/m0[0:k]-1/mf[0:k])))/np.sum(ve*(1/m0-1/mf))
    #
    return dmLdp

#
# Exemplo 8.4 do livro - Calculo das derivadas de carga util com respeito
# as massas estruturais e de propelente.
# TEWARI, A. Atmospheric and Space Flight Dynamics: 
# Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
# Utilizacao dos dados dos exemplos 8.2 e 8.3
## Dados
# Massa no inicio da queima de cada estagio
m0=np.array([1.094266769684211e+05, 5.284772960000001e+04, 19806.1060, 6046.9099])
# Massa de propelente consumida em cada estagio
mp = np.array([55000, 2.8978709948e+04, 1.279605244274179e+04, 4568.81457])
# Carga util - kg
mL = 1.124811598501497e+03;
# Velocidade de exaustao de cada estagio (estagios zero e 1 sao
# equivalentes)
ve = 1.0e+03*np.array([2.363318181818182,   2.844900000000000,   2.844900000000000,   4.463550000000001])
# Determina as massas estruturais de cada estagio
N=m0.size
ms=np.zeros(N)
ms[N-1]=m0[N-1]-mL-mp[N-1]
ms[1:N-1]=m0[1:N-1]-m0[2:N]-mp[1:N-1]
# Massas ao final da queima de cada estagio
mf=m0-mp;
# Mostra as massas
print('Massa da carga util (kg)')
print(mL)
print('Massas no inicio da queima de cada estagio')
print(m0)
print('Massas no final da queima de cada estagio')
print(mf)
print('Massas estruturais de cada estagio')
print(ms)
print('Massas de propelente de cada estagio')
print(mp)
## Calculo das derivadas de massa de carga util com respeito as massas estruturais
dmLds=der_mL_dms(mp,ms,mL,ve)
print('Vetor de derivadas da carga util com respeito as massas estruturais: ');
print(dmLds)
## Calculo das derivadas de massa de carga util com respeito as massas de propelente
dmLdp=der_mL_dmp(mp,ms,mL,ve)
print('Vetor de derivadas da carga util com respeito as massas de propelente: ');
print(dmLdp)
