# -*- coding: utf-8 -*-
"""
Created on Thu May 18 16:25:22 2023

@author: ANDRÉ LUIS DA SILVA
"""
#  
# Exemplo 8.2 do livro TEWARI, A. Atmospheric and Space Flight Dynamics: Modelling
# and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
#
import numpy as np
from scipy.optimize import fsolve
#
## Variaveis globais para passagem de parametros
#
## Constante
g=9.81;  # m/s^2
## Dados
Isp=np.array([290, 290, 455])  # s
sig=np.array([0.07, 0.07, 0.07])   # Razoes estruturais
alf=np.array([1, 1.2, 0.65])   # Razao de carga util de cada estagio normalizada pela do primeiro
mL=1000    # kg - Massa de carga util
Dv=13000   # m/s - Impulso de velocidade total do foguete
#
## Funcao objetivo para encontrar a razao de carga util do primeiro estagio
def obj_eq_fog(lam1):
    ## Calculo da funcao cujo zero deve ser encontrado
    # Soma dos incrementos de delta v em cada estagio
    sdv=np.sum(-g*Isp*np.log(sig+(1-sig)*lam1*alf))
    # A diferenca entre o delta v desejado e o delta v provido pelos 3 estagios
    # deve ser nula
    y=Dv-sdv
    return y
#
## Determinacao da razao de carga util do primeiro estagio
# Chute inicial
lam1=0.5
# Usa a funcao fzero
lam1= fsolve(obj_eq_fog,lam1)
## Determinacao da massa no inicio da queima de cada estagio
# Razao de carga util de cada estagio
lam=lam1*alf
# Razao de carga util total
lamT=lam[0]*lam[1]*lam[2]
# Massa no inicio da queima do primeiro estagio
m01=mL/lamT
# Massa no inicio da queima do segundo e terceiro estagios
m02=m01*lam[0]
m03=m02*lam[1]
## Determinacao da massa de propelente de cada estagio
# Massa estrutural e total de propelente em cada estagio
msp=np.empty(3)
msp[0]=m01-m02;msp[1]=m02-m03;msp[2]=m03-mL
# Massa estutural em cada estagio
ms=msp*sig
# Massa de propelente em cada estagio
mp=msp-ms
# Massa de propelente total
mpT=np.sum(mp)
# Percentual de massa de propelente com respeito a massa total do foguete
pmp=100*mpT/m01
# Massa estrutural total
msT=np.sum(ms)
# Percentual de massa estrutural com respeito a massa total do foguete
pms=100*msT/m01
## Saida de resultados
print('Razao de carga util de cada estagio: ');
print('lamba_1 = ',lam[0]);print('lamba_2 = ',lam[1]);print('lamba_3 = ',lam[2])
print('Razao de carga util total: lamb_T = ',lamT)
print('Massa no inicio da queima de cada estagio: ')
print('m01 =',m01,'kg');print('m02 = ',m02,'kg');print('m03 =',m03,'kg')
print('Massa estrutural e de propelente em cada estagio: ')
print('Primeiro estagio: m_s1+m_p1 = ',msp[0],'kg')
print('Segundo estagio: m_s2+m_p2 = ',msp[1],'kg')
print('Terceiro estagio: m_s3+m_p3 = ',msp[2],'kg')
print('Massa de propelente em cada estagio: ');
print('Primeiro estagio: m_p1 = ',mp[0],'kg')
print('Segundo estagio: m_p2 = ',mp[1],'kg')
print('Terceiro estagio: m_p3 = ',mp[2],'kg')
print('Massa total de propelente: m_pT = ',mpT,'kg')
print('Massa de propelente com respeito a massa total do foguete: ',pmp,'%')
print('Massa estrutural em cada estagio: ')
print('Primeiro estagio: m_s1 = ',ms[0],'kg')
print('Segundo estagio: m_s2 = ',ms[1],'kg')
print('Terceiro estagio: m_s3 = ',ms[2],'kg')
print('Massa estrutural total: m_sT = ',msT,'kg')
print('Massa estrutural com respeito a massa total do foguete: ',pms,'%')