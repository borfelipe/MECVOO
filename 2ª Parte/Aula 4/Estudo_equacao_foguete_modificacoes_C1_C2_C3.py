# -*- coding: utf-8 -*-
"""
Created on Mon May 22 16:24:16 2023

@author: taura
"""

## Aplicação prática de cálculo de impulso de velocidade total
# Modificações das configurações de foguetes apresentados na referência
# Pedro L. K. da Cás, Carlos A. G. Veras, Olexiy Shynkarenko, and Rodrigo Leonardi. A
# brazilian space launch system for the small satellite market. Aerospace, 6(11), 2019
#       1° estágio| 2º estágio | 3° estágio
# C-1-n:  4xS50     | 1xS50      | 1xS44
# C-2-n:  4xS50     | 1xS50      | 1xRD843
# C-3-n:  4xS50     | 1xS50      | 4xRD843
#
import numpy as np
#
## Constantes
g=9.81; # m/s**2
## Dados
# Massa de carga útil - Supor um CubeSat 8U (8*1,33+2,36 kg)
mL=13; # kg
# Massa estutural
ms_c1=np.array([4650*4/3, 1367, 166.5])   # kg - C-1-n
ms_c2=np.array([4650*4/3, 1367, 161.9])   # kg - C-2-n
ms_c3=np.array([4650*4/3, 1367, 228.7])   # kg - C-3-n
# Massa de propelente
mp_c1=np.array([33157*4/3, 11058, 813])   # kg - C-1-n
mp_c2=np.array([33157*4/3, 11058, 609])   # kg - C-2-n
mp_c3=np.array([33157*4/3, 11058, 811])   # kg - C-3-n
# Impulso específico ideal
Isp_c1=np.array([251, 271, 270]) # s - C-1-n
Isp_c2=np.array([251, 271, 315]) # s - C-2-n
Isp_c3=np.array([251, 271, 315]) # s - C-3-n
## Cálculos
# Razões estruturais
sigma_c1=ms_c1/(ms_c1+mp_c1)
sigma_c2=ms_c2/(ms_c2+mp_c2)
sigma_c3=ms_c3/(ms_c3+mp_c3)
# Massa total no inicio da queima de cada estagio
m0_c1=np.empty(3);m0_c2=np.empty(3);m0_c3=np.empty(3)
# C-1-n
m0_c1[0]=np.sum(ms_c1)+np.sum(mp_c1)+mL
m0_c1[1]=m0_c1[0]-ms_c1[0]-mp_c1[0]
m0_c1[2]=m0_c1[1]-ms_c1[1]-mp_c1[1]
# C-2-n
m0_c2[0]=sum(ms_c2)+sum(mp_c2)+mL
m0_c2[1]=m0_c2[0]-ms_c2[0]-mp_c2[0]
m0_c2[2]=m0_c2[1]-ms_c2[1]-mp_c2[1]
# C-3-n
m0_c3[0]=sum(ms_c3)+sum(mp_c3)+mL
m0_c3[1]=m0_c3[0]-ms_c3[0]-mp_c3[0]
m0_c3[2]=m0_c3[1]-ms_c3[1]-mp_c3[1]
#
# Razões de carga útil
lamb_c1=np.empty(3);lamb_c2=np.empty(3);lamb_c3=np.empty(3)
# C-1-n
lamb_c1[0]=m0_c1[1]/m0_c1[0];
lamb_c1[1]=m0_c1[2]/m0_c1[1];
lamb_c1[2]=mL/m0_c1[2];
# C-2-n
lamb_c2[0]=m0_c2[1]/m0_c2[0];
lamb_c2[1]=m0_c2[2]/m0_c2[1];
lamb_c2[2]=mL/m0_c2[2];
# C-3-n
lamb_c3[0]=m0_c3[1]/m0_c3[0];
lamb_c3[1]=m0_c3[2]/m0_c3[1];
lamb_c3[2]=mL/m0_c3[2];
# Razão de carga útil total
# C-1-n
lambL_c1=np.prod(lamb_c1);
# C-2-n
lambL_c2=np.prod(lamb_c2);
# C-3-n
lambL_c3=np.prod(lamb_c3);

# Impulso de velocidade
# C-1-n
Dv_c1=-np.sum(g*Isp_c1*np.log(sigma_c1+(1-sigma_c1)*lamb_c1))
# C-2-n
Dv_c2=-np.sum(g*Isp_c2*np.log(sigma_c2+(1-sigma_c2)*lamb_c2))
# C-3-n
Dv_c3=-np.sum(g*Isp_c3*np.log(sigma_c3+(1-sigma_c3)*lamb_c3))
## Resultados
print('*********************************************');
print('C-1-n');
print('*********************************************');
print('Massas iniciais antes da queima de cada estagio - kg');
print(m0_c1);
print('Razoes estruturais');
print(sigma_c1);
print('Razoes de carga util');
print(lamb_c1);
print('Razao de carga útil total');
print(lambL_c1);
print('Impulso de velocidade total - km/s');
print(Dv_c1/1e3);
print('*********************************************');
print('C-2-n');
print('*********************************************');
print('Massas iniciais antes da queima de cada estagio - kg');
print(m0_c2);
print('Razoes estruturais');
print(sigma_c2);
print('Razoes de carga util');
print(lamb_c2);
print('Razao de carga útil total');
print(lambL_c2);
print('Impulso de velocidade total - km/s');
print(Dv_c2/1e3);
print('*********************************************');
print('C-3-n');
print('*********************************************');
print('Massas iniciais antes da queima de cada estagio - kg');
print(m0_c3);
print('Razoes estruturais');
print(sigma_c3);
print('Razoes de carga util');
print(lamb_c3);
print('Razao de carga útil total');
print(lambL_c3);
print('Impulso de velocidade total - km/s');
print(Dv_c3/1e3);