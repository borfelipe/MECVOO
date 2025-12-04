import numpy as np
# Exemplo 8.3
# TEWARI, A. Atmospheric and Space Flight Dynamics: 
# Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
## Dados do exemplo 8.3
# Massa de propelente dos boosters
mpb=30000;  # kg
# Impulso específico dos boosters
Ispb=200;   # s
# Razao estrutural dos boosters
sigb=0.05;
# Massa de propelente do primeiro estagio do veiculo nucleo queimada
# durante o uso dos boosters
mp=np.zeros(4)
mp[0]=25000; # kg
## Resultados do exemplo 8.2
# Massa total antes no inicio da queima de cada estagio do veiculo nucleo
m0=np.zeros(4)
m0[1]=77847.7296; # kg
m0[2]=19806.1060;   # kg
m0[3]=6046.9099;  # kg
# Razao estrutural de cada estagio do veiculo nucleo
sigk=0.07;
sig=np.zeros(4)
sig[2]=sigk;sig[3]=sigk;    # Nao muda nos estagios 2 e 3
# Razao de carga util do segundo e terceiro estagios do veiculo nucleo
lam=np.zeros(4)
lam[2]=0.305305336; # Nao muda
lam[3]=0.165373724; # Nao muda
# Impulsos especificos
Isp=np.zeros(4)
Isp[1]=290;   # s
Isp[2]=290;   # s
Isp[3]=455;   # s
# Massa de propelente do terceiro estagio
mp[3]=4693.62617; # 
## Calculos
print('Massa estrutural do primeiro estagio do veiculo nucleo');
ms=np.zeros(4)
ms[1]=sigk*(m0[1]-m0[2]);
print(ms[1]);
print('Massa de propelente do primeiro estagio do veiculo nucleo');
mp[1]=m0[1]-m0[2]-ms[1];
print(mp[1]);
print('Massa estrutural dos boosters');
msb=(sigb/(1-sigb))*mpb;
print(msb);
print('Massa inicial do foguete com os boosters');
m0[0]=m0[1]+mpb+msb;
print(m0[0]);
print('Razao estrutural do estagio zero');
print(msb);print(ms[1]);print(mpb);print(mp[0])
sig[0]=(msb+ms[1])/(msb+ms[1]+mpb+mp[0]);
print(sig[0]);
print('Razao de carga util do estagio zero');
lam[0]=(m0[1]-mp[0])/m0[0];
print(lam[0]);
print('Razao estrutural do primeiro estagio modificado');
sig[1]=ms[1]/(ms[1]+mp[1]-mp[0]);
print(sig[1]);
print('Razao de carga util do primeiro estagio modificado');
lam[1]=m0[2]/(m0[1]-mp[0]);
print(lam[1]);
print('Velocidade de exaustao media do estagio zero');
ve=np.zeros(4)
ve[1]=9.81*Isp[1];
veb=9.81*Ispb;
ve[0]=(mpb*veb+mp[0]*ve[1])/(mpb+mp[0]);
print(ve[0]);
print('Impulso de velocidade total');
ve[2]=Isp[2]*9.81;ve[3]=Isp[3]*9.81;
Dv=-np.sum(ve*np.log(sig+(1-sig)*lam))
print(Dv);
##
print('Analise de desempenho versus eficiencia');
print('Acrescimo de impulso de velocidade (m/s)');
DDv=Dv-13000;
print(DDv);
print('Variacao percentual do desempenho');
des=100*DDv/13000;
print(des);
print('Nova razao de carga util total');
lamT=np.prod(lam);
print(lamT);
print('Variacao da razao de carga util total');
DlamT=lamT-0.0128;
print(DlamT);
print('Variacao percentual da eficiencia');
ef=100*DlamT/0.0128;
print(ef);
## Novas razoes estrutural e de carga util do terceiro estagio para um
# aumento de carga util mantendo o impulso de velocidade constante
k=-(np.sum(ve[0:3]*np.log(sig[0:3]+(1-sig[0:3])*lam[0:3]))+13000)/ve[3];
c1=np.exp(k);
ms[3]=m0[3]-mp[3]-1000;
c2=ms[3]/m0[3];
print('Nova razao estrutural do 3º estagio para Dv constante');
sig[3]=c2/(1-c1+c2);
print(sig[3]);
print('Nova razao de carga util do 3º estagio para Dv constante');
lam[3]=c1-c2;
print(lam[3]);
print('Nova massa de carga util (kg)');
mL=lam[3]*m0[3];
print(mL);
print('Aumento de carga util (kg)');
print(mL-1000);
print('Aumento percentual de carga util');
print(100*(mL-1000)/1000);
print('Nova massa de propelente do terceiro estagio (kg)');
Dm=mL-1000;
mp[3]=mp[3]-(Dm);   # Troca propelente do 3º estagio por carga util
print(mp[3]);
print('Nova razao de carga util total');
lamTn=np.prod(lam);
print(lamTn);
print('Variacao da razao de carga util total');
DlamT=lamTn-0.0128;
print(DlamT);
print('Variacao percentual na eficiencia');
ef=100*DlamT/0.0128;
print(ef);