from numpy import *

G = 6.67430e-11  # Constante gravitacional em m^3 kg^-1 s^-2
M = 5.972e24     # Massa da Terra em kg
R = 6371e3       # Raio médio da Terra em metros
h = 400e3        # Altitude da órbita em metros
margem = 1.5e3   # Margem de velocidade em m/s

raio_orbita = R + h
velocidade_orbital = sqrt(G * M / raio_orbita)
impulso_total = velocidade_orbital + margem

print(f"Altitude da órbita: {h/1000} km")
print(f"Velocidade orbital calculada: {velocidade_orbital/1000:.2f} km/s")
print(f"Impulso de velocidade total necessário (com margem): {impulso_total/1000:.2f} km/s")