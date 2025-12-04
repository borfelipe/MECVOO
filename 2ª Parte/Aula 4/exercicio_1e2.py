from numpy import *

# Constantes
g = 9.81
mL = 5
mp_original = array([677, 898])
ms_original = array([284, 320])
Isp_original = array([260.6, 261.1])
sigma_3, Isp_3 = 0.21, 315
mp3_range = linspace(115, 312, 100)


def calcular_dv(mp_stages, ms_stages, Isp_stages, mp3):
    """Calcula impulso de velocidade para massa de propelente mp3"""
    mp = append(mp_stages, mp3)
    ms = append(ms_stages, (sigma_3 / (1 - sigma_3)) * mp3)
    Isp = append(Isp_stages, Isp_3)

    sigma = ms / (ms + mp)
    m0 = array([sum(ms) + sum(mp) + mL, 0, 0])
    m0[1] = m0[0] - ms[0] - mp[0]
    m0[2] = m0[1] - ms[1] - mp[1]

    lamb = array([m0[1] / m0[0], m0[2] / m0[1], mL / m0[2]])
    return -sum(g * Isp * log(sigma + (1 - sigma) * lamb))


# Exercício 1: Ordem original
Dv_ex1 = array([calcular_dv(mp_original, ms_original, Isp_original, mp3)
                for mp3 in mp3_range])
idx1 = argmax(Dv_ex1)

print("=" * 60)
print("EXERCÍCIO 1: Ordem Original")
print(f"  ΔV máximo: {Dv_ex1[idx1] / 1e3:.4f} km/s")
print(f"  Massa propelente 3° estágio: {mp3_range[idx1]:.2f} kg")
print(f"  Massa estrutural 3° estágio: {(sigma_3 / (1 - sigma_3)) * mp3_range[idx1]:.2f} kg")

# Exercício 2: Ordem invertida
Dv_ex2 = array([calcular_dv(mp_original[::-1], ms_original[::-1],
                            Isp_original[::-1], mp3) for mp3 in mp3_range])
idx2 = argmax(Dv_ex2)

print("\n" + "=" * 60)
print("EXERCÍCIO 2: Ordem Invertida")
print(f"  ΔV máximo: {Dv_ex2[idx2] / 1e3:.4f} km/s")
print(f"  Massa propelente 3° estágio: {mp3_range[idx2]:.2f} kg")
print(f"  Massa estrutural 3° estágio: {(sigma_3 / (1 - sigma_3)) * mp3_range[idx2]:.2f} kg")

print("\n" + "=" * 60)
print("COMPARAÇÃO")
print(f"  Diferença ΔV: {(Dv_ex2[idx2] - Dv_ex1[idx1]):.2f} m/s")
print(f"  Variação: {((Dv_ex2[idx2] / Dv_ex1[idx1] - 1) * 100):.2f}%")