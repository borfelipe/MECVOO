from scipy.optimize import minimize, NonlinearConstraint
from numpy import *


def otimizar_foguete(mL, Isp1, Isp2, sig1, sig2, Dv_km):
    """Otimiza distribuição de massa entre estágios"""
    g0 = 9.80665
    ve1, ve2 = Isp1 * g0, Isp2 * g0
    Dv = Dv_km * 1000
    beta = array([1, ve2 / ve1])
    sig = array([sig1, sig2])

    # Maximizar produto das razões de carga útil
    obj = lambda lam: -prod(lam)

    # Restrição: equação de foguete
    eq = lambda lam: -Dv - sum(ve1 * beta * log(sig + (1 - sig) * lam))

    res = minimize(obj, [0.3, 0.3],
                   bounds=[(0.01, 0.999)] * 2,
                   constraints=NonlinearConstraint(eq, 0, 0))

    lam = res.x
    lamT = -res.fun
    m01 = mL / lamT

    return {
        'λ₁': lam[0], 'λ₂': lam[1], 'λT': lamT,
        'α₂': lam[1] / lam[0],
        'm_inicial': m01,
        'm_estágio_1': m01 * (1 - lam[0]),
        'm_estágio_2': m01 * lam[0] - mL,
        'Δv_verif': -sum(ve1 * beta * log(sig + (1 - sig) * lam)) / 1000
    }


# Parâmetros
mL, Isp1, Isp2, sig1, sig2 = 350, 200, 350, 0.07, 0.05

# Resultados
for caso, Dv in [('a', 9.5), ('b', 10.5)]:
    print(f"\n{'=' * 60}\nCASO ({caso}): Δv = {Dv} km/s\n{'=' * 60}")
    r = otimizar_foguete(mL, Isp1, Isp2, sig1, sig2, Dv)
    print(f"λ₁ = {r['λ₁']:.6f}  |  λ₂ = {r['λ₂']:.6f}  |  λT = {r['λT']:.6f}")
    print(f"α₂ = {r['α₂']:.6f}")
    print(f"\nMassa inicial: {r['m_inicial']:.2f} kg")
    print(f"Massa estágio 1: {r['m_estágio_1']:.2f} kg")
    print(f"Massa estágio 2: {r['m_estágio_2']:.2f} kg")
    print(f"Verificação Δv: {r['Δv_verif']:.4f} km/s")