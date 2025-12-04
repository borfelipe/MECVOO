from numpy import *

m_e = [20000, 8000, 3200] # Massas de cada estágio (estrutura + propelente) em kg
sigma = [0.10, 0.11, 0.09] # Razões estruturais
I_sp = [270, 350, 455] # Impulsos específicos em segundos
m_L = 1000 # Massa da carga útil em kg

# (i)
m_0 = sum(m_e, initial = m_L)

# (ii)
m_01 = sum(m_e[0:], initial = m_L)
m_02 = sum(m_e[1:], initial = m_L)
m_03 = sum(m_e[2:], initial = m_L)

# (iv)
m_p = [(1-sigma[0])*m_e[0], (1-sigma[1])*m_e[1], (1-sigma[2])*m_e[2]]

# (iii)
m_s = [m_e[0] - m_p[0], m_e[1] - m_p[1], m_e[2] - m_p[2]]

# (v)
lambda_k = [m_L/m_01, m_L/m_02, m_L/m_03]

# (vi)
lambda_T = prod(lambda_k)

# (vii)
v_ek = [impulso * 9.80665 for impulso in I_sp]

v_ek_arr = array(v_ek)
sigma_arr = array(sigma)
lambda_k_arr = array(lambda_k)

Delta_v = -sum(v_ek_arr * log(sigma_arr + ((1-sigma_arr)*sigma_arr*lambda_k_arr)))


# (i) Massa Total
print(f"\n(i) Massa total do foguete (m₀): {m_0:.2f} kg")

# (ii) Massas Parciais
print("\n(ii) Massas parciais antes da queima de cada estágio:")
print(f"     - Estágio 1 (m₀₁): {m_01:.2f} kg")
print(f"     - Estágio 2 (m₀₂): {m_02:.2f} kg")
print(f"     - Estágio 3 (m₀₃): {m_03:.2f} kg")

# (iii) Massa de Estrutura
print("\n(iii) Massa de estrutura de cada estágio:")
print(f"     - Estágio 1 (mₛ₁): {m_s[0]:.2f} kg")
print(f"     - Estágio 2 (mₛ₂): {m_s[1]:.2f} kg")
print(f"     - Estágio 3 (mₛ₃): {m_s[2]:.2f} kg")

# (iv) Massa de Propelente
print("\n(iv) Massa de propelente de cada estágio:")
print(f"     - Estágio 1 (mₚ₁): {m_p[0]:.2f} kg")
print(f"     - Estágio 2 (mₚ₂): {m_p[1]:.2f} kg")
print(f"     - Estágio 3 (mₚ₃): {m_p[2]:.2f} kg")

# (v) Razões de Carga Útil
print("\n(v) Razões de carga útil de cada estágio:")
print(f"     - Estágio 1 (λ₁): {lambda_k[0]:.4f}")
print(f"     - Estágio 2 (λ₂): {lambda_k[1]:.4f}")
print(f"     - Estágio 3 (λ₃): {lambda_k[2]:.4f}")

# (vi) Razão de Carga Útil Total
print(f"\n(vi) Razão de carga útil total (λT = produto de λₖ): {lambda_T:.8f}")

# (vii) Impulso de Velocidade Total (Delta-v)
print(f"\n(vii) Impulso de velocidade total (Δv): {Delta_v:.2f} m/s")
print(f"      (ou {Delta_v/1000:.2f} km/s)")
print("---------------------------------------------------------")
