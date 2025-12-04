from numpy import *

lambda_k = full(3, 1/3)
sigma_k = [0.15, 0.1, 0.07] # Razões estruturais
I_sp = [270, 350, 455] # Impulsos específicos em segundos
m_L = [1e3, 2e3] # Massa da carga útil em kg

for i in range(0,2):
  # (i)
  lambda_T = prod(lambda_k)
  
  
  # (ii)
  v_ek = [impulso * 9.80665 for impulso in I_sp]
  v_ek_arr = array(v_ek)
  sigma_arr = array(sigma_k)
  lambda_k_arr = array(lambda_k)
  
  Delta_v = -sum(v_ek_arr * log(sigma_arr + ((1-sigma_arr)*sigma_arr*lambda_k_arr)))
  
  # (iii)
  m_0 = m_L[i]/lambda_T
  
  # (iv)
  m_01 = m_0
  m_02 = lambda_k[0]*m_01
  m_03 = lambda_k[1]*m_02
  
  # (vi)
  m_e = [m_01-m_02, m_02-m_03, m_03-m_L[i]]
  m_p = [(1-sigma_k[0])*m_e[0], (1-sigma_k[1])*m_e[1], (1-sigma_k[2])*m_e[2]]
  
  # (v)
  m_s = [m_e[0] - m_p[0], m_e[1] - m_p[1], m_e[2] - m_p[2]]
  
  # --- Seção de Impressão dos Resultados ---
  print(f"\n\n=======================================================")
  print(f"=== Carga Útil de: {m_L[i]:.2f} kg ===")
  print(f"=======================================================")

  # (i) Razão de Carga Útil Total
  print(f"\n(i) Razão de carga útil total (λT): {lambda_T:.6f}")

  # (ii) Impulso de Velocidade
  print(f"\n(ii) Impulso de velocidade (Δv): {Delta_v:.2f} m/s (ou {Delta_v/1000:.2f} km/s)")

  # (iii) Massa Total
  print(f"\n(iii) Massa total inicial (m₀): {m_0:.2f} kg")

  # (iv) Massas Iniciais por Estágio
  print("\n(iv) Massas iniciais antes da queima de cada estágio:")
  print(f"     - m₀₁ (início do estágio 1): {m_01:.2f} kg")
  print(f"     - m₀₂ (início do estágio 2): {m_02:.2f} kg")
  print(f"     - m₀₃ (início do estágio 3): {m_03:.2f} kg")
  
  # (v) Massa de Estrutura por Estágio
  print("\n(v) Massa de estrutura (mₛ) de cada estágio:")
  print(f"     - Estágio 1: {m_s[0]:.2f} kg")
  print(f"     - Estágio 2: {m_s[1]:.2f} kg")
  print(f"     - Estágio 3: {m_s[2]:.2f} kg")

  # (vi) Massa Total e de Propelente por Estágio
  print("\n(vi) Massa total (mₑ) e de propelente (mₚ) de cada estágio:")
  print(f"     - Estágio 1: mₑ = {m_e[0]:.2f} kg | mₚ = {m_p[0]:.2f} kg")
  print(f"     - Estágio 2: mₑ = {m_e[1]:.2f} kg | mₚ = {m_p[1]:.2f} kg")
  print(f"     - Estágio 3: mₑ = {m_e[2]:.2f} kg | mₚ = {m_p[2]:.2f} kg")
