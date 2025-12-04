from numpy import *

v = 9000  # m/s - velocidade relativa
phi = 5  # graus - ângulo de elevação
A = 70  # graus - azimute
h = 200e3  # m - altitude (200 km)
delta = 0  # graus - latitude
l = -35  # graus - longitude
R_e = 6378137  # m - raio equatorial terrestre

phi_rad = radians(phi)
A_rad = radians(A)
delta_rad = radians(delta)
l_rad = radians(l)

# Distância radial
r = R_e + h

# (a) Taxa de variação da distância radial
r_dot = v * sin(phi_rad)
print(f"\n(a) Taxa de variação da distância radial:")
print(f"    ṙ = {r_dot:.4f} m/s")

# (b) Taxa de variação da longitude
l_dot_rad = (v * cos(phi_rad) * sin(A_rad)) / (r * cos(delta_rad))
l_dot_deg = degrees(l_dot_rad)
print(f"\n(b) Taxa de variação da longitude:")
print(f"    l̇ = {l_dot_rad:.6e} rad/s")
print(f"    l̇ = {l_dot_deg:.6f} °/s")

# (c) Taxa de variação da latitude
delta_dot_rad = (v * cos(phi_rad) * cos(A_rad)) / r
delta_dot_deg = degrees(delta_dot_rad)
print(f"\n(c) Taxa de variação da latitude:")
print(f"    δ̇ = {delta_dot_rad:.6e} rad/s")
print(f"    δ̇ = {delta_dot_deg:.6f} °/s")

# (d) Vetor velocidade angular do LVLH com respeito ao PCPF (no LVLH)
omega_lvlh_pcpf = array([
    l_dot_rad * sin(delta_rad),
    -delta_dot_rad,
    l_dot_rad * cos(delta_rad)
])
print(f"\n(d) Vetor velocidade angular ω_LVLH/PCPF (no LVLH):")
print(f"    ω = [{omega_lvlh_pcpf[0]:.6e}, {omega_lvlh_pcpf[1]:.6e}, {omega_lvlh_pcpf[2]:.6e}]ᵀ rad/s")

# (e) Vetor velocidade relativa do foguete (no LVLH)
v_lvlh = array([
    v * sin(phi_rad),
    v * cos(phi_rad) * sin(A_rad),
    v * cos(phi_rad) * cos(A_rad)
])
print(f"\n(e) Vetor velocidade relativa do foguete (no LVLH):")
print(f"    v = [{v_lvlh[0]:.4f}, {v_lvlh[1]:.4f}, {v_lvlh[2]:.4f}]ᵀ m/s")

# (f) Vetor velocidade relativa do foguete (no PCPF)
# Matriz de transformação LVLH para PCPF
C_lvlh_pcpf = array([
    [cos(delta_rad) * cos(l_rad), -sin(l_rad), -cos(l_rad) * sin(delta_rad)],
    [cos(delta_rad) * sin(l_rad),  cos(l_rad), -sin(l_rad) * sin(delta_rad)],
    [sin(delta_rad),                  0,              cos(delta_rad)]
])
v_pcpf = C_lvlh_pcpf @ v_lvlh
print(f"\n(f) Vetor velocidade relativa do foguete (no PCPF):")
print(f"    v_PCPF = [{v_pcpf[0]:.4f}, {v_pcpf[1]:.4f}, {v_pcpf[2]:.4f}]ᵀ m/s")

# (g) Vetor posição com respeito ao PCPF (no PCPF)
r_pcpf = array([
    r * cos(delta_rad) * cos(l_rad),
    r * cos(delta_rad) * sin(l_rad),
    r * sin(delta_rad)
])
print(f"\n(g) Vetor posição com respeito ao PCPF (no PCPF):")
print(f"    r_PCPF = [{r_pcpf[0]:.4f}, {r_pcpf[1]:.4f}, {r_pcpf[2]:.4f}]ᵀ m")
print(f"    ||r|| = {linalg.norm(r_pcpf):.4f} m")