from numpy import *

R_e = 6378137  # Raio equatorial terrestre [m]
omega_e = 7.292107e-5  # Velocidade de rotação da Terra [rad/s]
mu_t = 3.986e14  # Parâmetro gravitacional padrão da Terra [m³/s²]
v = 9500  # Velocidade relativa [m/s]
phi = deg2rad(-23)  # Ângulo de elevação [rad]
A = deg2rad(-50)  # Azimute [rad]
h = 300000  # Altitude [m]
delta = deg2rad(36)  # Latitude [rad]
l = deg2rad(-126)  # Longitude [rad]

# Forças
D = 0  # Arrasto desprezível
f_y = 0  # Força lateral desprezível
L = 0  # Sustentação desprezível
f_T = 0  # Propulsor não atuante
epsilon = 0  # Ângulo de defleção vertical
mu_angle = 0  # Ângulo de defleção lateral

m = 1  # [kg] Massa (não fornecida, mas cancelará nas equações) - valor unitário pois não afeta o resultado

r = R_e + h # Distância radial

# Aceleração da gravidade (Terra esférica)
g_c = mu_t / r**2  # Componente centrípeta
g_delta = 0  # Componente latitudinal (Terra esférica)

# (a) Derivadas da magnitude, azimute e elevação da velocidade relativa
print("\n" + "=" * 70)
print("(a) Derivadas da velocidade relativa")
print("=" * 70)

# Equação 23: Derivada da magnitude da velocidade
v_dot = (1/m) * (-D + f_T*cos(epsilon)*cos(mu_angle) - m*g_c*sin(phi) +
         m*g_delta*cos(A)*cos(phi)) - \
        r*omega_e**2*cos(delta)*(cos(A)*sin(delta)*cos(phi) -
                                     cos(delta)*sin(phi))

# Equação 24: Derivada do azimute
A_dot = (1/(m*v*cos(phi))) * (-f_y + f_T*sin(mu_angle) - m*g_delta*sin(A)) - \
        (1/(v*cos(phi))) * 2*v*omega_e*(cos(A)*cos(delta)*sin(phi) -
                                            sin(delta)*cos(phi)) + \
        (1/(v*cos(phi))) * (r*omega_e**2*sin(A)*sin(delta)*cos(delta) +
                                (v**2/r)*sin(A)*tan(delta)*cos(phi)**2)

# Equação 25: Derivada do ângulo de elevação
phi_dot = (1/(m*v)) * (L + f_T*cos(mu_angle)*sin(epsilon) - m*g_c*cos(phi) -
           m*g_delta*cos(A)*sin(phi)) + \
          (1/v) * (2*v*omega_e*sin(A)*cos(delta) + (v**2/r)*cos(phi) +
                   r*omega_e**2*cos(delta)*(cos(A)*sin(delta)*sin(phi) +
                                                cos(delta)*cos(phi)))

print(f"\ndv/dt = {v_dot:.6f} m/s²")
print(f"dA/dt = {rad2deg(A_dot):.8f} °/s = {A_dot:.10f} rad/s")
print(f"dφ/dt = {rad2deg(phi_dot):.8f} °/s = {phi_dot:.10f} rad/s")

# (b) Vetor aceleração inercial no SRV
print("\n" + "=" * 70)
print("(b) Vetor aceleração inercial no SRV")
print("=" * 70)

# Equações 19
a_xv = v_dot + r*omega_e**2*cos(delta)*(cos(A)*sin(delta)*cos(phi) -
                                             cos(delta)*sin(phi))

a_yv = v*cos(phi)*A_dot + 2*v*omega_e*(cos(A)*cos(delta)*sin(phi) -
                                            sin(delta)*cos(phi)) - \
       r*omega_e**2*sin(A)*sin(delta)*cos(delta) - \
       (v**2/r)*sin(A)*tan(delta)*cos(phi)**2

a_zv = -v*phi_dot + r*omega_e**2*cos(delta)*(cos(A)*sin(delta)*sin(phi) +
                                                  cos(delta)*cos(phi)) + \
       2*v*omega_e*sin(A)*cos(delta) + (v**2/r)*cos(phi)

a_srv = array([a_xv, a_yv, a_zv])

print(f"\nAceleração inercial no SRV:")
print(f"a_xv = {a_xv:.6f} m/s²")
print(f"a_yv = {a_yv:.6f} m/s²")
print(f"a_zv = {a_zv:.6f} m/s²")
print(f"\nMagnitude: {linalg.norm(a_srv):.6f} m/s²")

# (c) Vetor aceleração inercial no LVLH
print("\n" + "=" * 70)
print("(c) Vetor aceleração inercial no LVLH")
print("=" * 70)

# Matriz de transformação C_lvlh_srv (Equação 1)
C_lvlh_srv = array([
    [sin(phi), cos(phi)*sin(A), cos(phi)*cos(A)],
    [0, cos(A), -sin(A)],
    [-cos(phi), sin(phi)*sin(A), sin(phi)*cos(A)]
])

# Transformação inversa: SRV -> LVLH
C_srv_lvlh = C_lvlh_srv.T

# Aceleração no LVLH
a_lvlh = C_srv_lvlh @ a_srv

print(f"\nAceleração inercial no LVLH:")
print(f"a_x = {a_lvlh[0]:.6f} m/s² (direção radial)")
print(f"a_y = {a_lvlh[1]:.6f} m/s² (direção leste)")
print(f"a_z = {a_lvlh[2]:.6f} m/s² (direção norte)")
print(f"\nMagnitude: {linalg.norm(a_lvlh):.6f} m/s²")