import numpy as np
import matplotlib.pyplot as plt

# Parámetros del sistema
mass = 30000  # kg
nominalPower = 100e3  # W
n_tr = 0.9    # eficiencia
Af = 9.5   # área frontal m^2
rho = 1.225   # kg/m³
CD0 = 0.7     # coeficiente base
minDist = 60  # m
energia_combustible = 34.2e6  # J/L (34.2 MJ/L)
gravity = 9.81

# Parámetros del viento oscilatorio
amplitude = 10.5 # WindSpeed Amplitude m/s
omega_hz  = .01# Times per second
omega = omega_hz*2*np.pi 

print(omega,'rad/s')

t_0, t_f, dt = 0, 3600, 0.01
tiempo = np.arange(t_0, t_f, dt)

N = len(tiempo)
P_truck1 = np.zeros(N)
P_truck2 = np.zeros(N)
P_truck3 = np.zeros(N)
fuel_truck1 = np.zeros(N)
fuel_truck2 = np.zeros(N)
fuel_truck3 = np.zeros(N)
x1 = np.zeros(N)
v1 = np.zeros(N)
x2 = np.zeros(N)
v2 = np.zeros(N)
x3 = np.zeros(N)
v3 = np.zeros(N)

v1[0] = 0.1
x1[0] = 0
v2[0] = 0.1
x2[0] = -minDist
v3[0] = 0.1
x3[0] = -2 * minDist

phi1 = 0  # Fase fija para Camión 1
phi2 = np.pi / 2  # Fase fija para Camión 2 (90 grados)
phi3 = np.pi  # Fase fija para Camión 3 (180 grados)

def cdfunction(distance, area, cdNominal):
    if distance <= 0.5 * np.sqrt(area):
        cd_ad = cdNominal * (1 - 0.3)  
    elif distance >= 1.0 * np.sqrt(area):
        cd_ad = cdNominal * (1 - 0.2)  
    else:
        cd_ad = cdNominal  
    return cd_ad

def control_power(distance, minDist, maxPower):
    if distance < minDist:
        return 0
    else:
        return maxPower
    
def compute_acceleration(P, v, F_drag, mass, n_tr):
    return ((P * n_tr / v) - F_drag-(0.005*gravity*mass)) / mass

for i in range(N - 1):
    # Wind speed (only one value for now)
    v_wind1 = amplitude * np.sin(omega * tiempo[i])

    # --- Truck 1 ---
    P1 = nominalPower
    F_drag1 = 0.5 * rho * (v1[i] + v_wind1)**2 * CD0 * Af
    a1 = lambda v: compute_acceleration(P1, v, 0.5 * rho * (v + v_wind1)**2 * CD0 * Af, mass, n_tr)

    k1_v1 = dt * a1(v1[i])
    k1_x1 = dt * v1[i]

    k2_v1 = dt * a1(v1[i] + 0.5 * k1_v1)
    k2_x1 = dt * (v1[i] + 0.5 * k1_v1)

    k3_v1 = dt * a1(v1[i] + 0.5 * k2_v1)
    k3_x1 = dt * (v1[i] + 0.5 * k2_v1)

    k4_v1 = dt * a1(v1[i] + k3_v1)
    k4_x1 = dt * (v1[i] + k3_v1)

    v1[i+1] = v1[i] + (k1_v1 + 2*k2_v1 + 2*k3_v1 + k4_v1) / 6
    x1[i+1] = x1[i] + (k1_x1 + 2*k2_x1 + 2*k3_x1 + k4_x1) / 6
    fuel_truck1[i] = P1 / (n_tr * energia_combustible )
    P_truck1[i] = P1

    # --- Truck 2 ---
    delta_x12 = abs(x2[i] - x1[i])
    P2 = control_power(delta_x12, minDist, nominalPower)
    CD2 = cdfunction(delta_x12, Af, CD0)
    a2 = lambda v: compute_acceleration(P2, v, 0.5 * rho * (v + v_wind1)**2 * CD2 * Af, mass, n_tr)

    k1_v2 = dt * a2(v2[i])
    k1_x2 = dt * v2[i]

    k2_v2 = dt * a2(v2[i] + 0.5 * k1_v2)
    k2_x2 = dt * (v2[i] + 0.5 * k1_v2)

    k3_v2 = dt * a2(v2[i] + 0.5 * k2_v2)
    k3_x2 = dt * (v2[i] + 0.5 * k2_v2)

    k4_v2 = dt * a2(v2[i] + k3_v2)
    k4_x2 = dt * (v2[i] + k3_v2)

    v2[i+1] = v2[i] + (k1_v2 + 2*k2_v2 + 2*k3_v2 + k4_v2) / 6
    x2[i+1] = x2[i] + (k1_x2 + 2*k2_x2 + 2*k3_x2 + k4_x2) / 6
    fuel_truck2[i] = P2 / (n_tr * energia_combustible )
    P_truck2[i] = P2

    # --- Truck 3 ---
    delta_x23 = abs(x3[i] - x2[i])
    P3 = control_power(delta_x23, minDist, nominalPower)
    CD3 = cdfunction(delta_x23, Af, CD0)
    a3 = lambda v: compute_acceleration(P3, v, 0.5 * rho * (v + v_wind1)**2 * CD3 * Af, mass, n_tr)

    k1_v3 = dt * a3(v3[i])
    k1_x3 = dt * v3[i]

    k2_v3 = dt * a3(v3[i] + 0.5 * k1_v3)
    k2_x3 = dt * (v3[i] + 0.5 * k1_v3)

    k3_v3 = dt * a3(v3[i] + 0.5 * k2_v3)
    k3_x3 = dt * (v3[i] + 0.5 * k2_v3)

    k4_v3 = dt * a3(v3[i] + k3_v3)
    k4_x3 = dt * (v3[i] + k3_v3)

    v3[i+1] = v3[i] + (k1_v3 + 2*k2_v3 + 2*k3_v3 + k4_v3) / 6
    x3[i+1] = x3[i] + (k1_x3 + 2*k2_x3 + 2*k3_x3 + k4_x3) / 6
    fuel_truck3[i] = P3 / (n_tr * energia_combustible )
    P_truck3[i] = P3

# Distancias entre camiones
dist12 = x1 - x2
dist23 = x2 - x3

# # Calcular consumo
# promedio_consumo_truck1 = np.mean(fuel_truck1)
# promedio_consumo_truck2 = np.mean(fuel_truck2)
# promedio_consumo_truck3 = np.mean(fuel_truck3)

# consumo_total_truck1 = np.sum(fuel_truck1) * dt
# consumo_total_truck2 = np.sum(fuel_truck2) * dt
# consumo_total_truck3 = np.sum(fuel_truck3) * dt

# print(f"Promedio de consumo de combustible Camion 1: {promedio_consumo_truck1:.8f} L/s")
# print(f"Promedio de consumo de combustible Camion 2: {promedio_consumo_truck2:.8f} L/s")
# print(f"Promedio de consumo de combustible Camion 3: {promedio_consumo_truck3:.8f} L/s")
# print(f"Consumo total de combustible Camion 1: {consumo_total_truck1:.2f} L")
# print(f"Consumo total de combustible Camion 2: {consumo_total_truck2:.2f} L")
# print(f"Consumo total de combustible Camion 3: {consumo_total_truck3:.2f} L")



# Cálculo del consumo total
consumo_total_truck1 = np.sum(fuel_truck1) * dt
consumo_total_truck2 = np.sum(fuel_truck2) * dt
consumo_total_truck3 = np.sum(fuel_truck3) * dt

print(f"Consumo total de combustible Camion 1: {consumo_total_truck1:.4f} L")
print(f"Consumo total de combustible Camion 2: {consumo_total_truck2:.4f} L")
print(f"Consumo total de combustible Camion 3: {consumo_total_truck3:.4f} L")
# # Crear figura con subplots
# fig, axs = plt.subplots(3, 2, figsize=(12, 10))
# fig.suptitle(f'Viento en Contra Sinusoidal con Vel. {amplitude} m/s (Veces por segundo: {omega_hz})')

# # Subplot 1: Position evolution
# axs[0, 0].plot(tiempo, x1/1000, label="Camión 1")
# axs[0, 0].plot(tiempo, x2/1000, ':', label="Camión 2")
# axs[0, 0].plot(tiempo, x3/1000, '--', label="Camión 3")
# axs[0, 0].set_xlabel("Tiempo (s)")
# axs[0, 0].set_ylabel("Posición (km)")
# axs[0, 0].set_title("Evolución de la posición (Viento Oscilatorio)")
# axs[0, 0].legend()

# # Subplot 2: Velocity evolution
# axs[0, 1].plot(tiempo, v1, label="Velocidad Camión 1")
# axs[0, 1].plot(tiempo, v2, ':', label="Velocidad Camión 2")
# axs[0, 1].plot(tiempo, v3, '--', label="Velocidad Camión 3")
# axs[0, 1].set_xlabel("Tiempo (s)")
# axs[0, 1].set_ylabel("Velocidad (m/s)")
# axs[0, 1].set_title("Evolución de la velocidad (Viento Oscilatorio)")
# axs[0, 1].legend()

# # Subplot 3: Power comparison
# axs[1, 0].plot(tiempo, P_truck1, label="Potencia Camión 1")
# axs[1, 0].plot(tiempo, P_truck2, ':', label="Potencia Camión 2")
# axs[1, 0].plot(tiempo, P_truck3, '--', label="Potencia Camión 3")
# axs[1, 0].set_xlabel("Tiempo (s)")
# axs[1, 0].set_ylabel("Potencia (W)")
# axs[1, 0].set_title("Comparación de la potencia (Viento Oscilatorio)")
# axs[1, 0].grid(True)
# axs[1, 0].legend()

# # Subplot 4: Fuel consumption comparison (instantaneous)
# axs[1, 1].plot(tiempo, fuel_truck1, label="Consumo Camión 1 (L/s)")
# axs[1, 1].plot(tiempo, fuel_truck2, ':', label="Consumo Camión 2 (L/s)")
# axs[1, 1].plot(tiempo, fuel_truck3, '--', label="Consumo Camión 3 (L/s)")
# axs[1, 1].set_xlabel("Tiempo (s)")
# axs[1, 1].set_ylabel("Consumo de combustible (L/s)")
# axs[1, 1].set_title("Comparación de consumo (Viento Oscilatorio)")
# axs[1, 1].grid(True)
# axs[1, 1].legend()

# # Subplot 5: Distance between trucks
# axs[2, 0].plot(tiempo, dist12, label="Distancia Camión 1 - Camión 2")
# axs[2, 0].plot(tiempo, dist23, ':', label="Distancia Camión 2 - Camión 3")
# axs[2, 0].set_xlabel("Tiempo (s)")
# axs[2, 0].set_ylabel("Distancia (m)")
# axs[2, 0].set_title("Distancias entre camiones (Viento Oscilatorio)")
# axs[2, 0].grid(True)
# axs[2, 0].legend()

# # Subplot 6: Total fuel consumption (bar plot)
# trucks = ['Camión 1', 'Camión 2', 'Camión 3']
# totales = [consumo_total_truck1, consumo_total_truck2, consumo_total_truck3]
# axs[2, 1].bar(trucks, totales, color=['blue', 'orange', 'green'])
# axs[2, 1].set_ylabel("Combustible total (L)")
# axs[2, 1].set_title("Consumo total de combustible")
# axs[2, 1].grid(axis='y', linestyle='--', alpha=0.6)

# # Ajustar el diseño
# plt.tight_layout()

plt.figure(figsize=(10, 4))
plt.plot(tiempo, x1/1000, label="Camión 1")
plt.plot(tiempo, x2/1000, ':', label="Camión 2")
plt.plot(tiempo, x3/1000, '--', label="Camión 3")
plt.xlabel("Tiempo (s)")
plt.ylabel("Posición (km)")
plt.title(f"Evolución de la posición (Viento Oscilatorio {omega_hz} Veces por segundo)")
plt.legend()
plt.grid(True)
plt.tight_layout()


plt.figure(figsize=(10, 4))
plt.plot(tiempo, v1, label="Velocidad Camión 1")
plt.plot(tiempo, v2, ':', label="Velocidad Camión 2")
plt.plot(tiempo, v3, '--', label="Velocidad Camión 3")
plt.xlabel("Tiempo (s)")
plt.ylabel("Velocidad (m/s)")
plt.title(f"Evolución de la velocidad (Viento Oscilatorio {omega_hz} Veces por segundo)")
plt.legend()
plt.grid(True)
plt.tight_layout()

# plt.figure(figsize=(10, 4))
# plt.plot(tiempo, P_truck1, label="Potencia Camión 1")
# plt.plot(tiempo, P_truck2, ':', label="Potencia Camión 2")
# plt.plot(tiempo, P_truck3, '--', label="Potencia Camión 3")
# plt.xlabel("Tiempo (s)")
# plt.ylabel("Potencia (W)")
# plt.title("Comparación de la potencia (Viento Oscilatorio)")
# plt.legend()
# plt.grid(True)
# plt.tight_layout()

plt.figure(figsize=(10, 4))
plt.plot(tiempo, fuel_truck1, label="Consumo Camión 1 (L/s)")
plt.plot(tiempo, fuel_truck2, ':', label="Consumo Camión 2 (L/s)")
plt.plot(tiempo, fuel_truck3, '--', label="Consumo Camión 3 (L/s)")
plt.xlabel("Tiempo (s)")
plt.ylabel("Consumo de combustible (L/s)")
plt.title(f"Comparación de consumo (Viento Oscilatorio {omega_hz} Veces por segundo)")
plt.legend()
plt.grid(True)
plt.tight_layout()



plt.figure(2)
plt.plot(tiempo,dist12,'--r',label='Dist12')
plt.plot(tiempo,dist23,'--g',label='Dist23')
plt.grid()
plt.title('Distance Differences per Time')
plt.xlabel('Time /s')
plt.ylabel('Distance /m')
plt.legend()


plt.figure(figsize=(10, 8))
trucks = ['Camión 1', 'Camión 2', 'Camión 3']
totales = [consumo_total_truck1, consumo_total_truck2, consumo_total_truck3]
colors = ['blue', 'orange', 'green']

bars = plt.bar(trucks, totales, color=colors)

# Añadir el valor numérico encima de cada barra
for bar in bars:
    height = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width() / 2,  # posición x
        height,                            # posición y
        f'{height:.2f}',                   # texto
        ha='center', va='bottom', fontsize=10
    )

plt.ylabel("Combustible total (L)")
plt.title(f"Consumo total de combustible (Viento Oscilatorio {omega_hz} veces por segundo)")
plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.tight_layout()


# Subplot 3: Power comparison
figpower, axspower = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

# Plot for Truck 1
axspower[0].plot(tiempo, P_truck1, '-b', label="Potencia Camión 1")
axspower[0].set_ylabel("Potencia [W]")
axspower[0].legend()
axspower[0].grid(True)

# Plot for Truck 2
axspower[1].plot(tiempo, P_truck2, '-g', label="Potencia Camión 2")
axspower[1].set_ylabel("Potencia [W]")
axspower[1].legend()
axspower[1].grid(True)

# Plot for Truck 3
axspower[2].plot(tiempo, P_truck3, '-r', label="Potencia Camión 3")
axspower[2].set_xlabel("Tiempo [s]")
axspower[2].set_ylabel("Potencia [W]")
axspower[2].legend()
axspower[2].grid(True)

# Add an overall title and adjust layout
figpower.suptitle("Truck Power Comparative Respect to Time (Oscillating Wind {omega_hz} Times per Second)")
plt.tight_layout(rect=[0, 0, 1, 0.96])


plt.show()


