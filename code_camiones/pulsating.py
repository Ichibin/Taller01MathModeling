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

# Parámetros del viento oscilatorio
wind_speed = 10.5 # WindSpeed Amplitude m/s
pulse_dt = 25 #sec
pulse_width = .05

tk1 = np.arange(0, t_f, pulse_dt)  
tk2 = np.arange(10, t_f + 10, pulse_dt)  
tk3 = np.arange(20, t_f + 20, pulse_dt)  

def gaussian_pulse(t, tk, pulse_width):
    return wind_speed*np.sum(np.exp(-((t - tk) ** 2) / (2 * pulse_width ** 2)) / (pulse_width * np.sqrt(2 * np.pi)))

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
    
    
for i in range(N-1):
    t = tiempo[i]

    # Pulsación de viento en t
    wind_pulse = gaussian_pulse(t, tk1, pulse_width)
    wind_effect = wind_speed * wind_pulse

    # Estados actuales
    x1i, v1i = x1[i], v1[i]
    x2i, v2i = x2[i], v2[i]
    x3i, v3i = x3[i], v3[i]

    def compute_derivatives(x1i, v1i, x2i, v2i, x3i, v3i):
        # Velocidades relativas
        vrel1 = v1i + wind_effect
        vrel2 = v2i + wind_effect
        vrel3 = v3i + wind_effect

        # Distancias
        dx12 = abs(x2i - x1i)
        dx23 = abs(x3i - x2i)

        # Potencias
        P1 = nominalPower
        P2 = control_power(dx12, minDist, nominalPower)
        P3 = control_power(dx23, minDist, nominalPower)

        # Coeficientes de arrastre
        cd2 = cdfunction(dx12, Af, CD0)
        cd3 = cdfunction(dx23, Af, CD0)

        # Fuerzas de arrastre
        drag1 = 0.5 * rho * vrel1**2 * CD0 * Af
        drag2 = 0.5 * rho * vrel2**2 * cd2 * Af
        drag3 = 0.5 * rho * vrel3**2 * cd3 * Af

        # Aceleraciones
        a1 = ((P1 * n_tr / max(v1i, 1e-3)) - drag1) / mass
        a2 = ((P2 * n_tr / max(v2i, 1e-3)) - drag2) / mass
        a3 = ((P3 * n_tr / max(v3i, 1e-3)) - drag3) / mass

        return a1, a2, a3, P1, P2, P3

    # k1
    a1_k1, a2_k1, a3_k1, P1, P2, P3 = compute_derivatives(x1i, v1i, x2i, v2i, x3i, v3i)

    # k2
    x1_k2 = x1i + 0.5 * dt * v1i
    v1_k2 = v1i + 0.5 * dt * a1_k1
    x2_k2 = x2i + 0.5 * dt * v2i
    v2_k2 = v2i + 0.5 * dt * a2_k1
    x3_k2 = x3i + 0.5 * dt * v3i
    v3_k2 = v3i + 0.5 * dt * a3_k1
    a1_k2, a2_k2, a3_k2, _, _, _ = compute_derivatives(x1_k2, v1_k2, x2_k2, v2_k2, x3_k2, v3_k2)

    # k3
    x1_k3 = x1i + 0.5 * dt * v1_k2
    v1_k3 = v1i + 0.5 * dt * a1_k2
    x2_k3 = x2i + 0.5 * dt * v2_k2
    v2_k3 = v2i + 0.5 * dt * a2_k2
    x3_k3 = x3i + 0.5 * dt * v3_k2
    v3_k3 = v3i + 0.5 * dt * a3_k2
    a1_k3, a2_k3, a3_k3, _, _, _ = compute_derivatives(x1_k3, v1_k3, x2_k3, v2_k3, x3_k3, v3_k3)

    # k4
    x1_k4 = x1i + dt * v1_k3
    v1_k4 = v1i + dt * a1_k3
    x2_k4 = x2i + dt * v2_k3
    v2_k4 = v2i + dt * a2_k3
    x3_k4 = x3i + dt * v3_k3
    v3_k4 = v3i + dt * a3_k3
    a1_k4, a2_k4, a3_k4, _, _, _ = compute_derivatives(x1_k4, v1_k4, x2_k4, v2_k4, x3_k4, v3_k4)

    # Actualizar estados
    x1[i+1] = x1i + dt/6.0 * (v1i + 2*v1_k2 + 2*v1_k3 + v1_k4)
    v1[i+1] = v1i + dt/6.0 * (a1_k1 + 2*a1_k2 + 2*a1_k3 + a1_k4)

    x2[i+1] = x2i + dt/6.0 * (v2i + 2*v2_k2 + 2*v2_k3 + v2_k4)
    v2[i+1] = v2i + dt/6.0 * (a2_k1 + 2*a2_k2 + 2*a2_k3 + a2_k4)

    x3[i+1] = x3i + dt/6.0 * (v3i + 2*v3_k2 + 2*v3_k3 + v3_k4)
    v3[i+1] = v3i + dt/6.0 * (a3_k1 + 2*a3_k2 + 2*a3_k3 + a3_k4)

    # Registrar potencias y consumo
    P_truck1[i] = P1
    P_truck2[i] = P2
    P_truck3[i] = P3

    fuel_truck1[i] = P1 / (n_tr * energia_combustible)
    fuel_truck2[i] = P2 / (n_tr * energia_combustible)
    fuel_truck3[i] = P3 / (n_tr * energia_combustible)

# Distancias entre camiones
dist12 = x1 - x2
dist23 = x2 - x3



# Cálculo del consumo total
consumo_total_truck1 = np.sum(fuel_truck1) * dt
consumo_total_truck2 = np.sum(fuel_truck2) * dt
consumo_total_truck3 = np.sum(fuel_truck3) * dt

plt.figure(figsize=(10, 4))
plt.plot(tiempo, x1/1000, label="Camión 1")
plt.plot(tiempo, x2/1000, ':', label="Camión 2")
plt.plot(tiempo, x3/1000, '--', label="Camión 3")
plt.xlabel("Tiempo (s)")
plt.ylabel("Posición (km)")
plt.title(f"Evolución de la posición ")
plt.legend()
plt.grid(True)
plt.tight_layout()


plt.figure(figsize=(10, 4))
plt.plot(tiempo, v1, label="Velocidad Camión 1")
plt.plot(tiempo, v2, ':', label="Velocidad Camión 2")
plt.plot(tiempo, v3, '--', label="Velocidad Camión 3")
plt.xlabel("Tiempo (s)")
plt.ylabel("Velocidad (m/s)")
plt.title(f"Evolución de la velocidad ")
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
plt.title(f"Comparación de consumo")
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
plt.title(f"Consumo total de combustible ")
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