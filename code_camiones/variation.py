import numpy as np
import matplotlib.pyplot as plt
import os

def simulate_trucks(mass, nominalPower, CD0):
    # Parámetros del sistema
    n_tr = 0.9    # eficiencia
    Af = 9.5      # área frontal m^2
    rho = 1.225   # kg/m³
    minDist = 60  # m
    energia_combustible = 34.2e6  # J/L (34.2 MJ/L)
    gravity = 9.81

    # Parámetros del viento oscilatorio
    amplitude = 10.5  # WindSpeed Amplitude m/s
    omega_hz = 0.01   # Times per second
    omega = omega_hz * 2 * np.pi

    # Tiempo de simulación
    t_0, t_f, dt = 0, 500, 0.01
    tiempo = np.arange(t_0, t_f, dt)
    N = len(tiempo)

    # Inicialización de arrays
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

    # Condiciones iniciales
    v1[0] = 0.1
    x1[0] = 0
    v2[0] = 0.1
    x2[0] = -minDist
    v3[0] = 0.1
    x3[0] = -2 * minDist

    # Funciones auxiliares
    def cdfunction(distance, area, cdNominal):
        if distance <= 0.5 * np.sqrt(area):
            cd_ad = cdNominal * (1 - 0.3)
        elif distance >= 1.0 * np.sqrt(area):
            cd_ad = cdNominal * (1 - 0.2)
        else:
            cd_ad = cdNominal
        return cd_ad

    def control_power(distance, minDist, maxPower):
        return 0 if distance < minDist else maxPower

    def compute_acceleration(P, v, F_drag, mass, n_tr):
        return ((P * n_tr / v) - F_drag - (0.005 * gravity * mass)) / mass

    # Simulación
    for i in range(N - 1):
        # Velocidad del viento
        v_wind1 = amplitude * np.sin(omega * tiempo[i])

        # --- Camión 1 ---
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
        fuel_truck1[i] = P1 / (n_tr * energia_combustible)
        P_truck1[i] = P1

        # --- Camión 2 ---
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
        fuel_truck2[i] = P2 / (n_tr * energia_combustible)
        P_truck2[i] = P2

        # --- Camión 3 ---
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
        fuel_truck3[i] = P3 / (n_tr * energia_combustible)
        P_truck3[i] = P3

    # Distancias entre camiones
    dist12 = x1 - x2
    dist23 = x2 - x3

    # Cálculo del consumo total
    consumo_total_truck1 = np.sum(fuel_truck1) * dt
    consumo_total_truck2 = np.sum(fuel_truck2) * dt
    consumo_total_truck3 = np.sum(fuel_truck3) * dt

    # Imprimir resultados
    # print(f"Consumo total de combustible Camion 1: {consumo_total_truck1:.4f} L")
    # print(f"Consumo total de combustible Camion 2: {consumo_total_truck2:.4f} L")
    # print(f"Consumo total de combustible Camion 3: {consumo_total_truck3:.4f} L")

    # Retornar resultados
    return {
        'consumo_total_truck1': consumo_total_truck1,
        'consumo_total_truck2': consumo_total_truck2,
        'consumo_total_truck3': consumo_total_truck3,
        'x1': x1, 'v1': v1,
        'x2': x2, 'v2': v2,
        'x3': x3, 'v3': v3,
        'P_truck1': P_truck1, 'P_truck2': P_truck2, 'P_truck3': P_truck3,
        'fuel_truck1': fuel_truck1, 'fuel_truck2': fuel_truck2, 'fuel_truck3': fuel_truck3,
        'dist12': dist12, 'dist23': dist23,
        'tiempo': tiempo,
        'omega_hz': omega_hz
    }

def plot_truck_results(results, output_dir):
    # Crear el directorio si no existe
    os.makedirs(output_dir, exist_ok=True)

    tiempo = results['tiempo']
    omega_hz = results['omega_hz']
    x1, x2, x3 = results['x1'], results['x2'], results['x3']
    v1, v2, v3 = results['v1'], results['v2'], results['v3']
    P_truck1, P_truck2, P_truck3 = results['P_truck1'], results['P_truck2'], results['P_truck3']
    fuel_truck1, fuel_truck2, fuel_truck3 = results['fuel_truck1'], results['fuel_truck2'], results['fuel_truck3']
    dist12, dist23 = results['dist12'], results['dist23']
    consumo_total_truck1 = results['consumo_total_truck1']
    consumo_total_truck2 = results['consumo_total_truck2']
    consumo_total_truck3 = results['consumo_total_truck3']

    # Posición
    plt.figure(figsize=(10, 4))
    plt.plot(tiempo, x1 / 1000, label="Camión 1")
    plt.plot(tiempo, x2 / 1000, ':', label="Camión 2")
    plt.plot(tiempo, x3 / 1000, '--', label="Camión 3")
    plt.xlabel("Tiempo (s)")
    plt.ylabel("Posición (km)")
    plt.title(f"Evolución de la posición (Viento Oscilatorio {omega_hz} Veces por segundo)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'position_evolution.png'))

    # Velocidad
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
    plt.savefig(os.path.join(output_dir, 'velocity_evolution.png'))

    # Consumo de combustible
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
    plt.savefig(os.path.join(output_dir, 'fuel_consumption.png'))

    # Distancias entre camiones
    plt.figure(figsize=(10, 4))
    plt.plot(tiempo, dist12, '--r', label='Dist12')
    plt.plot(tiempo, dist23, '--g', label='Dist23')
    plt.grid()
    plt.title('Diferencias de distancia por tiempo')
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Distancia (m)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'distance_differences.png'))

    # Consumo total (barras)
    plt.figure(figsize=(10, 8))
    trucks = ['Camión 1', 'Camión 2', 'Camión 3']
    totales = [consumo_total_truck1, consumo_total_truck2, consumo_total_truck3]
    colors = ['blue', 'orange', 'green']
    bars = plt.bar(trucks, totales, color=colors)
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, height, f'{height:.2f}', ha='center', va='bottom', fontsize=10)
    plt.ylabel("Combustible total (L)")
    plt.title(f"Consumo total de combustible (Viento Oscilatorio {omega_hz} veces por segundo)")
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'total_fuel_consumption.png'))

    # Comparación de potencia
    figpower, axspower = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    axspower[0].plot(tiempo, P_truck1, '-b', label="Potencia Camión 1")
    axspower[0].set_ylabel("Potencia [W]")
    axspower[0].legend()
    axspower[0].grid(True)
    axspower[1].plot(tiempo, P_truck2, '-g', label="Potencia Camión 2")
    axspower[1].set_ylabel("Potencia [W]")
    axspower[1].legend()
    axspower[1].grid(True)
    axspower[2].plot(tiempo, P_truck3, '-r', label="Potencia Camión 3")
    axspower[2].set_xlabel("Tiempo [s]")
    axspower[2].set_ylabel("Potencia [W]")
    axspower[2].legend()
    axspower[2].grid(True)
    figpower.suptitle(f"Comparación de potencia (Viento Oscilatorio {omega_hz} Veces por segundo)")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(os.path.join(output_dir, 'power_comparison.png'))
    
    plt.close('all')
    print(f'Succesfully Printed into {output_dir}')

results1 = simulate_trucks(mass=30000, nominalPower=200e3, CD0=0.7)
plot_truck_results(results1,output_dir='resultados_simulacion/power_200')

results2 = simulate_trucks(mass=30000, nominalPower=50e3, CD0=0.7)
plot_truck_results(results2,output_dir='resultados_simulacion/power_50')

results3 = simulate_trucks(mass=30000, nominalPower=100e3, CD0=0.7)
plot_truck_results(results2,output_dir='resultados_simulacion/power_100')



results4 = simulate_trucks(mass=30000, nominalPower=100e3, CD0=0.7)
plot_truck_results(results4,output_dir='resultados_simulacion/cd_0.7')

results5 = simulate_trucks(mass=30000, nominalPower=100e3, CD0=0.5)
plot_truck_results(results5,output_dir='resultados_simulacion/cd_0.5')

results6 = simulate_trucks(mass=30000, nominalPower=100e3, CD0=1.2)
plot_truck_results(results6,output_dir='resultados_simulacion/cd_1.2')



results7 = simulate_trucks(mass=50000, nominalPower=100e3, CD0=0.7)
plot_truck_results(results7,output_dir='resultados_simulacion/m_50ton')

results8 = simulate_trucks(mass=100000, nominalPower=100e3, CD0=0.7)
plot_truck_results(results8,output_dir='resultados_simulacion/m_100ton')

results9 = simulate_trucks(mass=15000, nominalPower=100e3, CD0=0.7)
plot_truck_results(results9,output_dir='resultados_simulacion/m_15ton')




# results2 = simulate_trucks(mass=30000, nominalPower=120e3, CD0=0.8)
# plot_truck_results(results2,output_dir='resultados_simulacion/omega_0.1Hz')