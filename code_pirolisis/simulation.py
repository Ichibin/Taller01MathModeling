
import numpy as np
import matplotlib.pyplot as plt
from plotterDef import plot_particle_analysis, plottingCp, plottingVolFlow
import os



def particleSimulation(volumetricFlow, Cp_p):  # Modelo Sencillo
    from crossArea import crossArea
    from physicalFunc import tempInterp, interpolationMachine, dragCoeff
    import numpy as np

    # --- Parameters ---
    volumetricFlow = volumetricFlow / (1000 * 60)  # liters/min → m³/s
    totalHeight = 120e-2  # m
    g = 9.81

    # Particle Parameters
    densityP = 810  # kg/m³
    r = (150e-6) / 2  # m
    A_p = 4 * np.pi * r**2  # Surface area (m²)
    m_p = densityP * (4/3) * np.pi * r**3  # Mass (kg)

    # --- Time Setup ---
    t0 = 0
    tf = 3
    dt = 0.0001
    time = np.arange(t0, tf + dt, dt)
    N = len(time)

    # --- Initialization ---
    y = np.zeros(N)
    v = np.zeros(N)
    Tp = np.zeros(N)
    

    rhoPlot = np.zeros(N)
    kPlot = np.zeros(N)
    prPlot = np.zeros(N)
    nuPlot = np.zeros(N)
    muPlot = np.zeros(N)
    velfluidPlot = np.zeros(N)
    hPlot = np.zeros(N)
    dragForcePlot = np.zeros(N)
    acceleration = np.zeros(N)
    reynoldsPlot = np.zeros(N)
    n2Plot = np.zeros(N)

    # --- Initial Conditions ---
    y[0] = 0
    v[0] = 0
    Tp[0] = 550  # °C
    # --- Time Marching ---
    for t in range(N - 1):
        
        
        velFluid = volumetricFlow / crossArea(y[t])
        nitrogenTemp = tempInterp(y[t])
        rho, k, pr, nu, mu = interpolationMachine(nitrogenTemp)
        reynolds = (abs(velFluid - v[t]) * 2 * r) / nu
        dragCoefficient = dragCoeff(reynolds)
        # dv_dt = (0.5 * rho * dragCoefficient * (np.pi*r**2) * (velFluid - v[t])**2) / m_p
        dv_dt = (0.5 * rho * dragCoefficient * (np.pi * r**2) * (velFluid - v[t])**2) / m_p - g
        h = (k / (2 * r)) * (2 + (0.4 * reynolds**0.5 + 0.06 * reynolds**(2/3)) * pr**0.4)

        # Update vectors
        rhoPlot[t] = rho
        kPlot[t] = k
        prPlot[t] = pr
        nuPlot[t] = nu
        muPlot[t] = mu
        velfluidPlot[t] = velFluid
        dragForcePlot[t] = 0.5 * rho * A_p * dragCoefficient * (velFluid - v[t])**2
        hPlot[t] = h
        acceleration[t] = dv_dt
        reynoldsPlot[t] = reynolds
        n2Plot[t] = nitrogenTemp

        # Update kinematics and temperature
        v[t + 1] = v[t] + dt * dv_dt
        y[t + 1] = y[t] + dt * v[t]
        Tp[t + 1] = Tp[t] + dt * (h * A_p / (m_p * Cp_p)) * (nitrogenTemp - Tp[t])
        
        
        if y[t] > totalHeight:
            y[t:] = y[t]
            v[t:] = v[t]
            Tp[t:] = Tp[t]
            reynoldsPlot[t:] = reynoldsPlot[t]
            rhoPlot[t:] = rhoPlot[t]
            kPlot[t:] = kPlot[t]
            prPlot[t:] = prPlot[t]
            nuPlot[t:] = nuPlot[t]
            muPlot[t:] = muPlot[t]
            velfluidPlot[t:] = velfluidPlot[t]
            dragForcePlot[t:] = dragForcePlot[t]
            hPlot[t:] = hPlot[t]
            acceleration[t:] = acceleration[t]
            n2Plot[t:] =n2Plot[t]
            break
        

    return time, y, v, reynoldsPlot, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, Tp, velfluidPlot, n2Plot



def particleSimulationRadiationandFlotation (volumetricFlow, Cp_p): #modelo refinado
    
    #importing functions \ interpolations
    from crossArea import crossArea
    from physicalFunc import tempInterp,interpolationMachine, dragCoeff
    

    ## Parameters
    volumetricFlow = volumetricFlow/(1000*60) #units convertion
    totalHeight = 120E-2 

    # points = np.arange(0,totalHeight+0.01,0.01)

    #Particle Parameteres
    densityP  = 810 # Density of the particle (kg/m^3)
    r = (150e-6) / 2  # Radius of the particle (m)
    emissivitySteel = 0.89
    stefanBoltzmann = 5.67e-8 # W/m^2.K^4
    g = 9.81 
    ## Parameters


    ##Prepping Vectors
    t0 = 0 ; tf = 3; dt = 0.0001

    time = np.arange(t0,tf+dt,dt)
    N = len(time)
    y = np.zeros(N)
    v = np.zeros(N)
    K= 273.15
    ##

    Tp = np.zeros(N)
    reynolds = np.zeros(N)
    # printReynolds = np.zeros(N)

    A_p = 4 * np.pi * r**2  # Surface area for convection
    m_p = densityP *  (4/3) * np.pi * r**3 # Mass of the particle


    

    ## Vectors to Plot##
    rhoPlot = np.zeros(N)
    kPlot = np.zeros(N)
    prPlot = np.zeros(N)
    nuPlot = np.zeros(N)
    muPlot = np.zeros(N)
    velfluidPlot = np.zeros(N)
    hPlot = np.zeros(N)
    dragForcePlot = np.zeros(N)
    acceleration = np.zeros(N)
    reynoldsPlot = np.zeros(N)
    n2Plot = np.zeros(N)
    
    ## intial Condtidions
    t = 0
    y[0]= 0 
    v[0] = 0
    Tp[0] = 550+K #Kelvin
    ## intial Condtidions

    for t in range(N - 1):
        velFluid = volumetricFlow / crossArea(y[t])
        nitrogenTemp = tempInterp(y[t])
        rho, k, pr, nu, mu = interpolationMachine(nitrogenTemp)
        reynolds = (abs(velFluid - v[t]) * 2 * r) / nu
        dragCoefficient = dragCoeff(reynolds)
        dv_dt = ((rho *(4/3) * np.pi * r**3* g)- (m_p * g) + (0.5 * rho * dragCoefficient * (np.pi * r**2) * (velFluid - v[t])**2)  # drag
        ) / m_p
        h = (k / (2 * r)) * (2 + (0.4 * reynolds**0.5 + 0.06 * reynolds**(2/3)) * pr**0.4)
        
        # Update vectors
        rhoPlot[t] = rho
        kPlot[t] = k
        prPlot[t] = pr
        nuPlot[t] = nu
        muPlot[t] = mu
        velfluidPlot[t] = velFluid
        dragForcePlot[t] = 0.5 * rho * A_p * dragCoefficient * (velFluid - v[t])**2
        hPlot[t] = h
        acceleration[t] = dv_dt
        reynoldsPlot[t] = reynolds
        n2Plot[t] =nitrogenTemp
        
        v[t + 1] = v[t] + dt * dv_dt
        y[t + 1] = y[t] + dt * v[t]
        T_gas_K = nitrogenTemp + K
        T_particle_K = Tp[t]

        conv = h * A_p * (T_gas_K - T_particle_K)
        rad = emissivitySteel * stefanBoltzmann * A_p * (T_particle_K**4 - T_gas_K**4)
        dTdt = (conv - rad) / (m_p * Cp_p)
        Tp[t + 1] = Tp[t] + dt * dTdt
        
        if y[t] > totalHeight:
            y[t:] = y[t]
            v[t:] = v[t]
            Tp[t:] = Tp[t]-K
            reynoldsPlot[t:] = reynoldsPlot[t]
            rhoPlot[t:] = rhoPlot[t]
            kPlot[t:] = kPlot[t]
            prPlot[t:] = prPlot[t]
            nuPlot[t:] = nuPlot[t]
            muPlot[t:] = muPlot[t]
            velfluidPlot[t:] = velfluidPlot[t]
            dragForcePlot[t:] = dragForcePlot[t]
            hPlot[t:] = hPlot[t]
            acceleration[t:] = acceleration[t]
            n2Plot[t:] = n2Plot[t]
            break
        
    # TpCelcius = Tp-K
    
    return time, y, v, reynoldsPlot, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, Tp, velfluidPlot, n2Plot


time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, Tp, velfluidPlot, n2Plot = particleSimulation(230.5,1.52*1000)
plot_particle_analysis(time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, Tp, velfluidPlot,n2Plot, 'simple_Sim')

time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, TpCelcius, velfluidPlot, n2Plot= particleSimulationRadiationandFlotation(230.5, 1.52 * 1000)
plot_particle_analysis(time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, TpCelcius, velfluidPlot,n2Plot,'rad_and_Flot')



Cp_values = np.array([1.52, 0.52, 2.52, 5.])
Cp_values = Cp_values*1000

for i in Cp_values:
    time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, Tp, velfluidPlot, n2Plot = particleSimulation(230.5,i)
    plottingCp(time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, Tp, velfluidPlot, n2Plot,i,'sensibility_noRad')
    time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, TpCelcius, velfluidPlot, n2Plot= particleSimulationRadiationandFlotation(230.5,i)
    plottingCp(time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, Tp, velfluidPlot, n2Plot,i,'sensibility_Rad_and_Float')
    
    

volFlows = np.array([230.5, 500., 50., 1000.])
for i in volFlows:
    time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, Tp, velfluidPlot, n2Plot = particleSimulation(i,1.52 *1000)
    plottingVolFlow(time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, Tp, velfluidPlot, n2Plot,i,'sensibility_noRad')
    time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, TpCelcius, velfluidPlot, n2Plot= particleSimulationRadiationandFlotation(i,1.52 *1000)
    plottingVolFlow(time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, Tp, velfluidPlot, n2Plot,i,'sensibility_Rad_and_Float')


print('All Done Sir.')