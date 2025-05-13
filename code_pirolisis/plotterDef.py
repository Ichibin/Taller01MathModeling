import matplotlib.pyplot as plt
import os
def plot_particle_analysis(time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, Tp, velfluidPlot,nitrogenTemp,dir):
    
    
    directory = f'./{dir}'
    os.makedirs(directory, exist_ok=True)
    
    # Plot 1: Particle Motion Analysis (Position and Speed vs Time)
    fig, ax = plt.subplots(1, 2, figsize=(8, 5))
    fig.suptitle(f'Particle Motion Analysis ')

    ax[0].set_title(f'Position-Time from Particle ')
    ax[0].plot(time, y, '-r')
    ax[0].set_xlabel('Time /s')
    ax[0].set_ylabel('Position /m')
    ax[0].grid(True)

    ax[1].set_title(f'Speed-Time from Particle ')
    ax[1].plot(time, v, '-b')
    ax[1].set_xlabel('Time /s')
    ax[1].set_ylabel('Speed /m')
    ax[1].grid(True)

    plt.tight_layout()
    fig.savefig(f'{directory}/Particle Motion Analysis.png', dpi=300)
    plt.close()

    # Plot 2: Reynolds vs Height
    plt.figure(figsize=(8, 5))
    plt.plot(y, reynolds)
    plt.title(f'Reynolds - Height ')
    plt.xlabel('Height /m')
    plt.ylabel('Reynolds /1')
    plt.grid()
    plt.savefig(f'{directory}/Reynolds_Height.png', dpi=300)
    plt.close()

    # Plot 3: Fluid Properties (2x2 grid)
    figProps = plt.figure(figsize=(10, 8))
    figProps.suptitle(f'Fluid Properties over Height ')
    axProps = [
        plt.subplot(2, 2, 1),
        plt.subplot(2, 2, 2),
        plt.subplot(2, 2, 3),
        plt.subplot(2, 2, 4),
    ]

    axProps[0].set_title(f'Density vs Height ')
    axProps[0].plot(y, rhoPlot, 'b-', label='Density (ρ)')
    axProps[0].set_ylabel('Density (kg/m³)')
    axProps[0].legend()
    axProps[0].grid(True)

    axProps[1].set_title(f'Thermal Conductivity vs Height ')
    axProps[1].plot(y, kPlot, 'r-', label='Thermal Conductivity (k)')
    axProps[1].set_ylabel('Thermal Conductivity (W/m·K)')
    axProps[1].legend()
    axProps[1].grid(True)

    axProps[2].set_title(f'Prandtl Number vs Height')
    axProps[2].plot(y, prPlot, 'g-', label='Prandtl Number (Pr)')
    axProps[2].set_ylabel('Prandtl Number')
    axProps[2].legend()
    axProps[2].grid(True)

    axProps[3].set_title(f'Kinematic Viscosity vs Height ')
    axProps[3].plot(y, nuPlot, 'm-', label='Kinematic Viscosity (ν)')
    axProps[3].set_ylabel('Kinematic Viscosity (m²/s)')
    axProps[3].legend()
    axProps[3].grid(True)

    x_limits = axProps[0].get_xlim()
    for ax in axProps:
        ax.set_xlabel('Height /m')
        ax.set_xlim(x_limits)

    figProps.subplots_adjust(top=0.9, hspace=0.4, wspace=0.3)
    figProps.savefig(f'{directory}/fluidProperties.png', dpi=300)
    plt.close()

    # # Plot 4: Particle Temperature vs Time
    # plt.figure(figsize=(8, 5))
    # plt.plot(time, TpCelcius, '-g')
    # plt.title(f'Particle Temperature - Time ')
    # plt.xlabel('Time /s')
    # plt.ylabel('Particle Temp (°C)')
    # plt.grid()
    # plt.savefig(f'{directory}/Particle Temperature.png', dpi=300)
    

    # Plot 4: Particle Temperature vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(time, Tp, '-g', label='Particle Temperature')  # Added label for clarity
    # Find the maximum temperature and corresponding time
    max_temp = max(Tp)
    max_time = time[Tp.argmax()]
    # Plot a marker at the maximum value
    plt.scatter(max_time, max_temp, color='red', s=100, marker='*', label=f'Max Temp: {max_temp:.8f}°C')
    plt.title(f'Particle Temperature - Time ')
    plt.xlabel('Time /s')
    plt.ylabel('Particle Temp (°C)')
    plt.grid()
    plt.legend()  # Add legend to show the label of the marker
    plt.savefig(f'{directory}/Particle Temperature.png', dpi=300)
    plt.close()

    # Plot 5: Fluid Speed vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(time, velfluidPlot)
    plt.title(f'Fluid Speed over Time ')
    plt.xlabel('Time /s')
    plt.ylabel('Fluid Speed m/s')
    plt.grid()
    plt.savefig(f'{directory}/Fluid Speed_Time.png', dpi=300)
    plt.close()

    # Plot 6: Comparison of Fluid and Particle Speed
    plt.figure(figsize=(8, 5))
    plt.plot(y, velfluidPlot, ':r', label='Fluid Speed')
    plt.plot(y, v, ':g', label='Particle Speed')
    plt.title(f'Comparison of Fluid and Particle Speed Over Height  ')
    plt.legend()
    plt.xlabel('Height /m')
    plt.ylabel('Speed m/s')
    plt.grid()
    plt.savefig(f'{directory}/Comparison of Fluid and Particle Speed Over Height.png', dpi=300)
    plt.close() 
    # Plot 7: Convection Coefficient vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(time, hPlot, ':r', label='Convection Coefficient')
    plt.title(f'Convection Coeff. over Time ')
    plt.xlabel('Time /s')
    plt.ylabel('Convection Coefficient W/m²·K')
    plt.grid()
    plt.savefig(f'{directory}/Convection Coeff.png', dpi=300)
    plt.close()
    
    # Plot 8: Drag Force
    plt.figure(figsize=(8, 5))
    plt.plot(time, dragForcePlot, ':r', label='DragForce')
    plt.title('DragForce over Time')
    plt.xlabel('Time /s')
    plt.ylabel('DragForce N')
    plt.grid()
    plt.savefig(f'{directory}/DragForce.png', dpi=300)
    plt.close()
    
    
    # Plot 8: Comparative Nitrogen and Particle Temp
    plt.figure(figsize=(8, 5))
    plt.plot(y, nitrogenTemp, ':r', label='N2 Temp')
    plt.plot(y, Tp, ':b', label='Particle Temp')
    plt.title('Comparative Nitrogen and Particle Temp')
    plt.xlabel('Height /m')
    plt.ylabel('Temperature Celcius')
    plt.grid()
    plt.legend()
    plt.savefig(f'{directory}/Comparative Nitrogen and Particle Temp.png', dpi=300)
    plt.close()
    
    print(f'Succesfully Printed into {directory}')
    plt.close('all')
    
    
def plottingCp(time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, Tp, velfluidPlot,nitrogenTemp,Cp_p,dir):
    
    # Create directory for saving plots
    directory = f'./{dir}/{Cp_p/1000}'
    os.makedirs(directory, exist_ok=True)

    # Plot 1: Particle Motion Analysis (Position and Speed vs Time)
    fig, ax = plt.subplots(1, 2, figsize=(8, 5))
    fig.suptitle(f'Particle Motion Analysis ($C_p$ = {Cp_p/1000:.2f} kJ/kg·K)')

    ax[0].set_title(f'Position-Time from Particle ')
    ax[0].plot(time, y, '-r')
    ax[0].set_xlabel('Time /s')
    ax[0].set_ylabel('Position /m')
    ax[0].grid(True)

    ax[1].set_title(f'Speed-Time from Particle ')
    ax[1].plot(time, v, '-b')
    ax[1].set_xlabel('Time /s')
    ax[1].set_ylabel('Speed /m')
    ax[1].grid(True)

    plt.tight_layout()
    fig.savefig(f'{directory}/Particle Motion Analysis.png', dpi=300)

    # Plot 2: Reynolds vs Height
    plt.figure(figsize=(8, 5))
    plt.plot(y, reynolds)
    plt.title(f'Reynolds - Height ($C_p$ = {Cp_p/1000:.2f} kJ/kg·K)')
    plt.xlabel('Height /m')
    plt.ylabel('Reynolds /1')
    plt.grid()
    plt.savefig(f'{directory}/Reynolds_Height.png', dpi=300)

    # Plot 3: Fluid Properties (2x2 grid)
    figProps = plt.figure(figsize=(10, 8))
    figProps.suptitle(f'Fluid Properties over Height ($C_p$ = {Cp_p/1000:.2f} kJ/kg·K)')
    axProps = [
        plt.subplot(2, 2, 1),
        plt.subplot(2, 2, 2),
        plt.subplot(2, 2, 3),
        plt.subplot(2, 2, 4),
    ]

    axProps[0].set_title(f'Density vs Height ')
    axProps[0].plot(y, rhoPlot, 'b-', label='Density (ρ)')
    axProps[0].set_ylabel('Density (kg/m³)')
    axProps[0].legend()
    axProps[0].grid(True)

    axProps[1].set_title(f'Thermal Conductivity vs Height ')
    axProps[1].plot(y, kPlot, 'r-', label='Thermal Conductivity (k)')
    axProps[1].set_ylabel('Thermal Conductivity (W/m·K)')
    axProps[1].legend()
    axProps[1].grid(True)

    axProps[2].set_title(f'Prandtl Number vs Height')
    axProps[2].plot(y, prPlot, 'g-', label='Prandtl Number (Pr)')
    axProps[2].set_ylabel('Prandtl Number')
    axProps[2].legend()
    axProps[2].grid(True)

    axProps[3].set_title(f'Kinematic Viscosity vs Height ')
    axProps[3].plot(y, nuPlot, 'm-', label='Kinematic Viscosity (ν)')
    axProps[3].set_ylabel('Kinematic Viscosity (m²/s)')
    axProps[3].legend()
    axProps[3].grid(True)

    x_limits = axProps[0].get_xlim()
    for ax in axProps:
        ax.set_xlabel('Height /m')
        ax.set_xlim(x_limits)

    figProps.subplots_adjust(top=0.9, hspace=0.4, wspace=0.3)
    figProps.savefig(f'{directory}/fluidProperties.png', dpi=300)

    # # Plot 4: Particle Temperature vs Time
    # plt.figure(figsize=(8, 5))
    # plt.plot(time, TpCelcius, '-g')
    # plt.title(f'Particle Temperature - Time ($C_p$ = {Cp_p:.2f} kJ/kg·K)')
    # plt.xlabel('Time /s')
    # plt.ylabel('Particle Temp (°C)')
    # plt.grid()
    # plt.savefig(f'{directory}/Particle Temperature.png', dpi=300)
    

    # Plot 4: Particle Temperature vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(time, Tp, '-g', label='Particle Temperature')  # Added label for clarity
    # Find the maximum temperature and corresponding time
    max_temp = max(Tp)
    max_time = time[Tp.argmax()]
    # Plot a marker at the maximum value
    plt.scatter(max_time, max_temp, color='red', s=100, marker='*', label=f'Max Temp: {max_temp:.8f}°C')
    plt.title(f'Particle Temperature - Time ($C_p$ = {Cp_p/1000:.2f} kJ/kg·K)')
    plt.xlabel('Time /s')
    plt.ylabel('Particle Temp (°C)')
    plt.grid()
    plt.legend()  # Add legend to show the label of the marker
    plt.savefig(f'{directory}/Particle Temperature.png', dpi=300)


    # Plot 5: Fluid Speed vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(time, velfluidPlot)
    plt.title(f'Fluid Speed over Time ($C_p$ = {Cp_p/1000:.2f} kJ/kg·K)')
    plt.xlabel('Time /s')
    plt.ylabel('Fluid Speed m/s')
    plt.grid()
    plt.savefig(f'{directory}/Fluid Speed_Time.png', dpi=300)

    # Plot 6: Comparison of Fluid and Particle Speed
    plt.figure(figsize=(8, 5))
    plt.plot(time, velfluidPlot, ':r', label='Fluid Speed')
    plt.plot(time, v, ':g', label='Particle Speed')
    plt.title(f'Comparison of Fluid and Particle Speed Over Time ($C_p$ = {Cp_p/1000:.2f} kJ/kg·K)')
    plt.legend()
    plt.xlabel('Time /s')
    plt.ylabel('Speed m/s')
    plt.grid()
    plt.savefig(f'{directory}/Comparison of Fluid and Particle Speed Over Time.png', dpi=300)

    # Plot 7: Convection Coefficient vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(time, hPlot, ':r', label='Convection Coefficient')
    plt.title(f'Convection Coeff. over Time ($C_p$ = {Cp_p/1000:.2f} kJ/kg·K)')
    plt.xlabel('Time /s')
    plt.ylabel('Convection Coefficient W/m²·K')
    plt.grid()
    plt.savefig(f'{directory}/Convection Coeff.png', dpi=300)
    
    # Plot 8: Drag Force
    plt.figure(figsize=(8, 5))
    plt.plot(time, dragForcePlot, ':r', label='DragForce')
    plt.title('DragForce over Time')
    plt.xlabel('Time /s')
    plt.ylabel('DragForce N')
    plt.grid()
    plt.savefig(f'{directory}/DragForce.png', dpi=300)
    
    # Plot 8: Comparative Nitrogen and Particle Temp
    plt.figure(figsize=(8, 5))
    plt.plot(y, nitrogenTemp, ':r', label='N2 Temp')
    plt.plot(y, Tp, ':b', label='Particle Temp')
    plt.title(f'Comparative Nitrogen and Particle Temp ($C_p$ = {Cp_p/1000:.2f} kJ/kg·K)')
    plt.xlabel('Height /m')
    plt.ylabel('Temperature (°C)')
    plt.grid()
    plt.legend()
    plt.savefig(f'{directory}/Comparative Nitrogen and Particle Temp.png', dpi=300)
    plt.close()
    
    print(f'Succesfully Printed into {directory}')
    plt.close('all')
    

def plottingVolFlow(time, y, v, reynolds, rhoPlot, kPlot, prPlot, nuPlot, dragForcePlot, hPlot, Tp, velfluidPlot,nitrogenTemp,volFlow,dir):
    
    # Create directory for saving plots
    directory = f'./{dir}/Vol_{volFlow}'
    os.makedirs(directory, exist_ok=True)
    # Plot 1: Particle Motion Analysis (Position and Speed vs Time)
    fig, ax = plt.subplots(1, 2, figsize=(8, 5))
    fig.suptitle(f'Particle Motion Analysis (Vol. Flow = {volFlow:.2f} L/min)')

    ax[0].set_title(f'Position-Time from Particle')
    ax[0].plot(time, y, '-r')
    ax[0].set_xlabel('Time /s')
    ax[0].set_ylabel('Position /m')
    ax[0].grid(True)

    ax[1].set_title(f'Speed-Time from Particle ')
    ax[1].plot(time, v, '-b')
    ax[1].set_xlabel('Time /s')
    ax[1].set_ylabel('Speed /m')
    ax[1].grid(True)

    plt.tight_layout()
    fig.savefig(f'{directory}/Particle Motion Analysis.png', dpi=300)

    # Plot 2: Reynolds vs Height
    plt.figure(figsize=(8, 5))
    plt.plot(y, reynolds)
    plt.title(f'Reynolds - Height (Vol. Flow = {volFlow:.2f} L/min)')
    plt.xlabel('Height /m')
    plt.ylabel('Reynolds /1')
    plt.grid()
    plt.savefig(f'{directory}/Reynolds_Height.png', dpi=300)

    # Plot 3: Fluid Properties (2x2 grid)
    figProps = plt.figure(figsize=(10, 8))
    figProps.suptitle(f'Fluid Properties over Height (Vol. Flow = {volFlow:.2f} L/min)')
    axProps = [
        plt.subplot(2, 2, 1),
        plt.subplot(2, 2, 2),
        plt.subplot(2, 2, 3),
        plt.subplot(2, 2, 4),
    ]

    axProps[0].set_title(f'Density vs Height ')
    axProps[0].plot(y, rhoPlot, 'b-', label='Density (ρ)')
    axProps[0].set_ylabel('Density (kg/m³)')
    axProps[0].legend()
    axProps[0].grid(True)

    axProps[1].set_title(f'Thermal Conductivity vs Height ')
    axProps[1].plot(y, kPlot, 'r-', label='Thermal Conductivity (k)')
    axProps[1].set_ylabel('Thermal Conductivity (W/m·K)')
    axProps[1].legend()
    axProps[1].grid(True)

    axProps[2].set_title(f'Prandtl Number vs Height ')
    axProps[2].plot(y, prPlot, 'g-', label='Prandtl Number (Pr)')
    axProps[2].set_ylabel('Prandtl Number')
    axProps[2].legend()
    axProps[2].grid(True)

    axProps[3].set_title(f'Kinematic Viscosity vs Height ')
    axProps[3].plot(y, nuPlot, 'm-', label='Kinematic Viscosity (ν)')
    axProps[3].set_ylabel('Kinematic Viscosity (m²/s)')
    axProps[3].legend()
    axProps[3].grid(True)

    x_limits = axProps[0].get_xlim()
    for ax in axProps:
        ax.set_xlabel('Height /m')
        ax.set_xlim(x_limits)

    figProps.subplots_adjust(top=0.9, hspace=0.4, wspace=0.3)
    figProps.savefig(f'{directory}/fluidProperties.png', dpi=300)

    # # Plot 4: Particle Temperature vs Time
    # plt.figure(figsize=(8, 5))
    # plt.plot(time, TpCelcius, '-g')
    # plt.title(f'Particle Temperature - Time (Vol. Flow = {volFlow:.2f} L/min)')
    # plt.xlabel('Time /s')
    # plt.ylabel('Particle Temp (°C)')
    # plt.grid()
    # plt.savefig(f'{directory}/Particle Temperature.png', dpi=300)
    # Plot 4: Particle Temperature vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(time, Tp, '-g', label='Particle Temperature')  # Added label for clarity
    # Find the maximum temperature and corresponding time
    max_temp = max(Tp)
    max_time = time[Tp.argmax()]
    # Plot a marker at the maximum value
    plt.scatter(max_time, max_temp, color='red', s=100, marker='*', label=f'Max Temp: {max_temp:.8f}°C')
    plt.title(f'Particle Temperature - Time (Vol. Flow = {volFlow:.2f} L/min)')
    plt.xlabel('Time /s')
    plt.ylabel('Particle Temp (°C)')
    plt.grid()
    plt.legend()  # Add legend to show the label of the marker
    plt.savefig(f'{directory}/Particle Temperature.png', dpi=300)

    # Plot 5: Fluid Speed vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(time, velfluidPlot)
    plt.title(f'Fluid Speed over Time (Vol. Flow = {volFlow:.2f} L/min)')
    plt.xlabel('Time /s')
    plt.ylabel('Fluid Speed m/s')
    plt.grid()
    plt.savefig(f'{directory}/Fluid Speed_Time.png', dpi=300)

    # Plot 6: Comparison of Fluid and Particle Speed
    plt.figure(figsize=(8, 5))
    plt.plot(time, velfluidPlot, ':r', label='Fluid Speed')
    plt.plot(time, v, ':g', label='Particle Speed')
    plt.title(f'Comparison of Fluid and Particle Speed Over Time (Vol. Flow = {volFlow:.2f} L/min)')
    plt.legend()
    plt.xlabel('Time /s')
    plt.ylabel('Speed m/s')
    plt.grid()
    plt.savefig(f'{directory}/Comparison of Fluid and Particle Speed Over Time.png', dpi=300)

    # Plot 7: Convection Coefficient vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(time, hPlot, ':r', label='Convection Coefficient')
    plt.title(f'Convection Coeff. over Time (Vol. Flow = {volFlow:.2f} L/min)')
    plt.xlabel('Time /s')
    plt.ylabel('Convection Coefficient W/m²·K')
    plt.grid()
    plt.savefig(f'{directory}/Convection Coeff.png', dpi=300)
    
    
    # Plot 8: Drag Force
    plt.figure(figsize=(8, 5))
    plt.plot(time, dragForcePlot, ':r', label='DragForce')
    plt.title('DragForce over Time')
    plt.xlabel('Time /s')
    plt.ylabel('DragForce N')
    plt.grid()
    plt.savefig(f'{directory}/DragForce.png', dpi=300)
    
    # Plot 8: Comparative Nitrogen and Particle Temp
    plt.figure(figsize=(8, 5))
    plt.plot(y, nitrogenTemp, ':r', label='N2 Temp')
    plt.plot(y, Tp, ':b', label='Particle Temp')
    plt.title(f'Comparative Nitrogen and Particle Temp ($C_p$ = {volFlow:.1f} L/min)')
    plt.xlabel('Height /m')
    plt.ylabel('Temperature (°C)')
    plt.grid()
    plt.legend()
    plt.savefig(f'{directory}/Comparative Nitrogen and Particle Temp.png', dpi=300)
    plt.close()
    
    print(f'Succesfully Printed into {directory}')
    plt.close('all')
