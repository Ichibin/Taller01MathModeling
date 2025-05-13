
def tempInterp(z_eval):
    # This function interpolates the temperatures from the gas in Kelvin, depending on the height 
    
    # Returns the temperature value at said height in Celcius
    import numpy as np  
    import matplotlib.pyplot as plt
    
    z_data = np.array([10e-2, 22e-2, 78e-2, 112e-2])
    T_data = np.array([550, 700, 1000, 850])
    # T_data = T_data+273.15
    # print('T-data',T_data)
    
    T_g = lambda z: np.interp(z, z_data, T_data)
    
    # Generate values for visualization
    
    # Compute temperature at specific height
    T_at_z = T_g(z_eval)
    # print(f"Temperatura a z = {z_eval} m: {T_at_z:.2f} [C]")
    return T_at_z


# def interpolationMachine(temp_eval):
#     import numpy as np      
#     #This function interpolates Nitrogen Values from Fluids Mechanics Cengel, in order to evaluate these properties as the height changes. 
    
#     #Recives as an Argument the Temperature in Celcius
#     #Returns the fluid ploperties at said Temp.
    
#     temps = np.array([100,150,200,300,400,500, 1000])  # in °C
#     density = np.array([.9149,.8068,.7215,.5956,.5072,0.4416,0.2681])
#     thermalConductivity = np.array([0.03090, 0.03418, 0.03727, 0.04309, 0.04848, 0.05358, 0.07938])
#     prandtl = np.array([.7056,.7025,.7025,.7078,.7153,.7215,.7022])
#     kinematicViscosity = np.array([2.289e-5,2.851e-5,3.457e-5,4.783e-5,6.242e-5,7.816e-5,1.713e-4])
#     dinamycViscosity = np.array([2.094e-5,2.300e-5,2.494e-5,2.849e-5,3.166e-5,3.451e-5,4.594e-5])
    
#     values = [density,thermalConductivity,prandtl,kinematicViscosity,dinamycViscosity]

#     polysVal = []
    
#     for prop in values:
#         polys = np.polyfit(temps,prop,3)
#         polysVal.append(np.polyval(polys,temp_eval))
        
        
#         #    rho,          k    ,         pr,      nu,         mu
#     return polysVal[0],polysVal[1],polysVal[2],polysVal[3],polysVal[4]


def dragCoeff (reynolds):
    Cd = (24/reynolds)*(1+.15*reynolds**.687)
    return Cd
    



# def dragCoeff(reynolds):
#     if reynolds < 3e-5:
#         return 0
#     Cd = (24 / reynolds) * (1 + 0.15 * reynolds**0.687)
#     return Cd



def interpolationMachine(temp_eval):
    from scipy.interpolate import interp1d
    import numpy as np
    
    temps = np.array([100,150,200,300,400,500,1000])
    props = [
        np.array([.9149,.8068,.7215,.5956,.5072,0.4416,0.2681]),  # density
        np.array([0.03090, 0.03418, 0.03727, 0.04309, 0.04848, 0.05358, 0.07938]),  # k
        np.array([.7056,.7025,.7025,.7078,.7153,.7215,.7022]),  # Pr
        np.array([2.289e-5,2.851e-5,3.457e-5,4.783e-5,6.242e-5,7.816e-5,1.713e-4]),  # ν
        np.array([2.094e-5,2.300e-5,2.494e-5,2.849e-5,3.166e-5,3.451e-5,4.594e-5])   # μ
    ]
    
    results = []
    for prop in props:
        interp_fn = interp1d(temps, prop, kind='cubic', fill_value='extrapolate')
        results.append(interp_fn(temp_eval))
    return tuple(results)
