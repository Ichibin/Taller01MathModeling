


def crossArea(x):
    
    # this calculates the cross sectional area per xs, from the air channel

    # Returns the value of said area in m^2

    import numpy as np 
    lowerD = 6.7E-2 # Lower diameter
    upperD = 8.1E-2 # Upper diameter
    area = 0
    if x <40E-2:
        area = np.pi*lowerD**2/4
    elif 40E-2 <= x <= 60E-2:
        z = x - 40E-2  # Distance into the transition region
        diameter_x = lowerD + (upperD - lowerD) * (z / (60E-2 - 40E-2))  # Correct scaling
        area = np.pi * (diameter_x / 2) ** 2  # Compute area
    elif x > 60E-2:
        area = np.pi*upperD**2/4
    else:
        raise ValueError(f"x value {x} is out of bounds (expected between 0 and 60 cm).")
    return area
