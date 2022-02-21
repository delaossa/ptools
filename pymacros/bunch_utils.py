import numpy as np
from scipy.constants import c


def gaussian_bunch(q_tot, n_part, gamma0, s_g, s_z, emit_x, s_x,
                   z_c=0., tf=0, x_c=0.):
    n_part = int(n_part)
    
    z = z_c + s_z * np.random.standard_normal(n_part)
    x = x_c + s_x * np.random.standard_normal(n_part)
    y = s_x * np.random.standard_normal(n_part)
        
    gamma = np.random.normal(gamma0, s_g, n_part)

    s_ux = emit_x / s_x
    ux = s_ux * np.random.standard_normal(n_part)
    uy = s_ux * np.random.standard_normal(n_part)

    uz = np.sqrt((gamma ** 2 - 1) - ux ** 2 - uy ** 2)

    if tf != 0.:
        x = x - ux * c * tf / gamma
        y = y - uy * c * tf / gamma
        z = z - uz * c * tf / gamma

    q = np.ones(n_part) * q_tot / n_part

    return x, y, z, ux, uy, uz, q


def flattop_bunch(q_tot, n_part, gamma0, s_g, length, s_z, emit_x, s_x,
                  z_c=0., tf=0, x_c=0.):
    n_part = int(n_part)

    norma = length + np.sqrt(2 * np.pi) * s_z
    n_plat = int(n_part * length / norma)
    n_gaus = int(n_part * np.sqrt(2 * np.pi) * s_z / norma)

    # Create flattop and gaussian profiles
    z_plat = np.random.uniform(0., length, n_plat)
    z_gaus = s_z * np.random.standard_normal(n_gaus)

    # Concatenate both profiles
    z = np.concatenate((z_gaus[np.where(z_gaus <= 0)],
                        z_plat,
                        z_gaus[np.where(z_gaus > 0)] + length))
       
    z = z - length / 2. + z_c  # shift to final position
    
    n_part = len(z)
    x = x_c + s_x * np.random.standard_normal(n_part)
    y = s_x * np.random.standard_normal(n_part)
        
    gamma = np.random.normal(gamma0, s_g, n_part)

    s_ux = emit_x / s_x
    ux = s_ux * np.random.standard_normal(n_part)
    uy = s_ux * np.random.standard_normal(n_part)

    uz = np.sqrt((gamma ** 2 - 1) - ux ** 2 - uy ** 2)

    if tf != 0.:
        x = x - ux * c * tf / gamma
        y = y - uy * c * tf / gamma
        z = z - uz * c * tf / gamma
    
    q = np.ones(n_part) * q_tot / n_part

    return x, y, z, ux, uy, uz, q
