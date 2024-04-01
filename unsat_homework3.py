#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import matplotlib.pyplot as plt
import scienceplots
import pandas as pd
from scipy.optimize import curve_fit
from scipy.misc import derivative
from sklearn.metrics import r2_score

plt.style.use(['science','notebook', 'grid'])


# In[4]:


df = pd.read_excel("D:\\Users\\Prabhat Paudyal\\Documents\\for_hw3_unsat.xlsx", sheet_name = 1)
deg_sat = df['sat'].values
suction = df['S'].values
vol_water = df['theta'].values
fredlund_fit = [99.98914818, 99.92426358, 99.88121032, 99.38067384, 95.69972697, 92.30738497, 67.29883859, 24.88601575, 19.93857605, 16.85829771, 14.48288079, 12.58985287, 11.00972887, 9.654388606,
8.47203427,
7.428882803,
6.50033786,
5.667375315,
4.914818388,
4.230376478,
3.604029561,
3.027569316,
2.49423963,
1.998450355,
1.532542075,
1.100072356,
0.693930955,
0.307701296,
0
]


# In[5]:


#plot the SWCC in terms of degree of saturation

plt.plot(suction, deg_sat, color = 'green')
plt.xscale('log')
plt.xlabel('Matric Suction (Ψ), kPa')
plt.ylabel('Degree of saturation (S), %')

plt.show()


# In[6]:


deg_sat_bc = df['brooks_corey (fit, S)'].values
deg_sat_vg = df['vangenuchten (sat)'].values

plt.plot(suction, deg_sat, label = 'Provided Data')
plt.plot(suction, deg_sat_bc, label = 'Brooks & Corey')
plt.plot(suction, deg_sat_vg, label = 'Vangenuchten & Mualem')
plt.xscale('log')
plt.xlabel('Matric Suction (Ψ), kPa')
plt.ylabel('Degree of saturation (S), %')
plt.legend()

plt.show()


# In[7]:


plt.plot(suction, deg_sat_vg, label = 'Vangenuchten & Mualem')
plt.xscale('log')
plt.xlabel('Matric Suction (Ψ), kPa')
plt.ylabel('Degree of saturation (S), %')

plt.show()


# In[8]:


plt.plot(suction, fredlund_fit, label='Torres (2011)')
plt.plot(suction, deg_sat_vg, label='van Genuchten (1980)')
plt.xscale('log')
plt.xlabel('Matric Suction (Ψ), kPa')
plt.ylabel('Degree of saturation (S), %')
plt.legend(loc='upper right')

plt.show()


# In[9]:


def Brooks_Corey_hydraulic_list(k_s, psi_b, psi_list, eta):
    k_w_list = []  # Initialize an empty list to store the calculated k_w values
    for psi in psi_list:  # Iterate over each psi value in the input list
        if psi < psi_b:
            k_w = k_s  # If psi is less than psi_b, k_w is equal to k_s
        else:
            k_w = k_s * np.power((psi_b / psi), eta)  # Calculate k_w based on the Brooks-Corey equation
        k_w_list.append(k_w)  # Append the calculated k_w to the list of k_w values
    return k_w_list  # Return the list of calculated k_w values


# In[13]:


def VG_hydraulic_list(k_s, a, psi_list, n, m):
    k_w_list = []  # Initialize an empty list to store the calculated k_w values
    for psi in psi_list:
        k_w = k_s * np.power(1 - np.power(a * psi, n-1) * np.power(1 + np.power(a * psi, n), -m), 2) / np.power(1 + np.power(a * psi, n), m/2)                               
        k_w_list.append(k_w)  # Append the calculated k_w to the list of k_w values
    return k_w_list  # Return the list of calculated k_w values


# In[20]:


k_w = [3.20523E-06,
2.15118E-06,
1.83897E-06,
1.06415E-06,
8.72794E-07,
8.32269E-07,
7.12546E-07,
4.13565E-07,
3.17743E-07,
2.46433E-07,
1.8855E-07,
1.43689E-07,
1.09111E-07,
8.25988E-08,
6.23556E-08,
4.69587E-08,
3.52868E-08,
2.64646E-08,
1.98137E-08,
1.48112E-08,
1.10564E-08,
8.2432E-09,
6.13891E-09,
4.56718E-09,
3.38801E-09,
2.51847E-09,
1.87184E-09,
1.38723E-09,
1.07868E-09
]

plt.plot(suction, k_w)
plt.xscale('log')
plt.xlabel('Matric Suction (Ψ), kPa')
plt.ylabel('Unsaturated Hydraulic Conductivity (k), m/s')
plt.show()


# In[15]:


k_w


# In[17]:


k_w = Brooks_Corey_hydraulic_list(1e-4, 124.183, suction, 2.288)

plt.plot(suction, k_w)
plt.xscale('log')
plt.xlabel('Matric Suction (Ψ), kPa')
plt.ylabel('Unsaturated Hydraulic Conductivity (k), m/s')
plt.show()


# In[18]:


k_w


# In[38]:


def brooks_and_corey(suction_value, lambda_, s_r, s_s, psi):
    # Ensure parameters make physical sense
    sat = np.zeros_like(suction_value)
    for i, s in enumerate(suction_value):
        if s > psi:            
    sat[i] = s_r + (s_s - s_r) * np.power((psi / s), lambda_)
        else:
            sat[i] = s_s  # Suction values below psi are considered fully saturated
    return sat


# In[18]:


def vanGenuchten(sat, suction_value, s_r, s_s, a, n):
    m = 1 - 1/n
    sat = s_r + (s_s - s_r) / pow((np.log(np.exp + np.pow((suction_value/a),n)),m))
    return sat


# In[48]:


initial_guess = [0.25, 0, 0.387, 0.1]

# Fit the model to your data
params_opt, params_cov = curve_fit(brooks_and_corey, suction, vol_water, p0=initial_guess, maxfev=5000)

# Use the optimized parameters to calculate the fitted values
deg_sat_fitted = brooks_and_corey(suction, *params_opt)

# Calculate R-squared to assess the fit
r_squared = r2_score(deg_sat, deg_sat_fitted)
print(f'Optimized parameters: {params_opt}')
print(f'R-squared: {r_squared}')


# In[27]:


# Calculate fitted values using the optimized parameters
deg_sat_fitted = brooks_and_corey(suction, *params_opt)

# Plot original data
plt.figure(figsize=(8, 6))
plt.scatter(suction, deg_sat, color='blue', label='Original Data')

# Plot fitted curve
# Since 'suction' values were used directly for fitting, ensure they are sorted for plotting
sorted_indices = np.argsort(suction)
plt.plot(suction[sorted_indices], deg_sat_fitted[sorted_indices], color='red', label='Fitted Curve')

# Enhance plot
plt.xscale('log')  # Set x-axis to logarithmic scale
plt.xlabel('Matric Suction (Ψ), kPa')
plt.ylabel('Degree of Saturation (S), %')
plt.title('Soil Water Characteristic Curve (SWCC)')
plt.legend()
plt.grid(True)

# Show plot
plt.show()


# In[ ]:




