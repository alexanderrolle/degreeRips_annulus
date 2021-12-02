#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 00:18:11 2021

@author: arolle
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from nu import nu

# computes the derivative of nu in the case 
# that the s-ball about (c,0) intersects both 
# C(O, R) and C(O, Q)
def nu_diff(c, r1, r2, s, sigma, tau):
    
    x1 = ((c ** 2) + (r1 ** 2) - (s ** 2)) / (2 * c)
    x2 = ((c ** 2) + (r2 ** 2) - (s ** 2)) / (2 * c)
    
    y1 = np.sqrt((r1 ** 2) - (x1 ** 2))
    y2 = np.sqrt((r2 ** 2) - (x2 ** 2))
    
    return -2 * ((sigma - tau) * y1 + tau * y2)

TOL = 0.00001

############################################################
## The case w > 0 ##########################################
############################################################

# choose inner radius r1, outer radius 2, and w
r1 = 0.4
r2 = 0.5
w = 25 / 500

sigma = w / (np.pi * (r1 ** 2))
tau = (1 - w) / ((np.pi * (r2 ** 2)) - (np.pi * (r1 ** 2)))

num_evaluation_points = 2000
X = np.linspace(0, 0.5, num_evaluation_points)

omega = np.zeros(shape=num_evaluation_points)
for i in range(num_evaluation_points):
    
    s = X[i]
    
    if s <= 0.5 * (r2 - r1):
        omega[i] = r1 + s
        
    if s > 0.5 * (r2 - r1) and s < 0.5 * (r2 + r1):
        
        radius_sq = 0.5 * ((r1 ** 2) + (r2 ** 2))
        guess = np.sqrt(radius_sq - (s ** 2))

        x = fsolve(func=nu_diff, x0=guess, args=(r1, r2, s, sigma, tau))
        
        omega[i] = x[0]
        
    if s >= 0.5 * (r2 + r1):
    
        omega[i] = 0
        
phi_0_s_vals = []
phi_0_vals = []
        
phi_1_s_vals = []
phi_1_vals = []

phi_1_fill_s_vals = []
phi_1_fill_vals_above = []
phi_1_fill_vals_below = []

phi_2_s_vals = []
phi_2_vals = []

phi_2_fill_s_vals = []
phi_2_fill_vals_above = []
phi_2_fill_vals_below = []

phi_infty_s_vals = []
phi_infty_vals = []

for i in range(num_evaluation_points):
    s = X[i]
    
    rho_1 = s / np.sqrt(3)
    rho_2 = s / (2 * np.sin(2 * np.pi / 5))
    rho_infty = s / 2
    
    if rho_1 < omega[i]:
        phi_1_val = nu(w, r1, r2, s, rho_1)
        
        # rho_1 < omega, so we should plot phi_0
        phi_0_val = nu(w, r1, r2, s, omega[i])
        phi_0_s_vals.append(s)
        phi_0_vals.append(phi_0_val)
        
        phi_1_fill_s_vals.append(s)
        phi_1_fill_vals_above.append(phi_0_val)
        phi_1_fill_vals_below.append(phi_1_val)
    
    if rho_1 >= omega[i]:
        phi_1_val = nu(w, r1, r2, s, omega[i])
        
    if rho_2 < omega[i]:
        phi_2_val = nu(w, r1, r2, s, rho_2)
    
    if rho_2 >= omega[i]:
        phi_2_val = nu(w, r1, r2, s, omega[i])
        
    if rho_infty < omega[i]:
        phi_infty_val = nu(w, r1, r2, s, rho_infty)
        
    if rho_infty >= omega[i]:
        phi_infty_val = nu(w, r1, r2, s, omega[i])
        
    if np.abs(phi_1_val - phi_2_val) > TOL:
        
        phi_1_s_vals.append(s)
        phi_1_vals.append(phi_1_val)
        
        phi_2_fill_s_vals.append(s)
        phi_2_fill_vals_above.append(phi_1_val)
        phi_2_fill_vals_below.append(phi_2_val)
        
    if np.abs(phi_2_val - phi_infty_val) > TOL:
        
        phi_2_s_vals.append(s)
        phi_2_vals.append(phi_2_val)
        
    phi_infty_s_vals.append(s)
    phi_infty_vals.append(phi_infty_val)
    
figsize = (10,10)
plt.figure(figsize=figsize)

plt.xlim(0, 0.5)
plt.ylim(0, 0.5)
plt.xticks([i * 0.05 for i in range(11)])
plt.yticks([i * 0.05 for i in range(11)])

plt.xlabel('Rips parameter')
plt.ylabel('degree parameter')

width=2
plt.plot(phi_0_s_vals, phi_0_vals, color='r', linewidth=width)
plt.plot(phi_1_s_vals, phi_1_vals, color='b', linewidth=width)
plt.plot(phi_2_s_vals, phi_2_vals, color='#FFA500', linewidth=width)
plt.plot(phi_infty_s_vals, phi_infty_vals, color='k', linewidth=width)

plt.fill_between(phi_1_fill_s_vals, phi_1_fill_vals_above, phi_1_fill_vals_below, color='b', alpha=0.15)
plt.fill_between(phi_2_fill_s_vals, phi_2_fill_vals_above, phi_2_fill_vals_below, color='#FFA500', alpha=0.15)

plt.show()


############################################################
## The case w = 0 ##########################################
############################################################

# when w=0 there is a simple formula for omega
def omega_w0(r1, r2, s):
    
    if 0 <= s and s <= 0.5 * (r2 - r1):
        return r1 + s
    
    if s > 0.5 * (r2 - r1) and s < 0.5 * (r2 + r1):
        
        radius_sq = 0.5 * ((r1 ** 2) + (r2 ** 2))
        c = np.sqrt(radius_sq - (s ** 2))
        
        if c >= r1:
            return c
        if c < r1:
            return r1
        
    if s >= 0.5 * (r2 + r1):
        return r1

TOL = 0.00001

r1 = 0.4
r2 = 0.5
w = 0

sigma = 0
tau = 1 / ((np.pi * (r2 ** 2)) - (np.pi * (r1 ** 2)))

num_evaluation_points = 4000
X = np.linspace(0, 1.0, num_evaluation_points)

omega = np.zeros(shape=num_evaluation_points)
for i in range(num_evaluation_points):
    s = X[i]
    omega[i] = omega_w0(r1, r2, s)

phi_0_s_vals = []
phi_0_vals = []
        
phi_1_s_vals = []
phi_1_vals = []

phi_1_fill_s_vals = []
phi_1_fill_vals_above = []
phi_1_fill_vals_below = []

phi_2_s_vals = []
phi_2_vals = []

phi_2_fill_s_vals = []
phi_2_fill_vals_above = []
phi_2_fill_vals_below = []

phi_3_s_vals = []
phi_3_vals = []

phi_3_fill_s_vals = []
phi_3_fill_vals_above = []
phi_3_fill_vals_below = []

phi_rest_fill_s_vals = []
phi_rest_fill_vals_above = []
phi_rest_fill_vals_below = []

phi_infty_s_vals_small = []
phi_infty_vals_small = []

phi_infty_s_vals_large = []
phi_infty_vals_large = []

for i in range(num_evaluation_points):
    s = X[i]
    
    phi_0_val = nu(w, r1, r2, s, omega[i])
    rho_1 = s / np.sqrt(3)
    rho_2 = s / (2 * np.sin(2 * np.pi / 5))
    rho_3 = s / (2 * np.sin(3 * np.pi / 7))
    rho_infty = s / 2
    
    # calculate phi_1
    if rho_1 <= r1:
        phi_1_val = 0
    
    if rho_1 > r1 and rho_1 < omega[i]:
        phi_1_val = nu(w, r1, r2, s, rho_1)
        
    if rho_1 > r1 and rho_1 >= omega[i]:
        phi_1_val = nu(w, r1, r2, s, omega[i])
    
    # calculate phi_2
    if rho_2 <= r1:
        phi_2_val = 0
    
    if rho_2 > r1 and rho_2 < omega[i]:
        phi_2_val = nu(w, r1, r2, s, rho_2)
        
    if rho_2 > r1 and rho_2 >= omega[i]:
        phi_2_val = nu(w, r1, r2, s, omega[i])
        
    # calculate phi_3
    if rho_3 <= r1:
        phi_3_val = 0
    
    if rho_3 > r1 and rho_3 < omega[i]:
        phi_3_val = nu(w, r1, r2, s, rho_3)
        
    if rho_3 > r1 and rho_3 >= omega[i]:
        phi_3_val = nu(w, r1, r2, s, omega[i])
        
    # calculate phi_infty
    if rho_infty <= r1:
        phi_infty_val = 0
    
    if rho_infty > r1 and rho_infty < omega[i]:
        phi_infty_val = nu(w, r1, r2, s, rho_infty)
        
    if rho_infty > r1 and rho_infty >= omega[i]:
        phi_infty_val = nu(w, r1, r2, s, omega[i])
    
    if np.abs(phi_0_val - phi_1_val) > TOL:
        
        phi_0_s_vals.append(s)
        phi_0_vals.append(phi_0_val)
        
        phi_1_fill_s_vals.append(s)
        phi_1_fill_vals_above.append(phi_0_val)
        phi_1_fill_vals_below.append(phi_1_val)
    
    if np.abs(phi_1_val - phi_2_val) > TOL:
        
        phi_1_s_vals.append(s)
        phi_1_vals.append(phi_1_val)
        
        phi_2_fill_s_vals.append(s)
        phi_2_fill_vals_above.append(phi_1_val)
        phi_2_fill_vals_below.append(phi_2_val)
        
    if np.abs(phi_2_val - phi_3_val) > TOL:
        
        phi_2_s_vals.append(s)
        phi_2_vals.append(phi_2_val)
        
        phi_3_fill_s_vals.append(s)
        phi_3_fill_vals_above.append(phi_2_val)
        phi_3_fill_vals_below.append(phi_3_val)
        
    if np.abs(phi_3_val - phi_infty_val) > TOL:
        
        phi_3_s_vals.append(s)
        phi_3_vals.append(phi_3_val)
        
        phi_rest_fill_s_vals.append(s)
        phi_rest_fill_vals_above.append(phi_3_val)
        phi_rest_fill_vals_below.append(phi_infty_val)
    
    if rho_infty <= r1:
    
        phi_infty_s_vals_small.append(s)
        phi_infty_vals_small.append(phi_infty_val)
        
    if rho_infty > r1:
    
        phi_infty_s_vals_large.append(s)
        phi_infty_vals_large.append(phi_infty_val)
    
figsize = (10,10)
plt.figure(figsize=figsize)

plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xticks([i * 0.1 for i in range(11)])
plt.yticks([i * 0.1 for i in range(11)])

plt.xlabel('Rips parameter')
plt.ylabel('degree parameter')

width=2
plt.plot(phi_0_s_vals, phi_0_vals, color='r', linewidth=width)
plt.plot(phi_1_s_vals, phi_1_vals, color='b', linewidth=width)
plt.plot(phi_2_s_vals, phi_2_vals, color='#FFA500', linewidth=width)
plt.plot(phi_3_s_vals, phi_3_vals, color=(0.8,0.8,0.8), linewidth=width)
plt.plot(phi_infty_s_vals_small, phi_infty_vals_small, color='k', linewidth=width)
plt.plot(phi_infty_s_vals_large, phi_infty_vals_large, color='k', linewidth=width)

plt.fill_between(phi_1_fill_s_vals, phi_1_fill_vals_above, phi_1_fill_vals_below, color='b', alpha=0.15)
plt.fill_between(phi_2_fill_s_vals, phi_2_fill_vals_above, phi_2_fill_vals_below, color='#FFA500', alpha=0.15)
plt.fill_between(phi_3_fill_s_vals, phi_3_fill_vals_above, phi_3_fill_vals_below, color='#FFDAB9', alpha=0.25)
plt.fill_between(phi_rest_fill_s_vals, phi_rest_fill_vals_above, phi_rest_fill_vals_below, color=(0.8,0.8,0.8), alpha=0.25)

plt.show()