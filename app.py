#!/bin/python3
# -*- coding: utf8 -*-

"""
@Author : Pierre VAUDRY

Release date: 24/12/2022

[Description]
	Orbital calculator is a calculator script to analyze an orbit and plot the ground track. 

[Class]
	{class of project}

[Function]
	See all fonction in ./src/lib/orbital_lib.py

[Other variable]:
	{variable}
"""
#############################################
# --------- Import module section --------- #
#############################################
#from math import abs, acos, sqrt, pow
import math
import numpy as np
import pandas as pd
from src.lib.orbital_lib import *



#############################################
# -------- Variable Global section -------- #
#############################################

#"v °": ['t (sec)', 'la', 'L0 °', 'L °']
v_degree = [v for v in range(-210, 300, 30)]
time_parsing, latitude, Longitude, Longitude_corr = [], [], [], []


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

omega = 270 #
e_ = 0.8320 #eccentricite e
a_ = 40708 #demi-petit axe a
inclinaison = 61 #inclinaison i
Rt = 6378 #Rayon Terre
MuT = 398600 #
la_ = -61 #
rp = 6708 #rayon periastre
ra = 6790 #rayon apoastre
v = -omega
vc_ = Vc(e_)
Lohm = 120
alpha = (360 / 86164)

condition_correction = [
    [-4 * math.pi + vc_, -2 * math.pi - vc_],
    [-2 * math.pi - vc_, -2 * math.pi + vc_],
    [-2 * math.pi + vc_, -vc_],
    [-vc_, vc_],
    [vc_, 2 * math.pi - vc_],
    [2 * math.pi - vc_, 2 * math.pi + vc_],
    [-4 * math.pi + vc_, -2 * math.pi - vc_]
                        ]
correction = [
    [[-3 * math.pi, -1],
        [-2 * math.pi, 1],
        [-math.pi, -1],
        [0,1],
        [math.pi, -1],
        [2 * math.pi, 1],
        [3 * math.pi, -1]],
    [[3 * math.pi, -1],
        [2 * math.pi, 1],
        [math.pi, -1],
        [0,1],
        [-math.pi, -1],
        [-2 * math.pi, 1],
        [-3 * math.pi, -1]],
                ]

correction_L0 = [
    [-omega - 90 - 2 * 180, -omega - 90 - 180, 2 * 180],
    [-omega - 90 - 180, -omega - 90, -180],
    [-omega - 90, -omega + 90, 0],
    [-omega + 90, -omega + 90 + 180, 180],
    [-omega + 90 + 180, -omega + 90 + 2 * 180, 2 * 180]
                ]

#############################################
# ------ Process / Function section ------- #
#############################################

if __name__ == '__main__':
    # Run the program
    # si pas donnee
    #eccentricity_ = eccentricity(ra, rp) 
    #a_axis = a(rp, ra)
    #b_axis = b(a_axis, eccentricity_)
    
    texte = """

Apoapsis altitude : Za = {} km
Periapsis altitude : Zb = {} km
Apoapsis radius : ra = {} km
Periapsis radius : rb = {} km
Orbital period (sec) : T = {} s
Orbital period (sexagesimal form) : T = {} 
Mean motion (mouvement moyen) : n = {} rd/s
Mean motion (mouvement moyen) : n = {} rd/day

    
    """
    Ra_ = Ra(e_, a_)
    Rp_ = Rp(e_, a_)
    Za_ = Za(Ra_, Rt)
    Zp_ = Zp(Rp_, Rt)
    T_ = T(a_, MuT)
    T_sexagesimal = T_ % 60
    rest = T_ // 60
    n_ = n(MuT, a_)
        
    print(texte.format(Za_, Zp_, Ra_, Rp_, T_, T_sexagesimal, n_, n_ * 3600))
                     
    if inclinaison < 90: #
        prograde = correction[0]
        for corr_1 in range(len(condition_correction)): 
            if condition_correction[corr_1][0] < math.radians(v) < condition_correction[corr_1][1]: 
                tp_ = Tp(a_, MuT, prograde[corr_1][0], prograde[corr_1][1], e_, math.radians(v)) 
                for v_iter in range(len(v_degree)): 
                    la_ = la(math.radians(omega), math.radians(inclinaison), math.radians(v_degree[v_iter])) 
                    latitude.append(round(math.degrees(la_), 2)) 
                    for corr_2 in range(len(condition_correction)):
                        if condition_correction[corr_2][0] < math.radians(v_degree[v_iter]) < condition_correction[corr_2][1]:
                            t_ = t(a_, MuT , prograde[corr_2][0], prograde[corr_2][1], e_, math.radians(v_degree[v_iter]), tp_)
                            time_parsing.append(round(t_, 2))
                    for L0_iter in range(len(correction_L0)):
                        if correction_L0[L0_iter][0] <= v_degree[v_iter] <= correction_L0[L0_iter][1]:
                            corr_11 = correction_L0[L0_iter][2]
                    L0_ = L0(math.radians(corr_11), math.radians(latitude[v_iter]), math.radians(inclinaison)) 
                    Longitude.append(round(math.degrees(L0_), 2))     
                    Long = L(Lohm, Longitude[v_iter], alpha, time_parsing[v_iter])    
                    Longitude_corr.append(round(Long, 2))      
        
    elif inclinaison > 90:
        retrograde = correction[1]
        for corr_1 in range(len(condition_correction)):
            if condition_correction[corr_1][0] < math.radians(v) < condition_correction[corr_1][1]:
                tp_ = Tp(a_, MuT, retrograde[corr_1][0], retrograde[corr_1][1], e_, math.radians(v)) 
                for v_iter in range(len(v_degree)):
                    la_ = la(math.radians(omega), math.radians(inclinaison), math.radians(v_degree[v_iter]))
                    latitude.append(round(math.degrees(la_), 2))
                    for corr_2 in range(len(condition_correction)):
                        if condition_correction[corr_2][0] < math.radians(v_degree[v_iter]) < condition_correction[corr_2][1]:
                            t_ = t(a_, MuT , retrograde[corr_2][0], retrograde[corr_2][1], e_, math.radians(v_degree[v_iter]), tp_)
                            time_parsing.append(round(t_, 2))
                    for L0_iter in range(len(correction_L0)):
                        if correction_L0[L0_iter][0] <= v_degree[v_iter] <= correction_L0[L0_iter][1]:
                            corr_11 = correction_L0[L0_iter][2]
                    L0_ = L0(corr_11, math.radians(latitude[v_iter]), math.radians(inclinaison)) 
                    Longitude.append(round(math.degrees(L0_), 2))  
                    Long = L(Lohm, L0_, alpha, time_parsing[v_iter])
                    Longitude_corr.append(round(Long, 2))
                            
    
    elif inclinaison == 90:
        pass
    
    print('################################################')
    
    data = {'v °': v_degree,'t (sec)': time_parsing, 'la': latitude, 'L0 °': Longitude , 'L °': Longitude_corr}
    df = pd.DataFrame(data) #, columns = v_degree
    #print(data)
    print(df)
    
    """import matplotlib.pyplot as plt
    plt.scatter(df['L °'], df['la'])
    plt.show()"""
    """img = plt.imread("./src/data/GroundTrack.jpg")
    fig, ax = plt.subplots(Longitude_corr, latitude)
    ax.imshow(img)"""
