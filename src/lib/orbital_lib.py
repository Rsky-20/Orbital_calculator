#!/bin/python3
# -*- coding: utf8 -*-

"""
@Author : Pierre VAUDRY

Release date: 24/12/2022

[Description]
	Orbital lib is a asset tool of function to help and solve for orbital analyse 

[Class]
	{class of project}

[Function]
	{fonction of project}

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

#############################################
# ------ Process / Function section ------- #
#############################################

def e(a, b):
    """Formulas of eccentricity

    Args:
        a : demi-grand axe
        b : demi-petit axe

    Returns:
        float: eccentricity
    """    
    return math.sqrt((a ** 2) - (b ** 2))/ a

def a(rp:float, ra:float):
    """Function to solve demi-petit axe

    Args:
        rp (float): _description_
        ra (float): _description_

    Returns:
        _type_: _description_
    """
    return (rp + ra) / 2

def b(a, e):
    """Function to solve demi-grand axe

    Args:
        a (float): _description_
        e (float): _description_

    Returns:
        _type_: _description_
    """
    return math.sqrt((a ** 2) * (1 - e ** 2))

def eccentricity(R_a, R_p):
    """Function to solve eccentricity

    Args:
        R_a (_type_): _description_
        R_p (_type_): _description_

    Returns:
        _type_: _description_
    """
    return (R_a - R_p) / (R_a + R_p)

def Ra(e, a):
    """Function to solve the radius at apoge

    Args:
        e (_type_): _description_
        a (_type_): _description_

    Returns:
        _type_: _description_
    """
    return a * (1 + e)

def Rp(e, a ):
    """Function to solve the radius at perige

    Args:
        e (_type_): _description_
        a (_type_): _description_

    Returns:
        _type_: _description_
    """
    return a * (1 - e)

def Za(Ra, R):
    """Function to solve apoge altitude % the earth ground

    Args:
        Ra (_type_): _description_
        R (_type_): _description_

    Returns:
        _type_: _description_
    """
    return Ra - R

def Zp(Rp, R):
    """Function to solve perige altitude % the earth ground

    Args:
        Rp (_type_): _description_
        R (_type_): _description_

    Returns:
        _type_: _description_
    """
    return Rp - R

def T(a, MuT):
    """Function to solve the orbital period

    Args:
        a (_type_): _description_
        MuT (_type_): _description_

    Returns:
        _type_: _description_
    """
    return 2 * math.pi * math.sqrt(math.pow(a, 3) / MuT)

def n(MuT, a_):
    """Function to solve the mean motion (jsp ce que c'est)

    Args:
        MuT (_type_): _description_
        a_ (_type_): _description_

    Returns:
        _type_: _description_
    """
    return math.sqrt(MuT / math.pow(a_, 2))

def Vc(e):
    """Function to solve the critical true anomaly

    Args:
        e (_type_): _description_

    Returns:
        _type_: _description_
    """
    return np.arccos(-e)

def t_tp(a, MuT,  e, v):
    """_summary_

    Args:
        a (_type_): _description_
        MuT (_type_): _description_
        e (_type_): _description_
        v (_type_): _description_

    Returns:
        _type_: _description_
    """
    return math.sqrt(math.pow(a, 3) / MuT) * ((math.asin((math.sqrt(1 - math.pow(e, 2)) * math.sin(v)) / (1 + (e * math.cos(v))))) - e * ((math.sqrt(1 - math.pow(e, 2)) * math.sin(v)) / (1 + e * math.cos(v))))

def t(a, MuT, correction, correction_signe, e, v, tp):
    """Function to solve the passing times

    Args:
        a (_type_): _description_
        MuT (_type_): _description_
        correction (_type_): _description_
        correction_signe (_type_): _description_
        e (_type_): _description_
        v (_type_): _description_
        tp (_type_): _description_

    Returns:
        _type_: _description_
    """
    return math.sqrt(math.pow(a, 3) / MuT) * (correction + correction_signe * (math.asin((math.sqrt(1 - math.pow(e, 2)) * math.sin(v)) / (1 + (e * math.cos(v))))) - e * ((math.sqrt(1 - math.pow(e, 2)) * math.sin(v)) / (1 + e * math.cos(v)))) + tp

def Tp(a, MuT, correction, correction_signe, e, v):
    """Function to solve the passing time to the periapsis

    Args:
        a (_type_): _description_
        MuT (_type_): _description_
        correction (_type_): _description_
        correction_signe (_type_): _description_
        e (_type_): _description_
        v (_type_): _description_

    Returns:
        _type_: _description_
    """
    return - math.sqrt(math.pow(a, 3) / MuT) * (correction + correction_signe * (math.asin((math.sqrt(1 - math.pow(e, 2)) * math.sin(v)) / (1 + (e * math.cos(v))))) - e * ((math.sqrt(1 - math.pow(e, 2)) * math.sin(v)) / (1 + e * math.cos(v))))

def L0(correction, la, i):
    """Function to solve the Longitude

    Args:
        correction (_type_): _description_
        la (_type_): _description_
        i (_type_): _description_

    Returns:
        _type_: _description_
    """
    return correction + math.asin(math.tan(la) / math.tan(i))

def L(Lohm, L0, alpha, t):
    """Function to solve the Longitude with earth rotation

    Args:
        Lohm (_type_): _description_
        L0 (_type_): _description_
        alpha (_type_): _description_
        t (_type_): _description_

    Returns:
        _type_: _description_
    """
    return Lohm + L0 - alpha * t

def la(w, i, v):
    """Function to solve the Latitude

    Args:
        w (_type_): _description_
        i (_type_): _description_
        v (_type_): _description_

    Returns:
        _type_: _description_
    """
    return math.asin(math.sin( w + v) * math.sin(i))

def V(MuT, a, r):
    """Function to solve the speed 

    Args:
        MuT (_type_): _description_
        a (_type_): _description_
        r (_type_): _description_

    Returns:
        _type_: _description_
    """
    return math.sqrt(2 * (-(MuT / (2 * a)) + (MuT / r)))

def r(a, e, v):
    """Function to solve 

    Args:
        a (_type_): _description_
        e (_type_): _description_
        v (_type_): _description_

    Returns:
        _type_: _description_
    """
    return (a * ( 1 - e ** 2)) / (1 + e * math.cos(v))

if __name__ == '__main__':
    # Run & Test the program
    e_ = 0.8320
    a_ = 40708
    i_ = 61
    Rt = 6378
    MuT = 398600
    la_ = -61
    rp = 6708
    ra = 6790
    eccentricity_ = eccentricity(ra, rp)
    a_axis = a(rp, ra)
    b_axis = b(a_axis, eccentricity_)
    
    
    """Ra_ = Ra(e_, a_)
    print(Ra_)
    Rp_ = Rp(e_, a_)
    print(Rp_)
    Za_ = Za(Ra_, Rt)
    print(Za_)
    Zp_ = Zp(Rp_, Rt)
    print(Zp_)
    T_ = T(a_, MuT)
    print(T_)
    n_ = n(MuT, a_)
    print(n_)
    tp_ = Tp(a_, -2 * math.pi, e_, math.radians(-270))
    print(tp_)
    
    
    L0_ = L0(90, la_, i_)"""
    
    VA = V(MuT,a_, r(a_,e_,Vc(e_)))
    print(VA)
    #print(L0_)
    la_ = la(math.radians(270), math.radians(61), math.radians(-180))
    print(math.degrees(la_))
    L0_ = L0(math.radians(-180), (la_), math.radians(i_))
    print(math.degrees(L0_))
    

    

    
    