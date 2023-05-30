import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import os
import pandas as pd
from scipy import interpolate

#Paramters - convert everything to SI units because this is way less confusing
R = 8.314 # J/(mole K) - universal gas constant
Mw = 18.015/1000.0 # kg/mole - molar mass of water, Weast et al. 1987
Ma = 28.96521/1000 # kg/mol - mass of dry air
Na = 6.0221414e23 # Avogadro's number, molecules/mole

def svpice(temp):
    # saturation vapor pressure over ice, after Murphy and Koop (2005), in Pa
    # T>110 K
    # temp - K, vp_ice - hPa
    a0 = 9.550426
    a1 = 5723.265
    a2 = 3.53068
    a3 = 0.00728332

    vpice=np.exp(a0-a1/temp+a2*np.log(temp)-a3*temp)

    return vpice
def Li(temp):
    ## J/kg (190 < T < 273), Beard and Pruppacher, 1971
    Li0 = 2.836*1e6 # K/kg

    ci = 2106 # J/kg/K
    cpv = 1885 # J/kg/K

    Li = Li0 - (ci-cpv)*(temp-273.15)

    return Li
def densityice(temp):
    # density of ice in g/cm3
    # Pruppacher and Klett, 1997
    p0=0.9167
    p1=1.75e-4
    p2=5.0e-7

    rho = p0-p1*(temp-273.15)-p2*(temp-273.15)**2

    return rho*1000 # convert to kg/m3
def densitydryair(temp,press):
    # Temp - K, Press - hPa
    # density of dry air, g/cm3
    # from the ideal gas law
    Rs = R/(Ma/1000.0) # J/(mole K)/(kg/mole) -> specific gas constant, J/(kg K)
    rho = (press*100)/(temp*Rs)*0.001 # g/cm3

    return rho*1000 # convert to kg/m3
def ka(temp):
    #heat conductivity of dry air, Beard and Pruppacher (1971)
    #from Max's paper - originally from P&K, 1978

    k0 = 4.3783e-3
    k1 = 7.1128e-5

    heatconductivity=k0+k1*temp # J/(m s K)

    return heatconductivity
def diffusivity(temp,press):
    # diffusivity for water molecules in air, Hall and Pruppacher, 1976
    Dw0=0.211   # cm2/s , Hall and Pruppacher, 1976
    p0=1013.25  #standard pressure, hPa
    T0=273.15   #standard temperature, K
    gam=1.94 #temperature coefficient

    Dw=Dw0*(p0/press)*(temp/T0)**(gam) # cm2/s

    return Dw*1e-4 # m2/s
def icemasspart(h2o_ppmv_vapor,h2o_ppmv_total,cn_welas_ice,tg,press):
    # average particle mass (g/particle)
    h2o_ppmv_ice=(h2o_ppmv_total-h2o_ppmv_vapor)/cn_welas_ice
    rhoa = densitydryair(tg,press) # kg/(m3)

    icedens=rhoa*h2o_ppmv_ice/(1e6)*Mw/Ma

    return icedens/1000.0 #kg/particle
