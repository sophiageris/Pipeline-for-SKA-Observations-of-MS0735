#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 23:24:43 2023

@author: sophiageris
"""

import astropy
from timeit import default_timer as timer
from matplotlib import pyplot
import scipy
import numpy as np
from scipy.integrate import quad
import math
from scipy import integrate
from scipy.special import kn
import scipy.special as sc
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

start = timer()
# constants
gamma = 0.3081 # shape parameter
alpha = 1.0510 # shape parameter
beta = 5.4905 # shape parameter
sigma_t = 6.98e-74 # Thomson cross section in Mpc^2
c = 9.72e-15 # speed of light in Mpc s^-1
m_e = 9.10938356e-31 # electron mass in kg
m_ec2 = 511 # keV
h = 4.13E-15 # eV/Hz
k = 8.617e-5 # eV/K
T_cmb = 2.7255 # K
H_0 = 70 # km s^-1 Mpc^-1
omega_lambda_0 = 0.7
c_km = 299792
radian = 3437.75 # 1radian = 3437.75arcmin
c_500 = 1.177
#j = 16020000 # conversion factor from joules to Janskys
c_m = 3e8 # speed of light in m/s

s = np.linspace(-10, 9, 1000)
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

#defining pixel grid
def create_pixgrid(num_pix, pix_size):
    ex = np.arange((-num_pix/2)*pix_size, ((((num_pix/2))*pix_size)), pix_size) # pix width in arcmin
    #print(ex)
    y = np.arange((-num_pix/2)*pix_size, ((((num_pix/2))*pix_size)), pix_size) # pix height in armin
    x_1, y_1 = np.meshgrid(ex, y)
    theta_proj = ((((x_1)**2)+((y_1)**2))**0.5) # projected radius of each pixel on the sky in arcmin
    return x_1, y_1, theta_proj, ex, y

def x(nu):
    return (h*nu)/(k*T_cmb)

def D_A(z): 
    #d_phys = (c_km*z)/(H_0*np.sqrt(omega_lambda_0))
    d_a = (cosmo.angular_diameter_distance(z).value) # in Mpc
    return (d_a/radian) # D_A in Mpc/arcmin

def rho_crit(z): 
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    rho = cosmo.critical_density(z) # rho_critical in g/cm^3
    return rho.to(u.solMass / u.Mpc**3) # rho_critical in solar masses/ Mpc^3

def r_500(z, M_500): # in Mpc
    return np.cbrt((3*M_500)/(4*np.pi*500*(rho_crit(z).value))) # M_500 in solMass and rho_crit in solMass/Mpc^3 therefore r_500 is in Mpc

def r_s(z, M_500): #in Mpc
    return r_500(z, M_500)/c_500

def theta_s(z, M_500):
    return r_s(z, M_500)/ D_A(z) # in arcmin

def unnorm_GNFW(r, r_s):
    x = r/r_s
    return (x**-gamma)*(1+x**alpha)**((gamma-beta)/alpha)

def inte(r, z, M_500):
    inte = (r**2)*unnorm_GNFW(r, r_s(z, M_500))
    return inte

def inte_total(z, M_500):
    return quad(inte, 0, np.inf, args=(z, M_500))[0]

def inte_500(z, M_500):
    return quad(inte, 0, r_500(z, M_500), args=(z, M_500))[0]

Y_int_ana = scipy.special.gamma((3-gamma)/alpha)*scipy.special.gamma((beta-3)/alpha)/alpha/scipy.special.gamma((beta-gamma)/alpha)*r_s(0.5, 3.25e14)**3
Y_int_num = inte_total(0.5, 3.25e14)
#print(Y_int_ana)
#print(Y_int_num)

def y_500_num(z, M500): #in arcmin^2 (scaling relation from Planck Collaboration: https://arxiv.org/pdf/1303.5080.pdf)
    Ez=cosmo.efunc(z)
    DA=cosmo.angular_diameter_distance(z).value/radian # in Mpc/arcmin
    DA2Y500=10**(-0.175)*(M500/6e14)**1.77 * Ez**(2./3) * 1e-4
    return DA2Y500/DA**2

def y_tot(z, M_500): # in arcmin^2
    return ((inte_total(z, M_500))/(inte_500(z, M_500)))*(y_500_num(z, M_500))

def beta_th(kT_e):
    return (m_ec2)/kT_e # unitless

def g_x(nu):
    return (((x(nu)**4)*np.exp(x(nu)))/((np.exp(x(nu))-1)**2))*((x(nu)*((np.exp(x(nu))+1)/(np.exp(x(nu))-1)))-4) # unitless

def Pressure(l, ytot, thetas, DA, theta_proj):
    P_norm = ytot/(4*np.pi*DA*thetas**3)*alpha*scipy.special.gamma((beta-gamma)/alpha)/scipy.special.gamma((3-gamma)/alpha)/scipy.special.gamma((beta-3)/alpha) # in Mpc^-1
    theta_los = l/DA # arcmin (when the integral is through the bubbles then the z and -z intercepts need to be in Mpc which is why in the bounds function below I multiplied them by da, to convert from arcmin to Mpc so that theta_los ends up being in arcmin as it should)
    theta = np.sqrt(theta_proj**2 + theta_los**2) # arcmin
    return P_norm*unnorm_GNFW(theta, thetas) # Mpc^-1

def Bounds(rb, xb, yb, zb, num_pix, pix_size, z):
    da = D_A(z)
    Z_plus = []
    Z_minus = []
    x_1, y_1, theta_proj, ex, y = create_pixgrid(num_pix, pix_size)
    y_los1 = np.zeros((len(ex),len(y)))
    for i in range(len(ex)):
        for j in range(len(y)):
            zed_plus = (np.sqrt((rb**2)-((x_1[i][j]-xb)**2)-((y_1[i][j]-yb)**2)) + zb)*da # in Mpc
            Z_plus.append(zed_plus)
    z_plus = np.reshape(np.array(Z_plus),y_los1.shape)
            
    for i in range(len(ex)):
        for j in range(len(y)):
            zed_minus = -(np.sqrt((rb**2)-((x_1[i][j]-xb)**2)-((y_1[i][j]-yb)**2)) + zb)*da # in Mpc
            Z_minus.append(zed_minus)
    z_minus = np.reshape(np.array(Z_minus),y_los1.shape)
    
    return z_plus, z_minus

'''def Bounds2(rb, x_1, xb2, y_1, yb2, zb2, start, stop, pix_size):
    Z2 = []
    x_1, y_1, theta_proj, ex, y = create_pixgrid(start, stop, pix_size)
    for i in range(len(ex)):
        for j in range(len(y)):
            z2 = np.sqrt((rb**2)-((x_1[i][j]-xb2)**2)-((y_1[i][j]-yb2)**2)) + zb2
            Z2.append(z2)
    return Z2'''

def smax(p):
    return 2*np.arcsinh(p)

def P_t_p(t, p):
    f1 = 3*np.abs(1-t)/(32*(p**6)*t)
    brac1 = 1+(10+(8*p**2)+(4*p**4))*t+t**2
    f2 = 3*(1+t)/(8*p**5)
    brac2 = (3+3*p**2+p**4)/(1+p**2)**0.5 - (3+2*p**2)/(2*p)*(2*np.arcsinh(p)-np.abs(np.log(t)))
        
    if np.abs(np.log(t))>smax(p):
        return 0
    else: 
        return -f1*brac1+f2*brac2
    
def blackbody(t,nu):
    return ((x(nu)/t)**3)/(np.exp(x(nu)/t)-1)
        
def thermaldist(p, kT_e):
    return (beta_th(kT_e)/kn(2, beta_th(kT_e)))*(p**2)*np.exp(-beta_th(kT_e)*((1+(p**2))**0.5))

def integrand(p, t, nu, kT_e):
    return P_t_p(t, p)*thermaldist(p, kT_e)*blackbody(t, nu)

def p_s(s, p):
    if np.abs(s) > smax(p):
        return 0
    else: 
        return ((((-1*(3*np.abs(1-np.exp(s)))/(32*(p**6)*np.exp(s)))*(1+(10+(8*(p**2))+(4*(p**4)))*(np.exp(s))+(np.exp(s)**2)))) + (((3*(1+np.exp(s)))/(8*(p**5)))*(((3+(3*(p**2))+(p**4))/((1+(p**2))**0.5))-(((3+(2*(p**2)))/(2*p))*(2*np.arcsinh(p)-np.abs(np.log(np.exp(s))))))))
  
def powerlaw(p, alpha2, p_1, p_2):
    if p<p_1 or p>p_2:
        return 0
    else:   
        return (alpha2-1)*p**-alpha2/(p_1**(1-alpha2) - p_2**(1-alpha2))

def integrand_nonthermal(p, s, alpha2, p_1, p_2, nu):
    if p<p_1 or p>p_2:
        return 0
    else:
        return p_s(s, p)*powerlaw(p, alpha2, p_1, p_2)*np.exp(s)*blackbody(np.exp(s), x(nu))

def ymap(nu, z, M_500, xb, yb, rb, num_pix, pix_size, lower, upper):
    x_1, y_1, theta_proj, ex, y = create_pixgrid(num_pix, pix_size)
    y_los = np.zeros((len(ex),len(y)))
    ytot =  y_tot(z, M_500)
    DA = D_A(z)
    
    for i in range(len(ex)):
        for j in range(len(y)):
            if not isinstance(lower, np.ndarray):
                I = quad(Pressure, lower, upper, args=(ytot, theta_s(z, M_500), DA, theta_proj[j,i]))[0]
                y_los[j][i] = I
            elif ((ex[i]-xb)**2)+((y[j]-yb)**2)<rb**2:
                I = quad(Pressure, lower[j,i], upper[j,i], args=(ytot, theta_s(z, M_500), DA, theta_proj[j,i]))[0]
                y_los[j][i] = I
            else:
                I = 0
                y_los[j][i] = I
    return y_los

def ymap_cluster(nu, z, M_500, xb, yb, rb, num_pix, pix_size): # ymap of the cluster whole cluster if there was no bubbles
    ymap_cl = ymap(nu, z, M_500, xb, yb, rb, num_pix, pix_size, -np.inf, np.inf)
    return ymap_cl


#ymap_cluster2 = ymap_cluster(14.11e9, 0.5, 3.25e14, 0.24, 0.67, 0.5, -total_range, total_range, 0.033333333)
#pyplot.imshow(ymap_cluster2)

def g_tilde_thermal(nu, kT_e):
    j_x_values = integrate.dblquad(integrand, 0, np.inf, 0, np.inf, args=(nu, kT_e))[0] 
    return ((j_x_values-blackbody(1, nu))*m_ec2)/(kT_e)

def my_betainc(a, b, x):
    if a>0 and b>0:
      return sc.betainc(a,b,x)*sc.beta(a,b)
    elif a>0:
      return (x**a)/a*sc.hyp2f1(a,1-b,a+1,x)
    else:
      return sc.beta(a,b)-((1-x)**b)*(x**a)/b*sc.hyp2f1(a,a+b,b+1,1-x)

def g_tilde_nonthermal(nu, p_1, p_2, alpha2):
    # integral over p
    inte = np.array([scipy.integrate.quad(integrand_nonthermal, p_1, p_2, args=(sj, alpha2, p_1, p_2, nu))[0] for sj in s])
    # integral over s
    j_x_values_nonthermal = scipy.integrate.trapz(inte, s)
    kT_top = m_ec2*(alpha2 - 1) 
    kT_bottom = 6*((p_1**(1-alpha2)) - (p_2**(1-alpha2))) 
    beta2 = (my_betainc((alpha2-2)/2, (3-alpha2)/2, 1/(1+(p_1**2)))) - (my_betainc((alpha2-2)/2, (3-alpha2)/2, 1/(1+(p_2**2))))
    kT_e_values = (kT_top/kT_bottom)*beta2
    return (j_x_values_nonthermal - blackbody(1, nu))*(m_ec2/(kT_e_values))

def f_thermal(nu, kT_e):
    return 1 - g_tilde_thermal(nu, kT_e)/g_x(nu)

def f_nonthermal(nu, p_1, p_2, alpha2):
    return 1 - g_tilde_nonthermal(nu, p_1, p_2, alpha2)/g_x(nu)

def ymap_bubbles(nu, z, M_500, xb, yb, zb, xb2, yb2, zb2, rb, supp_func, supp_args, num_pix, pix_size): #ymap of the bubbles ie. ymap is zero everywhere but in the bubble regions I'm calculating ymap by integrating over line of sight but when line of sight from -z to +z and am suppressing the signal by f (so leaving out the l.o.s parts -inf to -z and +z to +inf)
    Z_plus, Z_minus = Bounds(rb, xb, yb, zb, num_pix, pix_size, z)
    Z_plus2, Z_minus2 = Bounds(rb, xb2, yb2, zb2, num_pix, pix_size, z)
    #print(Z_plus, Z_minus)
    
    f = supp_func(*supp_args)
    print(f)

    bubble1 = f*ymap(nu, z, M_500, xb, yb, rb, num_pix, pix_size, Z_minus, Z_plus)
    bubble2 = f*ymap(nu, z, M_500, xb2, yb2, rb, num_pix, pix_size, Z_minus2, Z_plus2)
    bubbles = bubble1 + bubble2

    return bubbles

def ymap_cluster_bubbles(nu, z, M_500, xb, yb, zb, xb2, yb2, zb2, rb, supp_func, supp_args, num_pix, pix_size):
    y_bubbles = ymap_bubbles(nu, z, M_500, xb, yb, zb, xb2, yb2, zb2, rb, supp_func, supp_args, num_pix, pix_size)
    y_cluster = ymap_cluster(nu, z, M_500, xb, yb, rb, num_pix, pix_size)
    
    total_y = y_cluster - y_bubbles
    
    return total_y

y_map = ymap_cluster_bubbles(14.11e9, 0.216, 8e14, 0.24, 0.67, 0, -0.33, -0.82, 0, 0.5, f_thermal, [14.11e9, 1258.93 ], 512, 2/60)


def calculating_ymap_total(nu, z, M_500, xb, yb, zb, xb2, yb2, zb2, rb, supp_func, supp_args, num_pix, pix_size):
    
    ymap_whole = ymap_cluster_bubbles(nu, z, M_500, xb, yb, zb, xb2, yb2, zb2, rb, supp_func, supp_args, num_pix, pix_size)

    ymap_total = np.sum(ymap_whole)*((pix_size)**2)
    
    return ymap_total

#total = calculating_ymap_total(14.11e9, 0.216, 8e14, 0.24, 0.67, 0, -0.33, -0.82, 0, 0.5, f_thermal, [14.11e9, 1000], 512, 2/60)
#print(total)

def signal_MJy_sr(nu, z, M_500, xb, yb, zb, xb2, yb2, zb2, rb, supp_func, supp_args, num_pix, pix_size): # calculating the signal in MJy/sr by multiplying the cluster (without bubbles) ymap by spectral distortion g_x and subtracting bubble y map (suppressed integral from -z to +z) times modified spectral distortion g_tilde times the factor 270.33 which comes from Tony Mroczkowski [2019] sz review paper 
    return 270.33*g_x(nu)*(ymap_cluster(nu, z, M_500, xb, yb, rb, num_pix, pix_size) - ymap_bubbles(nu, z, M_500, xb, yb, zb, xb2, yb2, zb2, rb, supp_func, supp_args, num_pix, pix_size)) # in MJy/sr 

def signal_Jy_pix(nu, z, M_500, xb, yb, zb, xb2, yb2, zb2, rb, supp_func, supp_args, num_pix, pix_size): # to get the signal in Jy/pix I multiply by 1e6 to convert from MJy to Jy (pix_size/60))*(np.pi/180))**2 aqnd then multiply by pix_size (divided by 60 to convert to degrees I think) and then squared to get the area and I think the pi/180 to converting to sr or something but this is what Yvette did in the code she sent me called 'sim_MS0735.py' so I think it's right
    signal = signal_MJy_sr(nu, z, M_500, xb, yb, zb, xb2, yb2, zb2, rb, supp_func, supp_args, num_pix, pix_size)
    return signal*(1e6)*((((pix_size/60))*(np.pi/180))**2)

signal_output = signal_Jy_pix(14.11e9, 0.216, 8e14, 0.24, 0.67, 0, -0.33, -0.82, 0, 0.5, f_thermal, [14.11e9, 1258.93 ], 512, 2/60) # pix size in arcmin

pyplot.imshow(signal_output)
pyplot.colorbar()
pyplot.show()

from astropy.io import fits

f = fits.open('./MS0735_sims/MS0735_correct/ymap_template.fits')

f.info()

f[0].data = signal_output

#cl = ymap_cluster(14.11e9, 0.216, 8.4e14, 0.24, 0.67, 0.5, -10, 10, 0.078125)*g_x(14.11e9)
#print(cl)

#Print the corresponding values
f.info()

f[0].header['cdelt1'] = -0.000555556

f[0].header['cdelt2'] = 0.000555556

f[0].header['crval1'] = 15

f[0].header['crpix1'] = 257

f[0].header['crval2'] = -30

f[0].header['crpix2'] = 257

f.verify('fix')

f.writeto('./MS0735_sims/MS0735_correct/ms0735_thermal_Jy_pix_MUSTANG_upperkTe.fits', overwrite=True)

f.info()

end = timer()
#print(end - start)

end = timer()
#print(end - start)




