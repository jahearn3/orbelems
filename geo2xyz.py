# Geocentric Orbital Elements to State Vectors
#
# Author: Joseph A'Hearn
# Created 03/20/2017
# Significant edits made 06/23/2017
#
# This program converts initial orbital elements to initial Cartesian position and velocity components
#

import numpy as np
from constants_of_mercury6 import G, AU
from bodies import R_Sat, mu_Sat, J2, J4, J6
# ------------------------------------------------------------------------------------------------------------------

def geo2xyz(a=167537.5, e=0.000457806317, i=0.001319397343, mean_longitude=207, longit_perictr=157.983521, long_ascnode=319.6712727):
	# input parameters should be done in km and degrees
	a = a * 1.0E+03 	# converting to meters
	i = np.deg2rad(i)	# converting to radians
	mean_longitude = np.deg2rad(mean_longitude)  	
	longit_perictr = np.deg2rad(longit_perictr)  			
	long_ascnode   = np.deg2rad(long_ascnode)  				
	# ----------------------------------------------------------------------------------------------------------	
	# Extra computations to make the next lines look simpler and prettier
	b = R_Sat() / a 
	h = mu_Sat() / (a**3.0)
	# Computing the Frequencies (Equations 14-21 in Renner-Sicardy 2006)
	n        = np.sqrt(h) * (1.0 + (0.75 * (b**2.0) * J2()) - (0.9375 * (b**4.0) * J4()) + (1.09375 * (b**6.0) * J6()) - (0.28125 * (b**4.0) * (J2()**2.0)) + (0.703125  * (b**6.0) * J2() * J4()) + (0.2109375 * (b**6.0) * J2()**3.0) + (3.0 * (b**2.0) * J2() * (e**2.0)) - (12.0  * (b**2.0) * J2() * (i**2.0)))
	kappa    = np.sqrt(h) * (1.0 - (0.75 * (b**2.0) * J2()) + (2.8125 * (b**4.0) * J4()) - (5.46875 * (b**6.0) * J6()) - (0.28125 * (b**4.0) * (J2()**2.0)) + (2.109375  * (b**6.0) * J2() * J4()) - (0.2109375 * (b**6.0) * J2()**3.0) - (9.0 * (b**2.0) * J2() * (i**2.0)))
	nu       = np.sqrt(h) * (1.0 + (2.25 * (b**2.0) * J2()) - (4.6875 * (b**4.0) * J4()) + (7.65625 * (b**6.0) * J6()) - (2.53125 * (b**4.0) * (J2()**2.0)) + (10.546875 * (b**6.0) * J2() * J4()) + (5.6953125 * (b**6.0) * J2()**3.0) + (6.0 * (b**2.0) * J2() * (e**2.0)) - (12.75 * (b**2.0) * J2() * (i**2.0)))
	eta_sq   = h * (1.0 -  (2.0  * (b**2.0) * J2()) + (9.375  * (b**4.0) * J4()) - (21.875  * (b**6.0) * J6()))
	chi_sq   = h * (1.0 +  (7.5  * (b**2.0) * J2()) - (21.875 * (b**4.0) * J4()) + (45.9375 * (b**6.0) * J6()))
	alpha_1  = ((2.0 * nu) + kappa) / 3.0
	alpha_2  = ((2.0 * nu) - kappa)
	alpha_sq = alpha_1 * alpha_2
	# To Cylindrical Coordinates (Equations 2-7 in Renner-Sicardy 2006)
	r = a * (1.0 - (e * np.cos(mean_longitude - longit_perictr)) + ((e**2.0) * ((1.5 * (eta_sq/(kappa**2.0))) - 1.0 - (eta_sq/(2.0 * (kappa**2.0))) * np.cos(2.0*(mean_longitude - longit_perictr)))) + i**2.0 * ((0.75 * (chi_sq/(kappa**2.0))) - 1.0 + (chi_sq/(4.0 * alpha_sq)) * np.cos(2.0*(mean_longitude - long_ascnode)))) 
	L = mean_longitude + (2.0 * e * (n/kappa) * np.sin(mean_longitude - longit_perictr)) + ((e**2.0) * (0.75 + (eta_sq/(2.0 * (kappa**2.0)))) * (n/kappa) * np.sin(2.0*(mean_longitude - longit_perictr))) - ((i**2.0) * (chi_sq/(4.0 * alpha_sq)) * (n/nu) * np.sin(2.0*(mean_longitude - long_ascnode)))
	z = a * i * (np.sin(mean_longitude - long_ascnode) + (e * chi_sq/(2.0 * kappa * alpha_1)) * np.sin((2.0 * mean_longitude) - longit_perictr - long_ascnode) - (e * 1.5 * (chi_sq/(kappa * alpha_2)) * np.sin(longit_perictr - long_ascnode)))
	rdot = a * kappa * ((e * np.sin(mean_longitude - longit_perictr)) + ((e**2.0) * (eta_sq/(kappa**2.0)) * np.sin(2.0*(mean_longitude - longit_perictr))) - ((i**2.0) * (chi_sq/(2.0 * alpha_sq)) * (nu/kappa) * np.sin(2.0*(mean_longitude - long_ascnode))))
	Ldot = n * (1.0 + (2.0 * e * np.cos(mean_longitude - longit_perictr)) + (e**2.0 * (3.5 - (3.0 * eta_sq/(kappa**2.0)) - (kappa**2.0 / (2.0 * (n**2.0))) + ((1.5 + (eta_sq/(kappa**2.0))) * np.cos(2.0*(mean_longitude - longit_perictr))))) + ((i**2.0) * (2.0 - ((kappa**2.0)/(2.0 * (n**2.0))) - (1.5 * (chi_sq/(kappa**2.0))) - ((chi_sq/(2.0 * alpha_sq)) * np.cos(2.0*(mean_longitude - long_ascnode))))))
	zdot = a * i * nu * (np.cos(mean_longitude - long_ascnode) + (e * ((chi_sq * (kappa + nu)) / (2.0 * kappa * alpha_1 * nu)) * np.cos((2.0 * mean_longitude) - longit_perictr - long_ascnode)) + (e * ((1.5) * (chi_sq * (kappa - nu)) / (kappa * alpha_2 * nu)) * np.cos(longit_perictr - long_ascnode)))
	# To Cartesian Coordinates (Equations 8-13 in Renner-Sicardy 2006) N.B. z and zdot are the same in cylindrical and cartesian
	x = r * np.cos(L)
	y = r * np.sin(L)
	xdot = (rdot * np.cos(L)) - (r * Ldot * np.sin(L))
	ydot = (rdot * np.sin(L)) + (r * Ldot * np.cos(L))
	return x / AU(), y / AU(), z / AU(), xdot * 86400.0 / AU(), ydot * 86400.0 / AU(), zdot * 86400.0 / AU() # units of AU and AU/day
