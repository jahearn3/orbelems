# Cartesian Position Vectors to Geocentric Orbital Elements
#
# Author: Joseph A'Hearn
# Created 11/07/2016
# Significant edits were made 07/11/2017
#
# This program converts the cartesian position vectors to geocentric orbital elements
#   It is based on Maryame el Moutamid's Fortran code
#   which is based on Renner & Sicardy 2006

import numpy as np

import constants_of_mercury6 
import bodies as bd 
from useful import print_progress
from data_loader import aei_data
import plot_assistant as pa 
import matplotlib.pyplot as plt 

# ------------------------------------------------------------------------------------------------------------------
# Functions

def arctangent(x, y):
	if(x == 0):
		if(y > 0):
			L = np.pi / 2
		else:
			L = (3 / 2) * np.pi
	elif(x < 0):
		L = np.arctan(y / x) + np.pi
	elif(np.arctan(y / x) < 0):
		L = np.arctan(y / x) + (2 * np.pi)
	else:
		L = np.arctan(y / x)
	return L

def compute_cylindrical_state_vectors(x, y, z, u, v, w, AU): # Computing the cylindrical state vectors (Equations 22-29 in Renner-Sicardy 2006)
	s = np.sqrt(x**2 + y**2) * AU
	L = arctangent(x, y)
	z_m = z * AU
	s_dot =  (( u * np.cos(L)) + (v * np.sin(L)))      * AU / 8.64E+04
	L_dot = (((-u * np.sin(L)) + (v * np.cos(L))) / s) * AU / 8.64E+04
	z_dot = w * AU / 8.64E+04
	return (s, L, z_m, s_dot, L_dot, z_dot)

def initialize_orbital_elements(s): 	# Initial assumption: perfect circle in the equatorial plane
	return (s, 0, 0, 0, 0, 0) # These become (a, e, i, mean_longitude, longit_perictr, long_ascndnode)

def initialize_second_order():
	return (0, 0, 0, 0, 0, 0) # These become (s_C, L_C, z_C, s_dot_C, L_dot_C, z_dot_C)

def mean_motion(R, mu, J2, J4, J6, geometric_elements): # Equation 14 in Renner-Sicardy 2006
	b = (R/geometric_elements[0]) 
	n_1st_order = (1.0 + (0.75 * ((b)**2.0) * J2) - (0.9375 * ((b)**4.0) * J4) + (1.09375 * ((b)**6.0) * J6))
	n_2nd_order = (- (0.28125 * ((b)**4.0) * J2**2.0) + (0.703125 * ((b)**6.0) * J2 * J4) + (0.2109375 * ((b)**6.0) * J2**3.0) + (3.0 * ((b)**2.0) * J2 * geometric_elements[1]**2.0)  - (12.0 * ((b)**2.0) * J2 * geometric_elements[2]**2.0))
	n = ((mu/((geometric_elements[0])**3.0))**0.5) * (n_1st_order + n_2nd_order)
	return n

def square_frequencies(R, mu, J2, J4, J6, geometric_elements):
	b = R / geometric_elements[0] 
	eta_sq = (mu / ((geometric_elements[0])**3)) * (1 -  (2.0  * (b**2) * J2) + (9.375  * (b**4) * J4) - (21.875  * (b**6) * J6))
	chi_sq = (mu / ((geometric_elements[0])**3)) * (1 +  (7.5  * (b**2) * J2) - (21.875 * (b**4) * J4) + (45.9375 * (b**6) * J6))
	return (eta_sq, chi_sq)

def freq_nsq(R, mu, J2, J4, J6, s0):
	b = R / s0
	nsq = (mu / (s0**3)) * (1.0 + (1.5 * J2 * (b**2)) - (1.875 * J4 * (b**4)) + (2.1875 * J6 * (b**6)))
	return nsq

def horizontal_epicyclic(R, mu, J2, J4, J6, geometric_elements):
	b = R / geometric_elements[0] 
	k1 = 1 - ((3 / 4) * J2 * (b**2)) + ((45 / 16) * J4 * (b**4)) - ((175 / 32) * J6 * (b**6))
	k2 = (- ((9 / 32) * (J2**2) * (b**4))) + ((135 / 64) * J2 * J4 * (b**6)) - ((27 / 128) * (J2**3) * (b**6)) - (9 * (b**2) * J2 * (geometric_elements[2]**2))
	return np.sqrt(mu / (geometric_elements[0]**3)) * (k1 + k2)

def vertical_epicyclic(R, mu, J2, J4, J6, geometric_elements):
	b = R / geometric_elements[0]
	n1 = 1 + ((9 / 4) * (b**2) * J2) - ((75 / 16) * (b**4) * J4) + ((245 / 32) * (b**6) * J6)
	n2 = (- ((81 / 32) * (b**4) * (J2**2))) + ((675 / 64) * (b**6) * J2 * J4) + ((729 / 128) * (b**6) * (J2**3))
	n3 = (6 * (b**2) * J2 * (geometric_elements[1]**2)) - ((51 / 4) * (b**2) * J2 * (geometric_elements[1]**2))
	return np.sqrt(mu / (geometric_elements[0]**3)) * (n1 + n2 + n3)

def all_frequencies(R, mu, J2, J4, J6, geometric_elements):  # Computing the Frequencies (Equations 14-21 in Renner-Sicardy 2006)
	n = mean_motion(R, mu, J2, J4, J6, geometric_elements)
	sq_frequencies = square_frequencies(R, mu, J2, J4, J6, geometric_elements)
	kappa = horizontal_epicyclic(R, mu, J2, J4, J6, geometric_elements)
	nu = vertical_epicyclic(R, mu, J2, J4, J6, geometric_elements)
	eta_sq = sq_frequencies[0]
	chi_sq = sq_frequencies[1]
	alpha_1  = (((2 * nu) + kappa) / 3)
	alpha_2  =  ((2 * nu) - kappa)
	alpha_sq = alpha_1 * alpha_2
	return (n, kappa, nu, eta_sq, chi_sq, alpha_1, alpha_2, alpha_sq)

def compute_geometric_elements(cylindrical_state_vectors, second_order, frequencies): # Computing the new values for the geometric elements (Equations 42-47 in Renner-Sicardy 2006)
	s_C = second_order[0]
	L_C = second_order[1]
	z_C = second_order[2]
	s_dot_C = second_order[3]
	L_dot_C = second_order[4]
	z_dot_C = second_order[5]
	n        = frequencies[0]
	kappa    = frequencies[1]
	nu       = frequencies[2]
	eta_sq   = frequencies[3]
	chi_sq   = frequencies[4]
	alpha_1  = frequencies[5]
	alpha_2  = frequencies[6]
	alpha_sq = frequencies[7]
	s = cylindrical_state_vectors[0] 
	L = cylindrical_state_vectors[1]
	z = cylindrical_state_vectors[2]
	s_dot = cylindrical_state_vectors[3]
	L_dot = cylindrical_state_vectors[4]
	z_dot = cylindrical_state_vectors[5]
	a = (s - s_C) / (1 - ((L_dot - L_dot_C - n) / (2 * n)))
	e = np.sqrt(((L_dot - L_dot_C - n) / (2 * n))**2 + ((s_dot - s_dot_C)/(a * kappa))**2)
	i = np.sqrt(((z - z_C) / (a))**2 + ((z_dot - z_dot_C)/(a * nu))**2)
	mean_longitude = L - L_C - (2 * (n / kappa) * ((s_dot - s_dot_C) / (a * kappa)))

	if(e < 1.0E-13): # to prevent longitude of pericenter from messing things up when the orbit is nearly circular
		longit_perictr = 0
	else: 
		lpx = a * kappa * (1 - ((s - s_C) / a))
		lpy = s_dot - s_dot_C 
		longit_perictr = mean_longitude - arctangent(lpx, lpy)
	if((z_dot - z_dot_C) == 0):
		long_ascndnode = 0 # it's actually undefined in this case, but for computation it is by convention set to zero
	else: 
		lanx = z_dot - z_dot_C
		lany = nu * (z - z_C)
		long_ascndnode = mean_longitude - arctangent(lanx, lany)

	return (a, e, i, mean_longitude, longit_perictr, long_ascndnode)

def compute_second_order(geometric_elements, frequencies):
	a = geometric_elements[0]
	e = geometric_elements[1]
	i = geometric_elements[2]
	mean_longitude = geometric_elements[3]
	longit_perictr = geometric_elements[4]
	long_ascndnode = geometric_elements[5]
	n        = frequencies[0]
	kappa    = frequencies[1]
	nu       = frequencies[2]
	eta_sq   = frequencies[3]
	chi_sq   = frequencies[4]
	alpha_1  = frequencies[5]
	alpha_2  = frequencies[6]
	alpha_sq = frequencies[7]
	# Computing the second order terms (Equations 36-41 in Renner-Sicardy 2006) 
	s_C_1 = ((3 / 2) * (eta_sq / (kappa**2)) - 1 - ((eta_sq / (2 * kappa**2)) * np.cos(2 * (mean_longitude - longit_perictr))))
	s_C_2 = ((3 / 4) * (chi_sq / (kappa**2)) - 1 + ((chi_sq / (4 * alpha_sq)) * np.cos(2 * (mean_longitude - long_ascndnode))))
	s_C = (a * (e**2) * s_C_1) + (a * (i**2) * s_C_2)
	L_C_1 = ((3 / 4) + (eta_sq / (2 * (kappa**2)))) * (n / kappa) * np.sin(2 * (mean_longitude - longit_perictr))  
	L_C_2 =            (chi_sq / (4 * alpha_sq))    * (n / nu)    * np.sin(2 * (mean_longitude - long_ascndnode))
	L_C = ((e**2) * L_C_1) - ((i**2) * L_C_2)
	z_C_1 =           (chi_sq / (2 * kappa * alpha_1)) * np.sin((2 * mean_longitude) - longit_perictr - long_ascndnode) 
	z_C_2 = (3 / 2) * (chi_sq / (kappa * alpha_2))     * np.sin(longit_perictr - long_ascndnode)
	z_C = a * i * e * (z_C_1 - z_C_2)
	s_dot_C_1 = (eta_sq / kappa)               * np.sin(2 * (mean_longitude - longit_perictr))
	s_dot_C_2 = (chi_sq / (2 * alpha_sq)) * nu * np.sin(2 * (mean_longitude - long_ascndnode))
	s_dot_C = a * (((e**2) * s_dot_C_1) - ((i**2) * s_dot_C_2))
	L_dot_C_1 = (7 / 2) - (kappa**2 / (2 * (n**2))) - (   3    *  eta_sq / (kappa**2))  + (((3 / 2) + (eta_sq / (kappa**2))) * np.cos(2 * (mean_longitude - longit_perictr)))
	L_dot_C_2 =    2    - (kappa**2 / (2 * (n**2))) - ((3 / 2) * (chi_sq / (kappa**2))) - ( chi_sq  / (2 * alpha_sq))        * np.cos(2 * (mean_longitude - long_ascndnode))
	L_dot_C = ((e**2) * n * L_dot_C_1) + ((i**2) * n * L_dot_C_2)	
	z_dot_C_1 =           (chi_sq * (kappa + nu) / (2 * kappa * alpha_1)) * np.cos(((2 * mean_longitude) - longit_perictr - long_ascndnode)) 
	z_dot_C_2 = (3 / 2) * (chi_sq * (kappa - nu) / (    kappa * alpha_2)) * np.cos(longit_perictr - long_ascndnode)
	z_dot_C = a * i * e * (z_dot_C_1 + z_dot_C_2)
	return (s_C, L_C, z_C, s_dot_C, L_dot_C, z_dot_C)

def refine_semimajor_axis(R, mu, J2, J4, J6, geometric_elements, cylindrical_state_vectors, epsilon, hz):
	nsq = freq_nsq(R, mu, J2, J4, J6, cylindrical_state_vectors[0])
	n0 = np.sqrt(nsq)
	s0c = 0
	iterate = True
	itercount = 0
	while(iterate == True):
		s0 = np.sqrt(hz / n0)
		if(np.absolute(s0c - s0) < epsilon):
			iterate = False
		s0c = s0
		nsq = freq_nsq(R, mu, J2, J4, J6, s0)
		n0 = np.sqrt(nsq)
		itercount += 1
	a_refined = s0 * (1 + (geometric_elements[1]**2) + (geometric_elements[2]**2))
	return a_refined
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

def xyz2geo_0(x, y, z, u, v, w):
	AU = constants_of_mercury6.AU()
	R = bd.R()
	mu = bd.mu()
	J2 = bd.J2()
	J4 = bd.J4()
	J6 = bd.J6()
	cylindrical_state_vectors = compute_cylindrical_state_vectors(x, y, z, u, v, w, AU)
	geometric_elements = initialize_orbital_elements(cylindrical_state_vectors[0])
	frequencies = all_frequencies(R, mu, J2, J4, J6, geometric_elements)
	second_order = initialize_second_order()
	a_prev = 0 
	iterate = True
	iteration = 0
	max_iterations = 150
	epsilon=1.0E-04
	while(iterate == True):
		geometric_elements = compute_geometric_elements(cylindrical_state_vectors, second_order, frequencies)
		second_order = compute_second_order(geometric_elements, frequencies)
		# breaks out if epsilon has been reached
		if(np.absolute(geometric_elements[0] - a_prev) < epsilon):
			iterate = False
		if(iteration + 1 == max_iterations):
			print(str(iteration) + " iterations complete: a = " + str(geometric_elements[0]) + ", delta_a = " + str(geometric_elements[0] - a_prev)) 
			iterate = False
		a_prev = geometric_elements[0]
		frequencies = all_frequencies(R, mu, J2, J4, J6, geometric_elements)
		iteration += 1
	return geometric_elements[0], geometric_elements[1], np.rad2deg(geometric_elements[2]), np.rad2deg(geometric_elements[3]), np.rad2deg(geometric_elements[4]), np.rad2deg(geometric_elements[5])


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Main
def xyz2geo(body_list=bd.full_body_list(), increment=1, epsilon=1.0E-04, use_L_z_for_a=False, max_iterations=150, plots=False):
	AU = constants_of_mercury6.AU()
	R = bd.R()
	mu = bd.mu()
	J2 = bd.J2()
	J4 = bd.J4()
	J6 = bd.J6()
	
	for body in body_list:
		with open(body + ".geo", "a") as myfile:
			myfile.write("\n\t\t\t\t\t\t\t  " + body + "\n\n")
			myfile.write("\tTime (days)\t\ta\t\t\t\te\t\t\t\t\ti\t\t\t\t\tlambda\t\t\tomega bar\t\tnode\n")
			t, x, y, z, u, v, w = aei_data(body)
			hz = ((x * v) - (y * u)) * (AU**2) / 8.64E+04
			j = 0
			if(plots == True):
				ml = []
				nn = []
				kk = []
				uu = []
				ee = []
				cc = []
				a1 = []
				a2 = []
				sc = []
				lc = []
				zc = []
				sd = []
				ld = []
				zd = []
			while(j < len(t)):
				cylindrical_state_vectors = compute_cylindrical_state_vectors(x[j], y[j], z[j], u[j], v[j], w[j], AU)
				geometric_elements = initialize_orbital_elements(cylindrical_state_vectors[0])
				frequencies = all_frequencies(R, mu, J2, J4, J6, geometric_elements)
				second_order = initialize_second_order()
				a_prev = 0 
				iterate = True
				iteration = 0
				if(plots == True):
					ite = []
					a = []
					e = []
					i = []
					m = []
					p = []
					l = []
				while(iterate == True):
					geometric_elements = compute_geometric_elements(cylindrical_state_vectors, second_order, frequencies)
					if(geometric_elements[0] < 0):
						print('Negative semimajor axis encountered on iteration ' + str(iteration)) 
						iterate = False
						break
					if(geometric_elements[1] > 1):
						print('Eccentricity > 1 encountered on iteration ' + str(iteration)) 
						iterate = False
						break
					second_order = compute_second_order(geometric_elements, frequencies)
					# breaks out if epsilon has been reached
					if(np.absolute(geometric_elements[0] - a_prev) < epsilon):
						iterate = False
					# breaks out after 150 iterations even if epsilon has not been reached
					if(iteration + 1 == max_iterations):
						print(str(iteration) + " iterations complete: a = " + str(geometric_elements[0]) + ", delta_a = " + str(geometric_elements[0] - a_prev)) 
						iterate = False
					a_prev = geometric_elements[0]
					frequencies = all_frequencies(R, mu, J2, J4, J6, geometric_elements)
					if(plots == True):
						ite.append(iteration)
						a.append(geometric_elements[0])
						e.append(geometric_elements[1])
						i.append(geometric_elements[2])
						m.append(geometric_elements[3])
						p.append(geometric_elements[4])
						l.append(geometric_elements[5])
					if(plots == True):
						if(iteration == 0):
							nn.append(np.rad2deg(frequencies[0]) * 8.64E+04)
							kk.append(np.rad2deg(frequencies[1]) * 8.64E+04)
							uu.append(np.rad2deg(frequencies[2]) * 8.64E+04)
							ee.append(np.rad2deg(np.sqrt(frequencies[3])) * 8.64E+04)
							cc.append(np.rad2deg(np.sqrt(frequencies[4])) * 8.64E+04)
							a1.append(np.rad2deg(frequencies[5]) * 8.64E+04)
							a2.append(np.rad2deg(frequencies[6]) * 8.64E+04)
							sc.append(second_order[0])
							lc.append(second_order[1])
							zc.append(second_order[2])
							sd.append(second_order[3])
							ld.append(second_order[4])
							zd.append(second_order[5])
					iteration += 1
				if(plots == True):
					ml.append(np.rad2deg(geometric_elements[3] - geometric_elements[4]) % 360)
					

				#-------------------------------------------------------------------------------------------------------------------
				myfile.write("\t\t%.3f" % t[j])  
				
				# semi-major axis
				if (use_L_z_for_a == True):
					a_refined = refine_semimajor_axis(R, mu, J2, J4, J6, geometric_elements, cylindrical_state_vectors, epsilon, hz[j])
					myfile.write("\t\t%.3f" % a_refined)
				else:
					myfile.write("\t\t%.3f" % geometric_elements[0]) 				
				# other orbital elements
				myfile.write("\t%.15f" % geometric_elements[1]) 								# eccentricity
				myfile.write("\t%.15f" % (np.rad2deg(geometric_elements[2]) % 360)) 			# inclination
				myfile.write("\t%.10f" % (np.rad2deg(geometric_elements[3]) % 360))				# mean longitude
				myfile.write("\t%.10f" % (np.rad2deg(geometric_elements[4]) % 360))				# longitude of pericenter
				myfile.write("\t%.10f" % (np.rad2deg(geometric_elements[5]) % 360) + "\n") 		# longitude of the ascending node

				
				if(plots == True):
					#if((body == 'MIMAS') and (j < 20)):
					if((iteration == max_iterations) or (j == 1)):
						fig, ax1, ax2, ax3, ax4, ax5, ax6 = pa.initialize32plot()
						fig, ax1 = pa.title_and_axes(fig, ax1, 'Semimajor axis', 'iteration', 'a',                     ite[0], ite[-1], 0.9999 * np.amin(a[-5:-1]), 1.0001 * np.amax(a[-5:-1]))
						fig, ax2 = pa.title_and_axes(fig, ax2, 'Eccentricity',   'iteration', 'e',                     ite[0], ite[-1])
						fig, ax3 = pa.title_and_axes(fig, ax3, 'Inclination',    'iteration', 'i (deg)',               ite[0], ite[-1])
						fig, ax4 = pa.title_and_axes(fig, ax4, 'Mean longitude', 'iteration', r'$\lambda$' + ' (deg)', ite[0], ite[-1])
						fig, ax5 = pa.title_and_axes(fig, ax5, 'Pericenter',     'iteration', r'$\varpi$' + ' (deg)',  ite[0], ite[-1])
						fig, ax6 = pa.title_and_axes(fig, ax6, 'Ascending node', 'iteration', r'$\Omega$' + ' (deg)',  ite[0], ite[-1])
						ax1.scatter(ite, a, s=5)
						ax2.scatter(ite, e, s=5)
						ax3.scatter(ite, np.rad2deg(i) % 360, s=5)
						ax4.scatter(ite, np.rad2deg(m) % 360, s=5)
						ax5.scatter(ite, np.rad2deg(p) % 360, s=5)
						ax6.scatter(ite, np.rad2deg(l) % 360, s=5)
						ax1.scatter(ite[-1], a[-1], s=10, label='main iteration')
						#ax1.scatter(ite[-1], a_refined, s=10, label='refinement')						
						ax2.scatter(ite[-1], e[-1], s=10)						
						ax3.scatter(ite[-1], np.rad2deg(i[-1]) % 360, s=10)
						ax4.scatter(ite[-1], np.rad2deg(m[-1]) % 360, s=10)
						ax5.scatter(ite[-1], np.rad2deg(p[-1]) % 360, s=10)
						ax6.scatter(ite[-1], np.rad2deg(l[-1]) % 360, s=10)
						ax1.legend(loc=0)
						plt.tight_layout()
						pa.save_and_clear32plot(fig, ax1, ax2, ax3, ax4, ax5, ax6, body + '_' + str(j) + '_' + str(iteration))

						

				if(t[j] == 0):
					print("Beginning conversion to geometric elements for " + str(body) + ".") 
	
				j += increment # moving on to the next data point
				print_progress(j, len(t))
		if(plots == True):
			fig, ax1, ax2, ax3, ax4, ax5, ax6 = pa.initialize32plot()
			fig, ax1 = pa.title_and_axes(fig, ax1, 'Mean motion',                    r'$\lambda - \varpi$', 'n'          + '(deg/day)', 0, 360)
			fig, ax2 = pa.title_and_axes(fig, ax2, 'Horizontal epicyclic frequency', r'$\lambda - \varpi$', r'$\kappa$'  + '(deg/day)', 0, 360)
			fig, ax3 = pa.title_and_axes(fig, ax3, 'Vertical epicyclic frequency',   r'$\lambda - \varpi$', r'$\nu$'     + '(deg/day)', 0, 360)
			fig, ax4 = pa.title_and_axes(fig, ax4, '',                               r'$\lambda - \varpi$', r'$\eta$'    + '(deg/day)', 0, 360)
			fig, ax5 = pa.title_and_axes(fig, ax5, '',                               r'$\lambda - \varpi$', r'$\chi$'    + '(deg/day)', 0, 360)
			fig, ax6 = pa.title_and_axes(fig, ax6, '',                               r'$\lambda - \varpi$', r'$\alpha$'  + '(deg/day)', 0, 360)
			ax1.scatter(ml, nn)                              
			ax2.scatter(ml, kk)
			ax3.scatter(ml, uu)
			ax4.scatter(ml, ee)
			ax5.scatter(ml, cc)
			ax6.scatter(ml, a1)
			ax6.scatter(ml, a2)
			plt.tight_layout()
			pa.save_and_clear32plot(fig, ax1, ax2, ax3, ax4, ax5, ax6, body + '_frequencies')

		if(plots == True):
			fig, ax1, ax2, ax3, ax4, ax5, ax6 = pa.initialize32plot()
			fig, ax1 = pa.title_and_axes(fig, ax1, '', r'$\lambda - \varpi$', r'$r_C$',       0, 360)
			fig, ax2 = pa.title_and_axes(fig, ax2, '', r'$\lambda - \varpi$', r'$\dot{r}_C$', 0, 360) 
			fig, ax3 = pa.title_and_axes(fig, ax3, '', r'$\lambda - \varpi$', r'$L_C$',       0, 360)
			fig, ax4 = pa.title_and_axes(fig, ax4, '', r'$\lambda - \varpi$', r'$\dot{L}_C$', 0, 360)
			fig, ax5 = pa.title_and_axes(fig, ax5, '', r'$\lambda - \varpi$', r'$z_C$'      , 0, 360)
			fig, ax6 = pa.title_and_axes(fig, ax6, '', r'$\lambda - \varpi$', r'$\dot{z}_C$', 0, 360)
			ax1.scatter(ml, sc)                              
			ax2.scatter(ml, sd)
			ax3.scatter(ml, lc)
			ax4.scatter(ml, ld)
			ax5.scatter(ml, zc)
			ax6.scatter(ml, zd)
			plt.tight_layout()
			pa.save_and_clear32plot(fig, ax1, ax2, ax3, ax4, ax5, ax6, body + '_2nd_order_terms')
	print("    Conversion to geometric orbital elements complete.")
		
def unitcheck():
	print('AU: ' + str(constants_of_mercury6.AU())) 
	print('R: ' + str(bd.R())) 
	print('mu: ' + str(bd.mu()))
	print('J2: ' + str(bd.J2()))
	print('J4: ' + str(bd.J4()))
	print('J6: ' + str(bd.J6()))



#xyz2geo()
#unitcheck()