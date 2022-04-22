# Bodies
#
# Author: Joseph A'Hearn
# Created 06/22/2017
#
# This program provides functions to get info about the bodies from the mercury6 input files
#   

import numpy as np 
import random
from constants_of_mercury6 import G, AU, M_Sun

def central_body():
	with open('param.in') as f:
		content = f.readlines()
	content = [x.strip() for x in content]  
	J = [] 
	for i in range(29,32):
		J.append(content[i][13:])
	return float(content[27][30:]) * AU(), float(content[28][23:]) * M_Sun(), float(J[0]), float(J[1]), float(J[2])

def R():
	with open('param.in') as f:
		content = f.readlines()
	content = [x.strip() for x in content]
	return float(content[27][30:]) * AU()

def M():
	with open('param.in') as f:
		content = f.readlines()
	content = [x.strip() for x in content]
	return float(content[28][23:]) * M_Sun()

def mu():
	return G() * M()
def J2():
	return central_body()[2]
def J4():
	return central_body()[3]
def J6():
	return central_body()[4]

def full_body_list(filename='big.in'):
	body_list = []
	with open(filename) as bf:
		content = bf.readlines()
	content = [x.strip() for x in content]  
	for i in range(6, len(content)):	
		if((i - 6) % 4 == 0): 
			body_list.append(content[i][:8].strip(' '))
	body_list = list(filter(None, body_list)) # to remove empty strings
	return body_list

def full_mass_list(filename='big.in'):
	mass_list = []
	with open(filename) as bf:
		content = bf.readlines()
	content = [x.strip() for x in content]
	for i in range(6, len(content)):	
		if((i - 6) % 4 == 0):	
			line_with_mass = content[i][13:].split(' ')[0].translate({ord(j):'e' for j in 'd'})
			if(line_with_mass != ''): # avoid trying to convert to float if it's an empty string
				mass = float(line_with_mass)					
				mass_list.append(mass * M_Sun())
	return mass_list

def sort_data(filename='big.in'):
	body_lines = []
	x = []
	y = []
	z = []
	u = []
	v = []
	w = []
	with open(filename) as bf:
		content = bf.readlines()
	content = [x.strip() for x in content] # strips the "\n" from the end of each string 
	for i in range(6, len(content)):	
		if((i - 6) % 4 == 0): 			# singles out the lines containing the names of the bodies 
			body_lines.append(content[i]) # reads the whole line (including mass, etc.)
		if((i - 6) % 4 == 1):
			x.append(float(content[i].split()[0]))
			y.append(float(content[i].split()[1]))
			z.append(float(content[i].split()[2]))
		if((i - 6) % 4 == 2):
			u.append(float(content[i].split()[0]))
			v.append(float(content[i].split()[1]))
			w.append(float(content[i].split()[2]))
		# if there is need to change the spin, a few more lines of code must be added here
		#if((i - 6) % 4 == 3):
	return body_lines, x, y, z, u, v, w

def initial_state_vector(body):
	body_list = full_body_list()
	body_lines, x, y, z, u, v, w = sort_data()
	for i in range(len(body_list)):
		if(body_list[i] == body):
			return x[i], y[i], z[i], u[i], v[i], w[i]
		else:
			print(body + ' not found.')

#print(G())
#print(M())
#print(mu())
#
#print(str(((mu() / ((2 * np.pi / 86400)**2))**(1/3)) / 1000) + ' km')
