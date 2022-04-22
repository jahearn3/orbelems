# Data Loader
#
# Author: Joseph A'Hearn
# Created 06/23/2017
#
# This program loads data from an .aei or a .geo file 
#   

import numpy as np 

def aei_data(body):
	data = np.loadtxt(body + '.aei', skiprows=4)
	t = data[:,0] 						 
	x = data[:,1] 
	y = data[:,2] 
	z = data[:,3] 
	u = data[:,4] 
	v = data[:,5] 
	w = data[:,6]
	return t, x, y, z, u, v, w

def initial_aei_data(body):
	data = np.loadtxt(body + '.aei', skiprows=4)					 
	x = data[:,1] 
	y = data[:,2] 
	z = data[:,3] 
	u = data[:,4] 
	v = data[:,5] 
	w = data[:,6]
	return x[0], y[0], z[0], u[0], v[0], w[0]

def aei_data_at_specific_t(body, t0):
	data = np.loadtxt(body + '.aei', skiprows=4)					 
	x = data[:,1] 
	y = data[:,2] 
	z = data[:,3] 
	u = data[:,4] 
	v = data[:,5] 
	w = data[:,6]
	try:
		return x[t0], y[t0], z[t0], u[t0], v[t0], w[t0]
	except IndexError:
		return 0, 0, 0, 0, 0, 0
	#return x[t0], y[t0], z[t0], u[t0], v[t0], w[t0]

def final_aei_data(body):
	data = np.loadtxt(body + '.aei', skiprows=4)					 
	x = data[:,1] 
	y = data[:,2] 
	z = data[:,3] 
	u = data[:,4] 
	v = data[:,5] 
	w = data[:,6]
	return x[-1], y[-1], z[-1], u[-1], v[-1], w[-1]

def geo_data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	t = data[:,0] 						 
	a = data[:,1] 
	e = data[:,2] 
	i = data[:,3] 
	m = data[:,4] 
	p = data[:,5] 
	n = data[:,6]
	return t, a, e, i, m, p, n 

def exact_cer014data(file):
	data = np.loadtxt(file, skiprows=4)
	t = data[:,0]
	a = data[:,1]
	m = data[:,4]
	return t, a, m

def exact_cer14data(file):
	data = np.loadtxt(file, skiprows=4)
	a = data[:,1]
	m = data[:,4]
	return a, m

def initial_exact_cer14data(file):
	data = np.loadtxt(file, skiprows=4)
	a = data[:,1]
	m = data[:,4]
	return a[0], m[0]

def time_data(body):
	tdata = np.loadtxt(body + '.aei', skiprows=4)
	t = tdata[:,0]
	return t

def geotime_data(body):
	tdata = np.loadtxt(body + '.geo', skiprows=4)
	t = tdata[:,0]
	return t

def xy_data(body):
	data = np.loadtxt(body + '.aei', skiprows=4)					 
	x = data[:,1] 
	y = data[:,2] 
	return x, y

	return x, y
def xyz_data(body):
	data = np.loadtxt(body + '.aei', skiprows=4)					 
	x = data[:,1] 
	y = data[:,2] 
	z = data[:,3] 
	return x, y, z

def uv_data(body):
	data = np.loadtxt(body + '.aei', skiprows=4)					 
	u = data[:,4] 
	v = data[:,5] 
	return u, v

def xyuv_data(body):
	data = np.loadtxt(body + '.aei', skiprows=4)					 
	x = data[:,1] 
	y = data[:,2] 
	u = data[:,4] 
	v = data[:,5] 
	return x, y, u, v

def uvw_data(body):
	data = np.loadtxt(body + '.aei', skiprows=4)					 
	u = data[:,4] 
	v = data[:,5] 
	w = data[:,6] 
	return u, v, w

def radial_extrema(body):
	data = np.loadtxt(body + '.aei', skiprows=4)
	x = data[:,1] 
	y = data[:,2] 
	z = data[:,3] 
	r = np.sqrt(x**2 + y**2 + z**2) * 1.4959787E+11 # conversion to m
	return np.amin(r), np.amax(r)

def geo1mean(body):
	data = np.loadtxt(body + '.geo', skiprows=4)			 
	a = data[:,1] 
	return np.mean(a)   

def geo1data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	a = data[:,1] 
	return a

def geo12data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	a = data[:,1] 
	e = data[:,2] 
	return a, e

def geo2data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	e = data[:,2] 
	return e

def geo04data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	t = data[:,0] 
	m = data[:,4] 
	return t, m 

def geo045data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	t = data[:,0] 
	m = data[:,4] 
	p = data[:,5] 
	return t, m, p 

def geo124data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	a = data[:,1] 
	e = data[:,2] 
	m = data[:,4] 
	return a, e, m

def geo1245data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	a = data[:,1] 
	e = data[:,2] 
	m = data[:,4] 
	p = data[:,5] 
	return a, e, m, p 

def geo12356data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	a = data[:,1] 
	e = data[:,2] 
	i = data[:,3] 
	p = data[:,5] 
	n = data[:,6] 
	return a, e, i, p, n 

def geo125data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	a = data[:,1] 
	e = data[:,2] 
	p = data[:,5] 
	return a, e, p 

def geo23data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	e = data[:,2] 
	i = data[:,3] 
	return e, i

def geo25data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	e = data[:,2] 
	p = data[:,5] 
	return e, p

def geo14data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	a = data[:,1] 
	m = data[:,4] 
	return a, m

def initialgeo14data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	a = data[:,1] 
	m = data[:,4] 
	return a[0], m[0]

def geo4data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	m = data[:,4] 
	return m

def geo45data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	m = data[:,4] 
	p = data[:,5] 
	return m, p

def geo5data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	p = data[:,5] 
	return p

def geo456data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	m = data[:,4] 
	p = data[:,5] 
	n = data[:,6]
	return m, p, n

def geo56data(body):
	data = np.loadtxt(body + '.geo', skiprows=4)
	p = data[:,5] 
	n = data[:,6]
	return p, n

def ael_data(body):
	data = np.loadtxt(body + '.aei', skiprows=4)					 
	a = data[:,1] 
	e = data[:,2] 
	l = data[:,3] 
	return a, e, l 

def simulation_years():
	with open('param.in') as f:
		content = f.readlines()
	content = [x.strip() for x in content]  
	start_time = float(content[6][19:])
	stop_time  = float(content[7][19:])
	return (stop_time - start_time) / 365.25

