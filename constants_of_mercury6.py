# Constants of mercury6
#
# Author: Joseph A'Hearn
# Created 06/22/2017
#
# 
#   

def constants(cgs=False):
	with open('mercury.inc') as f:
		content = f.readlines()
	content = [x.strip() for x in content] # strips the "\n" from the end of each string 
	line_with_k2 = content[40][16:].translate({ord(i):None for i in ')'}) 
	k2_cgs  = float(line_with_k2.translate({ord(i):'e' for i in 'd'})) 
	AU_cm   = float(content[41][16:].translate({ord(i):None for i in ')'})) 
	M_Sun_g = float(content[42][18:].translate({ord(i):None for i in ')'})) 
	if(cgs == True):
		return k2_cgs, AU_cm, M_Sun_g
	else:
		AU_m = AU_cm / 100
		M_Sun_kg = M_Sun_g / 1000
		k2_SI = k2_cgs * (AU_m)**3 / ((86400)**2 * M_Sun_kg)
		return k2_SI, AU_m, M_Sun_kg

# the value of G that should be used for mercury6 conversions: 6.671984222963699e-11 
# accepted value today: 6.67259E-11

def G():
	return constants()[0]
def AU():
	return constants()[1]
def M_Sun():
	return constants()[2]