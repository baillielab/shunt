#!/usr/bin/python
# -*- coding: UTF-8 -*-


import cgitb, string, math, cgi

#-##############--------------------------------------------------

from math import *
import numpy as np
from scipy import optimize, integrate

# redefined in es function:
Q = 5 # l/min
VO2 = 0.25 # l/min
maxOER = 0.8
RER = 0.8
DPG = 0.00465 # assumed, doesn't make much difference anyway

# constant throughout
MCHC = 340

#----------------------------------------------------------------------------
# Name:		Input Constants
#----------------------------------------------------------------------------
# Fractional water contents
Wpl = 0.94 # plasma
Wrbc = 0.65 # red blood cells
# Constants from Dash and Bassingthwaight 2010
n0 = 1.7 # unitless
K_1 = 7.43e-7 # M
K_prime_1 = 1.35*10**-3 # unitless
K_2prime_1 = 5.5*10**-4 # M
K_2 = 2.95*10**-5 # unitless
K_2prime_2 = 1*10**-6 # M
K_prime_2 = K_2*K_2prime_2**-1 # M**-1
K_3 = 2.51*10**-5 # unitless
K_2prime_3 = 1*10**-6 # M
K_prime_3 = K_3*K_2prime_3**-1 # M**-1
#K_2prime_4 = 202123 # M**-1
K_2prime_5 = 2.63*10**-8 # M
K_2prime_6 = 1.91*10**-8 # M
# Constants from Wagner & Pruss 1993
Temp_critical = 647.096 # K
Pres_critical = 22.064e3 # kPa
a1 = -7.85951783
a2 = 1.84408259
a3 = -11.7866497
a4 = 22.6807411
a5 = -15.9618719
a6 = 1.80122502
# Standard temperature and pressure
R = 8.3145 # J.K**-1.mol**-1
STP_T = 273.15 # K
STP_P = 101.325 # kPa
# Standard bicarbonate
StdBicarb = 24.5e-3 # M

def calculateglobalvariables(current_temp, thisHb):
	global alphaO2, alphaCO2, PH2O, Wbl, HbMol, Hct, Hb, Temp
	alphaO2 = alphaO2_func(current_temp, Wpl) # M/kPa
	alphaCO2 = alphaCO2_func(current_temp, Wpl) # M/kPa
	tau = 1 - current_temp*Temp_critical**-1 # fraction
	PH2O = Pres_critical*e**(Temp_critical*current_temp**-1*(a1*tau+a2*tau**1.5+a3*tau**3+a4*tau**3.5+a5*tau**4+a6*tau**7)) # kPa
	HbMol = thisHb*64458**-1 # M
	Hct  = thisHb/(MCHC*10)
	Hb = thisHb
	Wbl = (1-Hct)*Wpl+Hct*Wrbc # fraction
	Temp = current_temp

#------------------------------------------------------------------------
# Name:		 alpha O2 / CO2
# Source:	   Dash & Bassingthwaight 2010
#------------------------------------------------------------------------
def alphaO2_func(Temp, Wpl):
	return 0.1333**-1*(1.37-0.0137*(Temp-310.15)+0.00058*(Temp-310.15)**2)*\
		(1e-6*Wpl**-1) # M/kPa
#----------------------------------------------------------------------------
# Name:		 alphaCO2
# Source:	   Kelman 1967
#----------------------------------------------------------------------------
def alphaCO2_func(Temp, Wpl):
	return 0.1333**-1*(3.07-0.057*(Temp-310.15)+0.002*(Temp-310.15)**2)*\
		(1e-5*Wpl**-1) # M/kPa
#----------------------------------------------------------------------------
# Name:	 Henderson-Hasselbalch Equation
#----------------------------------------------------------------------------
def HH_1(PnCO2, pH):
	'''Calculate CO2 as bicarbonate in SOLUTION'''
	return (K_1*alphaCO2*PnCO2)*(10**-pH)**-1 # M
#----------------------------------------------------------------------------
# Name:		 p50
# Source:	   Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def P50(pH, PnCO2, DPG, Temp):
	p50_pH = -2.84*(pH+log10(0.69)-7.24)+1.18*(pH+log10(0.69)-7.24)**2
	p50_PnCO2 = 4.82e-2*(PnCO2-5.332)+3.64e-5*(PnCO2-5.332)**2
	p50_DPG = 1.06e2*(DPG-4.65e-3)-2.62e3*(DPG-4.65e-3)**2
	p50_Temp = 1.99e-1*(Temp-310.15)+5.78e-3*(Temp-310.15)**2+9.33e-05*(Temp-310.15)**3
	return 3.57+p50_PnCO2+p50_pH+p50_DPG+p50_Temp # kPa
#----------------------------------------------------------------------------
# Name:		 SnO2
# Source:	   Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def SnO2_0(PnO2, p50_SnO2):
	'''Calculate O2 SATURATION from PARTIAL PRESSURE''' # fraction
	return ((PnO2*p50_SnO2**-1)**(1+n0))/\
		(1+((PnO2*p50_SnO2**-1)**(1+n0)))	
def SnO2_1(PnO2, pH, PnCO2, DPG, Temp): # Equation B.3
	'''Calculate O2 SATURATION from PARTIAL PRESSURE''' # fraction
	return ((PnO2*P50(pH, PnCO2, DPG, Temp)**-1)**(1+n0))/\
		(1+((PnO2*P50(pH, PnCO2, DPG, Temp)**-1)**(1+n0)))
def SnO2_2_null(Sats, CnO2, p50_SnO2):
	'''returns 0 '''
	return Wbl*alphaO2*p50_SnO2*(Sats*(1-Sats)**-1)**((1+n0)**-1) + \
		(4*HbMol)*Sats - (CnO2*(R*STP_T*STP_P**-1*1e2)**-1)
def SnO2_2(CnO2, p50_SnO2):
	'''Calculate O2 SATURATION from CONTENT''' # fraction
	args = (CnO2, p50_SnO2)
	return optimize.brentq(SnO2_2_null,1e-15,(1-1e-15),args=args,rtol=0.000001)

#----------------------------------------------------------------------------
def SnO2_1_null(p50, PnO2, SnO2):
	'''returns zero when saturation is found'''
	estimated_SnO2 = SnO2_0(PnO2, p50)
	return SnO2 - estimated_SnO2

def infer_p50(p1,s1):
	'''Calculate O2 SATURATION from CONTENT''' # fraction
	args = (p1, s1)
	return optimize.brentq(SnO2_1_null,1e-3,12,args=args,rtol=0.000001)

#----------------------------------------------------------------------------
# Name:		 Blood O2 content
# Source:	   Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def CnO2_0(PnO2, p50_CnO2):
	'''Calculate O2 CONTENT from PARTIAL PRESSURE'''
	return (Wbl*p50_CnO2*(SnO2_0(PnO2, p50_CnO2)*(1\
	-SnO2_0(PnO2, p50_CnO2))**-1)**((1+n0)**-1)*alphaO2\
	+4*HbMol*SnO2_0(PnO2, p50_CnO2))*(R*STP_T*STP_P**-1*1e2)
	# ml of O2 per 100ml blood STP
def CnO2_1(PnO2, pH, PnCO2, DPG, Temp):
	'''Calculate O2 CONTENT from PARTIAL PRESSURE'''
	return (Wbl*P50(pH,PnCO2,DPG,Temp)*(SnO2_1(PnO2,pH,PnCO2,DPG,Temp)*(1\
	-SnO2_1(PnO2,pH,PnCO2,DPG,Temp))**-1)**((1+n0)**-1)*alphaO2\
	+4*HbMol*SnO2_1(PnO2,pH,PnCO2,DPG,Temp))*(R*STP_T*STP_P**-1*1e2)
	# ml of O2 per 100ml blood STP
#----------------------------------------------------------------------------
# Name:		 Blood O2 partial pressure
# Source:	   Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def PnO2_1(SnO2, p50_SnO2):
	'''Calculate O2 PARTIAL PRESSURE from SATURATION'''
	return p50_SnO2*(SnO2*(1-SnO2)**-1)**(((1+n0))**-1) # kPa
def PnO2_2(CnO2,p50_SnO2):
	'''Calculate O2 PARTIAL PRESSURE from CONTENT'''
	SnO2 = SnO2_2(CnO2, p50_SnO2)
	return p50_SnO2*(SnO2*(1-SnO2)**-1)**(((1+n0))**-1) # kPa
#----------------------------------------------------------------------------
# Name:		 CnCO2
# Source:	   Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def Kratio(PnCO2, pH, Temp, Wpl):
	Hrbc = 10**-(pH+log(0.69,10)) # M
	CO2 = PnCO2*alphaCO2_func(Temp, Wpl) # M
	return (K_prime_2*CO2*(1+K_2prime_2*Hrbc**-1)+(1+Hrbc*K_2prime_5**-1))*\
	(K_prime_3*CO2*(1+K_2prime_3*Hrbc**-1)+(1+Hrbc*K_2prime_6**-1))**-1 # unitless
def CnCO2_Dissolved(PnCO2):
	'''Calculate CO2 dissolved in WHOLE BLOOD'''
	return Wbl*alphaCO2*PnCO2 # M
def CnCO2_Bicarb(pH, PnCO2):
	'''Calculate CO2 as bicarbonate in WHOLE BLOOD'''
	return ((1-Hct)*Wpl+Hct*Wrbc*0.69)*HH_1(PnCO2, pH) # M
def CnCO2_HbBound(PnCO2, PnO2, pH):
	'''Calculate Hb-CO2 from PARTIAL PRESSURE'''
	K_prime_4 = (alphaO2*PnO2)**n0*Kratio(PnCO2,pH,Temp,Wpl)*(P50(pH, PnCO2, DPG, Temp)*alphaO2)**-(1+n0)
	SHbCO2 = ( \
	K_prime_2*alphaCO2*PnCO2*(1+K_2prime_2*(10**-(pH+log(0.69,10)))**-1) + \
	K_prime_3*alphaCO2*PnCO2*(1+K_2prime_3*(10**-(pH+log(0.69,10)))**-1)*K_prime_4*alphaO2*PnO2
	)*( \
	K_prime_2*alphaCO2*PnCO2*(1+K_2prime_2*(10**-(pH+log(0.69,10)))**-1) + \
	(1+(10**-(pH+log(0.69,10)))*K_2prime_5**-1) + \
	K_prime_4*alphaO2*PnO2*(K_prime_3*alphaCO2*PnCO2*(1+K_2prime_3*(10**-(pH+log(0.69,10)))**-1) + \
	(1+(10**-(pH+log(0.69,10)))*K_2prime_6**-1)) \
	)**-1
	return 4*HbMol*SHbCO2 # M
def CnCO2_1(pH, PnCO2, PnO2):
	'''Calculate CO2 CONTENT from PARTIAL PRESSURE'''
	return (CnCO2_HbBound(PnCO2,PnO2,pH)+CnCO2_Bicarb(pH,PnCO2)\
		+CnCO2_Dissolved(PnCO2))*(R*STP_T*STP_P**-1*1e2) # ml CO2 per 100ml blood STP
def PnCO2_null(PnCO2, CnCO2, pH, CnO2):
	'''returns 0'''
	p50_PnCO2_null = P50(pH, PnCO2, DPG, Temp)
	return CnCO2_1(pH, PnCO2, PnO2_2(CnO2, p50_PnCO2_null))-CnCO2 # null
def PnCO2_1(CnCO2, pH, CnO2, PnCO2estimate):
	'''Calculate CO2 PARTIAL PRESSURE from CONTENTS'''
	args = CnCO2,pH,CnO2
	PnCO2 = optimize.newton(PnCO2_null,PnCO2estimate,args=args,tol=0.01)
	return PnCO2 # kPa
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

def alvgas(f, c, r=RER):
	ag = f*95.05885 - c/r + f*c*(1-r)/r
	return max(ag, 0.5)

def es(fio2, pao2, paco2, pH, thisHb=8, thisTemp=309.65, thisVO2=0.25, thisQ=5, thismaxOER=0.8, thisRER=0.8, thisDPG=0.00465, mode="quiet"):
	''' standalone function to calculate effective shunt from ABG '''
	''' inputs are in kPa pH g/dl Celsius l/min lo2/min'''
	thisVO2 = thisVO2*1000 # convert VO2 to mls/min here
	calculateglobalvariables(thisTemp, thisHb)
	p50 = P50(pH, paco2, thisDPG, thisTemp)
	cao2 = CnO2_0(pao2, p50)
	pAo2 = alvgas(fio2,paco2,thisRER)
	cco2 = CnO2_0(pAo2, p50)
	caco2 = CnCO2_1(pH, paco2, pao2)
	DO2 = thisQ*cao2*10
	thisVO2 = min(thisVO2, DO2*thismaxOER)
	cvo2 = ((DO2-thisVO2)/thisQ)/10
	cvco2 = ((thisQ*caco2*10 + thisVO2 * thisRER)/thisQ)/10
	actualOER = (cao2-cvo2)/cao2
	qsqt = (cao2 - cco2) / (cvo2 - cco2)
	# sanity checks
	if False:
		print ("p50", p50)
		print ("cao2", cao2)
		print ("pAo2", pAo2)
		print ("cco2", cco2)
		print ("caco2", caco2)
		print ("DO2", DO2)
		print ("thisVO2", thisVO2)
		print ("cvo2", cvo2)
		print ("cvco2", cvco2)
		print ("actualOER", actualOER)
		print ("qsqt", qsqt)
	if cao2 > cco2+0.1:
		return fail('cao2 > cco2', mode=mode)
	elif qsqt < 0:
		return fail('qsqt < 0', mode=mode)
	elif qsqt > 1:
		return fail('qsqt > 1', mode=mode)
	return qsqt

def fail(message, mode="quiet"):
	if mode=="verbose":
		return "Failed {}".format(message)
	elif mode=="quiet":
		return np.nan

#-------------------------------------
def get_local_inputs():
	import traceback
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-fio2',			default=0.21,		type=float,	help='fraction')
	parser.add_argument('-pao2',			default=13.3,		type=float,	help='kPa')
	parser.add_argument('-paco2',			default=5.3,		type=float,	help='kPa')
	parser.add_argument('-pH',				default=7.4,		type=float,	help='pH')
	#-----
	parser.add_argument('-gasunit',			default='kPa',		type=str, 	help='kPa or mmHg')
	parser.add_argument('-acidunit',		default='pH',		help='pH or H')
	#-----
	parser.add_argument('-Hb',				default=80,			type=float,	help='g/l')
	parser.add_argument('-Temp',	 		default=309.65,		type=float,	help='K') # 36.5 C = 309.65 K
	parser.add_argument('-VO2',				default=0.25,		type=float,	help='l/min')
	parser.add_argument('-Q',				default=6.5,		type=float,	help='l/min')
	parser.add_argument('-maxOER',			default=0.8,		type=float,	help='fraction')
	parser.add_argument('-DPG',				default=0.00465,	type=float,	help='M')
	parser.add_argument('-RER',				default=0.8,		type=float,	help='fraction')
	parser.set_defaults()
	args = parser.parse_args()
	localvars = vars(args)
	return localvars


def get_online_inputs(thesevars):  # detect the source of input variables automatically
	form = cgi.FieldStorage()
	onlinevars = {}
	online = False
	for v in thesevars:
		try:
			float(thesevars[v])
		except: 
			onlinevars[v] = form.getvalue(v)
			continue # this doesn't have to be a float. all others must be floatable.
		try:
			onlinevars[v] = float(form.getvalue(v))
			online = True 
		except:
			onlinevars[v] = thesevars[v] # default value
	if online:
		print ("Access-Control-Allow-Origin: *")
		print ("Content-Type: text/plain;charset=utf-8")
		print
	return onlinevars
	
def getinputs():
	variables = get_local_inputs()
	try:
		variables = get_online_inputs(variables)
	except:
		pass
	return variables


def setunits(value,unit):
	# Mass
	if unit == 'kg':
		return value
	elif unit == 'stone':
		return value*6.35029318
	elif unit == 'lb':
		return value*0.45359237
	# Length
	elif unit == 'm':
		return value
	elif unit == 'feet':
		return value*0.3048
	elif unit == 'ft':
		return value*0.3048
	# Temperature
	elif unit == 'deg C':
		return value+273.15
	elif unit == 'deg F':
		return (5*(value-32)*9**-1)+273.15
	elif unit =='K':
		return value
	# Concentration
	elif unit == 'mmol/l':
		return value*1e-3
	elif unit == 'mol/l':
		return value
	elif unit == 'mEq/l':
		return value*1e-3
	elif unit == 'Eq/l':
		return value
	# Haemoglobin
	elif unit == 'g/dl':
		return value*10
	elif unit == 'g/l':
		return value
	# Unitless
	elif unit == 'fraction':
		return value
	elif unit == '%':
		return value*1e-2
	elif unit == 'unitless':
		return value
	# Volume
	elif unit == 'ml':
		return value*1e-3
	elif unit =='l':
		return value
	# Rate
	elif unit == 'l/min':
		return value
	elif unit == 'ml/min':
		return value*10**-3
	elif unit == 'bpm':
		return value
	elif unit == 'mlO2/min/kPa':
		return value
	# Percentages
	elif unit == "%":
		return float(value)/100
	elif unit == "percentage":
		return float(value)/100
	elif unit == "fraction":
		return value
	# Error
	else:
		return 'Unit conversion error - unit not recognised'

#-------------------------------------
if __name__ == "__main__":
	import cgi
	import argparse
	inputs = getinputs()
	#print (inputs)
	shunt = es(
		fio2 = inputs['fio2'],
		pao2 = inputs['pao2'],
		paco2 = inputs['paco2'],
		pH = inputs['pH'],
		thisHb = inputs['Hb'],
		thisTemp = inputs['Temp'],
		thisVO2 = inputs['VO2'],
		thisQ = inputs['Q'],
		thismaxOER = inputs['maxOER'],
		thisRER = inputs['RER'],
		thisDPG = inputs['DPG'],
		)
	print (shunt)
















