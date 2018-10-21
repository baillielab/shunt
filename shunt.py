#!/opt/local/bin/python
# -*- coding: UTF-8 -*-
# Authors: JK Baillie, A Bretherick

import numpy as np
from math import *
from scipy import optimize, integrate

# redefined in es function:
Q = 5 # l/min
VO2 = 250 # mlO2/min
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

def calculate_shunt(abg, thisVO2):
	DO2= Q*abg["cao2"]*10
	thisVO2 = min(thisVO2, DO2*maxOER) 
	cvo2 = ((DO2-thisVO2)/Q)/10
	cvco2 = ((Q*abg["caco2"]*10 + thisVO2*RER)/Q)/10
	OER = (abg["cao2"]-cvo2)/abg["cao2"]
	qsqt = (abg["cao2"] - abg["cco2"]) / (cvo2 - abg["cco2"])
	return cvo2, cvco2, OER, qsqt, thisVO2

def predict_from_shunt(abg, thisVO2, newfio2, newpaco2, newp50):
	# ? should include new CO2, new Hb, new p50
	cvo2, cvco2, OER, qsqt, thisVO2	= calculate_shunt(abg, thisVO2)
	Qt = Q 
	Qs = Q * qsqt  #mls/min
	newpAo2 = alvgas(newfio2, newpaco2) # could possibly correct for p50 in capillaries
	calculateglobalvariables(abg["Temp"], abg["Hb"])
	newcco2 = CnO2_0(newpAo2, newp50)
	newcao2 = newcco2 - (abg["cco2"]-abg["cao2"])
	newpao2 = PnO2_2(newcao2, newp50)
	expected_error = -1.84*(abg["fio2"]- newfio2)+0.26
	outdic = {
		"adiff": abs (abg["cao2"] - newcao2),
		"s_newcao2": newcao2,
		"s_newpao2": newpao2,
		"sc_newpao2": newpao2-expected_error,
		"cvo2": cvo2,
		"cvco2": cvco2,
		"OER": OER,
		"qsqt": qsqt,
		"Qs": Qs,
		"Qt": Qt,
		}
	return outdic

def getvo2(abg1, abg2):
	diffs={}
	for v in range(100,710,10):
		diffs[v] = abs(predict_from_shunt(abg1, v, abg2["fio2"])["newpao2"] - abg2["cao2"])
	return sorted(iter(list(diffs.items())), key=lambda k_v: (k_v[1],k_v[0]), reverse=True)[0][0]

def predict_from_aa(abg, newfio2):
	newpAo2 = alvgas(newfio2, abg["paco2"])
	return newpAo2 - abg["Aa"]

def predict_from_pf(abg, newfio2):
	return newfio2 * abg["PF"]

def predict_from_pf_corrected(abg, newfio2):
	prediction = newfio2 * abg["PF"]
	expected_error = -5.6*(abg["fio2"]-newfio2)-0.07
	return prediction - expected_error

def predict_from_pa(abg, newfio2):
	newpAo2 = alvgas(newfio2, abg["paco2"])
	return newpAo2 * abg["PA"]


#-------------------------------------
def es(fio2, pao2, paco2, pH, Hb=8, tempC=36.5, thisVO2=250, thisQ=5, thismaxOER=0.8, thisRER=0.8, thisDPG=0.00465):
	''' standalone function to calculate effective shunt from ABG '''
	''' inputs are in kPa pH g/dl Celsius l/min mlo2/min'''
	thisTemp = tempC + 273.15
	calculateglobalvariables(thisTemp, Hb)
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
	return qsqt

#-------------------------------------
def getinputs():  # detect the source of input variables automatically
	try:
		form = cgi.FieldStorage()
		variables =  {
			'fio2':[float(form.getvalue("fio2")),'%'],
			'Hb':[float(form.getvalue("Hb")),form.getvalue("Hb_unit")],
			'Temp':[float(form.getvalue("Temp")),form.getvalue("Temp_unit")],
			'VO2':[float(form.getvalue("VO2")),'ml/min'],
			'Q':[float(form.getvalue("Q")),form.getvalue("Q_unit")],
			'maxOER':[float(form.getvalue("MaxOER")),'fraction'],
			'RER':[float(form.getvalue("RER")),'fraction'],
			'DPG':[float(form.getvalue("DPG")),form.getvalue("DPG_unit")],
		}
		# if we get this far, we must be online, so send headers
		print("Access-Control-Allow-Origin: *")
		print("Content-Type: text/plain;charset=utf-8")
		print()
		# and activate error handling if required
		debugging = False
		if debugging:
			import cgitb
			cgitb.enable()
	except:
		import traceback
		import argparse
		parser = argparse.ArgumentParser()
		parser.add_argument('-fio2',			default=0.21,		type=float,	help='fraction')
		parser.add_argument('-Hb',				default=150,		type=float,	help='g/l')
		parser.add_argument('-Temp',	 		default=309.65,		type=float,	help='K') # 36.5 C = 309.65 K
		parser.add_argument('-VO2',				default=0.25,		type=float,	help='l/min')
		parser.add_argument('-Q',				default=6.5,		type=float,	help='l/min')
		parser.add_argument('-maxOER',			default=0.8,		type=float,	help='fraction')
		parser.add_argument('-DPG',				default=0.00465,	type=float,	help='M')
		parser.add_argument('-RER',				default=0.8,		type=float,	help='fraction')
		parser.set_defaults()
		args = parser.parse_args()
		variables = vars(args)
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
	# Error
	else:
		return 'Unit conversion error'
#-------------------------------------
if __name__ == "__main__":
	import cgi
	import argparse
	inputs = getinputs()
	print (inputs)














