'''
Calculate the number of photons, number of background photons, and noise equivalent angle for a natural guide star
By Matthew Freeman and Paolo Turri

Equations 65-68 from section 3B of
Clare, R. et al (2006). Adaptive optics sky coverage modelling for extremely large telescopes. Applied Optics 45, 35 (8964-8978)
'''

import numpy as np
import math

wfs_list = ['LBWFS', 'LGSWFS', 'TRICK-H', 'TRICK-K']
wfs = 'TRICK-K'
m = 7.3						#magnitude of NGS


def nea_photons(m,wfs):

	D = 11.8125 					#telescope diameter
	side = 0.5625					#side lenghth of square subaperture
	time = 1/800 					#Time interval (s)
	# r_0 = 0.178 					#Fried parameter. Median 17.8cm at 0.5 microns, from KAON303
									#Median 22.0 cm at 0.5 micron from Raven paper (Yoshito 2017)
									#Carlos' script uses 0.16*cos(parm.atm.zenithAngle )^(3/5); % coherence lenght in meters at 0.5microns
									
	# wavelength = guide star imaging wavelength
	# ps = pixel scale (arcsec).  
	# sigma_e = rms detector read noise per pixel
	# theta_beta = spot size on detector (radians)
	# pix_per_ap = pixels per subaperture, for noise calculation

	if wfs == 'LBWFS':
		band = "R"
		wavelength = 0.641e-6
		ps = 1.5										#1.5 for WFS and low bandwithth WFS from Blake's config
		sigma_e = 0										#from Carlos' config
		theta_beta = 0.49 *(math.pi/180)/(60*60)		#from KAON 1303 Table 16
		throughput = 0.03								#from KAON 1303 Table 8
		pix_per_ap = 4 									#quadcell
	elif wfs == 'LGSWFS':
		# band = "" 
		wavelength = 0.589e-6
		ps = 1.2										#from Carlos' config
		sigma_e = 0.5									#from Carlos' config
		theta_beta = 1.5 *(math.pi/180)/(60*60) 		#from KAON 1303 Table 20
		throughput = 1 									# KAON 1303 Table 8 states 0.36, but Np=1000 is already measured on the detector.	
		pix_per_ap = 4 									#quadcell
	elif wfs == 'TRICK-H':
		band = "H"
		wavelength = 1.63e-6
		ps = 0.06										#from Carlos' config
		sigma_e = 4										#from Carlos' config
		theta_beta = 0.055 *(math.pi/180)/(60*60)		#Using OSIRIS FWHM from KAON 1303 Table 13 (as suggested by Peter)
		throughput = 0.56								#from KAON 1303 Table 8
		pix_per_ap = 4 									#ROI reduces from 16x16 to 2x2 as residual is reduced.
	elif wfs == 'TRICK-K':
		band = "K"
		wavelength = 2.19e-6
		ps = 0.04										#from Carlos' config
		sigma_e = 4										#from Carlos' config
		theta_beta = 0.074 *(math.pi/180)/(60*60)		#Scaling the K band 0.055 by 2.19/1.63 (wavelength ratio)
		throughput = 0.62								#from KAON 1303 Table 8
		pix_per_ap = 4 									#ROI decreases from 16x16 to 2x2 as residual reduces
	else:
		print("Other wfs:", wfs)


	if wfs == 'LGSWFS':
		Np, Nb = 1000, 3	
	else:
		Np, Nb = n_photons(side,time,m,band,ps,throughput)			#number of photons and background photons (instead of eqs 69 and 70)
	# A_sa = side**2 												#area of supaperture
	# N_sa = A_sa/(math.pi*(D/2)**2)								#total number of subapertures for the NGS WFS
	# theta_beta = wavelength/(4*r_0*0.4258) 						#Effective spot size of the subaperture NGS assuming seeing limited image. (eq 67)
	# theta_beta = 3*math.pi*wavelength*np.sqrt(N_sa)/(16*D)		#Effective spot size of the subaperture NGS assuming a diffraction limited core. (eq 68)
	SNR = Np/np.sqrt(Np+pix_per_ap*Nb+pix_per_ap*sigma_e**2)		#signal to noise ratio of a single subaperture (eq 66)
	sigma_theta = theta_beta/SNR  *(180/math.pi) *60*60*1000		#noise equivalent angle in milliarcseconds (eq 65)
	print("\n SNR:                                               {0}".format(
	    SNR))
	print("\n NEA:                                               {0}".format(
	    sigma_theta))


"""Calculate the number of photons from a star and background incident on a
square area in a given time interval.

Author: Paolo Turri

Bibliography:
[1] Bessel et al. (1998)
[2] Mann & von Braun (2015)
[3] https://www.cfht.hawaii.edu/Instruments/ObservatoryManual/
    CFHT_ObservatoryManual_%28Sec_2%29.html
"""

def n_photons(side, time, m ,band, ps, throughput):

	# Parameters
	# side = 0.5625  # Side of the square area (m)
	# time = 1/800  # Time interval (s)
	# m = 7.3 # Apparent magnitude (Vega system)
	# band = "K"  # Band name ("U", "B", "V", "R", "I", "J", "H", "K")
	# ps = 0.04  # Pixel scale (arcsec px^-1)

	# Fixed parameters
	c = 2.99792458e8  # Speed of light (m s^-1)
	h = 6.6260755e-27  # Plank constant (erg s)
	bands = {'name': ["U", "B", "V", "R", "I", "J", "H", "K"],
	         'lambd': [0.366, 0.438, 0.545, 0.641, 0.798, 1.22, 1.63, 2.19],
	         'delta_lambd': [0.0665, 0.1037, 0.0909, 0.1479, 0.1042, 0.3268, 0.2607,
	                         0.5569],
	         'phi_erg': [417.5, 632, 363.1, 217.7, 112.6, 31.47, 11.38, 3.961],
	         'bkg_m': [21.6, 22.3, 21.1, 20.3, 19.2, 14.8, 13.4, 12.6]}
	# Bands' names, effective wavelengths (um), equivalent widths (um), fluxes
	# (10^-11 erg s^-1 cm^-2 A^-1) and background magnitudes (arcsec^-2) [1, 2, 3].

	# Get band's data
	band_idx = np.where(np.array(bands['name']) == band)[0][0]
	lambd = float(bands['lambd'][band_idx])  # Band effective wavelength (um)
	delta_lamb = float(bands['delta_lambd'][band_idx])  # Band equivalent width (um)
	phi_erg = float(bands['phi_erg'][band_idx])  # Flux (erg s^-1 cm^-2 A^-1)
	bkg_m = float(bands['bkg_m'][band_idx])  # Background magnitude (arcsec^-2)

	# Calculations of the star
	f = c / (lambd * 1e-6)  # Band frequency (s^-1)
	phi_n = phi_erg * 1e-11 / (h * f)  # Numeric flux (s^-1 cm^-2 A^-1)
	n_ph_0 = phi_n * ((side * 1e2) ** 2) * time * delta_lamb * 1e4 * throughput  # Zeropoint (m = 0) number of photons on detector
	n_ph_star = n_ph_0 * (10 ** (-0.4 * m))  # Number of star photons

	# Calculations of the background
	n_ph_bkg = n_ph_0 * (10 ** (-0.4 * bkg_m)) * (ps ** 2)  # Number of background
	# photons (px^-1)

	# Show results
	print("\n Number of photons from the star:                   {0}".format(
	    n_ph_star))
	print("\n Number of photons per pixel from the background:   {0}".format(
	    n_ph_bkg))

	# return {"n_ph_star":n_ph_star, "n_ph_bkg":n_ph_bkg}
	return n_ph_star, n_ph_bkg

nea_photons(m,wfs)
