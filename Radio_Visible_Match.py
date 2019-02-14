import numpy as np
import astropy
from astropy.io import ascii,fits
import math
from time import gmtime, strftime
import sys
from sys import exit,argv
import csv
import os
from astropy import wcs
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from string import upper,lower
from numpy import recfromcsv
import matplotlib.pyplot as plt
from numpy.random import normal
import matplotlib.cm as cm
from astropy.coordinates import SkyCoord, search_around_sky, Angle
from astropy.nddata import Cutout2D
from astropy import units as u
import pylab
from astropy import stats
from urllib2 import urlopen
from astroquery.sdss import SDSS

'''
def PLOT_TRIPLE_MATCHES(Radio_Path, Visible_Path, Visible_Sources, Visible_Image, Radio_Sources1, Radio_Sources2):
	hdulist_Radio1  =  fits.open(Radio_Path + Radio_Sources1)
	cat_Radio1  =  hdulist_Radio1[1].data
	cols_Radio1  =  hdulist_Radio1[1].columns
	hdulist_Radio1.close()
	RA_Radio1 = cat_Radio1['RA']
	DEC_Radio1  =  cat_Radio1['DEC']

	hdulist_Radio2  =  fits.open(Radio_Path + Radio_Sources2)
	cat_Radio2  =  hdulist_Radio2[1].data
	cols_Radio2  =  hdulist_Radio2[1].columns
	hdulist_Radio2.close()
	RA_Radio2 = cat_Radio2['RA']
	DEC_Radio2  =  cat_Radio2['DEC']

	hdulist_Visible  =  fits.open(Visible_Path + Visible_Image)
	Visible_w  =  wcs.WCS(hdulist_Visible[0].header)
	image_data  =  hdulist_Visible[0].data
	hdulist_Visible.close()

	hdulist_Visible  = fits.open(Visible_Path + Visible_Sources)
	cat_Visible = hdulist_Visible[1].data
	Cols_Visible = hdulist_Visible[1].columns
	hdulist_Visible.close()
	RA_Visible = cat_Visible['Ra']
	DEC_Visible = cat_Visible['Dec']

	TILE_RA_LOWER_LIMIT, TILE_DEC_LOWER_LIMIT =  Visible_w.wcs_pix2world(1.,1.,1)
	TILE_RA_UPPER_LIMIT, TILE_DEC_UPPER_LIMIT =  Visible_w.wcs_pix2world(11000.,11000.,1)

	Visible_coord = SkyCoord(RA_Visible,DEC_Visible, unit='deg')
	Radio1_coord = SkyCoord(RA_Radio1,DEC_Radio1, unit='deg')
	Radio2_coord = SkyCoord(RA_Radio2,DEC_Radio2, unit='deg')


	match1, match2,sep1, dist1= search_around_sky(Visible_coord,Radio1_coord,Angle('00d00m2.5s'))

	RA_Radio = RA_Radio1[match2]
	DEC_Radio = DEC_Radio1[match2]
	Radio_coord = SkyCoord(RA_Radio,DEC_Radio, unit='deg')

	match1, match2,sep1, dist1= search_around_sky(Radio2_coord,Radio_coord,Angle('00d00m2.5s'))

	NEW_RA_RADIO1 = RA_Radio[match2]
	NEW_DEC_RADIO1 = DEC_Radio[match2]
	NEW_RA_RADIO2 = RA_Radio2[match1]
	NEW_DEC_RADIO2 = DEC_Radio2[match1]

	#print len(NEW_RA_RADIO2)
	#ax1 = plt.subplot(1,2,1)
	#ax1.plot(NEW_RA_RADIO1, NEW_DEC_RADIO1, linestyle="None", marker='o')
	#ax2 = plt.subplot(1,2,2)
	#ax2.plot(NEW_RA_RADIO2, NEW_DEC_RADIO2, linestyle="None", marker='o')
	#plt.show()
	#for i in range(len(NEW_RA_RADIO1)):
	#	if TILE_RA_LOWER_LIMIT > NEW_RA_RADIO1[i] > TILE_RA_UPPER_LIMIT and TILE_DEC_LOWER_LIMIT < NEW_DEC_RADIO1[i] < TILE_DEC_UPPER_LIMIT:
	#		

'''
def RADIO_CONTOURS(RA, DEC, FILE):

	filename = Radio_Path + 'FIRST_' + np.string_(FILE)
	Radio_Source_hdulist = fits.open(filename)
	Radio_Source_Image = Radio_Source_hdulist[0].data
	Radio_WCS = wcs.WCS(Radio_Source_hdulist[0].header)
	Radio_Source_hdulist.close()
		
	Radio_WCS_N = Radio_WCS.dropaxis(2)
	N_Radio_WCS = Radio_WCS_N.dropaxis(2)
	Radio_Pixel_Scale = wcs.utils.proj_plane_pixel_scales(N_Radio_WCS) * 3600.
	print 'Radio image pixel scale: ', Radio_Pixel_Scale, ' arcsec arcsec'
	Radio_Npixels_Needed = np.ceil(Desired_Image_Size/Radio_Pixel_Scale[0])
	print Radio_Npixels_Needed
	size = (Radio_Npixels_Needed,Radio_Npixels_Needed)     # pixels
	px, py = N_Radio_WCS.wcs_world2pix(RA, DEC, 1)
	position = (px, py)
		
	Cutout_Radio2 = Cutout2D(Radio_Source_Image, position, size, wcs=N_Radio_WCS)
	Contour_Color = 'b'
	return Cutout_Radio2, Contour_Color

def RADIO_CATALOGUES_COORDINATES (Cutout_Radio, Radio_Npixels_Needed, Radio_Path, Radio_Catalogue):


	hdulist_Radio  =  fits.open(Radio_Path + Radio_Catalogue)
	cat_Radio  =  hdulist_Radio[1].data
	cols_Radio  =  hdulist_Radio[1].columns
	hdulist_Radio.close()
	if Radio_Catalogue == 'Heywood_Stripe82(VLA)(1).fit' or Radio_Catalogue == 'Heywood_Stripe82(VLA)(2).fit':
		RA_Radio  =  cat_Radio['RAJ2000']
		DEC_Radio  =  cat_Radio['DEJ2000']
	else:
		RA_Radio  =  cat_Radio['RA']
		DEC_Radio  =  cat_Radio['DEC']
	#print RA_Radio,DEC_Radio
	

	TILE_RA_LOWER_LIMIT, TILE_DEC_LOWER_LIMIT =  Cutout_Radio.wcs.wcs_pix2world(1.,1.,1)
	TILE_RA_UPPER_LIMIT, TILE_DEC_UPPER_LIMIT =  Cutout_Radio.wcs.wcs_pix2world(Radio_Npixels_Needed,Radio_Npixels_Needed,1)

			
	X_MATCH_RA_RADIO = []
	X_MATCH_DEC_RADIO = []
		
	

	for n in range(len(RA_Radio)):

		#print TILE_RA_LOWER_LIMIT,"(",RA_Radio[n],')',TILE_RA_UPPER_LIMIT,'========>',TILE_DEC_LOWER_LIMIT,'(',DEC_Radio[n],')',TILE_DEC_UPPER_LIMIT
		if TILE_RA_LOWER_LIMIT > RA_Radio[n] > TILE_RA_UPPER_LIMIT and TILE_DEC_LOWER_LIMIT < DEC_Radio[n] < TILE_DEC_UPPER_LIMIT:
			
			pxx, pyy = Cutout_Radio.wcs.wcs_world2pix(RA_Radio[n], DEC_Radio[n], 1)
			X_MATCH_RA_RADIO.append(pxx)
			X_MATCH_DEC_RADIO.append(pyy)

	
	return X_MATCH_RA_RADIO, X_MATCH_DEC_RADIO

def RADIO_VISIBLE_FIELD_MATCH (Radio_Path,Radio_Sources, Visible_Path, Visible_Image, NAME):

	### Takes as entries a radio source catalog and a visible image (tile). 
	### Finds the radio sources within the tiles and returns its coordinates and a file name, along with a table with these information


	FLAG = 0
	hdulist_Radio  =  fits.open(Radio_Path + Radio_Sources)
	cat_Radio  =  hdulist_Radio[1].data
	cols_Radio  =  hdulist_Radio[1].columns
	hdulist_Radio.close()
	RA_Radio  =  cat_Radio['RA']
	DEC_Radio  =  cat_Radio['DEC']
	
	hdulist_Visible  =  fits.open(Visible_Path + Visible_Image)
	Visible_w  =  wcs.WCS(hdulist_Visible[0].header)
	image_data  =  hdulist_Visible[0].data
	hdulist_Visible.close()

	
	TILE_RA_LOWER_LIMIT, TILE_DEC_LOWER_LIMIT =  Visible_w.wcs_pix2world(1.,1.,1)
	TILE_RA_UPPER_LIMIT, TILE_DEC_UPPER_LIMIT =  Visible_w.wcs_pix2world(11000.,11000.,1)
	
	name  =  NAME+'_Radio_Visible_Field_Matches.csv'
	file  =  open(os.path.join(Main_Path,name), 'w')
	
	file.write('RA,DEC,RA_HMS,DEC_DMS,FILE')
	file.write('\n')
	a = 0

	RA_HMS = []
	DEC_DMS = []
	FILE = []
	MATCH_RA_RADIO = []
	MATCH_DEC_RADIO = []

	for i in range(len(RA_Radio)):
		if TILE_RA_LOWER_LIMIT > RA_Radio[i] > TILE_RA_UPPER_LIMIT and TILE_DEC_LOWER_LIMIT < DEC_Radio[i] < TILE_DEC_UPPER_LIMIT:
			name  = NAME+'_%s+%s.fits' %(RA_Radio[i], DEC_Radio[i])
			name = name.replace('+-','-')
			a += 1
			Coord = '%s,%s,%s,%s,%s' %(RA_Radio[i],DEC_Radio[i],Angle(RA_Radio[i],u.degree).to_string(unit = u.hour, sep = ':'),Angle(DEC_Radio[i],u.degree).to_string(unit = u.degree, sep = ':'),name)
		
			file.write(Coord)

			file.write('\n')



			RA_HMS.append(Angle(RA_Radio[i],u.degree).to_string(unit = u.hour, sep = ':'))
			DEC_DMS.append(Angle(DEC_Radio[i],u.degree).to_string(unit = u.degree, sep = ':'))
			FILE.append(name)
			MATCH_RA_RADIO.append(RA_Radio[i])
			MATCH_DEC_RADIO.append(DEC_Radio[i])

	file.close()
	return RA_HMS, DEC_DMS, FILE,MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, TILE_RA_LOWER_LIMIT, TILE_DEC_LOWER_LIMIT, TILE_RA_UPPER_LIMIT, TILE_DEC_UPPER_LIMIT

def RADIO_VISIBLE_SOURCE_MATCH (Radio_Path, Radio_Sources, Visible_Path, Visible_Sources, Visible_Image, NAME):

	FLAG = 1

	hdulist_Radio  =  fits.open(Radio_Path + Radio_Sources)
	cat_Radio  =  hdulist_Radio[1].data
	cols_Radio  =  hdulist_Radio[1].columns
	hdulist_Radio.close()
	RA_Radio  =  cat_Radio['RA']
	DEC_Radio  =  cat_Radio['DEC']
	
	if NAME == 'HODGE':
		FPEAK_Radio = cat_Radio['FLUX_1P4_GHZ']
		FINT_Radio = cat_Radio['INT_FLUX_1P4_GHZ']
	if NAME == 'FIRST':
		FPEAK_Radio = cat_Radio['FPEAK']
		FINT_Radio = cat_Radio['FINT']
	
	hdulist_Visible  = fits.open(Visible_Path + Visible_Sources)
	cat_Visible = hdulist_Visible[1].data
	Cols_Visible = hdulist_Visible[1].columns
	hdulist_Visible.close()
	RA_Visible = cat_Visible['Ra']
	DEC_Visible = cat_Visible['Dec']
	Photo_z_Visible = cat_Visible['zb']

	hdulist_Visible  =  fits.open(Visible_Path + Visible_Image)
	Visible_w  =  wcs.WCS(hdulist_Visible[0].header)
	image_data  =  hdulist_Visible[0].data

	hdulist_Visible.close()
	TILE_RA_LOWER_LIMIT, TILE_DEC_LOWER_LIMIT =  Visible_w.wcs_pix2world(1.,1.,1)
	TILE_RA_UPPER_LIMIT, TILE_DEC_UPPER_LIMIT =  Visible_w.wcs_pix2world(11000.,11000.,1)

	Visible_coord = SkyCoord(RA_Visible,DEC_Visible, unit='deg')
	Radio_coord = SkyCoord(RA_Radio,DEC_Radio, unit='deg')


	match1, match2,sep1, dist1= search_around_sky(Visible_coord,Radio_coord,Angle('00d00m2.5s'))

	NEW_RA_RADIO = RA_Radio[match2]
	NEW_DEC_RADIO = DEC_Radio[match2]
	NEW_FPEAK_RADIO = FPEAK_Radio[match2]
	NEW_FINT_RADIO = FINT_Radio[match2]
	NEW_Photo_z_Visible = Photo_z_Visible[match2]
	NEW_RA_VISIBLE = RA_Visible[match1]
	NEW_DEC_VISIBLE = DEC_Visible[match1]

	#print len(NEW_RA_RADIO)

	FILE = []
	MATCH_RA_RADIO = []
	MATCH_DEC_RADIO = []
	MATCH_FPEAK_RADIO = []
	MATCH_FINT_RADIO = []
	MATCH_Photo_z_Visible = []
	MATCH_RA_VISIBLE = []
	MATCH_DEC_VISIBLE = []

	name  =  NAME+'_Radio_Visible_Source_Matches.csv'
	file  =  open(os.path.join(Main_Path,name), 'w')
	
	file.write('RA,DEC,RA_HMS,DEC_DMS,FILE')
	file.write('\n')
	a = 0

	RA_HMS = []
	DEC_DMS = []
	FILE = []
	MATCH_RA_RADIO = []
	MATCH_DEC_RADIO = []

	for i in range(len(NEW_RA_RADIO)):
		if TILE_RA_LOWER_LIMIT > NEW_RA_RADIO[i] > TILE_RA_UPPER_LIMIT and TILE_DEC_LOWER_LIMIT < NEW_DEC_RADIO[i] < TILE_DEC_UPPER_LIMIT:
			name  = NAME+'_%s+%s.fits' %(NEW_RA_RADIO[i], NEW_DEC_RADIO[i])
			name = name.replace('+-','-')
			Coord = '%s,%s,%s,%s,%s' %(NEW_RA_RADIO[i],NEW_DEC_RADIO[i],Angle(NEW_RA_RADIO[i],u.degree).to_string(unit = u.hour, sep = ':'),Angle(NEW_DEC_RADIO[i],u.degree).to_string(unit = u.degree, sep = ':'),name)
		
			file.write(Coord)

			file.write('\n')



			RA_HMS.append(Angle(NEW_RA_RADIO[i],u.degree).to_string(unit = u.hour, sep = ':'))
			DEC_DMS.append(Angle(NEW_DEC_RADIO[i],u.degree).to_string(unit = u.degree, sep = ':'))
			FILE.append(name)
			MATCH_RA_RADIO.append(NEW_RA_RADIO[i])
			MATCH_DEC_RADIO.append(NEW_DEC_RADIO[i])
			MATCH_FPEAK_RADIO.append(NEW_FPEAK_RADIO[i])
			MATCH_FINT_RADIO.append(NEW_FINT_RADIO[i])
			MATCH_Photo_z_Visible.append(NEW_Photo_z_Visible[i])
			#MATCH_RA_VISIBLE.append(NEW_RA_VISIBLE[i])
			#MATCH_DEC_RADIO.append(NEW_DEC_VISIBLE[i])

			
	file.close()

	return FILE, MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, MATCH_FPEAK_RADIO, MATCH_FINT_RADIO, MATCH_Photo_z_Visible, RA_HMS, DEC_DMS



def SDSS_QUERY_SPEC_Z(TILE_RA_LOWER_LIMIT, TILE_DEC_LOWER_LIMIT, TILE_RA_UPPER_LIMIT, TILE_DEC_UPPER_LIMIT):

	query = "SELECT objID, ra, dec,z FROM SpecPhoto WHERE (ra between %s and %s) and (dec between %s and %s)" %(TILE_RA_UPPER_LIMIT,TILE_RA_LOWER_LIMIT,TILE_DEC_LOWER_LIMIT, TILE_DEC_UPPER_LIMIT)	
	data = SDSS.query_sql(query, data_release = 15)

	name  =  'SDSS_Tile_Sources.csv'
	file  =  open(os.path.join(Main_Path,name), 'w')

	file.write('specobjID,ra,dec,z')
	file.write('\n')

	specobjID = data['objID']
	ra = data['ra']
	dec = data['dec']
	z = data['z']
	for i in range(len(z)):
		info = '%s,%s,%s,%s' %(specobjID[i],ra[i],dec[i],z[i])
		file.write(info)
		file.write('\n')
	file.close()

	print data

#def DOWNLOAD_RADIO_SOURCES (MATCH_RA_RADIO, MATCH_DEC_RADIO, RA_HMS, DEC_DMS, FILE, NAME):
#
#	### Downloads the images of the the radio catalog (FIRST catalogue. Must be changed for other catalogues) sources that are located in the visible tile
#
#	for i in range(len(RA_HMS)):
#
#		N_RA = RA_HMS[i].replace(':','%20') + '%20'
#		#print N_RA
#		E_DEC = DEC_DMS[i].replace(':','%20')
#		D_DEC = '%2B' + E_DEC
#		N_DEC = D_DEC.replace('%2B-','%2D')
#		#print N_DEC
#		if NAME == 'HODGE':
#			LINK ='https://third.ucllnl.org/cgi-bin/stripe82image?RA=%s%s&Dec=&Equinox=J2000&ImageSize=4.5&MaxInt=3&Survey=stripe82&FITS=1&Download=1' %(N_RA,N_DEC)
#			name = NAME+'_%s+%s.fits' %(MATCH_RA_RADIO[i], MATCH_DEC_RADIO[i])
#			name = name.replace('+-','-')
#			if not os.path.isfile(Radio_Path+name):
#				data = urlopen(LINK)
#				print "FILE!"
#				with open(os.path.join(Radio_Path ,name), "wb") as output:
#				    output.write(data.read())
#			else:
#				print "NOPE"
#
#		if NAME == 'FIRST':
#			if Angle('22h06m00s').degree < MATCH_RA_RADIO[i] < Angle('23h21m00s').degree and Angle('-01d21m00s').degree < MATCH_DEC_RADIO[i] <  Angle('+01d21m00s').degree: 
#					
#				LINK1 ='https://third.ucllnl.org/cgi-bin/stripe82image?RA=%s%s&Dec=&Equinox=J2000&ImageSize=4.5&MaxInt=3&Survey=stripe82&FITS=1&Download=1' %(N_RA,N_DEC)
#				name = 'HODGE'+'_'+NAME+'_%s+%s.fits' %(MATCH_RA_RADIO[i], MATCH_DEC_RADIO[i])
#				name = name.replace('+-','-')
#				if not os.path.isfile(Radio_Path+name):
#					data = urlopen(LINK1)
#					print "FILE!"
#					with open(os.path.join(Radio_Path ,name), "wb") as output:
#					    output.write(data.read())
#				else:
#					print "NOPE"
#
#			elif Angle('00h39m00s').degree < MATCH_RA_RADIO[i] < Angle('02h02m21s').degree and Angle('-01d21m00s').degree < MATCH_DEC_RADIO[i] < Angle('+01d21m00s').degree:
#
#				LINK1 ='https://third.ucllnl.org/cgi-bin/stripe82image?RA=%s%s&Dec=&Equinox=J2000&ImageSize=4.5&MaxInt=3&Survey=stripe82&FITS=1&Download=1' %(N_RA,N_DEC)
#				name = 'HODGE'+'_'+NAME+'_%s+%s.fits' %(MATCH_RA_RADIO[i], MATCH_DEC_RADIO[i])
#				name = name.replace('+-','-')
#				if not os.path.isfile(Radio_Path+name):
#					data = urlopen(LINK1)
#					print "FILE!"
#					with open(os.path.join(Radio_Path ,name), "wb") as output:
#					    output.write(data.read())
#				else:
#					print "NOPE"
#
#			LINK ='https://third.ucllnl.org//cgi-bin/firstimage?RA=%s%s&amp;Dec=&amp;Equinox=J2000&amp;ImageSize=4.5&amp;MaxInt=10&amp;FITS=1&amp;Download=1' %(N_RA,N_DEC)
#			
#			name = NAME+'_%s+%s.fits' %(MATCH_RA_RADIO[i], MATCH_DEC_RADIO[i])
#			name = name.replace('+-','-')
#			if not os.path.isfile(Radio_Path+name):
#				data = urlopen(LINK)
#				print "FILE!"
#				with open(os.path.join(Radio_Path ,name), "wb") as output:
#				    output.write(data.read())
#			else:
#				print "NOPE"


def DOWNLOAD_RADIO_SOURCES (MATCH_RA_RADIO, MATCH_DEC_RADIO, RA_HMS, DEC_DMS, FILE, NAME):

	### Downloads the images of the the radio catalog (FIRST catalogue. Must be changed for other catalogues) sources that are located in the visible tile

	for i in range(len(RA_HMS)):

		N_RA = RA_HMS[i].replace(':','%20') + '%20'
		#print N_RA
		E_DEC = DEC_DMS[i].replace(':','%20')
		D_DEC = '%2B' + E_DEC
		N_DEC = D_DEC.replace('%2B-','%2D')
		#print N_DEC
		if NAME == 'HODGE':
			
			LINK ='https://third.ucllnl.org/cgi-bin/stripe82image?RA=%s%s&Dec=&Equinox=J2000&ImageSize=4.5&MaxInt=3&Survey=stripe82&FITS=1&Download=1' %(N_RA,N_DEC)
			name = NAME+'_%s+%s.fits' %(MATCH_RA_RADIO[i], MATCH_DEC_RADIO[i])
			name = name.replace('+-','-')
			if not os.path.isfile(Radio_Path+name):
				data = urlopen(LINK)
				print "FILE!"
				with open(os.path.join(Radio_Path ,name), "wb") as output:
				    output.write(data.read())
			else:
				print "NOPE"

			LINK1 ='https://third.ucllnl.org//cgi-bin/firstimage?RA=%s%s&amp;Dec=&amp;Equinox=J2000&amp;ImageSize=4.5&amp;MaxInt=10&amp;FITS=1&amp;Download=1' %(N_RA,N_DEC)
			name = 'FIRST_'+NAME+'_%s+%s.fits' %(MATCH_RA_RADIO[i], MATCH_DEC_RADIO[i])
			name = name.replace('+-','-')
			if not os.path.isfile(Radio_Path+name):
				data = urlopen(LINK1)
				print "FILE!"
				with open(os.path.join(Radio_Path ,name), "wb") as output:
				    output.write(data.read())
			else:
				print "NOPE"	
		if NAME == 'FIRST':

			if Angle('22h06m00s').degree < MATCH_RA_RADIO[i] < Angle('23h21m00s').degree and  Angle('-01d21m00s').degree < MATCH_DEC_RADIO[i] <  Angle('+01d21m00s').degree:
				print "no" 
			elif Angle('00h39m00s').degree < MATCH_RA_RADIO[i] < Angle('02h02m21s').degree and Angle('-01d21m00s').degree < MATCH_DEC_RADIO[i] < Angle('+01d21m00s').degree:
				print "no"
			else:

				LINK1 ='https://third.ucllnl.org//cgi-bin/firstimage?RA=%s%s&amp;Dec=&amp;Equinox=J2000&amp;ImageSize=4.5&amp;MaxInt=10&amp;FITS=1&amp;Download=1' %(N_RA,N_DEC)
				name = NAME+'_%s+%s.fits' %(MATCH_RA_RADIO[i], MATCH_DEC_RADIO[i])
				name = name.replace('+-','-')
				if not os.path.isfile(Radio_Path+name):
					data = urlopen(LINK1)
					print "FILE!"
					with open(os.path.join(Radio_Path ,name), "wb") as output:
					    output.write(data.read())
				else:
					print "NOPE"	
	
		
def PLOT_MATCHES (FILE, MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, FPEAK, FINT, Photo_z, NAME):
	
	### Plots Radio images and Visible fields side by side in squares of the same area of the sky

	Source_match_spec_z = 'SDSS_Tile_Sources.csv'
	data = ascii.read(Main_Path + Source_match_spec_z)
	Spec_RA = data['ra']
	Spec_DEC = data['dec']
	Spec_Z = data['z']
	Spec_ID = data['specobjID']


	if not os.path.isdir(Main_Path + 'PLOT_MATCHES/'):
		os.mkdir(Main_Path + 'PLOT_MATCHES/')

	###name  =  'SDSS_SPEC_Z_MATCH.csv'
	###file  =  open(Main_Path + 'PLOT_SOURCE_MATCHES/' + name, 'a')
	#file.write('ID,RA_Radio,DEC_Radio,RA_SDSS,DEC_SDSS,FPEAK,FINT,SPLUS_Photo_z,SDSS_Spec_redshift')
	#file.write('\n')
	
	
	for i in range(len(FILE)):
 		print i
		pname=np.string_(FILE[i]) + '.png'
		filename = Radio_Path + np.string_(FILE[i]) # here
		

		print 'plotting' + filename
		position = (MATCH_RA_RADIO[i], MATCH_DEC_RADIO[i])
		px, py = Visible_w.wcs_world2pix(MATCH_RA_RADIO[i], MATCH_DEC_RADIO[i], 1)
		position = (px, py)
		Visible_Pixel_Scale = wcs.utils.proj_plane_pixel_scales(Visible_w) * 3600.

		print 'Visible image pixel scale: ', Visible_Pixel_Scale, ' arcsec arcsec'
		Visible_Npixels_Needed = np.ceil(Desired_Image_Size/Visible_Pixel_Scale[0])
		size = (Visible_Npixels_Needed,Visible_Npixels_Needed)     # pixels
		Cutout_Visible = Cutout2D(image_data, position, size, wcs = Visible_w)
			

		if np.mean(Cutout_Visible.data > 0.0) and Cutout_Visible.data[Visible_Npixels_Needed/2][Visible_Npixels_Needed/2] != 0:
			
			Sigma_Visible = stats.sigma_clip(Cutout_Visible.data,sigma=3,iters=5)
			Sigma_Visible = np.std(Sigma_Visible)
			
			Mean_Visible = np.median(Cutout_Visible.data)
			
			#plt.contourf(Cutout_Visible.data, cmap='binary', alpha=0.5)
			
			if FLAG == 1:
				
				Coord_info = ' RA = %s  \n DEC = %s ' %(MATCH_RA_RADIO[i], MATCH_DEC_RADIO[i])
				plt.text(0.5, 0.5, Coord_info)
				#print 'angle',Angle('23h21m00s').degree
				#exit(0)
				
			
			
			Radio_Source_hdulist = fits.open(filename)
			Radio_Source_Image = Radio_Source_hdulist[0].data
			Radio_WCS = wcs.WCS(Radio_Source_hdulist[0].header)
			Radio_Source_hdulist.close()

			

			Radio_WCS_N = Radio_WCS.dropaxis(2)
			N_Radio_WCS = Radio_WCS_N.dropaxis(2)

			Radio_Pixel_Scale = wcs.utils.proj_plane_pixel_scales(N_Radio_WCS) * 3600.
			print 'Radio image pixel scale: ', Radio_Pixel_Scale, ' arcsec arcsec'
			Radio_Npixels_Needed = np.ceil(Desired_Image_Size/Radio_Pixel_Scale[0])
			print Radio_Npixels_Needed
			size = (Radio_Npixels_Needed,Radio_Npixels_Needed)     # pixels
			px, py = N_Radio_WCS.wcs_world2pix(MATCH_RA_RADIO[i], MATCH_DEC_RADIO[i], 1)
			position = (px, py)
			
			Cutout_Radio = Cutout2D(Radio_Source_Image, position, size, wcs=N_Radio_WCS)

			Sigma_Radio = np.std(Cutout_Radio.data)
			Mean_Radio = np.median(Cutout_Radio.data)


			if NAME == 'HODGE':
				Cutout_Radio2, Contour_Color = RADIO_CONTOURS(MATCH_RA_RADIO[i], MATCH_DEC_RADIO[i], FILE[i])
				Color = 'g'
				ax1 = plt.subplot(2,2,1, projection = Cutout_Visible.wcs)
				ax2 = plt.subplot(2,2,3, projection = Cutout_Visible.wcs)
				ax3 = plt.subplot(2,2,2, projection=Cutout_Radio.wcs) 
				ax4 =plt.subplot(2,2,4, projection=Cutout_Radio2.wcs)

				ax1.set_yticklabels([])
				ax1.set_xticklabels([])
				ax2.set_yticklabels([])
				ax2.set_xticklabels([])
				ax3.set_yticklabels([])
				ax3.set_xticklabels([])
				ax4.set_yticklabels([])
				ax4.set_xticklabels([])

				ax2.contour(Cutout_Radio2.data, levels = [0.001,0.002,0.004,0.008,0.016], colors=Contour_Color, alpha=0.5, linewidths = 0.8, transform=ax2.get_transform(Cutout_Radio2.wcs))
				
				ax4.imshow(Cutout_Radio2.data,origin='lower',cmap='binary',vmin=Mean_Radio-Sigma_Radio, vmax=Mean_Radio+10*Sigma_Radio, zorder=0)
			else:
				Color = 'b'
				ax1 = plt.subplot(2,1,1, projection = Cutout_Visible.wcs)
				ax2 = plt.subplot(2,1,3, projection = Cutout_Visible.wcs)
				ax3 = plt.subplot(2,1,2, projection=Cutout_Radio.wcs) 
				ax1.set_yticklabels([])
				ax1.set_xticklabels([])
				ax2.set_yticklabels([])
				ax2.set_xticklabels([])
				ax3.set_yticklabels([])
				ax3.set_xticklabels([])
			if FLAG == 1:
				
				Coord_info = ' RA = %s  \n DEC = %s ' %(MATCH_RA_RADIO[i], MATCH_DEC_RADIO[i])
				ax2.text(0.5, 0.5, Coord_info, fontsize = 8)

			X_MATCH_RA_HODGE,X_MATCH_DEC_HODGE = RADIO_CATALOGUES_COORDINATES(Cutout_Visible, Visible_Npixels_Needed, Radio_Path, HODGE_Sources)
			X_MATCH_RA_FIRST,X_MATCH_DEC_FIRST = RADIO_CATALOGUES_COORDINATES(Cutout_Visible, Visible_Npixels_Needed, Radio_Path, FIRST_Sources)
			X_MATCH_RA_HEYWOOD_E, X_MATCH_DEC_HEYWOOD_E = RADIO_CATALOGUES_COORDINATES(Cutout_Visible, Visible_Npixels_Needed, Radio_Path, HEYWOOD_Sources_Eastern)
			X_MATCH_RA_HEYWOOD_W, X_MATCH_DEC_HEYWOOD_W = RADIO_CATALOGUES_COORDINATES(Cutout_Visible, Visible_Npixels_Needed, Radio_Path, HEYWOOD_Sources_Western)
			#print "---------------------------------",X_MATCH_RA_HODGE,X_MATCH_DEC_HODGE


			
			ax1.set_title('Visible Field')
			ax1.imshow(Cutout_Visible.data,vmin = Mean_Visible-Sigma_Visible, vmax=Mean_Visible+5*Sigma_Visible,origin='lower',cmap='binary')
		
			
			ax1.plot(X_MATCH_RA_FIRST,X_MATCH_DEC_FIRST, linestyle = 'None', marker='+', color='r', label="FIRST", alpha=0.5)
			ax1.plot(X_MATCH_RA_HODGE,X_MATCH_DEC_HODGE, linestyle = 'None', marker='.', color='b', label="HODGE+11", alpha=0.5)
			ax1.plot(X_MATCH_RA_HEYWOOD_E, X_MATCH_DEC_HEYWOOD_E, linestyle = 'None', marker='x', color='g', label="HEYWOOD+16E", alpha=0.5)
			ax1.plot(X_MATCH_RA_HEYWOOD_W, X_MATCH_DEC_HEYWOOD_W, linestyle = 'None', marker='x', color='m', label="HEYWOOD+16W", alpha=0.5)
			ax1.legend(fontsize=4)
			ax2.imshow(Cutout_Visible.data,vmin = Mean_Visible-Sigma_Visible, vmax=Mean_Visible+5*Sigma_Visible,origin='lower',cmap='binary')
			ax2.contour(Cutout_Radio.data, levels = [0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128, 0.256], colors= Color, alpha=0.5, linewidths = 0.8, transform=ax2.get_transform(Cutout_Radio.wcs))
			#ax2 = plt.subplot(1,2,2, projection=Cutout_Radio.wcs) 
			ax3.set_title('Radio Source')
			ax3.imshow(Cutout_Radio.data,origin='lower',cmap='binary',vmin=Mean_Radio-Sigma_Radio, vmax=Mean_Radio+10*Sigma_Radio, zorder=0)
			#ax4.imshow(Cutout_Radio2.data,origin='lower',cmap='binary',vmin=Mean_Radio-Sigma_Radio, vmax=Mean_Radio+10*Sigma_Radio, zorder=0)
			#plt.subplots_adjust(top=0.969,bottom=0.031,left=0.138,right=0.95,hspace=0.2,wspace=0.516)
			
			#plt.yticks([])
			#plt.xticks([])
			
			#plt.tight_layout()
			#plt.legend(fontsize=6)
			if FLAG == 1:


				#name  =  'SDSS_SPEC_Z_MATCH.csv'
				#file  =  open(os.path.join(Main_Path + 'PLOT_SOURCE_MATCHES/',name), 'w')
				#file.write('ID,RA_Radio,DEC_Radio,RA_SDSS,DEC_SDSS,FPEAK,FINT,SPLUS_Photo_z,SDSS_Spec_redshift')
				#file.write('\n')

				Flux_info = ' FPEAK = %s mJy \n FINT = %s mJy \n Photo-z = %s' %(FPEAK[i], FINT[i], Photo_z[i])
				for f in range(len(Spec_Z)):
					
					if np.around(MATCH_RA_RADIO[i],decimals=4) == np.around(Spec_RA[f],decimals=4) and np.around(MATCH_DEC_RADIO[i],decimals=4) == np.around(Spec_DEC[f],decimals=4) :
						 
						Spec_redshift=Spec_Z[f]
						Flux_info = ' ID = %s \n FPEAK = %s mJy \n FINT = %s mJy \n Photo-z = %s \n Spec-z = %s' %(Spec_ID[f],FPEAK[i], FINT[i], Photo_z[i], Spec_redshift)
						table_info = '%s,%s,%s,%s,%s,%s,%s,%s,%s' %(Spec_ID[f],MATCH_RA_RADIO[i],MATCH_DEC_RADIO[i],Spec_RA[f],Spec_DEC[f], FPEAK[i], FINT[i], Photo_z[i], Spec_redshift)
						
						
						###file.write(table_info)
						###file.write('\n')
						break
					#	#plt.text(0.5, 0.5, Flux_info)
					#	
					#else:
					#	Spec_redshift = None	
					#	#Flux_info = ' FPEAK = %s mJy \n FINT = %s mJy \n Photo-z = %s \n Spec-z = %s' %(FPEAK[i], FINT[i], Photo_z[i], Spec_redshift)
					#file.write(table_info)
					#file.write('\n')
				#print Flux_info
				#print table_info
				ax3.text(0.5, 0.5, Flux_info, fontsize = 8)
				pname = '*MATCH*' + np.string_(FILE[i]) + '.png'
				plt.savefig(os.path.join(Main_Path + 'PLOT_SOURCE_MATCHES/',pname))
			
			else:	
				plt.savefig(os.path.join(Main_Path + 'PLOT_MATCHES/',pname)) 
			#plt.show()

	
def main():

	global Desired_Image_Size, Main_Path, Radio_Path, FIRST_Sources,HODGE_Sources, HEYWOOD_Sources_Eastern, HEYWOOD_Sources_Western, Visible_Path, Visible_Sources, Visible_Image, Tiles_List, Tiles_Table


	Name  =  'rodrigo'
	
	if Name  ==  'rodrigo':
		Main_Path  =  '/home/rodrigo/Documentos/Project/Images_Match/'
	if Name  ==  'roderik':
		Main_Path  =  '' #enter your path here

	Desired_Image_Size = 120.
	
	Radio_Path  =  Main_Path + 'Radio_Path/'
	FIRST_Sources  =  'first_14dec17.fits'
	HODGE_Sources = 'Hodge_Stripe82(VLA).fits'
	HEYWOOD_Sources_Eastern = 'Heywood_Stripe82(VLA)(1).fit'
	HEYWOOD_Sources_Western = 'Heywood_Stripe82(VLA)(2).fit'

	Visible_Path  =  Main_Path + 'Visible_Path/'
	Visible_Sources  = 'SPLUS_STRIPE82_master_catalogue_edr_march2018.fits'

	Visible_Sources_Catalogue  = fits.open(Visible_Path + Visible_Sources)
	Tiles_List  =  'Tiles_List.csv' # .csv here
	Tiles_Table  =  ascii.read(Visible_Path + Tiles_List)
	Tiles  =  Tiles_Table['Tiles']


	print "==========\nEnter the number of the pieces of code you want to run:\n 1. Find Radio sources in the Visible area\n 2. Find Radio and Visible cross-matches\n 3. Query for sources  spectroscopic redshift\n 4. Download Radio images (better resolution images will be also downloaded whenever available)\n 5. Plot Radio images and Visible field\n 6.Plot Radio and Visible matches\n=========="
	_input = raw_input('numbers:')
	_radio_cat_input = raw_input('Enter the radio catalog you want to use:')
	if not 'FIRST' in _radio_cat_input  and not 'HODGE' in _radio_cat_input :
		print "ERROR!"

	#print _input
	###if not os.path.isdir(Main_Path + 'PLOT_SOURCE_MATCHES/'):
	###	os.mkdir(Main_Path + 'PLOT_SOURCE_MATCHES/')
###
	###name  =  'SDSS_SPEC_Z_MATCH.csv'
	###file  =  open(os.path.join(Main_Path + 'PLOT_SOURCE_MATCHES/',name), 'w')
	###file.write('ID,RA_Radio,DEC_Radio,RA_SDSS,DEC_SDSS,FPEAK,FINT,SPLUS_Photo_z,SDSS_Spec_redshift')
	###file.write('\n')

	for Tile in range(len(Tiles)):
		Visible_Image  =  Tiles[Tile]
		'''
		RA_HMS, DEC_DMS, FILE,MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, TILE_RA_LOWER_LIMIT, TILE_DEC_LOWER_LIMIT, TILE_RA_UPPER_LIMIT, TILE_DEC_UPPER_LIMIT = RADIO_VISIBLE_FIELD_MATCH(Radio_Path,FIRST_Sources, Visible_Path, Visible_Image, 'FIRST')
		DOWNLOAD_RADIO_SOURCES(MATCH_RA_RADIO, MATCH_DEC_RADIO, RA_HMS, DEC_DMS, FILE, 'FIRST')
		RA_HMS, DEC_DMS, FILE,MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, TILE_RA_LOWER_LIMIT, TILE_DEC_LOWER_LIMIT, TILE_RA_UPPER_LIMIT, TILE_DEC_UPPER_LIMIT = RADIO_VISIBLE_FIELD_MATCH(Radio_Path,HODGE_Sources, Visible_Path, Visible_Image, 'HODGE')
		#
		
		#print TILE_RA_UPPER_LIMIT,TILE_RA_LOWER_LIMIT,TILE_DEC_UPPER_LIMIT,TILE_DEC_LOWER_LIMIT
		#SDSS_QUERY_SPEC_Z(TILE_RA_LOWER_LIMIT, TILE_DEC_LOWER_LIMIT, TILE_RA_UPPER_LIMIT, TILE_DEC_UPPER_LIMIT)
		
		DOWNLOAD_RADIO_SOURCES(MATCH_RA_RADIO, MATCH_DEC_RADIO, RA_HMS, DEC_DMS, FILE, 'HODGE')
		
		#PLOT_MATCHES(FILE, MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, 0, 0, 0)
		'''

		#FILE, MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, MATCH_FPEAK_RADIO, MATCH_FINT_RADIO, MATCH_Photo_z_Visible, RA_HMS, DEC_DMS = RADIO_VISIBLE_SOURCE_MATCH (Radio_Path, FIRST_Sources, Visible_Path, Visible_Sources, Visible_Image, 'FIRST')
		#PLOT_MATCHES(FILE, MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, MATCH_FPEAK_RADIO, MATCH_FINT_RADIO, MATCH_Photo_z_Visible)
		#FILE, MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, MATCH_FPEAK_RADIO, MATCH_FINT_RADIO, MATCH_Photo_z_Visible, RA_HMS, DEC_DMS = RADIO_VISIBLE_SOURCE_MATCH (Radio_Path, HODGE_Sources, Visible_Path, Visible_Sources, Visible_Image, 'HODGE')
		#DOWNLOAD_RADIO_SOURCES(MATCH_RA_RADIO, MATCH_DEC_RADIO, RA_HMS, DEC_DMS, FILE, 'FIRST')
		#
		#PLOT_MATCHES(FILE, MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, MATCH_FPEAK_RADIO, MATCH_FINT_RADIO, MATCH_Photo_z_Visible)
		#PLOT_TRIPLE_MATCHES(Radio_Path, Visible_Path, Visible_Sources, Visible_Image, FIRST_Sources, HODGE_Sources)
		if '1' in _input:
			RA_HMS, DEC_DMS, FILE,MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, TILE_RA_LOWER_LIMIT, TILE_DEC_LOWER_LIMIT, TILE_RA_UPPER_LIMIT, TILE_DEC_UPPER_LIMIT = RADIO_VISIBLE_FIELD_MATCH(Radio_Path,FIRST_Sources, Visible_Path, Visible_Image, _radio_cat_input)
		if '2' in _input:
			if _radio_cat_input == 'HODGE':
				FILE, MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, MATCH_FPEAK_RADIO, MATCH_FINT_RADIO, MATCH_Photo_z_Visible, RA_HMS, DEC_DMS = RADIO_VISIBLE_SOURCE_MATCH (Radio_Path, HODGE_Sources, Visible_Path, Visible_Sources, Visible_Image, _radio_cat_input)
			else:
				FILE, MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, MATCH_FPEAK_RADIO, MATCH_FINT_RADIO, MATCH_Photo_z_Visible, RA_HMS, DEC_DMS = RADIO_VISIBLE_SOURCE_MATCH (Radio_Path, FIRST_Sources, Visible_Path, Visible_Sources, Visible_Image, _radio_cat_input)
		if '3' in _input:	
			SDSS_QUERY_SPEC_Z(TILE_RA_LOWER_LIMIT, TILE_DEC_LOWER_LIMIT, TILE_RA_UPPER_LIMIT, TILE_DEC_UPPER_LIMIT)
		if '4' in _input:
			DOWNLOAD_RADIO_SOURCES(MATCH_RA_RADIO, MATCH_DEC_RADIO, RA_HMS, DEC_DMS, FILE, _radio_cat_input)
		if '5' in _input:
			PLOT_MATCHES(FILE, MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, MATCH_FPEAK_RADIO, MATCH_FINT_RADIO, MATCH_Photo_z_Visible, _radio_cat_input)
		if '6' in _input:
			PLOT_MATCHES(FILE, MATCH_RA_RADIO, MATCH_DEC_RADIO, image_data, Visible_w, FLAG, MATCH_FPEAK_RADIO, MATCH_FINT_RADIO, MATCH_Photo_z_Visible, _radio_cat_input)
	#data = ascii.read(Main_Path + 'PLOT_SOURCE_MATCHES/'+'SDSS_SPEC_Z_MATCH.csv')
	#Spec_z = data['SDSS_Spec_redshift']
	#Photo_z = data['SPLUS_Photo_z']
	#plt.plot(Spec_z, Photo_z, linestyle = 'None', marker='o', color='k')
	#plt.xlabel('Spectroscopic redshift (SDSS)')
	#plt.ylabel('Photometric redshift (S-PLUS)')
	#plt.show()
if __name__  ==  '__main__':
     main()