#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 14:42:42 2021

@author: tomandrews_h
"""

from astropy.io import fits
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import CircularAnnulus
from astropy.table import Table
#%%
"""
Mask edges to avoid edge effects
"""

hdulist = fits.open('A1_mosaic.fits')

x_min = 263
y_min = 243

x_max = 2350
y_max = 4400

#%%
"""
Identify brightest stars/galaxies
"""

def max_brightness(image_input, mask_image_ones):
    image_high = image_input * mask_image_ones
    ind = np.unravel_index(np.argmax(image_high, axis=None), image.shape)
    return ind
        
"""
Move 0s to matrix mask_image_ones to mask off particular region
"""
def clean_image_new(y1, y2, x1, x2):
    for y in range(y1, y2):
        for x in range(x1,x2):
            if image[y][x] >= 0:
                mask_image_ones[y][x] = 0
            else:
                continue
""" Make a 2d gaussian plot

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    
    central_maximum is a scaling so intensity is
    approximately equal to that of a galaxy
    """    

def twodGaussian(size, fwhm ,central_max, center):
    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return central_max*np.exp(-((x-x0)**2 + (y-y0)**2) / fwhm**2) + 3418
"""
image = twodGaussian(size=30, fwhm = 4, central_max = 10000, center = None)
mask_image_ones = sp.ones(cropped_image.shape, dtype = int)
x_min = 0 
y_min = 0
x_max = 30
y_max = 30
"""

"""
Determines source magnitude from total pixel count of object
"""
def source_magnitude(pc):
    return -2.5*sp.log10(pc) + 25.3

"""
Determing the effective radius of the galaxy
"""
def object_find(x_centre, y_centre, mask_image_ones):
    radius_list = sp.arange(1,51,1)
    radius_plot = []
    sum_list = []
    bkg_list =[]
    bkg_diff = -10
    for radius in radius_list:
        if abs(bkg_diff) < 1 or (bkg_diff*(-1))<= 0:
            break
        else:
            position_centre = (x_centre,y_centre)
            aperture = CircularAperture(position_centre, r = radius)
            aperture.plot(color='white', lw=2)
            
            aperture2 = CircularAperture(position_centre, r = radius+1)
            
            annulus_aperture = CircularAnnulus(position_centre, r_in = radius, r_out = radius+1)
            annulus_aperture.plot(color='red', lw=2)
            
            apers = [aperture, annulus_aperture]
            mask_image_zeros = 1- mask_image_ones
        
            a = 1 - mask_image_ones
            mask_global = a.astype(bool)
            
     
            phot_table = aperture_photometry(image, apers,mask = mask_global)
            mask1 = aperture.to_mask(method = 'center')
            mask = mask1
            image1 = mask.to_image(shape = image.shape)
            
            aperture_area = aperture.area - np.sum(image1*np.copy(mask_image_zeros))
           
            
            mask2 = aperture2.to_mask(method = 'center')
            mask = mask2
            image2 = mask.to_image(shape= image.shape)
            annulus_area = (aperture2.area - np.sum(image2*np.copy(mask_image_zeros)))- aperture_area
            
            
            bkg_mean = phot_table['aperture_sum_1']/annulus_area
            bkg_list.append(bkg_mean)
            bkg_sum = bkg_mean*aperture_area
            final_sum = phot_table['aperture_sum_0'] - bkg_sum
            phot_table['residual_aperture_sum'] = final_sum
            sum_list.append(final_sum)
            radius_plot.append(radius)
           
            
            if radius == 1:
                continue
            #elif bkg_list[radius-1] < 3418:
                #break
            else:
                bkg_diff = bkg_list[radius-1] -bkg_list[radius-2]
                continue
                
            continue
           
    masks = aperture.to_mask(method ='center')
  
    mask = masks
  
    image_l = mask.to_image(shape = image.shape)
    b  = 1 - image_l
    
    mask_image_ones = mask_image_ones*b
    
    if len(radius_plot) <3:
        radius = None
        
    return radius, final_sum, mask_image_ones

#%%
"""
This will perform aperture analysis on either the next brightest object 
or a manually selected object by entering x_centre and y_centre
"""

def catalogue_object(radius, final_sum, mask_image_ones, y_centre, x_centre):
    
    if radius != None and final_sum>0:
        x.append(x_centre)
        y.append(y_centre)
        pxc.append(final_sum)
        eff_r.append(radius)
        error.append(sp.sqrt(final_sum))
        mag.append(source_magnitude(final_sum))
        t = Table([x,y,pxc, eff_r, error, mag], names = ('x center', 'y center', 'pixel count', 'effective radius', 'error', 'source magnitude'))
        y_centre, x_centre = max_brightness(image, mask_image_ones)
    
    return t, y_centre, x_centre
#%%%
"""
Testing our functions on a syntehtic 2d gaussian object
"""
image = twodGaussian(size=30, fwhm = 4, central_max = 10000, center = None)
mask_image_ones = sp.ones(image.shape, dtype = int)
x_min = 0 
y_min = 0
x_max = 30
y_max = 30

y_centre, x_centre = max_brightness(image, mask_image_ones)
radius, final_sum, mask_image_ones = object_find(x_centre, y_centre, mask_image_ones)
x= []
y = []
pxc =[]
eff_r = []
error = []
mag = []
catalogue_object(radius, final_sum, mask_image_ones, y_centre, x_centre)
#%%%
"""
Displaying function working with two adjacent objects
"""
image = hdulist[0].data
x_min = 118
y_min = 120
x_max = 2469
y_max = 4512
mask_image_ones = sp.ones(image.shape, dtype = int)
x_centre = 353
y_centre = 3037
radius, final_sum, mask_image_ones = object_find(x_centre, y_centre, mask_image_ones)
x= []
y = []
pxc =[]
eff_r = []
error = []
mag = []
catalogue_object(radius, final_sum, mask_image_ones)

#repeat for elliptical

#%%
"""
Considers whole image, eacgh tim next secytion is run, the next brightest
objects is catalogued
"""
image = hdulist[0].data
x_min = 118
y_min = 120
x_max = 2469
y_max = 4512
mask_image_ones = sp.ones(image.shape, dtype = int)
y_centre, x_centre = max_brightness(image, mask_image_ones)
x= []
y = []
pxc =[]
eff_r = []
error = []
mag = []
#%%
radius, final_sum, mask_image_ones = object_find(x_centre, y_centre, mask_image_ones)
catalogue_object(radius, final_sum, mask_image_ones)

f2 = plt.figure()
f2.add_subplot(111)
plt.title('Masked Image')
plt.spy(mask_image_ones, origin = 'lower')

#%%
"""
Clean up messy/ bleeding objects with updated approach
Block off region instead of individual pixels
"""

image = hdulist[0].data
mask_image_ones = sp.ones(image.shape, dtype = int)

clean_image_new(0, y_min, 0, 2570)
clean_image_new(y_max, 4611, 0, 2570)
clean_image_new(0, 4611, 0, x_min)
clean_image_new(0, 4611, x_max, 2570)


def circle_mask(position_centre, radius, mask_image, image_input):
    aperture = CircularAperture(position_centre, r = radius)
    masks = aperture.to_mask(method ='center')
    
    mask = masks
    big_circle = mask.to_image(shape = image.shape)

    central_circle_mask = 1 - big_circle
    
    mask_image_ones = mask_image*central_circle_mask
    return mask_image_ones


mask_image_ones = circle_mask(position_centre = (1434, 3199), radius = 290., mask_image = mask_image_ones,image_input = image)

for y in range(106, 450):
    for x in range(1014, 1692):
        if hdulist[0].data[y][x] > 6000:
            mask_image_ones[y][x] = 0
        else:
            continue

for y in range(449, 489):
    for x in range(1370,1470):
        if hdulist[0].data[y][x] > 6000:
            mask_image_ones[y][x] = 0
        else:
            continue
        
clean_image_new(2221,2363, 860, 950)
clean_image_new(3199, 3421, 724, 836)
clean_image_new(0, 4571, 1422, 1454)
clean_image_new(2700,2846, 930, 1024)
#clean_image_new(1917, 1960, 652, 700)
#clean_image_new(1475,1517,610,662)
#clean_image_new(2277, 2328, 422, 472)
clean_image_new(4052, 4140, 513, 609)
#clean_image_new(558, 602, 1752, 1800)
#mask_image_ones = circle_mask(position_centre = (191, 3924), radius = 27., mask_image = mask_image_ones, image_input = image)
#clean_image_new(2242, 2282, 684, 728)
clean_image_new(1398, 1454, 2056, 2138)
clean_image_new(3996, 4070, 1429,1495)
#clean_image_new(1634,1674, 947, 989)
clean_image_new(3704,3806, 2094, 2175)
#mask_image_ones = circle_mask(position_centre = (1312, 4400), radius = 22, mask_image = mask_image_ones, image_input = image)
clean_image_new(2276,2338, 2097, 2165)
#mask_image_ones = circle_mask(position_centre = (1365, 4330), radius = 25, mask_image = mask_image_ones, image_input = image)
mask_image_ones = circle_mask(position_centre = (2466,3410), radius = 34., mask_image = mask_image_ones, image_input = image)
#clean_image_new(2966,2990,1406, 1430)
#clean_image_new(3350,3392, 875, 915)
#clean_image_new(3828,3864, 2253, 2303)
#mask_image_ones = circle_mask(position_centre = (208,844), radius = 21., mask_image = mask_image_ones, image_input = image)
#mask_image_ones = circle_mask(position_centre = (577, 4315), radius = 20., mask_image = mask_image_ones, image_input = image)


#plt.clf()
plt.title('Mask Image Initial')
plt.spy(mask_image_ones, origin = 'lower')
plt.savefig('clean_image_df4')
#plt.show()
#%%%
"""
Plots histogram fior background noise of entire image
"""
hdulist = fits.open('A1_mosaic.fits')

data_list = []

for y in range(y_min, y_max):
    for x in range(x_min,x_max):
        if hdulist[0].data[y][x]*mask_image_ones[y][x] < 35000 and hdulist[0].data[y][x]*mask_image_ones[y][x] > 0:
            data_list.append(hdulist[0].data[y][x])
        else:
            continue
 
frequency,count, a = plt.hist(data_list, bins = 15000, color = 'darkkhaki')
 
bin_middle = []
for i in range(1,len(count)):
    x = (count[i]+count[i-1])/2 # calculates centre of each bin
    bin_middle.append(x)
 
def gaus(x,b,c,d):
    return d*sp.exp(-(x-b)**2/(2*c**2))

popt,pcov = curve_fit(gaus,bin_middle,frequency,p0=[3500, 10, 3000000])
g = sp.arange(min(bin_middle),max(bin_middle),0.5)

plt.plot(g,gaus(g, *popt), label="Fit", color='darkslategrey')
plt.vlines(popt[0], 0, popt[2], linestyles='dashed', label='Mean = {:.0f}'.format(popt[0]), color='darkslategrey')
plt.xlabel('Pixel Count', size='15')
plt.ylabel('Frequency of Count', size='15')
plt.title('Upper bound to remove')
plt.grid()
plt.legend(fontsize = 'large')

print(popt[0], popt[1], popt[2]) 
print(sp.sqrt(sp.diag(pcov)))         
plt.savefig('Global background with mask')
plt.show()
#%%
"""
Calculates centre and determines radius and brightness
This code runs over enire image until criteria for min brightness is met
"""

x= []   #creates columns for table t
y = []
pxc =[]
eff_r = []
error = []
mag = []

y_centre, x_centre = max_brightness(image, mask_image_ones)
i=1
while image[y_centre][x_centre] > 3464.6:
    radius, final_sum, mask_image_ones = object_find(x_centre, y_centre, mask_image_ones)
    if radius != None and final_sum >0:   
        print('object no. {i}  location {location}  pc {pc}'.format(i =i, location = (x_centre, y_centre), pc = image[y_centre][x_centre]))
        x.append(x_centre)
        y.append(y_centre)
        pxc.append(final_sum)
        eff_r.append(radius)
        error.append(sp.sqrt(((2.5/(sp.log(10)*final_sum))*sp.sqrt(final_sum))**2 + (0.02**2)))
        mag.append(source_magnitude(final_sum))
        t = Table([x,y,pxc, eff_r, error, mag], names = ('x center', 'y center', 'pixel count', 'effective radius', 'error', 'source magnitude'))
        y_centre, x_centre = max_brightness(image, mask_image_ones)
        i += 1
        continue
    else:
        y_centre, x_centre = max_brightness(image, mask_image_ones)
        continue

f2 = plt.figure()
f2.add_subplot(111)
plt.title('Final Mask Image')
plt.spy(mask_image_ones, origin = 'lower')
plt.savefig('final_mask_df4')

#%%%
from matplotlib.lines import Line2D

def logN(set_mag, mag):
    count = 0
    for m in mag:
        if m < set_mag:
            count += 1
            print(m, 'm')
            continue
        else:
            continue
    
    return count

data = []

for mag_input in range(0, int(max(mag))+2):
    count =  logN(mag_input, mag)
    if count == 0:
        continue
    else:
        print(count, 'count')
        data.append([mag_input, sp.log10(count)])

m_list =[]
logN_list =[]
m_listall = []
logN_listall =[]
for m, logN in data:
    if m<=17:   
        m_list.append(m)
        logN_list.append(logN)
        m_listall.append(m)
        logN_listall.append(logN)
        continue
    else:
        m_listall.append(m)
        logN_listall.append(logN)
        continue

#plt.scatter(*zip(*data))
f3 = plt.figure()
#f3.add_subplots(111)

fit, cov = np.polyfit(m_list, logN_list, 1, cov = True)
poly = sp.poly1d(fit)

#line_1 = Line2D(m_list,poly(m_list), color = 'darkred', linestyle = '--',label = 'Linear Fit')
#f3.add_line(line_1)
plt.plot(m_list,poly(m_list), color = 'darkred', linestyle = '--',label = 'gradient: {fit:.2f} +/- {error:.2f}'.format(fit = fit[0], error =sp.sqrt(cov[0][0])))
plt.plot(m_listall, logN_listall, color = 'plum',marker = 'x', linestyle ='None')
plt.xlabel('Magnitude')
plt.ylabel('Log N(magnitude)')
plt.title ('Log N(m) v m')

print('gradient: ', fit[0], ' +/-', sp.sqrt(cov[0][0]))
print('y_intercept: ', fit[1], ' +/-', sp.sqrt(cov[1][1]))
plt.legend()
plt.savefig('logN_plot_df4')
plt.show()
#%%
import pandas as pd

df_phototable = pd.DataFrame({'x': x,
                   'y': y,
                   'pxc':   pxc,
                   'eff_r': eff_r,
                   'error':error,
                   'mag':mag})
df_magdata = pd.DataFrame({'m_list':m_list,
                   'logN_list': logN_list})
df_magdata_all = pd.DataFrame({'m_listall':m_listall,
                   'logN_listall':logN_listall})
    
df_phototable.to_csv('phototable_df4.csv', index=False)
df_magdata.to_csv('magdata_df4.csv', index = False)
df_magdata_all.to_csv('magdata_all_df4.csv', index = False)
df_readme = pd.DataFrame({'data set 4': 'min pc 3464.6 with only saturated objects masked; see image "saturated mask locations"',
                          'x':  (x_min,x_max),
                          'y':  (y_min,y_max),
                          })
df_readme.to_csv('readme_df4.csv', index = False)

#%%
"""
crop image after analysis done
"""

cropped_image_x1 = sp.delete(hdulist[0].data, slice(x_max,2571),1)
cropped_image_x2 = sp.delete(cropped_image_x1, slice(0,x_min),1)

cropped_image_y1 = sp.delete(cropped_image_x2, slice(y_max, 4612), 0)
cropped_image = sp.delete(cropped_image_y1, slice(0,y_min), 0)

masked_image_x1 = sp.delete(mask_image_ones, slice(x_max,2571),1)
masked_image_x2 = sp.delete(masked_image_x1, slice(0,x_min),1)

masked_image_y1 = sp.delete(masked_image_x2, slice(y_max, 4612), 0)
mask_image_ones_cropped = sp.delete(masked_image_y1, slice(0,y_min), 0)













