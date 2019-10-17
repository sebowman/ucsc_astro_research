'''This script is designed to scrape photometry data for a single object based on its image coordinates from image data processed by the DOLPHOT stellar photometry package; and average the magnitude values from the observations at those coordinates, grouped by date and filter.  

The required arguments are the name of the main DOLPHOT output file {output_filename} and the X, Y image coordinates of the desired object. (Example command line input to run script: 'python dolphot_retrieval.py iPTF13bvn-WFPC2 438 520').  This script must be run in a directory containing the output files from DOLPHOT.  The paths to the filters.dat files for the relevant instruments must be updated in lines 124-128.

Photometry data from the object detected by DOLPHOT and having the shortest distance to the given image coordinates, within some radius, is selected and used for the analysis in the remainder of the script.  This radius can be changed by updating line 33.

The script imports the packages 'astropy.io.fits', 'astropy.io.ascii', and 'numpy', the class 'astropy.table.Table', and the modules 'sys' and 'os'.

The output will be a CSV file titled "phot_plot_data_{output_filename}.csv", containing the average magnitude and average magnitude uncertainty for each observation date and filter.
'''


# import dependencies
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import numpy as np
import sys
import os

# read in dolphot output file
output_filename = sys.argv[1]

# select data rows containing X,Y coordinates within a certain radius of desired coordinates
src = []

with open (output_filename) as outfile:
    for line in outfile:
        line = line.strip()
        columns = line.split()
        X = float(columns[2])
        Y = float(columns[3])
        r = (X-float(sys.argv[2]))**2 + (Y-float(sys.argv[3]))**2
        if r <= 2:
            src.append(columns)

# select the closest object detected by DOLPHOT to the desired coordinates for further analysis, 
# or report error if none
if len(src) == 0:
    print('There are no sources in the selected range.')
    sys.exit()

src = Table(rows=src)
    
if len(src) > 1:
    src1 = np.array(src)
    dist = []
    for i in src1:
        dist.append(np.sqrt((float(sys.argv[2])-float(i['col2']))**2 + (float(sys.argv[3])-float(i['col3']))**2))
    src = Table(src[dist.index(min(dist))])

print(src[0][2],src[0][3])

# read in columns file and make it a table object
cols = open(str(output_filename)+str('.columns'), 'r')
columns = []
for line in cols:
    line = line.strip()
    columns.append(str(line))
cols.close()

# get magnitude values from selected data
magnitudes = []
for i in columns:
    if '(' in i:
        if 'VEGAMAG' in i:
            magnitudes.append(i.split('. '))
magnitudes = Table(magnitudes)
mag = []
for i in magnitudes[0]:
    mag.append(src[0][str('col')+str(int(i)-1)])

# get magnitude uncertainty values from selected data
mag_error = []
for i in columns: 
    if '(' in i:
        if 'Magnitude uncertainty' in i:
            mag_error.append(i.split('. '))
mag_error = Table(mag_error)
mag_err = []
for i in mag_error[0]:
    mag_err.append(src[0][str('col')+str(int(i)-1)])
    
# get instrument filters of selected images
fltr1 = []
for i in magnitudes[1]:
    fltr1.append(i.split('('))
fltr1 = Table(fltr1)
fltr2 = []
for i in fltr1[1]:
    fltr2.append(i.split(', '))
fltr2 = Table(fltr2)
fltr = []
for i in fltr2[0]:
    fltr.append(i)
    
# get .fits file names of selected images
filenames1 = []
for i in columns: 
    if '(' in i: 
        if 'VEGAMAG' in i:
            filenames1.append(i.split(', '))
filenames1 = Table(filenames1)
filenames2 = []
for i in filenames1[1]:
    filenames2.append(i.split(' ('))
filenames2 = Table(filenames2)
filenames = []
for i in filenames2[0]:
    filenames.append(i)

# get observation dates of selected images
date_obs = []
for i in filenames: 
    hdul = fits.open(str(i)+str('.fits'))
    hdr = hdul[0].header
    date_obs.append(hdr['DATE-OBS'])

# get zero points of selected data, based on filters.dat files for each relevant instrument
hdul_1 = fits.open(str(filenames[0])+str('.fits'))
hdr_1 = hdul_1[0].header
instr = hdr_1['INSTRUME']

if instr == 'WFC3':
    fdata = open('/home/marley/dolphot/wfc3/data/filters.dat', 'r')
if instr == 'ACS':
    fdata = open('/home/marley/dolphot/acs/data/filters.dat', 'r')
if instr == 'WFPC2':
    fdata = open('/home/marley/dolphot/wfpc2/data/filters.dat', 'r')

dat = []
for line in fdata:
    line = line.strip()
    dat.append(str(line))

zpt_vals = []
for k in fltr:
    if k in dat:
        zpt_vals.append(dat[dat.index(k)+2])
        
zpts1 = []
if ' ' in zpt_vals[1]:
    for i in zpt_vals: 
        zpts1.append(i.split(' '))
    zpts1 = Table(zpts1)
    zpts2 = []
    if zpts1[0][1] == '-1':
        for i in zpts1[1]:
            zpts2.append(i.split(' '))
    else:
        for i in zpts1[0]:
            zpts2.append(i.split(' '))
    zpts2 = Table(zpts2)
    zpts = []
    for i in zpts2[0]:
        zpts.append(i)
else: 
    zpts = []
    for i in zpt_vals:
        zpts.append(i)
    
# get measured count values from selected data
cts1 = []
for i in columns: 
    if '(' in i:
        if 'Measured counts,' in i:
            cts1.append(i.split('. '))
cts1 = Table(cts1)
cts_index = []
for i in cts1[0]:
    cts_index.append(int(i)) 
cts = []
for i in cts_index:
    cts.append(src[0][str('col')+str(int(i)-1)])
    
# get exposure time values from selected data
fltr3 = []
for i in fltr2[1]:
    fltr3.append(i.split(' '))
fltr3 = Table(fltr3)
expt = []
for i in fltr3[0]:
    expt.append(i)

# write to a file the observation date, filter, measured counts, exposure time, zero point,
# normalized count rate, and normalized count rate uncertainty values
import csv
with open('phot_data.csv', mode='w+') as csv_file:
    fieldnames = ['date_obs', 'filter', 'cts', 'expt', 'mag_err', 'zpts']
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

    writer.writeheader()
    for i in range(0,len(zpts)):
        writer.writerow({'date_obs': date_obs[i], 'filter': fltr[i], 'cts': cts[i], 'expt': expt[i], 'mag_err': mag_err[i], 'zpts': zpts[i]})

# remove non-detections
phot_data = ascii.read('phot_data.csv')

for i, row in enumerate(phot_data):
    if row['cts'] < 0:
	for i, row in enumerate(phot_data):
            if row['cts'] < 0: 
		phot_data.remove_row(i)

# separate the data from the previously written file into groups by filter and observation date
grps = []
for i in set(phot_data['filter']):
    with open('grps1.csv', mode ='w+') as grps_i:
        writer = csv.writer(grps_i)
        for row in phot_data:
            if row[1] == i:
                writer.writerow(row)
    grps_i = ascii.read('grps1.csv')
    for k in set(grps_i['col1']):
        with open('grps2.csv', mode ='w+') as grps_k:
            writer = csv.writer(grps_k)
            for a in grps_i:
                if a[0] == k:
                    writer.writerow(a)
        grps_k = ascii.read('grps2.csv')
        grps.append(grps_k)

# calculate avgerage magnitude and magnitude uncertainty per group from measured counts, 
# exposure time, zero points, normalized counts and normalized counts uncertainty
def f_to_m(f,zpt):
    mag = -2.5*np.log10(f) + zpt
    return mag

def df_to_dm(df,f):
    dm = (1.086/f)*df
    return dm

def dm_to_df(dm,f):
    df = (dm*f)/(1.086)
    return df

f_avg = []
sigma_f_avg = []
m_avg = []
sigma_m_avg = []
dates = []
filters = []
for i in grps: 
    f_av = (1./len(i))*sum(i['col3']/i['col4'])
    f_avg.append(f_av)
    sigma_f = dm_to_df(i['col5'],i['col3']/i['col4'])
    sigma_f_av = (1./len(i))*sum(sigma_f)
    sigma_f_avg.append(sigma_f_av)
    zpt = i[0]['col6']
    m_avg.append(f_to_m(f_av,zpt))
    sigma_m_avg.append(df_to_dm(sigma_f_av,f_av))
    dates.append(i[0]['col1'])
    filters.append(i[0]['col2'])

# write to a file the final calculated average magnitude and magnitude uncertainty
# of our target for each observation date and filter

with open('phot_plot_data_' + output_filename + '.csv', mode='w+') as csv_file:
    fieldnames = ['date_obs', 'filter', 'm_avg', 'sigma_m_avg']
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

    writer.writeheader()
    for i in range(0,len(grps)):
        writer.writerow({'date_obs': dates[i], 'filter': filters[i], 'm_avg': m_avg[i], 'sigma_m_avg': sigma_m_avg[i]})

# remove intermediate files
os.remove('phot_data.csv')
os.remove('grps1.csv')
os.remove('grps2.csv')
