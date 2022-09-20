#!/bin/env python3
'''
READER for RIMFAX.

Simple reader for RIMFAX Calibrated radar products
Duncan A. Young
UTIG.
'''


import csv
import numpy as np
import os
import yaml
import xml.etree.ElementTree as ET
from matplotlib import pyplot as plt

ns = {'pds': 'http://pds.nasa.gov/pds4/pds/v1', 
    'mars2020':"http://pds.nasa.gov/pds4/mission/mars2020/v1",
    "geom": "http://pds.nasa.gov/pds4/geom/v1",
    "proc": "http://pds.nasa.gov/pds4/proc/v1",
    "msn_surface": "http://pds.nasa.gov/pds4/msn_surface/v1",
    "xsi": "http://www.w3.org/2001/XMLSchema-instance"
    }

def read_data_csv(path):
    with open(path) as csvfile:
        rimfaxreader = csv.reader(csvfile)
        for row in rimfaxreader:
            print(row)

def read_metadata(path):
    tree  = ET.parse(path)
    root = tree.getroot()
    for thing in root.findall('mars2020:RIMFAX_Parameters',ns):
        print(thing)
        print(thing.tag,thing.text)


def read_data_np(path):

    name=path.split('/')[-1].split('.')[0]

    rimfax=np.genfromtxt(path,delimiter=',',names=True,dtype=None)
    valid_records = []
    cal_records = []
    cal_trace = [0] * 1410

    lat=rimfax['ant_lat']
    lon=rimfax['ant_lon']
    sample_time_interval=rimfax['sample_time_increment']

    for i, record_type in enumerate(rimfax['record_type']):
        if record_type == 0:
            valid_records.append(i)
        if record_type == 8:
            cal_records.append(i)

    radargram = np.zeros((1410,len(valid_records)))
    raw_radargram = np.zeros((1410,len(valid_records)))
    
    for sample in range(1,1410,1):
        s = f's{sample:04d}'
        j = 0
        for k,value in enumerate(rimfax[s]):
            if k < 5:
                cal_trace[0] = cal_trace[0] + value
                j = j + 1

    for l,value in enumerate(cal_trace):
        cal_trace[l] = value/len(cal_records)

    times=[]
    distance=[0]
    for sample in range(1,1410,1):
        s = f's{sample:04d}'
        j = 0
        for k,value in enumerate(rimfax[s]):
            if k in valid_records:
                if not times:
                    time0=rimfax['sclk'][k]
                times.append(rimfax['sclk'][k]-time0)
                if j < 5:
                    cal_trace[sample] = cal_trace[sample] + value
                else:
                    raw_radargram[sample][j] = value 
                    dt = sample_time_interval[k]
                    radargram[sample][j] = value - (cal_trace[sample]/5)
                j = j + 1
                distance.append(distance[j-1]+rimfax['sounding_group_spacing'][k]/100)

    #print(times)
    filtered_radargram = np.zeros(radargram.shape)
    noise = np.zeros(radargram.shape)
    window=100
    left_limit = window/2 + 1
    right_limit = radargram.shape[1] - left_limit
    print(left_limit,right_limit)
    print(radargram.shape)
    for m,sample in enumerate(radargram):
        for n,trace in enumerate(sample):
            if ( n > left_limit ) and (  n < right_limit ):
                noise[m][n] =  np.mean( sample[n-int(window/2):n+int(window/2)] )
            elif (  n >= right_limit ):
                noise[m][n] = noise[m][int(right_limit-1)]
            else:
                noise[m][n] = 0


    filtered_radargram = radargram - noise
            
    extents=(min(distance),max(distance),radargram.shape[0]*dt,0)

    top=0.25e6
    bottom=-0.25e6

    print(np.min(raw_radargram))
    print(np.max(raw_radargram))
    plt.subplot(3,1,1)
    plt.imshow(raw_radargram,aspect='auto',cmap='gray',vmax=top,vmin=bottom,extent=extents)
    plt.title(f'Raw data')
    plt.ylabel('Fast Time (nsec)')
    plt.xlabel('Distance (m)')
    plt.subplot(3,1,2)
    plt.imshow(radargram,aspect='auto',cmap='gray',vmax=top,vmin=bottom,extent=extents)
    plt.title(f'Filter based on 1st 5 records')
    plt.xlabel('Distance (m)')
    plt.subplot(3,1,3)
    plt.imshow(filtered_radargram,aspect='auto',cmap='gray',vmax=top,vmin=bottom,extent=extents)
    plt.title(f'+ {window} record filter')
    plt.xlabel('Distance (m)')
    plt.tight_layout()
    plt.savefig(f'{name}.png')

def plot_data(number):

    path='/disk/kea/SDS/orig/supl/xlob-rimfax/RIMFAX/urn-nasa-pds-mars2020_rimfax/data_calibrated/2021'
    data_file=f'rimfax_calibrated_{number:05d}.csv'
    print(f'Doing {data_file}')
    read_data_np(os.path.join(path,data_file))
 
plot_data(31)
plot_data(48)
plot_data(86)
