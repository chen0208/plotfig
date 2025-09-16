#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 09:49:36 2025

@author: cjq
"""

###################################################
M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
G = 6.674e-8 # 万有引力常数
###################################################
import numpy as np
from unyt import unyt_array
import yt
import matplotlib.pyplot as plt
import os
import sys
import struct
import math
from dump_ic import ic_file


#set unit
unit_base = {
    "UnitLength_in_cm": 1,
    "UnitMass_in_g": 1,
    "UnitVelocity_in_cm_per_s": 1,
}

#load data
fname1 = '/home/cjq/gadget4_yt/BSS_star/output_1M50/snapshot_100.hdf5'
fname2 = '/media/cjq/zzq/cjq/BSS/output_2M100/snapshot_100.hdf5'

bbox_lim = 1e30  # kpc

bbox = [[-bbox_lim, bbox_lim], [-bbox_lim, bbox_lim], [-bbox_lim, bbox_lim]]

ds1 = yt.load(fname1, unit_base=unit_base, bounding_box=bbox)
ad1 = ds1.all_data()
ds2 = yt.load(fname2, unit_base=unit_base, bounding_box=bbox)
ad2 = ds2.all_data()

co1 = ad1["PartType0","Coordinates"]
mass1 = ad1["PartType0","Masses"]
velocities1 = ad1["PartType0","Velocities"]
u1 = ad1["PartType0","InternalEnergy"]

co1_c = ad1["PartType1","Coordinates"]
mass1_c = ad1["PartType1","Masses"]
velocities1_c = ad1["PartType1","Velocities"]

co2 = ad2["PartType0","Coordinates"]
mass2 = ad2["PartType0","Masses"]
velocities2 = ad2["PartType0","Velocities"]
u2 = ad2["PartType0","InternalEnergy"]

co2_c = ad2["PartType1","Coordinates"]
mass2_c = ad2["PartType1","Masses"]
velocities2_c = ad2["PartType1","Velocities"]

#caculate total number of particles
total_number_of_particles_gas = len(co1)+len(co2)
total_number_of_particles_halo = len(co1_c)+len(co2_c)

def binary_initial_conditions(M1,M2,a,e):
    M=M1+M2
    r_peri = a * (1+e)
    v_peri = (G*M*(1-e)/(a*(1+e)))**0.5
    
    r1 = r_peri * M2/M
    r2 = r_peri * M1/M
    v1 = v_peri * M2/M
    v2 = v_peri * M1/M
    
    return r1,v1,r2,v2

def calculate_angular_velocity(m1, m2, r):
	numerator = G * (m1 + m2)
	denominator = r ** 3
	angular_velocity = math.sqrt(numerator / denominator)
	return angular_velocity
	
def angular_to_linear_velocity(angular_velocity, radius):
    return np.cross(angular_velocity, radius)
    
def co_v_trans(r,w,co,v):
    P = 2*np.pi/w
    co = np.array(co)
    co[:,0] += r
    v = angular_to_linear_velocity([0,0,w],co)
    return co,v


def e_r1_r2(m1,m2,a):
    q = m2/m1
    r1 = q*a/(1+q)
    r2 = a/(1+q)
    return r1,r2

m1=np.sum(mass1+mass1_c)
m2=np.sum(mass2+mass2_c)

a = 4.3*R_to_solar
e = 0.75
r1,v1,r2,v2 = binary_initial_conditions(m1,m2,a,e)
print((r1+r2)/R_to_solar)
w = v1/r1

co1,velocities1 = co_v_trans(r1,w,co1,velocities1)
co2,velocities2 = co_v_trans(-r2,w,co2,velocities2)
co1_c,velocities1_c = co_v_trans(r1,w,co1_c,velocities1_c)
co2_c,velocities2_c = co_v_trans(-r2,w,co2_c,velocities2_c)

co_gas = np.concatenate((co1,co2),axis=0)
mass_gas = np.concatenate((mass1,mass2),axis=0)
velocities_gas = np.concatenate((velocities1,velocities2),axis=0)
u = np.concatenate((u1,u2),axis=0)
co_c = np.concatenate((co1_c,co2_c),axis=0)
mass_c = np.concatenate((mass1_c,mass2_c),axis=0)
velocities_c = np.concatenate((velocities1_c,velocities2_c),axis=0)
co = np.concatenate((co_gas,co_c),axis=0)
mass = np.concatenate((mass_gas,mass_c),axis=0)
velocities = np.concatenate((velocities_gas,velocities_c),axis=0)


npart=[total_number_of_particles_gas,total_number_of_particles_halo,0,0,0,0]
ic_file(npart,total_number_of_particles_gas+total_number_of_particles_halo,co,mass,velocities,u,"binary21_4.dat")
