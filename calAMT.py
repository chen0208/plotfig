#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 20:32:15 2025

@author: cjq
"""
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import yt
import os
import sys
sys.path.append('/home/cjq/MESA/0.25M_0.99M')
import plot_index
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator

fname = '/media/cjq/zzq/cjq/potential_BSS/test217/snapshot_000.hdf5'
def cal_total_angular_mometom(fname):
    unit_base = {
        "UnitLength_in_cm": 1,
        "UnitMass_in_g": 1,
        "UnitVelocity_in_cm_per_s": 1,
    }
    
    bbox_lim = 3.5e11  # kpc
    
    bbox = [[-bbox_lim, bbox_lim], [-bbox_lim, bbox_lim], [-bbox_lim, bbox_lim]]
    
    ds = yt.load(fname, unit_base=unit_base, bounding_box=bbox)
    ds.index
    
    ad = ds.all_data()
    
    coordinates = np.array(ad["PartType0", "Coordinates"])
    M = np.array(ad["PartType0", "Masses"])
    V = np.array(ad["PartType0", "Velocities"])
    
    total_angular_momentom = np.sum(np.cross(coordinates,M[0]*V))
    return total_angular_momentom


files = sorted(glob('snapshot_*.hdf5'))
for fn in files:
    total_angular_momentom = cal_total_angular_mometom(fn)
    with open('AMT_total.txt','a+') as f:
        f.write(str(total_angular_momentom)+'\n')
        
