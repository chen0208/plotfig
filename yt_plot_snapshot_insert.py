#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:23:26 2022

@author: cjq
"""

import numpy as np

import yt
import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', default=None)
    args = parser.parse_args()

    
    
fname = args.f
#fname = 'snapshot_000.hdf5'

unit_base = {
    "UnitLength_in_cm": 1,
    "UnitMass_in_g": 1,
    "UnitVelocity_in_cm_per_s": 1,
}

bbox_lim = 1e12  # kpc

bbox = [[-bbox_lim, bbox_lim], [-bbox_lim, bbox_lim], [-bbox_lim, bbox_lim]]

ds = yt.load(fname, unit_base=unit_base, bounding_box=bbox)
ds.index

ad = ds.all_data()

#px = yt.ProjectionPlot(ds, "x", ("gas", "density"))
#px.save()


sorted(ds.field_list)
print(ds.field_list)

density = ad["PartType0", "density"]
coordinates = ad["PartType0", "Coordinates"]
center = [0,0,0]



Sx=yt.SlicePlot(ds, "z", ("PartType0", "density"))
Sx.hide_colorbar()
Sx.annotate_timestamp(corner="upper_left",time_unit='s',text_args={"color":"white"},redshift=True,draw_inset_box=True)
Sx.annotate_scale(corner="upper_right")
Sx.save()
px = yt.ProjectionPlot(ds, "z", ("PartType0", "density"), center=center)
plot = px.plots["PartType0","density"]
colorbar = plot.cb
px._setup_plots()


