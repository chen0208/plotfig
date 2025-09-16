import numpy as np
import yt
import argparse
from glob import glob
import matplotlib.pyplot as plt
import imageio 
###################################################
M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
G = 6.674e-8 # 万有引力常数
###################################################
unit_base = {
    "UnitLength_in_cm": 1,
    "UnitMass_in_g": 1,
    "UnitVelocity_in_cm_per_s": 1,
}


bbox_lim = 2e30  # kpc
bbox = [[-bbox_lim, bbox_lim], [-bbox_lim, bbox_lim], [-bbox_lim, bbox_lim]]

@yt.particle_filter(requires=["density"],filtered_type="PartType0")
def high_density(pfilter, data):
    return data[pfilter.filtered_type,"density"] >1e-6

def yt_plot_slice(fname):


    ds = yt.load(fname, unit_base=unit_base, bounding_box=bbox)
    ds.add_particle_filter("high_density")
    
    ad = ds.all_data()
    coordinates = np.array(ad["PartType1", "Coordinates"])
    x_star1 = coordinates[0,0]
    x_star2 = coordinates[1,0]
    y_star1 = coordinates[0,1]
    y_star2 = coordinates[1,1]
    
    x_c = np.sqrt(x_star1**2+x_star2**2)
    y_c = np.sqrt(y_star1**2+y_star2**2)
    x1 = round(x_c/R_to_solar,2)
    y1 = round(y_c/R_to_solar,2)
    
    Sx=yt.SlicePlot(ds, "z", ("PartType0", "density"),width=(10.0,'Rsun'),center=[0,0,0])
#    Sx=yt.SlicePlot(ds, "z", ("PartType0", "density"),width=(10.0,'Rsun'),center=[x1*R_to_solar,y1*R_to_solar,0])
    colorbar = Sx.plots[("PartType0", "density")].cb
    ticks_values = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2]
    colorbar.set_ticks(ticks_values)
    ticks_label = [r"$10^{-6}$",r"$10^{-5}$",r"$10^{-4}$",r"$10^{-3}$",r"$10^{-2}$",r"$10^{-1}$",r"$10^{0}$",r"$10^{1}$",r"$10^{2}$"]
    colorbar.set_ticklabels(ticks_label)
    Sx.annotate_timestamp(corner="upper_left",time_unit='s',text_args={"color":"white"},redshift=True,draw_inset_box=True)
    Sx.annotate_scale(corner="upper_right")
    Sx.set_zlim(('PartType0', 'density'),1e-6,1e2)
    Sx.set_background_color(("PartType0", "density"))
    Sx.save()

yt_plot_slice('snapshot_066.hdf5')
files = sorted(glob('snapshot_*.hdf5'))
for fn in files:
    yt_plot_slice(fn)


