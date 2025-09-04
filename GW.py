#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 17:15:32 2025

@author: chenjingqi
"""

import legwork
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import legwork.source as source
import matplotlib.font_manager as fm
from matplotlib.ticker import AutoMinorLocator,MultipleLocator
font_path = '/Users/chenjingqi/Desktop/hybridCONe_plot/TimesNewRoman.ttf'
font_prop_legend = fm.FontProperties(fname=font_path,size=11)
font_prop_label = fm.FontProperties(fname=font_path,size=20)
font_prop_title = fm.FontProperties(fname=font_path,size=12.5)
font_prop_text = fm.FontProperties(fname=font_path,size=20)

#import legwork.visualisation as vis
data_binary = 'binary_history.data'
data = np.genfromtxt(data_binary, names=True, skip_header=5)
Age = np.log10(data['age'])
m1 = data['star_1_mass'] * u.Msun
m2 = data['star_2_mass'] * u.Msun
ecc = data['eccentricity']
dist = np.full(11298,553) * u.pc
f_orb = 2/(data['period_days']*24*60*60) * u.Hz
m_chrip = (m1*m2)**0.6/(m1+m2)**0.2
hc = 2.5e-20*(f_orb)**(7/6)*(m_chrip)**(5/3)*(15/0.553)
sources = source.Source(m_1=m1, m_2=m2, ecc=ecc, dist=dist,f_orb=f_orb)
snr = sources.get_snr(verbose=True)


fig, axs = plt.subplots(2, 2, figsize=(12,8), dpi=200)
fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.25,hspace=0.3)

axs[0][0].set_xlabel(xlabel=r'$\log_{10}(\rm t~/~{\rm yr})$',fontproperties=font_prop_label)
axs[0][0].set_ylabel(ylabel=r'$\rm f_{GW}(mHZ)$',fontproperties=font_prop_label)
axs[0][0].plot(Age,f_orb*1000)
axs[0][0].set_xlim(5.5,8.0)
axs[0][0].tick_params(axis='both',which='major',direction='in',width=1.2,length=10,top='on',right='on')
axs[0][0].tick_params(axis='both',which='minor',direction='in',width=1.0,length=5,top='on',right='on')
axs[0][0].xaxis.set_minor_locator(AutoMinorLocator(n=4))
axs[0][0].yaxis.set_minor_locator(AutoMinorLocator(n=4))
axs[0][0].text(7.6, 7.9, '(a)', fontproperties=font_prop_text, color='black')

axs[0][1].set_xlabel(xlabel=r'$\log_{10}(\rm t~/~{\rm yr})$',fontproperties=font_prop_label)
axs[0][1].set_ylabel(ylabel=r'$\rm M_{chirp}(\rm M_\odot)$',fontproperties=font_prop_label)
axs[0][1].plot(Age,m_chrip)
axs[0][1].set_xlim(5.5,8.0)
axs[0][1].set_ylim(0,0.5)
axs[0][1].tick_params(axis='both',which='major',direction='in',width=1.2,length=6,top='on',right='on')
axs[0][1].tick_params(axis='both',which='minor',direction='in',width=1.0,length=3,top='on',right='on')
axs[0][1].tick_params(axis='both',which='major',direction='in',width=1.2,length=10,top='on',right='on')
axs[0][1].tick_params(axis='both',which='minor',direction='in',width=1.0,length=5,top='on',right='on')
axs[0][1].xaxis.set_minor_locator(AutoMinorLocator(n=4))
axs[0][1].yaxis.set_minor_locator(AutoMinorLocator(n=4))
axs[0][1].text(7.6, 0.42, '(b)', fontproperties=font_prop_text, color='black')

axs[1][0].set_xlabel(xlabel=r'$\log_{10}(\rm t~/~{\rm yr})$',fontproperties=font_prop_label)
axs[1][0].set_ylabel(ylabel=r'$SNR$',fontproperties=font_prop_label)
axs[1][0].plot(Age,snr)
axs[1][0].set_xlim(5.5,8.0)
axs[1][0].set_yscale('log')
axs[1][0].tick_params(axis='both',which='major',direction='in',width=1.2,length=6,top='on',right='on')
axs[1][0].tick_params(axis='both',which='minor',direction='in',width=1.0,length=3,top='on',right='on')
axs[1][0].tick_params(axis='both',which='major',direction='in',width=1.2,length=10,top='on',right='on')
axs[1][0].tick_params(axis='both',which='minor',direction='in',width=1.0,length=5,top='on',right='on')
axs[1][0].xaxis.set_minor_locator(AutoMinorLocator(n=4))
axs[1][0].text(7.6, 800, '(c)', fontproperties=font_prop_text, color='black')
axs[1][0].text(5.75, 10, r'$LISA$', fontproperties=font_prop_text, color='red')
axs[1][0].axhline(y=7, color='red', linestyle='--', linewidth=2, alpha=0.7, label='y=7')


axs[1][1].set_xlabel(xlabel=r'$\rm f_{GW}(mHZ)$',fontproperties=font_prop_label)
axs[1][1].set_ylabel(ylabel=r'$\rm M_d(\rm M_\odot)$',fontproperties=font_prop_label)
axs[1][1].plot(f_orb*1000,m1)
axs[1][1].set_xlim(0,12)
axs[1][1].set_ylim(0,0.3)
axs[1][1].tick_params(axis='both',which='major',direction='in',width=1.2,length=6,top='on',right='on')
axs[1][1].tick_params(axis='both',which='minor',direction='in',width=1.0,length=3,top='on',right='on')
axs[1][1].tick_params(axis='both',which='major',direction='in',width=1.2,length=10,top='on',right='on')
axs[1][1].tick_params(axis='both',which='minor',direction='in',width=1.0,length=5,top='on',right='on')
axs[1][1].xaxis.set_minor_locator(AutoMinorLocator(n=4))
axs[1][1].yaxis.set_minor_locator(AutoMinorLocator(n=4))
axs[1][1].text(10.1, 0.25, '(d)', fontproperties=font_prop_text, color='black')
plt.savefig('GW.png',dpi=500,bbox_inches='tight')
plt.show()
