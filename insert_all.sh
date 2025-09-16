#!/usr/bin/env
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 16:28:14 2022

@author: cjq
"""
for i in {0,1,2}
do
	for j in {0,1,2,3,4,5,6,7,8,9}
	do
		for k in {0,1,2,3,4,5,6,7,8,9}
		do
			python yt_plot_snapshot_insert.py -f snapshot_$i$j$k.hdf5
		done
	done
done
	    
