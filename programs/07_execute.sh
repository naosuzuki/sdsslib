#!/bin/sh

#PBS -N 'GAIADR3'
 
#PBS -l nodes=1
 
#PBS -V 
 
#PBS -l walltime=48:00:00 
 
/work/nsuzuki/anaconda3/bin/python /work/nsuzuki/github/sdsslib/programs/06_sdssgaia_starquasar_plot_pipeline.py
 
 
