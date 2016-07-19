# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/raid1/lampe/.spyder2/.temp.py
"""

"""
Created on Mon Feb  9 12:27:06 2015

@author: fbeyer
"""

#from alff import create_alff
import sys
from RSFA import calc_RSFA
#import os
'''
Meta script to run lemon resting state preprocessing
---------------------------------------------------
Can run in two modes:
python run_alff.py s {subject_id}
python run_alff.py f {text file containing list of subjects}
'''
mode=sys.argv[1]
if mode == 's':
    subjects=[sys.argv[2]]
elif mode == 'f':
    with open(sys.argv[2], 'r') as f:
        subjects = [line.strip() for line in f]

for subject in subjects:
    print 'Running subject '+subject
    working_dir = '/scr/kennedy2/lampe/RSFA_ALFF/data/'+subject+'/'
    #os.makedirs(working_dir)
    data_dir = '/data/liem-1/LIFE/preprocessed/' + subject + '/'
    result_dir = '/scr/kennedy2/lampe/RSFA_ALFF/results/' +subject + '/'
    #os.makedirs(data_dir)
    ##change this depending on new and old freesurfer##
    #freesurfer_dir = '/scr/kennedy2/LIFE/freesurfer/subjects/' #  #' # ##
#    resting_dir = '/scr/kennedy2/lampe/RSFA_ALFF/results/'+subject+'/test/'
#    standard_brain = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'
#    standard_brain_resampled = '/scr/kennedy2/lampe/RSFA_ALFF/scripts/MNI_resampled.nii'
#    standard_brain_mask = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'
#    standard_brain_mask_resampled='/scr/kennedy2/lampe/RSFA_ALFF/scripts/MNI_resampled_brain_mask.nii'
    #os.makedirs(resting_dir)
    #os.makedirs(data_dir)
    #echo_space=0.00058 #in sec
    #te_diff=2.46 #in ms


    TR_ms = 2000  
    
    # PARAMETERS CUT OFF low RSFA
    low_hp_cutoff_freq = 0.01
    low_lp_cutoff_freq = 0.1

    # PARAMETERS CUT OFF high RSFA
    high_hp_cutoff_freq = 0.1
    high_lp_cutoff_freq = 0.25

    #epi_resolution = 3.0
    #TR=2.0
   # highpass=0.01
   # lowpass=0.1
    #vol_to_remove = 5
    #pe_dir = 'y-'
    #fwhm_smoothing = 6.0
    
    calc_RSFA(subject=subject,
              working_dir=working_dir,
              data_dir=data_dir, 
              result_dir=result_dir,
              TR_ms = TR_ms,
              low_hp_cutoff_freq =  low_hp_cutoff_freq,
              low_lp_cutoff_freq = low_lp_cutoff_freq,
              high_hp_cutoff_freq = high_hp_cutoff_freq,
              high_lp_cutoff_freq = high_lp_cutoff_freq)
                 
                 
#     highpass=highpass,
#    lowpass=lowpass,standard_brain = standard_brain, 
#    standard_brain_resampled = standard_brain_resampled, standard_brain_mask = standard_brain_mask,
#    standard_brain_mask_resampled = standard_brain_mask_resampled,
#    fwhm_smoothing = fwhm_smoothing)