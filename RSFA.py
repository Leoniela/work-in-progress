# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 15:17:28 2016

@author: lampe
"""

#from nipype.pipeline.engine import Node, Workflow
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl

import nipype.pipeline.engine as pe 
from nipype.interfaces.utility import Function
#from nipype.interfaces.nipy.preprocess import Trim
#from nipype.algorithms import misc

#from nilearn.image import resample_img
#import numpy
#import csv
#import pandas as pd


#def calc_RSFA(subject, working_dir,data_dir, out_dir,TR_ms, low_hp_cutoff_freq, low_lp_cutoff_freq, 
#              high_hp_cutoff_freq, high_lp_cutoff_freq,
#              standard_brain, standard_brain_resampled, standard_brain_mask, 
#              standard_brain_mask_resampled, fwhm_smoothing, epi_resolution,):
                  

def calc_RSFA(subject, 
              working_dir,
              data_dir, 
              result_dir,
              TR_ms, 
              low_hp_cutoff_freq, 
              low_lp_cutoff_freq, 
              high_hp_cutoff_freq, 
              high_lp_cutoff_freq):
    
                  
                  
                  
                  
    fsl.FSLCommand.set_default_output_type('NIFTI_GZ')
    
   # working_dir = '/scr/kennedy2/lampe/RSFA_ALFF/data/'+subject+'/'
    #os.makedirs(working_dir)
   # result_dir = '/scr/kennedy2/lampe/RSFA_ALFF/results/' +subject + '/'
  #  data_dir = '/data/liem-1/LIFE/preprocessed/'+subject+ '/'
        #os.makedirs(data_dir)
        ##change this depending on new and old freesurfer##
        #freesurfer_dir = '/scr/kennedy2/LIFE/freesurfer/subjects/' #  #' # ##
    #resting_dir = '/scr/kennedy2/lampe/RSFA_ALFF/results/'+subject+'/test/'
    #standard_brain = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'
    #standard_brain_resampled = '/scr/kennedy2/lampe/RSFA_ALFF/scripts/MNI_resampled.nii'
    #standard_brain_mask = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'
    #standard_brain_mask_resampled='/scr/kennedy2/lampe/RSFA_ALFF/scripts/MNI_resampled_brain_mask.nii'    
    
    
    
    
    mask_dir = '/scr/kennedy2/lampe/Alle_zusammengefasst/'

    TOADS_lesion_def = '_clone_transform_clone_reg_N3Corrected1_mask_cp_strip_durastripped_N3Corrected_clone_lToads_lesions_seg_def.nii.gz'
    TOADS_lesion =  '_clone_transform_clone_reg_N3Corrected1_mask_cp_strip_durastripped_N3Corrected_clone_lToads_lesions_seg.nii.gz'
    TOADS_full_segmentation = '_clone_transform_clone_reg_N3Corrected1_mask_cp_strip_durastripped_N3Corrected_clone_lToads_seg.nii.gz'
    
    templates={
    #'func': 'func/EPI_t2.nii',
    #'fmap_fullwarp' : 'unwarp/B0_ph.nii',
    #'mat_moco': '',
    #'transform_ts': 'lemon_resting/transform_timeseries/merge/rest2anat.nii.gz'
    
    'MNI_lesionmask' : mask_dir + 'deformed_lesions/' + subject + TOADS_lesion_def,
    'ind_lesionmask' : mask_dir + 'lesions/' + subject + TOADS_lesion,
    'ind_TOADS_segmentation' : mask_dir + 'lesions/' + subject + TOADS_full_segmentation,
    'ind_flair_TOADS' : mask_dir + 'flair_input/' + subject + '_spc_da-fl_irprep_sag_p2_iso_395_clone_transform_clone_reg_calc_N3Corrected_calc_filter.nii.gz',
    'ind_t1_TOADS' : mask_dir + 't1_input/' + subject + '*.nii.gz',
    'anat_head' : 'structural/T1.nii.gz', #either with mod or without
    'anat_brain' : 'structural/brain.nii.gz', #new version with brain_extraction from freesurfer  #T1_brain_brain.nii.gz',
    'brain_mask' : 'structural/T1_brain_mask.nii.gz', #T1_brain_brain_mask.nii.gz',
    'ants_affine': 'structural/transforms2mni/transform0GenericAffine.mat',
    'ants_warp':   'structural/transforms2mni/transform1Warp.nii.gz',
    'transform_ts': 'resting_state/coregister/rest_coregistered_nativespace.nii.gz' #nur coregistered, 
                                                                                #McFlirted (motion corrected, alle volumes aligned), 
                                                                                #3x3x3 downgesampelt, 5 volumes removed, kein slice time correction
                                                                                #kein denoising gemacht
        }
        
    RSFA = pe.Workflow(name='RSFA')
    RSFA.base_dir = working_dir
    RSFA.config['execution']['crashdump_dir'] = RSFA.base_dir + "/crash_files"
    
    selectfiles = pe.Node(nio.SelectFiles(templates,
                                          base_directory=data_dir), # data_dir = '/data/liem-1/LIFE/preprocessed/' +subject + '/' set in Metascript runRSFA
                                          name="Selectfiles")
                                          
   #      TR_ms = 2000

# PARAMETERS CUT OFF low RSFA
#low_hp_cutoff_freq = 0.01
#low_lp_cutoff_freq = 0.1

# PARAMETERS CUT OFF high RSFA
#high_hp_cutoff_freq = 0.1
#high_lp_cutoff_freq = 0.25


#    inputnode = pe.Node(interface=util.IdentityInterface(fields=['subject_id',
#                                                              'TR_ms',
#                                                              'low_lp_cutoff_freq',
#                                                              'low_hp_cutoff_freq',
#                                                              'high_lp_cutoff_freq',
#                                                              'high_hp_cutoff_freq']),
#                     name='inputnode')            

#this step is redundant, because Franz has already discarded the first 5 volumes

#Workflow to delete first five volume
#trim = pe.Node(Trim(), name ='Trim')
#trim.inputs.begin_index = 5 # remove 5 first volumes

#RSFA.connect(selectfiles, 'transform_ts', trim, 'in_file')                     

    # TR CONVERSION
    def get_TR_in_sec_fct(TR_ms):
        return TR_ms / 1000.0
        
        
        
    get_TR_in_sec = pe.Node(util.Function(input_names=['TR_ms'],
                                       output_names=['TR_sec'],
                                       function=get_TR_in_sec_fct),
                         name='get_TR_in_sec')
    get_TR_in_sec.inputs.TR_ms = TR_ms

    
    def calc_bp_sigma_fct(TR_sec, cutoff_freq):
        sigma = 1. / (2 * TR_sec * cutoff_freq)
        return sigma
        
        
   # calculate low and highpass sigmas     
        
    calc_low_lp_sigma = pe.Node(util.Function(input_names=['TR_sec', 'cutoff_freq'],
                                       output_names=['sigma'],
                                       function=calc_bp_sigma_fct), name='calc_low_lp_sigma')
    calc_low_lp_sigma.inputs.cutoff_freq = low_lp_cutoff_freq
                                       
    RSFA.connect(get_TR_in_sec, 'TR_sec', calc_low_lp_sigma, 'TR_sec')
    
    
    
    calc_low_hp_sigma = pe.Node(util.Function(input_names=['TR_sec', 'cutoff_freq'],
                                       output_names=['sigma'],
                                       function=calc_bp_sigma_fct), name='calc_low_hp_sigma')
    calc_low_hp_sigma.inputs.cutoff_freq = low_hp_cutoff_freq
    
    RSFA.connect(get_TR_in_sec, 'TR_sec', calc_low_hp_sigma, 'TR_sec')
    
    
    calc_high_lp_sigma = pe.Node(util.Function(input_names=['TR_sec', 'cutoff_freq'],
                                       output_names=['sigma'],
                                       function=calc_bp_sigma_fct), name='calc_high_lp_sigma')
    calc_high_lp_sigma.inputs.cutoff_freq = high_lp_cutoff_freq
    
    RSFA.connect(get_TR_in_sec, 'TR_sec', calc_high_lp_sigma, 'TR_sec')
    
    
    
    
    calc_high_hp_sigma = pe.Node(util.Function(input_names=['TR_sec', 'cutoff_freq'],
                                       output_names=['sigma'],
                                       function=calc_bp_sigma_fct), name='calc_high_hp_sigma')
    calc_high_hp_sigma.inputs.cutoff_freq = high_hp_cutoff_freq
    
    RSFA.connect(get_TR_in_sec, 'TR_sec', calc_high_hp_sigma, 'TR_sec')
                  
                  # bp filters
                  
    tempfilter_low = pe.Node(fsl.maths.TemporalFilter(), name= 'Low_filter')
    RSFA.connect(selectfiles, 'transform_ts', tempfilter_low, 'in_file')
    RSFA.connect(calc_low_lp_sigma, 'sigma', tempfilter_low, 'lowpass_sigma')
    RSFA.connect(calc_low_hp_sigma, 'sigma', tempfilter_low, 'highpass_sigma')
    
    tempfilter_high = pe.Node(fsl.maths.TemporalFilter(), 'High_filter')
    RSFA.connect(selectfiles, 'transform_ts', tempfilter_high, 'in_file')
    RSFA.connect(calc_high_lp_sigma, 'sigma', tempfilter_high, 'lowpass_sigma')
    RSFA.connect(calc_high_hp_sigma, 'sigma', tempfilter_high, 'highpass_sigma')
    
    stdev_low = pe.Node(fsl.maths.StdImage(), name= 'Low_std')
    RSFA.connect(tempfilter_low, 'out_file', stdev_low, 'in_file')
    
    stdev_high = pe.Node(fsl.maths.StdImage(), name= 'High_std')
    RSFA.connect(tempfilter_high, 'out_file', stdev_high, 'in_file')
    
    # register RSFA to t1
    
    flirt_low = pe.Node(interface=fsl.FLIRT(), name= 'Flirt_low')
    flirt_low.inputs.cost_func = 'mutualinfo'
    flirt_low.inputs.dof = 6
    flirt_low.inputs.out_matrix_file = "trans_matrix.mat"
    #flirt_low.inputs.reference = t1_resampled
    RSFA.connect(stdev_low, 'out_file', flirt_low, 'in_file')
    RSFA.connect(selectfiles, 'ind_t1_TOADS', flirt_low, 'reference')
    
    flirt_high = pe.Node(interface=fsl.FLIRT(), name= 'Flirt_high')
    flirt_high.inputs.cost_func = 'mutualinfo'
    flirt_high.inputs.dof = 6
    flirt_high.inputs.out_matrix_file = "trans_matrix.mat"
    #flirt_high.inputs.reference = t1_resampled
    RSFA.connect(stdev_high, 'out_file', flirt_high, 'in_file')
    RSFA.connect(selectfiles, 'ind_t1_TOADS', flirt_high, 'reference')
    
    applyxfm_low = pe.Node(interface=fsl.ApplyXfm(), name= 'Transfer_low')
    applyxfm_low.inputs.interp = 'nearestneighbour'
    applyxfm_low.inputs.apply_xfm = True
    applyxfm_low.inputs.terminal_output = 'file'
    #applyxfm_low.inputs.reference = temp_t1
    RSFA.connect(stdev_low, 'out_file', applyxfm_low, 'in_file')
    RSFA.connect(flirt_low, 'out_matrix_file', applyxfm_low, 'in_matrix_file')
    RSFA.connect(selectfiles, 'ind_t1_TOADS', applyxfm_low, 'reference')
    
    applyxfm_high = pe.Node(interface=fsl.ApplyXfm(), name= 'Transfer_high')
    applyxfm_high.inputs.interp = 'nearestneighbour'
    applyxfm_high.inputs.apply_xfm = True
    applyxfm_high.inputs.terminal_output = 'file'
    #applyxfm_high.inputs.reference = temp_t1
    RSFA.connect(stdev_high, 'out_file', applyxfm_high, 'in_file')
    RSFA.connect(flirt_high, 'out_matrix_file', applyxfm_high, 'in_matrix_file')
    RSFA.connect(selectfiles, 'ind_t1_TOADS', applyxfm_high, 'reference')
    
    # gray matter 
    gm_mask = pe.Node(fsl.ImageMaths(op_string = '-thr 14.5 -uthr 15.5 -bin'), name= 'GM_mask')
    RSFA.connect(selectfiles, 'ind_TOADS_segmentation', gm_mask, 'in_file')
    
    gm_low = pe.Node(fsl.maths.ApplyMask(), name= 'GM_low')
    RSFA.connect(applyxfm_low, 'out_file', gm_low, 'in_file')
    RSFA.connect(gm_mask, 'out_file', gm_low, 'mask_file')
    
    gm_high = pe.Node(fsl.maths.ApplyMask(), name= 'GM_high')
    RSFA.connect(applyxfm_high, 'out_file', gm_high, 'in_file')
    RSFA.connect(gm_mask, 'out_file', gm_high, 'mask_file')
    
    gm_low_mean = pe.Node(interface=fsl.ImageStats(op_string = '-M'), name = 'GM_low_mean')
    RSFA.connect(gm_low, 'out_file', gm_low_mean, 'in_file')
    
    gm_high_mean = pe.Node(interface=fsl.ImageStats(op_string = '-M'), name = 'GM_high_mean')
    RSFA.connect(gm_high, 'out_file', gm_high_mean, 'in_file')
    
    #white matter 
    wm_mask = pe.Node(fsl.ImageMaths(op_string = '-thr 24.5 -uthr 25.5 -bin'), name= 'WM_mask')
    RSFA.connect(selectfiles, 'ind_TOADS_segmentation', wm_mask, 'in_file')
    
    wm_low = pe.Node(fsl.maths.ApplyMask(), name= 'WM_low')
    RSFA.connect(applyxfm_low, 'out_file', wm_low, 'in_file')
    RSFA.connect(wm_mask, 'out_file', wm_low, 'mask_file')
    
    wm_high = pe.Node(fsl.maths.ApplyMask(), name= 'WM_high')
    RSFA.connect(applyxfm_high, 'out_file', wm_high, 'in_file')
    RSFA.connect(wm_mask, 'out_file', wm_high, 'mask_file')
    
    wm_low_mean = pe.Node(interface=fsl.ImageStats(op_string = '-M'), name = 'WM_low_mean')
    RSFA.connect(wm_low, 'out_file', wm_low_mean, 'in_file')
    
    wm_high_mean = pe.Node(interface=fsl.ImageStats(op_string = '-M'), name = 'WM_high_mean')
    RSFA.connect(wm_high, 'out_file', wm_high_mean, 'in_file')
    
    #WMH
    wmh_low = pe.Node(fsl.maths.ApplyMask(), name= 'WMH_low')
    RSFA.connect(applyxfm_low, 'out_file', wmh_low, 'in_file')
    RSFA.connect(selectfiles, 'ind_lesionmask', wmh_low, 'mask_file')
    
    wmh_high = pe.Node(fsl.maths.ApplyMask(), name= 'WMH_high')
    RSFA.connect(applyxfm_high, 'out_file', wmh_high, 'in_file')
    RSFA.connect(selectfiles, 'ind_lesionmask', wmh_high, 'mask_file')
    
    wmh_low_mean = pe.Node(interface=fsl.ImageStats(op_string = '-M'), name = 'WMH_low_mean')
    RSFA.connect(wmh_low, 'out_file', wmh_low_mean, 'in_file')
    
    wmh_high_mean = pe.Node(interface=fsl.ImageStats(op_string = '-M'), name = 'WMH_high_mean')
    RSFA.connect(wmh_high, 'out_file', wmh_high_mean, 'in_file')
    
    
    def store2file(SIC, 
               RSFA_l_wm, RSFA_h_wm , 
               RSFA_l_gm, RSFA_h_gm, 
               RSFA_l_wmh, RSFA_h_wmh):
        import pandas as pd
        import csv
    
        d = {'0subject' : [SIC],
                               '1RSFA_l_wm'  : [RSFA_l_wm],
                               '2RSFA_h_wm' : [RSFA_h_wm],
                               '3RSFA_l_gm' : [RSFA_l_gm],
                               '4RSFA_h_gm' : [RSFA_h_gm],
                               '5RSFA_l_wmh' : [RSFA_l_wmh],
                               '6RSFA_h_wmh' : [RSFA_h_wmh]
                            }
                            
        df = pd.DataFrame(d)
        
        with open('/scr/kennedy2/lampe/RSFA_ALFF/results/txt_files/RSFA_results.csv', 'a') as f:
            df.to_csv(f, header=False)
                            
                            
                            
    value_safe = pe.Node(name='RSFA_value_file', interface=Function(input_names=['SIC', 
                                                                                     'RSFA_l_wm', 'RSFA_h_wm',
                                                                                     'RSFA_l_gm', 'RSFA_h_gm',
                                                                                     'RSFA_l_wmh', 'RSFA_h_wmh'],
                                 output_names=["out_file"],
                                 function=store2file))
    
    value_safe.inputs.SIC = subject
    RSFA.connect(wm_low_mean, 'out_stat', value_safe, 'RSFA_l_wm')
    RSFA.connect(wm_high_mean, 'out_stat', value_safe, 'RSFA_h_wm')
    RSFA.connect(gm_low_mean, 'out_stat', value_safe, 'RSFA_l_gm')
    RSFA.connect(gm_high_mean, 'out_stat', value_safe, 'RSFA_h_gm')
    RSFA.connect(wmh_low_mean, 'out_stat', value_safe, 'RSFA_l_wmh')
    RSFA.connect(wmh_high_mean, 'out_stat', value_safe, 'RSFA_h_wmh')
    
    
    
    datasink = pe.Node(nio.DataSink(), name='sinker')
    datasink.inputs.base_directory = result_dir
    #datasink.remove_dest_dir = True
    #RSFA.connect(addrow, 'csv_file', datasink, 'container.txt_file')
    RSFA.connect(gm_low, 'out_file', datasink, 'container.low_freq_GM_mask')
    RSFA.connect(wmh_low, 'out_file', datasink, 'container.low_freq_lesion_mask')
    RSFA.connect(wm_low, 'out_file', datasink, 'container.low_freq_WM_mask')
    RSFA.connect(applyxfm_low, 'out_file', datasink, 'container.low_freq_RSFA')
    
    RSFA.connect(gm_high, 'out_file', datasink, 'container.high_freq_GM_mask')
    RSFA.connect(wmh_high, 'out_file', datasink, 'container.high_freq_lesion_mask')
    RSFA.connect(wm_high, 'out_file', datasink, 'container.high_freq_WM_mask')
    RSFA.connect(applyxfm_high, 'out_file', datasink, 'container.high_freq_RSFA')
    
    RSFA.run(plugin = 'CondorDAGMan') #plugin = "CondorDAGMan" 'linear'
