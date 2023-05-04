'''
Created on May 2, 2023

@author: dan
'''

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import os
from open_ephys.analysis import Session
import numpy as np
import pandas as pd
from usoe.usoe import *
from cmath import nan
#from usoe.usoe import get_oe_pulses, split_into_scripts, ScriptData

def build_csd_df(scriptdata, events, oe_trials):
    '''
    for csd, find rising edge events in the events list, for each trial. 
    There should be 6 in each trial.
    The resulting df has columns 
    t0,t1 = start and end time of trial
    i0-i5 = frame times indicating (in order) the onset of: W,B,W; the 4th is the end of the last W; the 5th can be ignored. 
    '''
    nda = np.full((len(oe_trials), 8), -1, dtype='float')
    for index,row in oe_trials.iterrows():
        ev_frames = events.loc[(events.state==1) & (events.timestamp > row.t0) & (events.timestamp < row.t1)]
        oe_frame_times = np.pad(ev_frames['timestamp'], (0,6-len(ev_frames)), 'constant', constant_values=(-1,-1))
        #dfstripped = dfwt.loc[(dfwt.w != 3) & (dfwt.w != 6)]
        
        #print('Found %d events for trial %d (iR = %d)' % (len(ev_frames), index, s.df.iloc[index][8]))
        #print(ev_frames['timestamp'])
        #print(np.shape(frame_times_relative))
        #print(np.shape(nda[index,1:7]))
        nda[index,0] = row['t0']
        nda[index,1:7] = oe_frame_times
        nda[index,7] = row['t1']
        
    return pd.DataFrame({'t0': nda[:,0], 'i0': nda[:,1], 'i1': nda[:,2], 'i2': nda[:,3], 'i3': nda[:,4], 'i4': nda[:,5], 'i5': nda[:,6], 't1': nda[:,7], 'iR': s.df[8]})
        


def build_tuning_df(scriptdata, events, oe_trials):
    '''
    for the tuning case, need the stim on time and trial end time, as well as success/failure.
    Find the times relative to the trial start time (in sp2 file), then add that relative difference 
    to the oe trial start time. This isn't perfect but will minimize error due to clock drifts (TODO)
    '''
    
    # each row in this array is a trial
    # the 5 columns represent the 5 times extracted from each tuning trial
    # shape of this array should be (ntrials, 5)
    sp2_times = np.array(scriptdata.df.iloc[:,1:6])
    sp2_times[sp2_times < 0] = nan
    
    #print('sp2_times', type(sp2_times), np.shape(sp2_times))
    #print('sp2_times[:,0]', type(sp2_times[:,0]), np.shape(sp2_times[:,0]))
    #sp2_trial_start_times = scriptdata.df.iloc[:,1]
    oe_trial_start_times = np.array(oe_trials['t0'])
    #print(type(oe_trial_start_times), np.shape(oe_trial_start_times))
    
    # subtract off the trial start time
    # oe_trial_start_times is made from an array, so it has shape (ntrials,)
    oe_times = (sp2_times[:,1:].transpose() - sp2_times[:,0] + oe_trial_start_times).transpose()
    oe_times = np.nan_to_num(oe_times, nan=-1)
    
    # build df
    #return pd.DataFrame({'v': scriptdata.df.iloc[:, 0], 't0': oe_trial_start_times, 's0': oe_times.iloc[:, 0], 's1': oe_times.iloc[:, 1], 's2': oe_times.iloc[:, 2], 's3': oe_times.iloc[:, 3], 't1': scriptdata.df.iloc[:,1]})
    #print('oe_times ',type(oe_times), np.shape(oe_times))
    return pd.DataFrame({'v': scriptdata.df.iloc[:, 0], 't0': oe_trial_start_times, 's0': oe_times[:, 0], 's1': oe_times[:, 1], 's2': oe_times[:, 2], 't1': oe_trials['t1']})
    #print(oe_times)
   
def print_hline():
    print('\n', 40 * '-', '\n') 

if __name__ == '__main__':
    
    
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--root-dir", dest="rootdir", required=True, help='join(root,intermediate) contains struture.oebin')
    parser.add_argument("-i", "--intermediate-dir", dest="intdir", default="", required=False, help="intermediate folder, join(root,intermediate) should be the folder with structure.oebin")
    parser.add_argument("-s", "--sp2-dir", dest="sp2dir", required=True, help='folder with spike2 extract files *.bdx, *.inx')
    parser.add_argument("-o", "--output-dir", dest="outdir", required=False, help='output dir, default: join(root,intermediate)')

    # Process arguments
    args = parser.parse_args()
    structure_folder = os.path.join(args.rootdir, args.intdir)
    if args.outdir == None:
        output_folder = structure_folder
    else:
        output_folder = args.outdir

    print_hline()
    
    # load oe event data
    print(f'Load open ephys events from {structure_folder}')
    br = BinaryRecording(structure_folder)
    br.load_events()
    ev = br.events

    # get pulses from open ephys digital in data, line 2
    df_pulses = get_oe_pulses(ev, 2)

    # split pulses into scripts
    df_scripts = split_into_scripts(df_pulses)
    #print(df_scripts)
    
    # load results from spike2 side
    script_results = ScriptData.load_sp2_results(args.sp2dir, ('inx', 'bdx'))
    print_hline()
    print('Loaded spike2 data for these scripts:')
    for s in script_results:
        print(str(s))

    # now look for AnyMatches
    print_hline()
    print('match scripts found in oe data with spike2 data...')
    tolerance = 0.1   # how close is a match, in sec   
    for index, row in df_scripts.iterrows():
        df_oe_trials = get_trials_from_pulses(df_pulses, t_start=row['t0'], t_end=row['t1'])
        print_hline()
        print('oe script type: %s, length %f, %d trials' % (row['type'], row['t1']-row['t0'], len(df_oe_trials)))
        #print(np.shape(df_trials))
        b_found = False
        for s in script_results:
            if abs(s.sp2_length() - (row['t1']-row['t0'])) < tolerance :
                output_file = os.path.join(output_folder, f'oe_{s.base}.csv')
                print('Found a match\n%s' % str(s))
                print('Write csv to %s' % (output_file))
                
                # Get the sp2 data converted to oe time
                if s.marker_width == 20:
                    df20 = build_csd_df(s, ev, df_oe_trials)
                    df20.to_csv(output_file, index=False)
                elif s.marker_width == 30:
                    df30 = build_tuning_df(s, ev, df_oe_trials)
                    df30.to_csv(output_file, index=False)
                elif s.marker_width == 25:
                    df25 = build_tuning_df(s, ev, df_oe_trials)
                    #print(df25)
                    df25.to_csv(output_file, index=False)
                    
                b_found = True
                break
        if not b_found:
            print('NO MATCH FOUND\b\n')   

    print_hline()
    print('Done.')
