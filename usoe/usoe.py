'''
Created on May 2, 2023

@author: dan
'''

from open_ephys.analysis import Session
from open_ephys.analysis.formats import BinaryRecording
import numpy as np
import pandas as pd
import os

allowed_pulse_widths = np.array([0.003,0.006,0.010,0.015,0.020,0.025,0.030,0.035,0.040,0.045,0.050,0.055,0.060]);
allowed_file_markers_dict = {15: 'BoxOfDonuts', 20: 'freeview', 25: 'BakersDozen', 30: 'Orientation', 35: 'Contrast', 40: 'Area', 45: 'SF', 50: 'TF', 55: 'XY'}


def load_oe_events(directory):
    '''
    load all events (from binary inputs and MessageCenter). 
    directory should contain 'structure.oebin', and at least the 'events' folder and its contents
    '''
    br = BinaryRecording(directory)
    br.load_events()
    return br.events
    
def get_oe_pulses(ev, line):
    '''
    Find pulses (up event followed by down event) on given line. 
    Returns start time and width rounded to closest allowed width.
    '''
    
    
    # just the event line of interest
    evp = ev[ev.line == line]
    
    # in-loop vars
    b1 = False  # set to True when an up transition is found
    t0 = -1
    
    # arrays for result (temp)
    t = []
    w = []

    for index, row in evp.iterrows():
        #print(index, row.timestamp, row.state)
        if (row.state == 1):
            b1 = True
            t0 = row.timestamp
        elif row.state == 0:
            if b1:
                imin = np.argmin(abs(allowed_pulse_widths - (row.timestamp - t0)))
                w.append(int(allowed_pulse_widths[imin]*1000))
                t.append(t0)
                b1 = False
                t0 = -1
                
    # result
    data={ 'w': w, 't': t}
    return pd.DataFrame(data)

def get_trials_from_pulses(dfwt, t_start=-1, t_end=-1):
    ''' 
    find trial start/end times using pulses with width 3/6ms for start/end. 
    '''
    
    
    # deal with start and end times. filter the input in two stages, the end
    # result being df3 - use that on subsequent steps
    if t_start < 0:
        df2 = dfwt
    else:
        df2 = dfwt[dfwt.t >= t_start]
    if t_end < 0:
        df3 = df2
    else:
        df3 = df2[df2.t <= t_end]
    bstarted = False
    t0 = -1
    t0s = []  # trial start/end times
    t1s = []
    for index, row in df3.iterrows():
        if row.w==3:
            if bstarted:
                t0s.append(t0)
                t1s.append(row.t)
                bstarted = True
                t0 = row.t
            else:
                t0 = row.t
                bstarted = True
        elif row.w == 6:
            if bstarted:
                t0s.append(t0)
                t1s.append(row.t)
                bstarted = False
                t0 = -1
    
    #print('There are %d trials' % (len(t0s)))
    #return np.array([t0s, t1s])
    return pd.DataFrame({'t0': t0s, 't1': t1s})

def split_into_scripts(dfwt):
    '''
    Given the width/time dataframe (returned from get_oe_pulses), split into spike2 scripts, as 
    defined by the start and end file markers. 
    '''
    if isinstance(dfwt, pd.DataFrame):

        # strip away all trial markers (w==3,6)
        dfstripped = dfwt.loc[(dfwt.w != 3) & (dfwt.w != 6)]

        bInScript = False
        w = 0
        t0 = -1
        stypes = []
        t0s = []
        t1s = []
        for index, row in dfstripped.iterrows():
            if row.w == 10:
                if bInScript:
                    stypes.append(w)
                    t0s.append(t0)
                    t1s.append(row.t)
                    bInScript = False
                else:
                    print("Warning: Found file end marker (w=10) without start.")
                    bInScript = False
            else:
                if bInScript:
                    stypes.append(w)
                    t0s.append(t0)
                    t1s.append(row.t)
                w = row.w
                t0 = row.t
                bInScript = True
        # in case we were inside a script when the end came
        if bInScript:
            stypes.append(w)
            t0s.append(t0)
            t1s.append(NaN)
        
        # make df with script results
        data = {'type': stypes, 't0': t0s, 't1': t1s}
        
        return pd.DataFrame(data)
    else:
        return None

def get_trial_start_offset(dfwt, t0, t1):             
    df3 = dfwt.loc[(dfwt.w == 3) & (dfwt.t > t0) & (dfwt.t < t1)]
    return np.array(df3['t']) - t0

class ScriptData:
    def __init__(self, txtfile):

        # save dir and filename
        (root, ext) = os.path.splitext(txtfile)
        self.ext = ext
        self.base = os.path.basename(root)
        self.directory = os.path.dirname(root)
        
        # get file sync marker time and width from first line
        with open(txtfile) as f:
            first_line = f.readline().strip('\n')
        fl = first_line.split(',')
        self.t_start = float(fl[0])
        self.t_end = float(fl[1])
        self.marker_width = int(fl[2])
        
        if allowed_file_markers_dict.get(self.marker_width) == None:
            print("File start marker width %d not found in list." % (self.marker_width))
        
        self.df = pd.read_csv(txtfile, header=None, skiprows=1)
    
    def __str__(self): 
        return f'Sp2 {self.base}({self.marker_width},{self.script_type()}), length={self.sp2_length()}, {len(self.df)} trials.'
          
    def script_type(self):
        if allowed_file_markers_dict.get(self.marker_width) == None:
            return 'Unknown'
        else:
            return allowed_file_markers_dict[self.marker_width]
 
    def sp2_length(self):
        return self.t_end - self.t_start
                
    def num_trials(self):
        return len(self.df)
    
    def get_trial_start_offset(self):
        return np.array(self.df[1]) - self.t_start
       
    @staticmethod
    def load_sp2_results(dirname, extensions):
        results = []
        for files in os.listdir(dirname):
            if files.endswith(extensions):
                results.append(ScriptData(os.path.join(dirname, files)))
        return results