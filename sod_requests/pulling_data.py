from obspy.fdsn import Client
from obspy import UTCDateTime
import matplotlib
import calendar
matplotlib.rcParams['backend'] = "TkAgg"
import matplotlib.pyplot as py
import numpy as np
import sys, os, os.path
import obspy
import obspy.signal
import scipy

search_start_time = UTCDateTime(2012,11,16,20,16,49)
search_end_time = UTCDateTime(2012,11,16,20,17,39)
search = client.get_waveforms("TA", "M44A", "--", "BHZ", search_start_time, search_end_time)



