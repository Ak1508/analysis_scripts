#!/usr/bin/python3.6

import sys, re, pickle
import numpy as np
import pandas as pd

if (sys.argv[1] != 'hms' and sys.argv[1] != 'shms') :
    print ('Usage: python debug1190.py hms (shms)')
    sys.exit(1)

if (sys.argv[1] == 'hms')  : 
    crno   = 3
    rdno   = 1
    spec   = 'h'
    dbg_fn = 'tst1190_main_debug_hr1577.txt'
    slno   = [2, 4, 5, 7,  8, 10, 13, 14, 16, 17]
if (sys.argv[1] == 'shms') : 
    crno   = 6
    rdno   = 2
    spec   = 'p'
    dbg_fn = 'tst1190_main_debug_pr2484.txt'
    slno   = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

td   = {}
for index, slot in enumerate(slno) :
    lines = []
    tmp_lines = []
    print ('Processing: %s slot = %d ...' % (sys.argv[1], slot))
    lines.clear()
    df = open(dbg_fn, 'rt')
    parsing = False
    for line in df:
        if (index+1 < len(slno)) :
            if   ('CodaDecode:: loading bank %d  %d   1190' % (crno, slno[index]) in line) : parsing = True
            elif ('CodaDecode:: loading bank %d  %d   1190' % (crno, slno[index+1]) in line) : parsing = False
            elif parsing : lines.append(line.strip('\n'))
        elif (index+1 == len(slno)) :
            if   ('CodaDecode:: loading bank %d  %d   1190' % (crno, slno[index]) in line) : parsing = True
            elif ('CodaDecode::Calling roc_decode %d' % rdno in line) : parsing = False
            elif parsing : lines.append(line.strip('\n'))
    td['slot-%d' % slot] = lines

tdc_hdr_dict = {}
for slot, line_list in td.items() :
    tdc_hdr_dict[slot] = [l for l in line_list if 'Caen1190Module:: 1190 TDC HEADER' in l]

chip_id = [0, 1, 2, 3]
chip_tdc_hdr_dict = {}
for slot, tdc_hdr_list in tdc_hdr_dict.items() :
    chip_tdc_hdr_dict[slot] = {}
    for index, cpid in enumerate(chip_id) :
        chip_tdc_hdr_dict[slot]['chip-%d' % chip_id[index]] = [tdc_hdr for tdc_hdr in tdc_hdr_list if 'chip id = %d' % chip_id[index] in tdc_hdr]

event_id_dict = {}
bunch_id_dict = {}
for slot, chip_tdc_hdr_list in chip_tdc_hdr_dict.items() :
    event_id_dict[slot] = {}
    bunch_id_dict[slot] = {}
    for chip, event_list in chip_tdc_hdr_dict[slot].items() : 
        event_id_list = []
        bunch_id_list = []
        event_id_dict[slot][chip] = []
        bunch_id_dict[slot][chip] = []
        for index, event in enumerate(event_list) :
            event_id_list.append(float(re.findall(r'\d+', event_list[index])[-2]))
            bunch_id_list.append(float(re.findall(r'\d+', event_list[index])[-1]))
        event_id_dict[slot][chip] = np.asarray(event_id_list)
        bunch_id_dict[slot][chip] = np.asarray(bunch_id_list)

edf = pd.DataFrame.from_dict(event_id_dict)
edf.to_pickle(spec+'edf.pkl')
bdf = pd.DataFrame.from_dict(bunch_id_dict)
bdf.to_pickle(spec+'bdf.pkl')
