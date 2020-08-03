#!/usr/bin/env python3.7

# Plot crude fatality ratio v.s. infected cases with day by day track.

# Usage:

# 1) Download data and save it to a file 'covid_cases.csv'.
#$   wget https://opendata.ecdc.europa.eu/covid19/casedistribution/csv -O covid_cases.csv
#    See also: https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide

# 2) Run this script. Plots will be saved in a directory called pic_tmp.
#$   python readcovid.py

# 3) make a video by ffmpeg (for example).
#$   ffmpeg -framerate 1/2 -start_number 0019 -i 'pic_tmp/pcum_d%04d.png' -r 25 -c:v libx264 -pix_fmt yuv420p out.mp4
# Or
#    ffmpeg -framerate 1 -pattern_type glob -i 'pic_tmp/*.png' -c:v libx264 -pix_fmt yuv420p out.mp4
# If audio is also need, add parameter '-i audio.wav'

# By https://github.com/bewantbe

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import csv

from numpy import nonzero, zeros, linspace, interp, mod, log10

# Load Chinese-English translation table of countries and territories.
import json
countrycode = json.loads(open('countrycode.json').read())
country_en2ch = {c['en']:c['cn'] for c in countrycode}

# Where and how the picture is saved.
def pic_output(st):
    fpath = 'pic_tmp/' + st + '.png'
    plt.savefig(fpath)

# Set plot defaults.
# See mpl.rcParams for full list.
mpl.rc('figure', figsize=(1600/240.0, 1200/240.0), dpi=240)
#mpl.rc('lines', markersize=2.0)
#mpl.rcParams['axes.formatter.limits'] = [-3,4]

font = {'family' : 'Adobe Fangsong Std', # 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 12}
mpl.rc('font', **font)

# Data file name.
fcovid = 'covid_cases.csv'
start_date = '2020-01-20'

st_date = []
li_epi = []
li_reg = []

# Load COVID infection data.
with open(fcovid, encoding='iso-8859-15') as f:
    reader = csv.reader(f)
    next(reader)   # skip first line (table header)
    for row in reader:
        st_date.append('%d-%02d-%02d'%(int(row[3]),int(row[2]),int(row[1])))
        li_epi.append([int(row[4]), int(row[5])])
        li_reg.append(row[6:10])

# Collect daily incremental data
# s_row is rows of: [date, [infect, dead], country]
s_row = [[st_date[k], np.array(li_epi[k]), li_reg[k][1]] for k in range(len(li_epi))]

# Get date with data
sorted_date = list(set(st_date))   # remove duplications
sorted_date.sort()

# Get list of all countries and territories.
dc_regime = {}  # indexed by regime abbreviation: [CN name, abbr, population]
for i in range(len(li_reg)):
    if li_reg[i][1] not in dc_regime:
        if li_reg[i][0].replace('_',' ') in country_en2ch.keys():
            li_reg[i][0] = country_en2ch[li_reg[i][0].replace('_',' ')]
        dc_regime[li_reg[i][1]] = [li_reg[i][0], li_reg[i][2], li_reg[i][3]]

# Compute cumulative [infect, dead] data from the incremental data s_row.
cum_by_day = []
cum_regime = {}
for j in range(len(sorted_date)):
    c_by_day = [l for l in s_row if l[0] == sorted_date[j]]  # TODO: can be improved.
    # Conceptually cum_regime.(geoCodes) += epinum
    for k in range(len(c_by_day)):
        if c_by_day[k][2] in cum_regime:
            cum_regime[c_by_day[k][2]] += c_by_day[k][1]
        else:
            cum_regime[c_by_day[k][2]] = c_by_day[k][1].copy()
    # deep copy data in cumulator cum_regime
    cum_by_day.append({k:cum_regime[k].copy() for k in cum_regime})

s_regime_name = list(dc_regime.keys())

# Get data matrix from list cum_by_day.
cum_mat_regime = zeros((len(dc_regime), len(sorted_date), 2), dtype='int64')
for j in range(len(sorted_date)):
    # cum_mat_regime[idx_geo, j, :] = cum_by_day[j]
    daily_geoid = cum_by_day[j].keys()
    for gi in daily_geoid:
        cum_mat_regime[s_regime_name.index(gi), j, :] = cum_by_day[j][gi]

print('cum_mat_regime.shape=', cum_mat_regime.shape)

"""
cum_mat_regime0 = cum_mat_regime
# cum_mat_regime, not correct due to skip data
cum_mat_regime = zeros((len(dc_regime), len(sorted_date), 2), dtype='int64')
for j in range(len(li_epi))[::-1]:  # assumed data in reversed order
    id_day = sorted_date.index(st_date[j])
    id_nat = s_regime_name.index(li_reg[j][1])
    if id_day == 0:
        cum_mat_regime[id_nat, id_day, 0] = li_epi[j][0]
        cum_mat_regime[id_nat, id_day, 1] = li_epi[j][1]
    else:
        cum_mat_regime[id_nat, id_day, 0] = li_epi[j][0] + cum_mat_regime[id_nat, id_day-1, 0]
        cum_mat_regime[id_nat, id_day, 1] = li_epi[j][1] + cum_mat_regime[id_nat, id_day-1, 1]
"""

trace_length = 14
off_day = 0  # shift between infected date and dead date (i.e. lag effect).
threshold_lowest_infected = 100

if 1:
    #plt.figure(100)
    fig, ax = plt.subplots()
    
    # trace line
    start_date_index = sorted_date.index(start_date) \
                         if start_date in sorted_date else 0
    for id_day in range(start_date_index,len(sorted_date)):
        print('id_day=', id_day)
        id_history_day = max(id_day - trace_length, 0)

        for id_nat in range(len(s_regime_name)):
            if cum_mat_regime[id_nat, id_day, 0] >= threshold_lowest_infected:
                infected = cum_mat_regime[id_nat, id_history_day:1+id_day, 0]
                infected_w1 = cum_mat_regime[id_nat, (id_history_day-off_day):(1+id_day-off_day), 0]
                dead     = cum_mat_regime[id_nat, id_history_day:1+id_day, 1]
                x = dead / infected_w1
                y = infected
                ax.plot(100*x, y, color=(0.8, 0.8, 0.8, 0.5), zorder=1)

        # construct colormap for dots
        n_cm = 7  # number of colors
        s_cm = mpl.cm.get_cmap('hsv')(np.array(range(n_cm))/n_cm) # base colormap
        cm3 = zeros((id_day-id_history_day+1, len(s_cm), 4))
        xi = linspace(0, 1, id_day-id_history_day+1)
        c0 = (0.8, 0.8, 0.8, 0.5)
        # the color of dots fade out to c0
        for id_c in range(s_cm.shape[0]):
            for cc in range(len(c0)):
                cm3[:,id_c,cc] = interp(xi, [0.0, 1.0], \
                  [c0[cc], s_cm[id_c, cc]] )

        # Marker size setting.
        s_mask_size = linspace(4, 20, id_day-id_history_day+1)
        s_mask_size[-1] = 40

        max_mortal_rate = 0.191
        
        # trace dot, from old day to latest
        for o_day in range(id_history_day, id_day+1):
            ids_regime, = nonzero(cum_mat_regime[:,o_day,0]>100)
            # draw only significant infected
            if len(ids_regime) == 0:
                continue
            infected = cum_mat_regime[ids_regime, o_day, 0]
            infected_w1 = cum_mat_regime[ids_regime, max(0, o_day - off_day), 0]
            dead     = cum_mat_regime[ids_regime, o_day, 1]
            x = dead / infected_w1
            y = infected
#            max_mortal_rate = min(max(np.max(x), max_mortal_rate), 1)
            cm_tmp = cm3[o_day-id_history_day, mod(ids_regime, n_cm), :]
            ax.scatter(100*x, y, s=s_mask_size[o_day-id_history_day], \
                       c=cm_tmp, zorder=2)

        # Print regime labels.
        
        # Fine tune the number of regime to show
        n_show_name = int(np.round(interp(id_day, [0, 70, 77, 85], [10,10,14,25] )))
        infected = cum_mat_regime[:, id_day, 0]
        infected_w1 = cum_mat_regime[:, max(0, id_day-off_day), 0]
        dead = cum_mat_regime[:, id_day, 1]
        # criteria(measure) for showing the regime name
        #interest_index = (log10(infected)*0.25)**2 + (10*(dead / infected_w1))**2 - 10000*(infected<100)
        interest_index = (log10(infected)*0.4)**2 + (10*(dead / infected_w1))**2 - 10000*(infected<100)
        #interest_index = (log10(infected)*0.4)**2 + (10*(dead / infected_w1))**2 - 10000*(infected<100)
        ids_interest = np.argsort(interest_index)[::-1]
        # remove NaN
        ids_interest = ids_interest[np.sum(np.isnan(interest_index)):] # remove Nan
        idx_high_interest = ids_interest[:min(len(ids_interest), n_show_name)]
        regime_interest = [s_regime_name[ii] for ii in idx_high_interest]
        regime_name_interest = [dc_regime[r][0] for r in regime_interest]
        for id_label in range(len(idx_high_interest)):
            x = cum_mat_regime[idx_high_interest[id_label], id_day, 1]
            y = cum_mat_regime[idx_high_interest[id_label], id_day, 0]
            y_w1 = cum_mat_regime[idx_high_interest[id_label], id_day-off_day, 0]
            #print('x=',x,' y=',y,' nat=',regime_name_interest[id_label])
            x = x / y_w1
            if y > 100:
                ax.text(100*x, y, regime_name_interest[id_label])

        ax.set_ylim([100, max(1e5, np.max(infected))])
        ax.set_yscale('log')
        ax.set_ylabel('感染人数', fontfamily='Adobe Fangsong Std', fontsize=16)
        ax.set_xlim([0, 100*max(0.1, 1.15*max_mortal_rate)])
        #ax.set_xlabel('crude fatality ratio', fontfamily='Times New Roman', fontsize=16)
        ax.set_xlabel('死亡(day)/感染(day)', fontfamily='Adobe Fangsong Std', fontsize=16)
        ax.xaxis.set_major_formatter(mpl.ticker.PercentFormatter())
        ax.set_title(sorted_date[id_day])
        #pic_output('pcum_%s'%(sorted_date[id_day]))
        pic_output('pcum_d%04d'%(id_day))
        ax.clear()
    
