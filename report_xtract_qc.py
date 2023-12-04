#!/usr/bin/env fslpython

# Report XTRACT QC metrics as pretty graphs
#
# Authors: Shaun Warrington
#
# Copyright (C) 2020 University of Oxford
# SHBASECOPYRIGHT

import sys,os,glob
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False

import jinja2
import seaborn as sns
import plotly
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from datetime import datetime

subject_list = sys.argv[1]
out          = sys.argv[2]
tracts       = sys.argv[3]
n_std        = float(sys.argv[4]) # number of standard deviations to consider an outlier
missing_list = sys.argv[5] # subjects missing the XTRACT folder
temp_loc     = sys.argv[6] # location of html report template

tracts=tracts.split(",")

# load subject list and remove any missing the xtract folder
subject_list = np.loadtxt(subject_list, unpack=False, dtype='str')

if os.path.isfile(missing_list):
    missing_list = np.loadtxt(missing_list, unpack=False, dtype='str')
    subject_list = np.setdiff1d(subject_list, missing_list)

ntracts = len(tracts)
nsubs = len(subject_list)

def get_stats(subject_list, tracts, filename): # list of subject, list of tracts, QC metric filename (e.g. waytotal, filecheck)
    metric = np.zeros((len(subject_list), len(tracts)))
    for idx,subpath in enumerate(subject_list):
        temp = np.loadtxt(os.path.join(subpath, 'xtract_qc', filename), dtype='int', usecols=(1))
        metric[idx,:] = temp
    return metric

filecheck = get_stats(subject_list, tracts, 'filecheck')
waytotal = pd.DataFrame(get_stats(subject_list, tracts, 'waytotal'), columns=tracts)
volume = pd.DataFrame(get_stats(subject_list, tracts, 'volume'), columns=tracts)

##########################################
# now report subject's missing results and give percent success for each tract and flag any subjects with all missing
# working with filecheck metric
sumcheck = np.sum(filecheck, axis=1)
sub_missing = np.argwhere(sumcheck < ntracts) # subjects missing data

sub_missing_txt = []
sub_missing_all_txt = []
for idx in sub_missing:
    sub = subject_list[idx]
    n_missing = int(ntracts - sumcheck[idx])
    if n_missing < ntracts:
        sub_missing_txt.append(f'{sub[0]} missing {n_missing} tracts')
    elif n_missing == ntracts:
        sub_missing_all_txt.append(f'{sub[0]} missing all tracts')

# get overall success
sumcheck_all = np.round((np.sum(sumcheck) / (ntracts * nsubs)) * 100, 3) # now overall percent success across subjects and tracts
sumcheck_all_txt = f'{sumcheck_all}% ({int(np.sum(sumcheck))} of {int(ntracts * nsubs)})'
tract_success = np.round((np.sum(filecheck, axis=0) / nsubs) * 100, 3) # percent of subjects successfully completed per tract

# report simple filecheck metric
ps_inc = 0
if sumcheck_all < 100:
    ps_inc = 1
    tract_success = pd.DataFrame({'tracts': tracts, 'percent_successful':tract_success})
    fig = px.line(tract_success, x=tracts, y="percent_successful",
        labels=dict(tracts="Tract", percent_successful="Percent successful (%)", )
        )
    fig.add_shape(type="line", x0=0, y0=100, x1=ntracts, y1=100, line=dict(color="orange", dash="dash"))
    fig.update_xaxes(title_text="Tract")
    fig.write_html(os.path.join(out, 'percent_success_fig.html'))

##########################################
# make metric distribution plots
def myround(x, base=5):
    y = base * round(x/base)
    if y < x:
        y = y + base
    return y

def plot_tractwise_metric(metric, metric_name, n_std):
    k = 4
    m = int(myround(ntracts)/k)
    fig, ax = plt.subplots(m, k, figsize=(16, m*1.75))
    n = 0
    for j in range(m):
        for i in range(k):
            if n < ntracts:
                # plot distribution
                x = metric.values[:,n]
                sns.kdeplot(x, ax=ax[j,i], fill=False, color='crimson')
                kdeline = ax[j,i].lines[0]
                xs = kdeline.get_xdata()
                ys = kdeline.get_ydata()
                middle = x.mean()
                median = np.median(x)
                sdev = x.std()
                left = middle - sdev
                right = middle + sdev
                ax[j,i].vlines(middle, 0, np.interp(middle, xs, ys), color='crimson', ls=':')
                ax[j,i].vlines(median, 0, np.interp(median, xs, ys), color='crimson', ls='--')
                ax[j,i].fill_between(xs, 0, ys, facecolor='crimson', alpha=0.5)
                ax[j,i].fill_between(xs, 0, ys, where=(left <= xs) & (xs <= right), interpolate=True, facecolor='crimson', alpha=0.5)
                # add markers for n_std standard deviations
                left = middle - n_std*sdev
                right = middle + n_std*sdev
                ax[j,i].axvline(left, color='k', ls='--')
                ax[j,i].axvline(right, color='k', ls='--')
                # set axis options
                ax[j,i].set_xlim(left=0)
                ax[j,i].set_title(tracts[n], size=18, y=1.05)
                ax[j,i].ticklabel_format(style='sci', scilimits=(0, 0))
            else:
                ax[j,i].axis('off')
            n += 1
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'{metric_name}_fig.jpg'))
    plt.close()

plot_tractwise_metric(waytotal, 'waytotal', n_std)
plot_tractwise_metric(volume, 'volume', n_std)


##########################################
# list of subjects with waytotal=0, but have tract files
# these subjects are likely to have failed due to bad registration
sub_ind_zerowt = np.argwhere((filecheck == 1) & (waytotal.values == 0))[:,0]
sub_ind_zerowt = np.unique(sub_ind_zerowt)

sub_zerowt_txt = []
for idx in sub_ind_zerowt:
    sub = subject_list[idx]
    nempty = np.where(waytotal.values[idx,:] == 0)[0].shape[0]
    sub_zerowt_txt.append(f'{sub} has empty tracts files for {nempty}')


# find outliers based on waytotal and volume
# add here a string of subjects where their volume/waytotal is >2*std away from mean
def get_outliers(metric):
    outliers = np.abs(metric-metric.mean()) >= (n_std*metric.std())
    outliers = np.logical_and(outliers, (metric > 0)) # subject > n_std away from mean but not 0
    outliers_summary = np.sum(outliers) # tract-wise outliers
    extreme_outliers = np.sum(outliers, axis=1) # now getting subject with more than 0.5*ntracts
    extreme_outliers = extreme_outliers >= np.rint(0.5*ntracts)
    return outliers, outliers_summary, extreme_outliers

# make outlier plots
fig = make_subplots(rows=2, cols=1,
    subplot_titles=('Waytotal',  'Volume'), shared_xaxes=True,
    specs=[[{"secondary_y": True}], [{"secondary_y": True}]])

# waytotal outlier plot
wt_outliers, outliers_summary, wt_extreme_outliers = get_outliers(waytotal)
outliers_summary = pd.DataFrame({'tracts': tracts, 'no_outliers':outliers_summary.values, 'perc':np.round(100*outliers_summary.values/nsubs,2)})
fig.add_trace(go.Scatter(x=outliers_summary['tracts'], y=outliers_summary['no_outliers'], name="Number", marker=dict(color='darkslategray')), row=1, col=1, secondary_y=False);

fig.add_trace(go.Scatter(x=outliers_summary['tracts'], y=outliers_summary['perc'], name="Percent", marker=dict(color='darkslategray')), row=1, col=1, secondary_y=True);

fig.update_layout(showlegend=False, hovermode="x unified");
fig.update_yaxes(title_text="# potential outliers", secondary_y=False);
fig.update_yaxes(title_text="% potential outliers", secondary_y=True);

# volume outlier plot
vol_outliers, outliers_summary, vol_extreme_outliers = get_outliers(volume)
outliers_summary = pd.DataFrame({'tracts': tracts, 'no_outliers':outliers_summary.values, 'perc':np.round(100*outliers_summary.values/nsubs,2)})
fig.add_trace(go.Scatter(x=outliers_summary['tracts'], y=outliers_summary['no_outliers'], name="Number", marker=dict(color='darkslategray')), row=2, col=1, secondary_y=False);

fig.add_trace(go.Scatter(x=outliers_summary['tracts'], y=outliers_summary['perc'], name="Percent", marker=dict(color='darkslategray')), row=2, col=1, secondary_y=True);

fig.update_layout(showlegend=False, hovermode="x unified")
fig.update_yaxes(title_text="# potential outliers", secondary_y=False)
fig.update_yaxes(title_text="% potential outliers", secondary_y=True)

fig.update_xaxes(showticklabels=False) # hide all the xticks
fig.update_xaxes(title_text="Tract", showticklabels=True, row=2, col=1, tickmode='linear')

fig.write_html(os.path.join(out, 'outlier_fig.html'))

# build outlier report txt/csv
# extreme outliers (those with more than 0.5*ntracts as outliers)
grot = np.argwhere(vol_extreme_outliers.values == 1)
grot = np.append(grot, np.argwhere(wt_extreme_outliers.values == 1))
grot = np.unique(grot)

sub_extreme_outliers_txt = []
if grot.shape[0] > 0:
    f = open(os.path.join(out, 'extreme_outliers.txt'), 'w')
    for idx in grot:
        sub_extreme_outliers_txt.append(f'{subject_list[idx]}\n')
        f.write(f'{subject_list[idx]}\n')
    f.close()

# if subject is n_std away from mean for waytotal/volume
# save outliers as .csv
sub_outliers_path = os.path.join(out, 'tract_outliers.csv')
wt_outliers.set_index(subject_list, inplace=True)
wt_outliers = wt_outliers.astype("int")
wt_outliers.to_csv(sub_outliers_path)
n_outliers = np.sum(np.sum(wt_outliers))

########################################## Create report
env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath=temp_loc))
template = env.get_template('template_qc_report.html')
filename = os.path.join(out, 'qc_report.html')
with open(filename, 'w') as fh:
    fh.write(template.render(
        time = datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
        ntracts = ntracts,
        nsubs = nsubs,
        n_std = n_std,
        sumcheck_all_txt = sumcheck_all_txt,
        sub_missing_txt = sub_missing_txt,
        sub_missing_all_txt = sub_missing_all_txt,
        sub_zerowt_txt = sub_zerowt_txt,
        sub_outliers_path = sub_outliers_path,
        n_outliers = n_outliers,
        sub_extreme_outliers_txt = sub_extreme_outliers_txt,
        ps_inc = ps_inc
        ))
