'''
The MIT License (MIT)
Copyright 2020 John Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

# import modules
import os
import numpy as np
import pandas as pd

from scipy import stats
# plotting modules
from array import array
import matplotlib as mpl
#mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

def pdf_open(filename):
    return PdfPages(filename)

def seaborn_style(style, despine = True):

    sns.set_style(style)
    if despine:
        sns.despine()

def make_scatter( ax, x_vals, y_vals, title=None,xlabel=None,ylabel=None, xlim=None,ylim=None,colors=None,do_stats=None,alphas=None,labels=None, legend=False, r2_pos=None, mark_origin=False,xticks=None,yticks=None):
# makes a scatter plot and calculates a linear regression line
# returns the r-squared value
    if len(do_stats)!=len(x_vals):
        do_stats=[True if i==0 else False for i in range(len(x_vals))]

    r_vals=[]
    did_stats = False
    for i in range(len(x_vals)):
        x_val=x_vals[i]
        y_val=y_vals[i]
        color=colors[i] if colors else "red"
        alpha=alphas[i] if alphas else 1
        label=labels[i] if labels else ''

        ax.scatter(x_val, y_val, c=color, alpha=alpha, edgecolors='black', label=label, zorder=i)

        if do_stats and do_stats[i]:
            slope, intercept, r_value, p_value, std_err = stats.linregress(x_val,y_val)
            did_stats = True
        # label with text
        # title
            tx = 0.1
            ty = 0.9
            if r2_pos:
                tx = r2_pos[0]
                ty = r2_pos[1]
            if p_value < 0.001:
                p_value = "< 0.001"
            else:
                p_value = str(round(p_value, 3))
            ax.text(tx , ty, r"$R^2$: " + str(round(r_value**2,4))+"    P-value: {}".format(p_value) ,
                        horizontalalignment="left",
                        verticalalignment="center",
                        transform=ax.transAxes, size='small')#,bbox=dict(facecolor='white', alpha=0.5,boxstyle='round'))
            # draw the line of best fit
            ax.plot(x_val, intercept + slope*x_val, color='black', label='Fit Line')
            r_vals.append(r_value**2) # technically here to support multiple fits but the output currently only supports the last fit
    if title is not None:
        ax.set_title(title)
    if xlabel is not None:
        ax.set_xlabel(xlabel, labelpad=3)
    if ylabel is not None:
        ax.set_ylabel(ylabel, labelpad=0)
    if xticks is not None:
        ax.set_xticks(xticks)
    if yticks is not None:
        ax.set_yticks(yticks)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if legend:
        ax.legend()
    if mark_origin:
        ax.axhline(0,color='gray',linestyle='-.',zorder=0, alpha=0.5)
        ax.axvline(0,color='gray',linestyle='-.',zorder=0, alpha=0.5)
    sns.despine()
    if did_stats:
        return {'r':r_value,'pval':p_value,'slope':slope,'intercept':intercept,'r2':r_vals[-1]}

def kdehist(ax, values, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, colors=None, alphas=None, labels=None, vertical=False, hist=None,kde=None,rug=None, show_axis=False, cut = 0, bins = None, norm_hist = False):

    for i in range(len(values)):
        data = values[i]
        color = colors[i] if colors else 'r'
        label = labels[i] if labels else ''
        alpha = alphas[i] if alphas else 0.25
        hist_bool = hist[i] if hist else True
        kde_bool = kde[i] if kde else True
        rug_bool = rug[i] if rug else False
        if len(data)>1:
            sns.distplot(data, color=color, ax=ax, hist=hist_bool, kde=kde_bool, rug=rug_bool, vertical=vertical, bins = bins, hist_kws = {'ec':'black'}, kde_kws = {'cut':cut}, norm_hist = norm_hist  )

    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    if not show_axis:
        ax.axis('off')


def make_hist( ax, values, title = "histogram", x_label = "x axis", y_label = "y axis", bins = 40, weights = None, data_label = None, color=None, histtype = 'bar', stacked = False, legend=False, alpha=1):

    #f, ax = plt.subplots()
    ax.hist(values, bins=bins, density=False, alpha=alpha, weights = weights, label=data_label, color = color, histtype = histtype, stacked = stacked, ec = 'black')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True)
    if legend:
        ax.legend(prop={'size': 10})

    #f.savefig(pdfPages, format='pdf', bbox_inches = 'tight')
    #plt.close(f)
    return ax

def make_bar(ax, values, title = "bar graph", x_label = "x axis", y_label = "y axis", x_ticks = []):

    ax.bar(np.arange(len(values)), values)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True)
    ax.set_xticklabels(x_ticks)
    return ax


def var_prop(ax, var_counts, title = 'proportion of reads' ):

    wt, single = 0, 0
    #if '0' in var_counts.index:
    wt = var_counts.loc[0,'count']
    #if '1' in var_counts.index:
    single = var_counts.loc[1,'count']
    rest = var_counts.loc[2:max(var_counts.index),'count'].sum()
    total = wt + single + rest

    data = [ x/total*100 for x in [wt, single, rest] ]
    #print(data)
    make_bar(ax, data, title = title, x_label = "number of mutants", y_label = "percentage of reads", x_ticks=[0,1,'>1'] )

    return ax

def heatmap(ax, df, title = None, xlabel=None, ylabel=None, cmap=None, vmin=None, vmax=None):

    #cmap = plt.cm.get_cmap("RdYlGn")
    #cmap.set_bad("#9aa6ba")
    sns.heatmap(df, ax=ax, vmin=vmin, vmax=vmax, square=True, cbar_kws={"shrink":0.5},cmap=cmap, linewidths= 0.1)

    return ax


def ECcurve(x, top, bottom, hill, ec50):
    y = bottom + (top-bottom)/(1+((ec50/x)**hill))
    #y = bottom + (top-bottom)/(1+(10**((log_ec50-x)*hill) ))
    return y

def x_from_y_ECcurve(y, top, bottom, hill, ec50):
    x = ec50 * ( ((top-bottom)/(y-bottom) -1) )**(-1/hill)
    #x = log_ec50 - (1/hill)*np.log10( ((top-bottom)/(y-bottom) -1) )
    return x

def addTextVertical(ax, data, h_pos = 0, v_pos=0.5, spacing = 1, size = "medium"):
    '''
    Wrapper function for printing a list as text. Uses the normal ax.text() function.
    Only propagates in the vertical direction. Prints from bottom up, just like the
    indexing on heat maps; this means the last value in the list will be printed at the
    top. Can be made to print downward by changing spacing to a negative
    number.
    '''
    for i, x in enumerate(data):
        ax.text( h_pos , i*spacing + v_pos, x,
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transData,
               size = size)

def addTextHorizontal(ax, data, v_pos = 0, h_pos = 0.5, spacing = 1, size = "medium", transform = 'data' ):
    '''
    Wrapper function for printing a list as text. Uses the normal ax.text() function.
    Only propagates in the horizontal direction. Prints from left to right, just like the
    indexing on heat maps; this means the last value in the list will be printed at the
    right. Can be made to print leftward by changing spacing to a negative
    number.
    '''
    if transform == 'data':
        transform_ax = ax.transData
    elif transform == 'axes':
        transform_ax = ax.transAxes
    for i, x in enumerate(data):
        ax.text(h_pos+i*spacing , v_pos, x,
                horizontalalignment="center",
                verticalalignment="center",
                transform=transform_ax,
               size = size)

# set the colormap and centre the colorbar
class MidpointNormalize(mpl.colors.Normalize):
    """Normalise the colorbar."""
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
