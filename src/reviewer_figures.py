#!/usr/bin/env python3
"""Independent grayscale scientific plots of source-defined Fig6B proportions and sharp strata.
Inputs must be existing derived *software-output* tables. No new biological inference.
"""
import csv, argparse
from pathlib import Path
from fractions import Fraction

def load_groups(p):
    with Path(p).open(newline='',encoding='utf-8') as f: data=list(csv.DictReader(f))
    d={x['group']:x for x in data}
    assert set(d)=={'Unperturbed','CAM_after_burst','CAM_throughout','rpoS_KO'}
    for x in data:
        assert int(x['Survive'])+int(x['No_Growth'])+int(x['Unknown'])==int(x['N_all'])
    return d

def load_strata(p):
    with Path(p).open(newline='',encoding='utf-8') as f: data=list(csv.DictReader(f))
    assert len(data)==8
    for x in data:
        lo,hi=Fraction(x['diff_low']),Fraction(x['diff_high'])
        assert lo<=hi and (lo>0)==(x['order_identified']=='True')
    return data

def make_plots(groups_csv,strata_csv,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    gg=load_groups(groups_csv)
    names=['Unperturbed','CAM_after_burst','CAM_throughout','rpoS_KO']
    lbls=['Unperturbed','CAM after 5 h','CAM throughout',r'$\Delta\mathit{rpoS}$']
    mid=[]; err=[]
    for n in names:
        d=gg[n];s=int(d['Survive']);n0=int(d['No_Growth']);nall=int(d['N_all'])
        p=s/(s+n0);mid.append(p);err.append((p*(1-p)/nall)**.5)
    fig,ax=plt.subplots(figsize=(7.4,4.5),layout='constrained')
    xx=list(range(4))
    ax.bar(xx,mid,width=.65,facecolor='.72',edgecolor='.25',linewidth=1.1)
    ax.errorbar(xx,mid,yerr=err,fmt='none',ecolor='.16',elinewidth=1.2,capsize=4)
    ax.set_xticks(xx,lbls)
    ax.set_ylim(0,1.075);ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_ylabel('Known-outcome regrowth fraction')
    # Paper captions supply the figure title; leave more vertical margin for data.

    ax.grid(axis='y',alpha=.15)
    for i,x in enumerate(mid):ax.text(i,x+.055,f'{x*100:.1f}%',ha='center',va='bottom',fontsize=9)
    f=out/'figure1_recovery_fraction.png';fig.savefig(f,dpi=200,bbox_inches='tight',pad_inches=0.2);plt.close(fig)
    st=sorted(load_strata(strata_csv),key=lambda r: (0 if r['group_A']=='Unperturbed' else 1,r['date'],r['STRESS']))
    # Two design families and dose-specific strata, deliberately not pooled as four-arm randomized experiment.
    fig,ax=plt.subplots(figsize=(9.0,5.25),layout='constrained')
    ys=list(range(8))[::-1]
    for i,r in enumerate(st):
        low=float(Fraction(r['diff_low']))*100; high=float(Fraction(r['diff_high']))*100
        y=ys[i];ax.hlines(y,low,high,linewidth=3,color='.35')
        ax.plot([low,high],[y,y],'|',color='.16',markersize=11,markeredgewidth=2)
        ax.plot((low+high)/2,y,'o',markersize=4,color='.12')
    ax.axvline(0,color='.32',linewidth=1.1,linestyle='--')
    ax.set_yticks(ys,[f"{r['date'][:4]}-{r['date'][4:6]}-{r['date'][6:]} | {r['STRESS'].replace('H2O2_','')}" for r in st])
    ax.set_xlim(-55,105);ax.set_xlabel('Sharp difference under arbitrary binary Unknown completion (percentage points)')
    # Caption is included in the manuscript, not within the plot image.

    ax.text(.02,.99,'First four: unperturbed − $\\Delta rpoS$\nLast four: CAM after 5 h − CAM throughout',transform=ax.transAxes,va='top',fontsize=9,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
    ax.grid(axis='x',alpha=.16)
    for side in ['top','right']:ax.spines[side].set_visible(False)
    f2=out/'figure5_matched_date_dose_bounds.png';fig.savefig(f2,dpi=200,bbox_inches='tight',pad_inches=0.2);plt.close(fig)
    return [str(f),str(f2)]

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--groups',required=True);p.add_argument('--strata',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();print('\n'.join(make_plots(a.groups,a.strata,a.out)))
if __name__=='__main__':main()
