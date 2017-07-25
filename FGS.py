import pandas as pd
import glob
import urllib
import numpy as np
import os
from progressbar import Percentage,ProgressBar,Bar,ETA,FileTransferSpeed
import progressbar
import tarfile
import datetime
from astropy.time import Time

def download(url,fileName,progress=True):
    '''
    Download a file from http either with or without a progress bar
    url = url to download
    fileName = filename to save to
    progress = display progress bar (True/False)
    '''
    
    if progress==True:
        widgets = ['Test: ', Percentage(), ' ', Bar(), ' ', ETA(), ' ', FileTransferSpeed()]
        pbar = ProgressBar(widgets=widgets)
        def dlProgress(count, blockSize, totalSize):
            if pbar.maxval is None:
                pbar.maxval = totalSize
                pbar.start()

            pbar.update(min(count*blockSize, totalSize))

        urllib.urlretrieve(url, fileName, reporthook=dlProgress)
        pbar.finish()
    else:
        urllib.urlretrieve(url, fileName)
        
    return 
 
def get_data(datadir='',mission='kepler',quarters=None):

    '''
    Obtain the FGS data for Kepler or K2 directly from MAST
    mission = kepler,k2
    '''
    if type(quarters)==int:
        quarters=[quarters]
    #Find the data for the right mission
    if mission=='kepler':
        prefix='kplr'
        pref='q'
    if mission=='k2':
        prefix='ktwo'
        pref='c'
    if (mission!='kepler')&(mission!='k2'):
        print "Choose 'kepler' or 'k2' for mission."
        return
    
    #Get ancilliary table
    if glob.glob(datadir+prefix+'-anc-2017013163000.txt')==[]:
        url='https://archive.stsci.edu/missions/'+mission+'/fgs/flc/'+prefix+'-anc-2017013163000.txt'
        fileName=datadir+prefix+'-anc-2017013163000.txt'
        download(url,fileName,progress=False)

    #Find all the unique stars
    stars=pd.read_csv(datadir+prefix+'-anc-2017013163000.txt',skiprows=[0,1,3])
    if quarters==None:
        quarters=np.unique(stars[stars.keys()[0]])
    for q in quarters:
        fname=''+prefix+'-flc-'+pref+str(q).zfill(2)+'-xx.txt'
        if glob.glob(datadir+fname)==[]:
            print 'Downloading Quarter '+str(q)
            url='https://archive.stsci.edu/missions/'+mission+'/fgs/flc/'+fname+'.tar.gz'
            fileName=datadir+fname+'.tar.gz'
            download(url,fileName)
            tar = tarfile.open(datadir+fname+'.tar.gz', "r:gz")
            tar.extractall(datadir+fname)
            tar.close()
            os.remove(datadir+fname+'.tar.gz')


def table(datadir='',return_loc=False):
    if glob.glob(datadir+'kplr-anc-2017013163000.txt')==[]:
        url='https://archive.stsci.edu/missions/kepler/fgs/flc/kplr-anc-2017013163000.txt'
        fileName=datadir+'kplr-anc-2017013163000.txt'
        try:
            download(url,fileName,progress=False)
        except:
            print 'Cannot download data'
        if glob.glob(datadir+'kplr-anc-2017013163000.txt')==[]:
            print 'No Kepler ancilliary data.'
            return None

    if glob.glob(datadir+'ktwo-anc-2017013163000.txt')==[]:
        url='https://archive.stsci.edu/missions/k2/fgs/flc/ktwo-anc-2017013163000.txt'
        fileName=datadir+'ktwo-anc-2017013163000.txt'
        try:
            download(url,fileName,progress=False)
        except:
            print 'Cannot download data'
        if glob.glob(datadir+'ktwo-anc-2017013163000.txt')==[]:
            print 'No K2 ancilliary data.'
            return None
        
    keys=['KEPLER_ID','RA','DEC','KEPMAG','mission','QUARTER','FGS_MODULE','STAR_INDEX']
    k2=pd.read_csv('ktwo-anc-2017013163000.txt',skiprows=[0,1,3])
#    k2=k2.loc[np.unique(k2['RA'],return_index=True)[1]].reset_index(drop=True)
    k2.columns=np.append('QUARTER',k2.keys()[1:])
    k2['mission']='k2'
    k2=k2[keys]

    kepler=pd.read_csv('kplr-anc-2017013163000.txt',skiprows=[0,1,3])
#    kepler=kepler.loc[np.unique(k2['RA'],return_index=True)[1]].reset_index(drop=True)
    kepler['mission']='kepler'
    kepler.columns=np.append('QUARTER',kepler.keys()[1:])
    kepler=kepler[keys]
    tab=k2.append(kepler).reset_index(drop=True)
    tab=tab[tab.KEPLER_ID!=-1].reset_index(drop=True)
    if return_loc==False:
        tab=tab.loc[np.unique(tab.RA,return_index=True)[1]].reset_index(drop=True)
        tab=tab.sort_values('mission')
        tab=tab[keys[0:-3]]

    return tab.reset_index(drop=True)


def fgs_lc(datadir,quarter,module,starno,div_trend=False,norm=True,raw=False,npoly=1,nsig=5.,mission='kepler'):

    '''
    Reads, cleans and returns an FGS light curve based on quarter, module, star number

    npoly = Order of polynomial to remove from flux
    nsig = sigma level to clip data
    raw = returns raw data with no cleaning
    norm = normalise data
    '''
    if type(quarter)==int:
        quarter=str(quarter).zfill(2)

    if type(module)==int:
        module=str(module).zfill(2)


    if mission=='kepler':
        prefix='kplr'
        pref='q'
    if mission=='k2':
        prefix='ktwo'
        pref='c'
        
    mods=np.asarray(['01','05','21','25'],dtype=str)
    fnames=glob.glob(datadir+prefix+'-flc-'+pref+quarter+'-xx.txt/*')
    if fnames==[]:
        try:
            get_data(datadir,mission=mission,quarters=[np.int(quarter)])
        except:
            print 'Cannot download data'
            return [],[],[],[]
        fnames=glob.glob(datadir+prefix+'-flc-'+pref+quarter+'-xx.txt/*')
        if fnames==[]:
            print 'No data'
            return None,None,None,None

    f=fnames[np.where(mods==module)[0][0]]
    df=pd.read_csv(f,skiprows=[0,2])
    if (int(starno)*3)+3>len(df.keys()):
        print 'No such star'
        return None,None,None,None


    #Clear out nans from the data
    #----------------------------
    keys=df.keys()[[0,(int(starno)*3)+1,(int(starno)*3)+2,(int(starno)*3)+3]]
    df=df[keys].dropna().reset_index(drop=True)

    
    dates=[datetime.datetime.strptime(d,'%m/%d/%y %H:%M:%S.%f') for d in df[df.keys()[0]]]
    t=Time(dates).jd
    pref=df.keys()[1][0:5]

    counts,col,row=np.asarray(df[keys[1]]),np.asarray(df[keys[2]]),np.asarray(df[keys[3]])
    if raw==True:
        return t,counts,col,row

    #Get rid of bad pointings
    #------------------------
    #If pointing drifts more than a few pixels it is likely wrong
    good=np.where((abs(col-np.median(col))<5.)&(abs(row-np.median(row))<5.))[0] 
    #Has to have at some values
    good=good[np.where(counts[good]>=250.)[0]]
    
    
    #Sigma to clip each quarter

    
    l=np.polyfit(t[good],col[good],npoly)
    line=np.polyval(l,t[good])
    good=good[np.where(np.abs(col[good]-line)<=nsig*np.std(col[good]-line))[0]]
    
    l=np.polyfit(t[good],row[good],npoly)
    line=np.polyval(l,t[good])
    good=good[np.where(np.abs(row[good]-line)<=nsig*np.std(row[good]-line))[0]]
     
    l=np.polyfit(t[good],counts[good],npoly)
    line=np.polyval(l,t[good])
    good=good[np.where(np.abs(counts[good]-line)<=nsig*np.std(counts[good]-line))[0]]
     
    if norm==True:
        counts/=np.nanmedian(counts[good])
    if div_trend==True:
        l=np.polyfit(t[good],counts[good],npoly)
        line=np.polyval(l,t)
        counts/=line
        

    return t[good],counts[good],col[good],row[good]



def gen_lc(datadir='',ID=None,norm=True,div_trend=False,quarters=None,raw=False,npoly=1,nsig=5):

    '''
    Generate light curve for FGS data. Reads data from multiple modules and stitches it together.

    '''
    tab=table(return_loc=True)
    tab=tab.sort_values('QUARTER')   
    loc=np.where(tab.KEPLER_ID==ID)[0]
    if len(loc)==0:
        print 'No such star'
        return
    if len(loc)>1:
        loc=[loc[0]]
    mission=np.asarray(tab.loc[loc,'mission'])[0]

    if mission=='kepler':
        prefix='kplr'
    if mission=='k2':
        prefix='ktwo'
    if (mission!='kepler')&(mission!='k2'):
        print "Choose 'kepler' or 'k2' for mission."
        return

    if type(quarters)==int:
        quarters=[quarters]
    #Get ancilliary table
    if glob.glob(datadir+prefix+'-anc-2017013163000.txt')==[]:
        url='https://archive.stsci.edu/missions/'+mission+'/fgs/flc/'+prefix+'-anc-2017013163000.txt'
        fileName=datadir+prefix+'-anc-2017013163000.txt'
        try:
            download(url,fileName,progress=False)
        except:
            print 'Cannot download data'
            return

    if ID==None:
        print 'Please provide an ID'
        return [],[],[],[]

    x,y,cols,rows=np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0)

    s_tab=tab[tab.KEPLER_ID==ID].reset_index(drop=True)
    s_qs=np.asarray(s_tab['QUARTER'],dtype=int)
    if quarters==None:
        quarters=s_qs
    for q in quarters:
        quarter=str(q).zfill(2)
        pos=np.where(s_qs==q)[0]
        if len(pos)==0:
            continue
        pos=pos[0]
        module=str(s_tab.loc[pos,'FGS_MODULE']).zfill(2)
        starno=s_tab.loc[pos,'STAR_INDEX']-1
        t,counts,col,row=fgs_lc(datadir,quarter,module,starno,raw=raw,norm=norm,mission=mission,npoly=npoly,nsig=nsig)
        x=np.append(x,t,axis=0)
        y=np.append(y,counts,axis=0)
        cols=np.append(cols,col,axis=0)
        rows=np.append(rows,row,axis=0)
    return x,y,cols,rows
