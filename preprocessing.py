# Functions for preprocessing data sets from IsoCloud4.

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import os
import pandas as pd
from scipy import interpolate

from icemicro import svpice

def interpdf(dftime,df,timeindex):
    # return a dataframe interpolated to the correct time
    # dftime - time index of df; df  - dataframe to interpolate; timeindex - time index to be interpolated to
    time1 = dftime.values

    interp_dict = {}
    for col in df.columns:
        currcol = df[col].values
        f = interpolate.interp1d(time1,currcol,bounds_error=False)
        interp_dict[col]=f(timeindex)
    dfinterp = pd.DataFrame.from_dict(interp_dict)
    dfinterp['Time_s']=timeindex

    return dfinterp.set_index('Time_s')
def loadaidadata(exp,verbose=False):
    # AutoMergeISO4/loadisoexp
    # instrument data directories
    #isodir="/Users/karalamb/ISOCLOUD/ISOCLOUD04/"
    isodir="ISOCLOUD04/"
    aerosoldir="Aerosol"
    apetdir="APeT"
    mbwdir="MBW"
    ptdir="pT"
    tdldir="TDL"
    tgvdir="tgv"
    welasdir="welas"
    welas2dir="welas2"
    chiwishdodir="ChiWisHDO"
    chiwish2odir="NewIsoCloudFits"
    bkscdir="Backscatter Corrections"

    # File name template strings
    aerosoltemp="ISO04_{:02.0f}_AIDA-Aerosol.csv"
    apettemp="ISO04_{:02.0f}_AIDA-APeT.csv"
    mbwtemp="ISO04_{:02.0f}_AIDA-MBW.csv"
    pttemp="ISO04_{:02.0f}_AIDA-pT.csv"
    tdltemp="ISO04_{:02.0f}_AIDA-TDL.csv"
    tgvtemp="ISO04_{:02.0f}_AIDA-tgv.csv"
    welastemp="ISO04_{:02.0f}_AIDA-welas.csv"
    welas2temp="ISO04_{:02.0f}_AIDA-welas2.csv"

    # back-scatter correction for ChiWis
    expdict = {"1": "20130311", "2": "20130311", "3": "20130311", "4": "20130311", "5": "20130311",
           "6": "20130313", "7": "20130313", "8": "20130313", "9": "20130313", "10": "20130313", "11": "20130313",
           "12": "20130314", "13": "20130314", "14": "20130314", "15": "20130314", "16": "20130314",
           "17": "20130314", "18": "20130315", "19": "20130315", "20": "20130315", "21": "20130315",
           "22": "20130315", "23": "20130315", "24": "20130318", "25": "20130318", "26": "20130318",
           "27": "20130318", "28": "20130319", "29": "20130319", "30": "20130319", "31": "20130319",
           "32": "20130319", "33": "20130319", "34": "20130320", "35": "20130320", "36": "20130320",
           "37": "20130320", "38": "20130320", "39": "20130321", "40": "20130321", "41": "20130321",
           "42": "20130321", "43": "20130321", "44": "20130322", "45": "20130322", "46": "20130322",
           "47": "20130322", "48": "20130322"}

    chiwishdotemp = "ISO04_{}_CHI-WIS.nas" #datestring,expnum
    chiwish2otemp = "{}_ISO04_EXP{:02.0f}.dat" #datestring,expnum
    bksctemp = "{} Exp {:02.0f} bksc ratio correction.nas" #datestring, expnum

    datestr = expdict[str(exp)]

    #IsoCloud Experiment list
    isoexplist = os.path.join(isodir,"exps.txt")

    # AIDA aerosol (CPC data) - cn_aida_dil dilution ratio during pumping
    aerosolfile = os.path.join(isodir,aerosoldir,aerosoltemp.format(exp))
    # AIDA housekeeping - gas temperatures - Kelvin
    tgvfile = os.path.join(isodir,tgvdir,tgvtemp.format(exp))
    # AIDA-TDL (either SP-APicT or APicT)
    tdlfile = os.path.join(isodir,tdldir,tdltemp.format(exp))
    # MBW
    mbwfile = os.path.join(isodir,mbwdir,mbwtemp.format(exp))
    # APeT data
    apetfile = os.path.join(isodir,apetdir,apettemp.format(exp))
    # AIDA housekeeping - time, p, tgas, twall
    ptfile = os.path.join(isodir,ptdir,pttemp.format(exp))
    # welas data - time, p, tgas, twall
    welasfile = os.path.join(isodir,welasdir,welastemp.format(exp))
    # welas2 data - time, p, tgas, twall
    welas2file = os.path.join(isodir,welas2dir,welas2temp.format(exp))
    # Chi-Wis data - Water Vapor Volumetric Mixing Ratio in ppmv: p(HDO)/p(ambient) * 10^6
    chiwishdofile = os.path.join(isodir,chiwishdodir,chiwishdotemp.format(datestr))
    # Chi-Wis data - Water Vapor Volumetric Mixing Ratio in ppmv: p(H2O)/p(ambient) * 10^6
    chiwish2ofile = os.path.join(isodir,chiwish2odir,chiwish2otemp.format(datestr,exp))
    # Back-scattering correction - relative extinction between Chi-Wis and SP-APicT
    bkscfile = os.path.join(isodir,bkscdir,bksctemp.format(datestr,exp))

    isoexpdf = pd.read_csv(isoexplist,header='infer',delim_whitespace=True)
    aerosoldf = pd.read_csv(aerosolfile,header='infer',skipinitialspace=True,skiprows=4)
    tgvdf = pd.read_csv(tgvfile,header='infer',skipinitialspace=True,skiprows=4)
    tdldf = pd.read_csv(tdlfile,header='infer',skipinitialspace=True,skiprows=4)
    mbwdf = pd.read_csv(mbwfile,header='infer',skipinitialspace=True,skiprows=4)
    apetdf = pd.read_csv(apetfile,header='infer',skipinitialspace=True,names=["time_apet","h2o_ppmv_apet"])
    if exp==5:
        ptdf = pd.read_csv(ptfile,header='infer',skipinitialspace=True,skiprows=0)
    else:
        ptdf = pd.read_csv(ptfile,header='infer',skipinitialspace=True,skiprows=4)
    welasdf = pd.read_csv(welasfile,header='infer',skipinitialspace=True,skiprows=4)
    welas2df = pd.read_csv(welas2file,header='infer',skipinitialspace=True,skiprows=4)
    chiwishdodf = pd.read_csv(chiwishdofile,header='infer')
    chiwish2odf = pd.read_csv(chiwish2ofile,delim_whitespace=True,skiprows=1,header='infer',
                              names=["Time (hrs)","P", "Temp","h2o_ppmv_chiwis"])
    if os.path.exists(bkscfile):
        bkscdf = pd.read_csv(bkscfile,delim_whitespace=True,header='infer',skiprows=2,
                             names=["Time_hrs", "backscatter_offset", "H2O_factor", "HDO_factor"])
    else:
        print("No backscatter correction")



    expdf = isoexpdf[isoexpdf['exps']==exp]
    expreftime = expdf['expreftime'].values[0]
    expstarttime = expdf['expstarttime'].values[0]
    expendtime = expdf['expendtime'].values[0]

    if verbose:
        print(datestr)
        print(expreftime,expstarttime,expendtime)

        print(aerosolfile)
        print(tgvfile)
        print(tdlfile)
        print(mbwfile)
        print(apetfile)
        print(ptfile)
        print(welasfile)
        print(welas2file)
        print(chiwishdofile)
        print(chiwish2ofile)
        print(bkscfile)

        # AutoMergeISO4/expansiontext
    aidadf = ptdf.set_index("Time_s")
    aidadf=aidadf.join(tgvdf.set_index("Time_s"),on="Time_s")
    aidadf=aidadf.join(aerosoldf.set_index("Time_s"),on="Time_s")
    aidadf=aidadf.join(tdldf.set_index("Time_s"),on="Time_s")
    aidadf=aidadf.join(mbwdf.set_index("Time_s"),on="Time_s")

    # Interpolate/upscale Welas to same time base
    welasdf = welasdf.set_index("Time_s")
    welasdf=interpdf(welasdf.index,welasdf,aidadf.index)
    welas2df = welas2df.set_index("Time_s")
    welas2df=interpdf(welas2df.index,welas2df,aidadf.index)

    aidadf=aidadf.join(welasdf,on="Time_s")
    aidadf=aidadf.join(welas2df,on="Time_s")

    # ApeT
    apetdf["time_apet"]=apetdf["time_apet"]-17.0
    apetdf = apetdf.set_index("time_apet")
    apetdf=interpdf(apetdf.index,apetdf,aidadf.index)
    aidadf=aidadf.join(apetdf,on="Time_s")

    # Chi-Wis + backscatter correction, interpolate to experiment time base
    chiwish2odf['Time_s']=chiwish2odf['Time (hrs)']*3600.0-expreftime
    chiwish2odf=chiwish2odf.drop(columns=['Time (hrs)'])
    chiwish2odf = chiwish2odf.set_index("Time_s")
    chiwish2odf=interpdf(chiwish2odf.index,chiwish2odf,aidadf.index)
    aidadf=aidadf.join(chiwish2odf,on="Time_s")

    if os.path.exists(bkscfile):
        bkscdf['Time_s']=bkscdf['Time_hrs']*3600-expreftime
        bkscdf = bkscdf.drop(columns=['Time_hrs'])
        bkscdf = bkscdf.set_index("Time_s")
        bkscdf=interpdf(bkscdf.index,bkscdf,aidadf.index)
        aidadf=aidadf.join(bkscdf,on="Time_s")
    else:
        bkscdata = np.ones((len(aidadf),3))
        bkscdf = pd.DataFrame(bkscdata,columns=['backscatter_offset','H2O_factor','HDO_factor'],index=aidadf.index)
        aidadf=aidadf.join(bkscdf)

    # Replace -9999 with nans
    for col in aidadf.columns:
        aidadf[col][aidadf[col]<-100]=np.nan
    #df.loc[df['First season'] > 1990, 'First Season'] = 1

    return aidadf
def plotisoexp(aidadf,mintime,maxtime):
    aidadfsel=aidadf.loc[(aidadf.index>mintime)]
    aidadfsel=aidadfsel.loc[(aidadfsel.index<maxtime)]
    aidadf = aidadfsel

    fig,axs = plt.subplots(8,1,sharex='col',figsize=(12,12))

    axs[0].plot(aidadf.index,aidadf["p_hPa"])
    axs[0].set_ylabel("Pressure (hPa)")
    axs[1].plot(aidadf.index,aidadf["Tg_K"])
    axs[1].plot(aidadf.index,aidadf["Tw_K"])
    axs[1].legend(["Tg_K","Tw_K"])
    axs[1].set_ylabel("Temperature (K)")
    axs[2].plot(aidadf.index,aidadf["Tgv1"])
    axs[2].plot(aidadf.index,aidadf["Tgv2"])
    axs[2].plot(aidadf.index,aidadf["Tgv3"])
    axs[2].plot(aidadf.index,aidadf["Tgv4"])
    axs[2].legend(["Tgv1","Tgv2","Tgv3","Tgv4"])
    axs[2].set_ylabel("Temperature (K)")

        # convert H2O measurements to ppmv
    axs[3].plot(aidadf.index,aidadf["pw_TDL_pa"]/aidadf["p_hPa"]*1e4)
    axs[3].plot(aidadf.index,aidadf["h2o_ppmv_chiwis"]*aidadf["H2O_factor"])
    axs[3].plot(aidadf.index,aidadf["pw_MBW_pa"]/aidadf["p_hPa"]*1e4)
    axs[3].plot(aidadf.index,aidadf["h2o_ppmv_apet"])
    axs[3].legend(["h2o_ppmv_TDL","h2o_ppmv_chiwis","h2o_ppmv_MBW","h2o_ppmv_apet"])
    axs[3].set_ylabel("H2O ppmv")

    # convert H2O measurements to RHi
    rhi_apet = (aidadf["h2o_ppmv_apet"])*aidadf["p_hPa"]/1e4/svpice(aidadf["Tg_K"])*100
    rhi_chiwis = (aidadf["h2o_ppmv_chiwis"]*aidadf["H2O_factor"])*aidadf["p_hPa"]/1e4/svpice(aidadf["Tg_K"])*100

    axs[4].plot(aidadf.index,aidadf["rhi_TDL"])
    axs[4].plot(aidadf.index,rhi_chiwis)
    axs[4].plot(aidadf.index,aidadf["rhi_MBW"])
    axs[4].plot(aidadf.index,rhi_apet)
    axs[4].legend(["rhi_TDL","rhi_chiwis","rhi_MBW","rhi_apet"])
    axs[4].set_ylabel("RHi")

    # Ice water content - ppmv
    axs[5].plot(aidadf.index,aidadf["h2o_ppmv_apet"]-aidadf["pw_TDL_pa"]/aidadf["p_hPa"]*1e4)
    axs[5].plot(aidadf.index,aidadf["h2o_ppmv_apet"]-aidadf["h2o_ppmv_chiwis"])
    axs[5].legend(["ice_ppmv_TDL","ice_ppmv_chiwis"])
    axs[5].set_ylabel("IWC (ppmv)")

    axs[6].plot(aidadf.index,aidadf["cn_welas"])
    axs[6].plot(aidadf.index,aidadf["cn_welas2"])
    axs[6].plot(aidadf.index,aidadf["cn_welas_ice"])
    axs[6].plot(aidadf.index,aidadf["cn_welas2_ice"])
    axs[6].legend(["cn_welas","cn_welas2","cn_welas_ice","cn_welas2_ice"])
    axs[6].set_ylabel("IN")

    axs[7].plot(aidadf.index,aidadf["fn_welas_ice"])
    axs[7].plot(aidadf.index,aidadf["fn_welas2_ice"])
    axs[7].legend(["fn_welas_ice","fn_welas2_ice"])
    axs[7].set_ylabel("Fraction Ice")

    plt.show()
class SSA(object):
    # Code from https://www.kaggle.com/jdarcy/introducing-ssa-for-time-series-decomposition
    __supported_types = (pd.Series, np.ndarray, list)

    def __init__(self, tseries, L, save_mem=True):
        """
        Decomposes the given time series with a singular-spectrum analysis. Assumes the values of the time series are
        recorded at equal intervals.

        Parameters
        ----------
        tseries : The original time series, in the form of a Pandas Series, NumPy array or list.
        L : The window length. Must be an integer 2 <= L <= N/2, where N is the length of the time series.
        save_mem : Conserve memory by not retaining the elementary matrices. Recommended for long time series with
            thousands of values. Defaults to True.

        Note: Even if an NumPy array or list is used for the initial time series, all time series returned will be
        in the form of a Pandas Series or DataFrame object.
        """

        # Tedious type-checking for the initial time series
        if not isinstance(tseries, self.__supported_types):
            raise TypeError("Unsupported time series object. Try Pandas Series, NumPy array or list.")

        # Checks to save us from ourselves
        self.N = len(tseries)
        if not 2 <= L <= self.N/2:
            raise ValueError("The window length must be in the interval [2, N/2].")

        self.L = L
        self.orig_TS = pd.Series(tseries)
        self.K = self.N - self.L + 1

        # Embed the time series in a trajectory matrix
        self.X = np.array([self.orig_TS.values[i:L+i] for i in range(0, self.K)]).T

        # Decompose the trajectory matrix
        self.U, self.Sigma, VT = np.linalg.svd(self.X)
        self.d = np.linalg.matrix_rank(self.X)

        self.TS_comps = np.zeros((self.N, self.d))

        if not save_mem:
            # Construct and save all the elementary matrices
            self.X_elem = np.array([ self.Sigma[i]*np.outer(self.U[:,i], VT[i,:]) for i in range(self.d) ])

            # Diagonally average the elementary matrices, store them as columns in array.
            for i in range(self.d):
                X_rev = self.X_elem[i, ::-1]
                self.TS_comps[:,i] = [X_rev.diagonal(j).mean() for j in range(-X_rev.shape[0]+1, X_rev.shape[1])]

            self.V = VT.T
        else:
            # Reconstruct the elementary matrices without storing them
            for i in range(self.d):
                X_elem = self.Sigma[i]*np.outer(self.U[:,i], VT[i,:])
                X_rev = X_elem[::-1]
                self.TS_comps[:,i] = [X_rev.diagonal(j).mean() for j in range(-X_rev.shape[0]+1, X_rev.shape[1])]

            self.X_elem = "Re-run with save_mem=False to retain the elementary matrices."

            # The V array may also be very large under these circumstances, so we won't keep it.
            self.V = "Re-run with save_mem=False to retain the V matrix."

        # Calculate the w-correlation matrix.
        self.calc_wcorr()

    def components_to_df(self, n=0):
        """
        Returns all the time series components in a single Pandas DataFrame object.
        """
        if n > 0:
            n = min(n, self.d)
        else:
            n = self.d

        # Create list of columns - call them F0, F1, F2, ...
        cols = ["F{}".format(i) for i in range(n)]
        return pd.DataFrame(self.TS_comps[:, :n], columns=cols, index=self.orig_TS.index)


    def reconstruct(self, indices):
        """
        Reconstructs the time series from its elementary components, using the given indices. Returns a Pandas Series
        object with the reconstructed time series.

        Parameters
        ----------
        indices: An integer, list of integers or slice(n,m) object, representing the elementary components to sum.
        """
        if isinstance(indices, int): indices = [indices]

        ts_vals = self.TS_comps[:,indices].sum(axis=1)
        return pd.Series(ts_vals, index=self.orig_TS.index)

    def calc_wcorr(self):
        """
        Calculates the w-correlation matrix for the time series.
        """

        # Calculate the weights
        w = np.array(list(np.arange(self.L)+1) + [self.L]*(self.K-self.L-1) + list(np.arange(self.L)+1)[::-1])

        def w_inner(F_i, F_j):
            return w.dot(F_i*F_j)

        # Calculated weighted norms, ||F_i||_w, then invert.
        F_wnorms = np.array([w_inner(self.TS_comps[:,i], self.TS_comps[:,i]) for i in range(self.d)])
        F_wnorms = F_wnorms**-0.5

        # Calculate Wcorr.
        self.Wcorr = np.identity(self.d)
        for i in range(self.d):
            for j in range(i+1,self.d):
                self.Wcorr[i,j] = abs(w_inner(self.TS_comps[:,i], self.TS_comps[:,j]) * F_wnorms[i] * F_wnorms[j])
                self.Wcorr[j,i] = self.Wcorr[i,j]

    def plot_wcorr(self, min=None, max=None):
        """
        Plots the w-correlation matrix for the decomposed time series.
        """
        if min is None:
            min = 0
        if max is None:
            max = self.d

        if self.Wcorr is None:
            self.calc_wcorr()

        ax = plt.imshow(self.Wcorr)
        plt.xlabel(r"$\tilde{F}_i$")
        plt.ylabel(r"$\tilde{F}_j$")
        plt.colorbar(ax.colorbar, fraction=0.045)
        ax.colorbar.set_label("$W_{i,j}$")
        plt.clim(0,1)

        # For plotting purposes:
        if max == self.d:
            max_rnge = self.d-1
        else:
            max_rnge = max

        plt.xlim(min-0.5, max_rnge+0.5)
        plt.ylim(max_rnge+0.5, min-0.5)
