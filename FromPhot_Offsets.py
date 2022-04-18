def sed_from_comps(wsed, flsedA, flsedB, ratio, slope):
    """Partendo dal flusso delle due componenti (A, disco di accrescimento, e B, toro), calcola la loro SED
    e la SED totale."""
    lfsed = np.log10(ratio*flsedB/0.59 + (1.-ratio)*flsedA/0.41 * (wsed/2200.)**slope )
    return lfsed
#
def phot_from_sed_AB(wo, stdev, blo, bup, wsed, flsedA, flsedB, zf1, ratio, norm, slope, IGM_tau):
 lwsed = np.log10(wsed)
 lfsed = sed_from_comps(wsed, flsedA, flsedB, ratio, slope)
 lfcorr = lfsed + norm
 wzsed = wsed*zf1
 lfm = [0.] * len(wo)
 for kk in range (0, len(wsed)) :
     if wsed[kk] < 1300. :         
        iw = int(round(wsed[kk]))
        z10 = int(round(10.*(zf1-1.)))
        key = (iw, z10)
        lfcorr[kk] -= IGM_tau[key]
        if IGM_tau[key]<0:
            print("TAU NOT OK")
     else :
        break
 for ll in range (0, len(wo)):
      wps = 0.
      for j in range (0, len(wzsed)):
        if blo[ll] < wzsed[j] and bup[ll] > wzsed[j] :
            prob = 2.- abs(wzsed[j]-wo[ll])/stdev[ll]
            if prob > 0. :
                wps += prob
                lfm[ll] += lfcorr[j]*prob
        elif wzsed[j] > bup[ll] :
            break
      if wps>0:
        lfm[ll] = lfm[ll]/wps
 return lfm
#
def residual(params, wsed, flsedA, flsedB, w, stdev, blo, bup, f, err, IGM_tau):
 zf1 = params['zf1'].value
 ratio = params['ratio'].value
 norm = params['norm'].value
 slope = params['slope'].value
 npar = 4
 emagmin = 0.1
 lw = np.log10(w)
 le = 0.4*err
 lf = [0.] * len(f)
 lfm = [0.] * len(lf)
 lfm = phot_from_sed_AB(w, stdev, blo, bup, wsed, flsedA, flsedB, zf1, ratio, norm, slope, IGM_tau)
 restot = 0.
 nphot = 0
 for ll in range (0, len(w)):
         sigma = emagmin if emagmin > le[ll] else le[ll]
         if f[ll]>0. :
             lf[ll] = -0.4*f[ll] -lw[ll] + 10.8
             res = (lf[ll]-lfm[ll])/sigma
             nphot += 1
         elif f[ll]<0. :
             lf[ll] = 0.4*f[ll] -lw[ll] + 10.8
             if lfm[ll] > lf[ll] :
                 res = (lf[ll]-lfm[ll])/sigma
                 nphot += 1
             else:
                 res = 0.
         restot += abs(res) 
 redux = nphot-1-npar if (nphot-1-npar) > 0 else 1
 restot = restot/redux/0.79788456
 return restot
#
def read_sed_file(sedfilename):
    wsed =  []
    flsed = []
    sed = open(sedfilename, 'r')
    for line in sed:
         line = line.strip()
         columns = line.split()
         wsed.append(float(columns[0]))
         flsed.append(float(columns[1]))
    sed.close()
    wsed = np.array(wsed)
    flsed = np.array(flsed)
    return wsed, flsed

def init_lock(l):
    global lock
    lock=l
    
def photometric_offset(phototab,offset_band):
    #lock.acquire()
    wsed, flsedA = read_sed_file('QSOsedA.dat')
    wsed, flsedB = read_sed_file('QSOsedB.dat')
    # CREATE IGM_TRANS E IGM_TAU DICTIONARIES
    a_file = open("taueff_NzS70.txt")
    IGM_trans = {}
    IGM_tau = {}
    for line in a_file:
        line = line.strip()
        value = line.split()
        iw = int(float(value[0]))
        for i in range(1, 71):
            key = (iw, i)
            IGM_trans[key] = float(value[i])
            IGM_tau[key] = -np.log10(IGM_trans[key])
    a_file.close()
    #lock.release()
    
    nrow = len(phototab)
    id = phototab['qid']
    z_spec = phototab['z_spec']
    nbands=41
    w = [None] * nbands
    fl = [None] * nbands
    err = [None] * nbands
    w = np.array([
    # GAIA
    5272.46, 6291.13, 7761.11,
    # 2MASS
    12390.58, 16487.19, 21634.04,
    # WISE
    33682.21, 46179.05, 120717.43, 221944.44,
    # skymapper
    3503.42, 3828.56, 5033.64, 6174.52, 7806.40, 8873.64,
    3503.42, 3828.56, 5033.64, 6174.52, 7806.40, 8873.64,
    # SDSS http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=SLOAN
    3565.05, 4700.33, 6174.48, 7533.63, 8781.69,
    3565.05, 4700.33, 6174.48, 7533.63, 8781.69,
    # PAN-STARRS
    4900.12, 6241.28, 7563.76, 8690.10, 9644.63,
    # VISTA VHS http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/filter-set
    #http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?&mode=search&search_text=vista
    12523.90, 16449.01, 21460.39,
    # GALEX
    2297
                     ])
    #BAND WIDTH 
    fwhm = np.array([
    #GAIA   
    2157.50, 4052.97, 2924.44,
    #2MASS
    1624.32, 2509.40, 2618.87,
    #WISE
    6626.42, 10422.66, 55055.71, 41016.83,
    #skymapper
    425.0, 318.8, 1477.5, 1523.3, 1201.6, 1112.8,       
    425.0, 318.8, 1477.5, 1523.3, 1201.6, 1112.8,           
    #SDSS
    582.28, 1262.68, 1149.52, 1238.95, 994.39,
    582.28, 1262.68, 1149.52, 1238.95, 994.39,
    #PanSTARRS                 
    1053.08, 1252.41, 1206.63, 997.71, 638.99,
    #VISTA VHS
    #1720.0, 2910.0, 3090.0, 
    1725.17, 2905.59, 3075.13,
    # GALEX
    732
               ])
    stdev = fwhm/2.355
    ##
    blo = np.array([
    #GAIA           3292.83, 3294.02, 6196.05,
                3390., 3880., 6370.,
    #2MASS          10806.47, 14787.38, 19543.69,
                11035., 15000., 19790.,                  
    #WISE           27540.97, 39633.26, 74430.44,
                27950., 40375., 75400., 197000.,
    #                195200.83
    #                201642.
    #SkyMapper      3067.30, 3550.45, 4102.78, 4925., 6929.17, 8159.09,
                3175., 3600., 4175., 5200., 7050., 8300.,
                3175., 3600., 4175., 5200., 7050., 8300.,
    #SDSS (trasmissione >1%, non >10%)
                3048.28, 3782.54, 5415.34, 6689.47, 7960.44,
                3048.28, 3782.54, 5415.34, 6689.47, 7960.44,
    #PanSTARRS      3943.40, 5386.23, 6778.45, 8028.00, 9100.50,
                3990., 5445.,  6850., 8115., 9150.,        
    #VISTA VHS      11000., 19500.,
    #                11035., 19790.,  
                11429.84, 14637.16, 19388.52,
    # GALEX
                1799.
                ])
    #BAND UPPER LIMIT 
    bup = np.array([
    #GAIA           6738.11, 10301.96, 10422.96,
                6660., 9900., 10150.,   
    #2MASS          14067.97, 18231.02, 23552.40,
                13510., 18040., 23525.,
    #WISE           38723.88, 53413.60, 172613.43,
                38225., 52700., 167500., 261800.,
    #SkyMapper      3866.93, 4216.67, 6570.00, 7231.82, 8647.37, 10679.17,
                3800., 4100., 6200., 7050., 8550., 10450.,
                3800., 4100., 6200., 7050., 8550., 10450.,
    #SDSS (trasmissione >1%, non >10%)
                4028.23, 5549.26, 6989.14, 8389.45, 10833.25,
                4028.23, 5549.26, 6989.14, 8389.45, 10833.25,
    #PanSTARRS      5593.27, 7035.65, 8304.37, 9346.00, 10838.50,
                5545., 6965., 8245., 9275., 10330., 
    #VISTA VHS      14000., 23500.,
    #                13510., 23525.,
                13668.39, 18340.58, 23661.88,
    # GALEX
                2805.
                ])
    ##
    z1fit = []
    z1fit =[-999. for x in range(nrow)]
    r_fit = []
    r_fit =[-999. for x in range(nrow)]
    s_fit = []
    s_fit =[-999. for x in range(nrow)]
    n_fit = []
    n_fit =[-999. for x in range(nrow)]
    chi = []
    chi =[-999. for x in range(nrow)]
    offset=np.full(nrow, np.nan)

    mag_names=['gaia_BP', 'gaia_G', 'gaia_RP',
               'mag_J', 'mag_H', 'mag_K',
               'unwise_w1', 'unwise_w2','wise_w3', 'wise_w4',
               'Q1_SkyM11_mag_u', 'Q1_SkyM11_mag_v', 'Q1_SkyM11_mag_g', 'Q1_SkyM11_mag_r', 'Q1_SkyM11_mag_i', 'Q1_SkyM11_mag_z',
               'Q1_SkyM3_mag_u', 'Q1_SkyM3_mag_v', 'Q1_SkyM3_mag_g', 'Q1_SkyM3_mag_r', 'Q1_SkyM3_mag_i', 'Q1_SkyM3_mag_z',
               'SDSS_DR14Q_mag_u', 'SDSS_DR14Q_mag_g', 'SDSS_DR14Q_mag_r', 'SDSS_DR14Q_mag_i', 'SDSS_DR14Q_mag_z',
               'SDSS_DR16Q_mag_u', 'SDSS_DR16Q_mag_g', 'SDSS_DR16Q_mag_r', 'SDSS_DR16Q_mag_i', 'SDSS_DR16Q_mag_z',
               'PanSTARRS1DR2_mag_g', 'PanSTARRS1DR2_mag_r', 'PanSTARRS1DR2_mag_i', 'PanSTARRS1DR2_mag_z', 'PanSTARRS1DR2_mag_Y',
               'jPetroMag', 'hPetroMag', 'ksPetroMag',
               'NUVmag']
    mag_idx=range(len(mag_names))
    mag_dict=dict(zip(mag_names,mag_idx))
    offset_idx=mag_dict[offset_band]
    ###
    #   BIG LOOP <<<<  BIG LOOP ##
    #   Grande loop su tutti gli oggetti nella tabella  #
    start = time.time()
    for i in range (nrow):
        # RIEMPIO L'ARRAY FL CON LE MAG AB - ATTENTO ALL'ORDINE
        # GAIA Michael Weiler A&A 617, A138 (2018)
        fl[0] = phototab['gaia_BP'][i]-25.362+25.3888
        fl[1] = phototab['gaia_G'][i]-25.6409+25.7455
        fl[2] = phototab['gaia_RP'][i]-24.7600+25.1185
        # 2MASS http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
        fl[3] = phototab['mag_J'][i]+0.91
        fl[4] = phototab['mag_H'][i]+1.39
        fl[5] = phototab['mag_K'][i]+1.85
        # WISE https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html
        fl[6] = phototab['unwise_w1'][i]
        fl[7] = phototab['unwise_w2'][i]
        fl[8] = phototab['wise_w3'][i]
        fl[9] = phototab['wise_w4'][i]+0.1
        # skymapper
        fl[10] = phototab['Q1_SkyM11_mag_u'][i]
        fl[11] = phototab['Q1_SkyM11_mag_v'][i]
        fl[12] = phototab['Q1_SkyM11_mag_g'][i]
        fl[13] = phototab['Q1_SkyM11_mag_r'][i]
        fl[14] = phototab['Q1_SkyM11_mag_i'][i]
        fl[15] = phototab['Q1_SkyM11_mag_z'][i]
        fl[16] = phototab['Q1_SkyM3_mag_u'][i]
        fl[17] = phototab['Q1_SkyM3_mag_v'][i]
        fl[18] = phototab['Q1_SkyM3_mag_g'][i]
        fl[19] = phototab['Q1_SkyM3_mag_r'][i]
        fl[20] = phototab['Q1_SkyM3_mag_i'][i]
        fl[21] = phototab['Q1_SkyM3_mag_z'][i]
        # SDSS
        fl[22] = phototab['SDSS_DR14Q_mag_u'][i]
        fl[23] = phototab['SDSS_DR14Q_mag_g'][i]
        fl[24] = phototab['SDSS_DR14Q_mag_r'][i]
        fl[25] = phototab['SDSS_DR14Q_mag_i'][i]
        fl[26] = phototab['SDSS_DR14Q_mag_z'][i]
        fl[27] = phototab['SDSS_DR16Q_mag_u'][i]
        fl[28] = phototab['SDSS_DR16Q_mag_g'][i]
        fl[29] = phototab['SDSS_DR16Q_mag_r'][i]
        fl[30] = phototab['SDSS_DR16Q_mag_i'][i]
        fl[31] = phototab['SDSS_DR16Q_mag_z'][i] 
        # PAN-STARRS
        fl[32] = phototab['PanSTARRS1DR2_mag_g'][i]
        fl[33] = phototab['PanSTARRS1DR2_mag_r'][i]
        fl[34] = phototab['PanSTARRS1DR2_mag_i'][i]
        fl[35] = phototab['PanSTARRS1DR2_mag_z'][i]
        fl[36] = phototab['PanSTARRS1DR2_mag_Y'][i]
        # VIS VHS http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/filter-set
        fl[37] = phototab['jPetroMag'][i]+0.916
        fl[38] = phototab['hPetroMag'][i]+1.366
        fl[39] = phototab['ksPetroMag'][i]+1.827
        # GALEX
        fl[40] = phototab['NUVmag'][i]
        #
        fl = np.array(fl)
        # IN FL CI SONO MAGNITUDINI AB
        lfl = np.array(fl)
        # ERRORS - ATTENTO ALL'ORDINE !!
        err[0]  = phototab['gaia_BP_err'][i]
        err[1]  = phototab['gaia_G_err'][i]
        err[2]  = phototab['gaia_RP_err'][i]
        err[3]  = phototab['mag_J_err'][i]
        err[4]  = phototab['mag_H_err'][i]
        err[5]  = phototab['mag_K_err'][i]
        err[6]  = phototab['unwise_w1_err'][i]
        err[7]  = phototab['unwise_w2_err'][i]
        err[8]  = phototab['wise_w3_err'][i]
        err[9]  = phototab['wise_w4_err'][i]
        err[10] = phototab['Q1_SkyM11_mag_u_err'][i]
        err[11] = phototab['Q1_SkyM11_mag_v_err'][i]
        err[12] = phototab['Q1_SkyM11_mag_g_err'][i]
        err[13] = phototab['Q1_SkyM11_mag_r_err'][i]
        err[14] = phototab['Q1_SkyM11_mag_i_err'][i]
        err[15] = phototab['Q1_SkyM11_mag_z_err'][i]
        err[16] = phototab['Q1_SkyM3_mag_u_err'][i]
        err[17] = phototab['Q1_SkyM3_mag_v_err'][i]
        err[18] = phototab['Q1_SkyM3_mag_g_err'][i]
        err[19] = phototab['Q1_SkyM3_mag_r_err'][i]
        err[20] = phototab['Q1_SkyM3_mag_i_err'][i]
        err[21] = phototab['Q1_SkyM3_mag_z_err'][i]
        err[22] = phototab['SDSS_DR14Q_mag_u_err'][i]
        err[23] = phototab['SDSS_DR14Q_mag_g_err'][i]
        err[24] = phototab['SDSS_DR14Q_mag_r_err'][i]
        err[25] = phototab['SDSS_DR14Q_mag_i_err'][i]
        err[26] = phototab['SDSS_DR14Q_mag_z_err'][i]
        err[27] = phototab['SDSS_DR16Q_mag_u_err'][i]
        err[28] = phototab['SDSS_DR16Q_mag_g_err'][i]
        err[29] = phototab['SDSS_DR16Q_mag_r_err'][i]
        err[30] = phototab['SDSS_DR16Q_mag_i_err'][i]
        err[31] = phototab['SDSS_DR16Q_mag_z_err'][i]
        err[32] = phototab['PanSTARRS1DR2_mag_g_err'][i]
        err[33] = phototab['PanSTARRS1DR2_mag_r_err'][i]
        err[34] = phototab['PanSTARRS1DR2_mag_i_err'][i]
        err[35] = phototab['PanSTARRS1DR2_mag_z_err'][i]
        err[36] = phototab['PanSTARRS1DR2_mag_Y_err'][i]
        err[37] = phototab['jPetroMagErr'][i]
        err[38] = phototab['hPetroMagErr'][i]
        err[39] = phototab['ksPetroMagErr'][i]
        err[40] = phototab['e_NUVmag'][i]
        # CHECK ERRORS - must be > emagmin.
        emagmin = 0.1
        for ee in range (0, len(err)):
            if not (err[ee] >emagmin):
                err[ee] = emagmin
        err = np.array(err)
        #
        #se l'oggetto non ha la magnitudine necessaria, lo salto
        if not np.isfinite(fl[offset_idx]) or abs(fl[offset_idx])>50:
            continue

        wo = []
        lflo = []
        erro = []
        stdo  = []
        bloo = []
        bupo = []
        for m in range (0, len(w)):
            if np.isfinite(lfl[m]) and abs(lfl[m])<50 and abs(lfl[m])>2 :
                wo.append(w[m])
                lflo.append(lfl[m])
                erro.append(err[m])
                stdo.append(stdev[m])
                bloo.append(blo[m])
                bupo.append(bup[m])
        wo = np.array(wo)
        lflo = np.array(lflo)
        erro = np.array(erro)
        stdo = np.array(stdo)
        bloo = np.array(bloo)
        bupo = np.array(bupo)
        
        #per il fit devo usare degli array senza la banda di cui devo calcolare l'offset
        wo2 = np.delete(wo, offset_idx)
        lflo2 = np.delete(lflo, offset_idx)
        erro2 = np.delete(erro, offset_idx)
        stdo2 = np.delete(stdo, offset_idx)
        bloo2 = np.delete(bloo, offset_idx)
        bupo2 = np.delete(bupo, offset_idx)
        
        #PERFORM FIT
        z_free=False
        if True:
            # PARAMETERS FOR THE MINIMIZATION  https://lmfit.github.io/lmfit-py/fitting.html
            params = Parameters()
            if z_free:
                #params.add('zf1', value=2.54, min=1.4, max=6.1)
                params.add('zf1', value=2.54, min=1.1, max=7.1)
            else:
                params.add('zf1', value=(z_spec[i]+1.), vary=False)
            params.add('ratio', value=0.59, min=0.1, max=0.9)
            params.add('norm', value=-0.7, min=-3., max=3.)
            params.add('slope', value=0., min=-0.5, max=0.5)
            out = minimize(residual, params, args=(wsed, flsedA, flsedB, wo2, stdo2, bloo2, bupo2, lflo2, erro2, IGM_tau),
                    nan_policy = 'omit', method='dual_annealing')

            z1fit[i] = out.params['zf1'].value
            r_fit[i] = out.params['ratio'].value
            n_fit[i] = out.params['norm'].value
            s_fit[i] = out.params['slope'].value
            chi[i] =  out.chisqr


        #  RECOMPUTE THE BEST FIT PHOTOMETRY
        #  PUT PHOTOMETRIC OBSERVATIONS IN LOG SCALE OF lambda * Flambda
        #  phot_from_sed(wo, stdev, blo, bup, wsed, lfsed, zf1, slope, norm, ampl, IGM_trans):
        lw = np.log10(wo)
        # DEVO METTERE GLI UPPER LIMITS
        uplims = [False]*len(wo)
        lf=[0.]*len(wo)
        for kk in range (0, len(wo)):
            if lflo[kk]>0. :
            #  VOGLIO PER COMODITÀ UNA SED CHE SIA 0. PER MAG_AB=17 @ 10^4 A
                lf[kk] = -0.4*lflo[kk] -lw[kk] + 10.8
            elif lflo[kk]<0. :
            # SE f<0. SI TRATTA DI UN UPPER LIMIT
                lf[kk] = 0.4*lflo[kk] -lw[kk] + 10.8
                uplims[kk] = True
        le = 0.4*erro
        # in lw,lf ci sono i DATI
        photfit = [0.] * len(wo)
        photfit = phot_from_sed_AB(wo, stdo, bloo, bupo, wsed, flsedA, flsedB, z1fit[i], r_fit[i], n_fit[i], s_fit[i], IGM_tau)
        # in photfit il fit ottimale
        #calcolo la differenza tra il fit e la magnitudine osservata
        offset[i] = (photfit[offset_idx]-lf[offset_idx])
    #
    z_fit = np.array(z1fit)-1.
    phototab_out = phototab
    col_a = Column(name='z_fit', data=z_fit)
    phototab_out.add_column(col_a)
    col_b = Column(name='ratio_fit', data=r_fit)
    phototab_out.add_column(col_b)
    col_c = Column(name='norm_fit', data=n_fit)
    phototab_out.add_column(col_c)
    col_d = Column(name='slope_fit', data=s_fit)
    phototab_out.add_column(col_d)
    col_e = Column(name='chisq', data=chi)
    phototab_out.add_column(col_e)
    col_f_name='offset_'+offset_band
    col_f = Column(name=col_f_name, data=offset)
    phototab_out.add_column(col_f)
    phototab_out.write('./Offsets/QSO_offset_{}_out.fits'.format(offset_band), overwrite=True)
    return phototab_out

###########################################################################################
import astropy
import math
from matplotlib.pyplot import *
import numpy as np
# import pylab
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from astropy.table import QTable, Column, MaskedColumn, vstack
from numpy.polynomial import Polynomial
from lmfit import minimize, Parameters, fit_report
import time
import multiprocessing as mp
from FromPhot_Offsets import photometric_offset, init_lock
import pandas as pd

##################   SET INITIAL OPTIONS   ###################
if __name__ == '__main__':
    z_selection='ALL'
    if z_selection=='ALL':
        phot_file='./QSO_Bright4.fits'
    elif z_selection=='TEST':
        phot_file='./TargetsTNG.fits'
    else:
        raise(Exception('Please select a subset'))    
    hdu_list = fits.open(phot_file, memmap=True)
    phototab_in = Table(hdu_list[1].data)
    hdu_list.close()
    
    mag_names=['gaia_BP', 'gaia_G', 'gaia_RP',
               'mag_J', 'mag_H', 'mag_K',
               'unwise_w1', 'unwise_w2','wise_w3', 'wise_w4',
               'Q1_SkyM11_mag_u', 'Q1_SkyM11_mag_v', 'Q1_SkyM11_mag_g', 'Q1_SkyM11_mag_r', 'Q1_SkyM11_mag_i', 'Q1_SkyM11_mag_z',
               'Q1_SkyM3_mag_u', 'Q1_SkyM3_mag_v', 'Q1_SkyM3_mag_g', 'Q1_SkyM3_mag_r', 'Q1_SkyM3_mag_i', 'Q1_SkyM3_mag_z',
               'SDSS_DR14Q_mag_u', 'SDSS_DR14Q_mag_g', 'SDSS_DR14Q_mag_r', 'SDSS_DR14Q_mag_i', 'SDSS_DR14Q_mag_z',
               'SDSS_DR16Q_mag_u', 'SDSS_DR16Q_mag_g', 'SDSS_DR16Q_mag_r', 'SDSS_DR16Q_mag_i', 'SDSS_DR16Q_mag_z',
               'PanSTARRS1DR2_mag_g', 'PanSTARRS1DR2_mag_r', 'PanSTARRS1DR2_mag_i', 'PanSTARRS1DR2_mag_z', 'PanSTARRS1DR2_mag_Y',
               'jPetroMag', 'hPetroMag', 'ksPetroMag',
               'NUVmag']    
    
    if True:
        input_list=[]
        for band in mag_names:
            input_list.append((phototab_in.copy(), band))
            
        proc_num=4
        l=mp.Lock()
        p = mp.Pool(processes = proc_num, initializer=init_lock, initargs=(l,))
        start = time.time()
        async_result = p.starmap_async(photometric_offset, input_list)
        p.close()
        p.join()
        print("Complete")
        end = time.time()
        print('total time (s)= ' + str(end-start))       
        
        phototabs_out = async_result.get()
        
