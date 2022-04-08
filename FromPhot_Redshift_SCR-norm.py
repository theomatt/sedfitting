def comp_mod_sed(wsed, flsedA, flsedB, ratio,slope):
 # MODELLIAMO LA SED COME DUE COMPONENTI IN FLUSSO
# L'ACCRETION DISK - OTTICO UV - SEDA
# l'IR BUMP- SEDB
# La SEDA vale, attorno a 1 micron, 0.41
# La SEDB vale, attrono a 1 micron, 0.59
# Le combino in modo che comunque attorno a 1 micron la somma valga 1)
 lfsed = np.log10(ratio*flsedB/0.59 + (1.-ratio)*flsedA/0.41 * (wsed/2200.)**slope )
 return lfsed
#
#
def phot_from_sed_AB(wo, stdev, blo, bup, wsed, flsedA, flsedB, zf1, ratio, norm, slope, IGM_tau):
 lwsed = np.log10(wsed)
 lfsed = comp_mod_sed(wsed, flsedA, flsedB, ratio, slope)
# renormalize
 lfcorr = lfsed + norm
# lfcorr sono i log dei flussi della SED - NB nel rest frame
# BRING TO OBSERVER'S FRAME
# wzsed sono le wavelengths osservate [redshiftate] NB in ordine crescente 
 wzsed = wsed*zf1
#
# lfm CONTERRÀ IL MODELLO DELLA FOTOMETRIA IN BASE ALLA SED, ALLO z, SLOPE E NORM
 lfm = [0.] * len(wo)
#
# IN QUESTA VERSIONE ANDIAMO FILTRO PER FILTRO
# DATO CHE I FILTRI POSSONO NON ESSERE ORDINATI IN WAVELENGTH E ANCHE OVERLAPPARE
# MENTRE LE WZSED SONO ORDINATE
#
# INNANZITUTTO CREIAMO LA SED ASSORBITA DALL'IGM
# DATO LO z
#
 for kk in range (0, len(wsed)) :
# IF NECESSARY CALCOLA L'ASSORBIMENTO
# CALCOLA IL TAU IGM - nel res frame >> wsed
     if wsed[kk] < 1300. :         
        iw = int(round(wsed[kk]))
        z10 = int(round(10.*(zf1-1.)))
        key = (iw, z10)
        lfcorr[kk] -= IGM_tau[key] #igmtau deve essere positivo
        if IGM_tau[key]<0:
            print("TAU NOT OK")
     else :
        break
#
# ORA PER OGNI FILTRO
 for ll in range (0, len(wo)):
# INIZIALIZZIAMO IL PESO
      wps = 0.
      for j in range (0, len(wzsed)):
        if blo[ll] < wzsed[j] and bup[ll] > wzsed[j] :
# SE SIAMO QUI IL PUNTO DELLA SED È DENTRO LA BANDA DEL FILTRO
# CALCOLO IL PESO - LINEARMENTE DATO CHE STO USANDO QUANTITÀ LOGARITMICHE
# E VOGLIO ESSERE RAPIDO
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
###################################

def residual(params, wsed, flsedA, flsedB, pixsed, w, stdev, blo, bup, f, err, IGM_tau):
#
### PARAMETRI INIZIALI
 plotview = False
 zf1 = params['zf1'].value
 ratio = params['ratio'].value
 norm = params['norm'].value
 slope = params['slope'].value
 npar = 4
##  MINIMUM ERROR IN THE PHOTOMETRY - E.G. DUE TO VARIABILITY
 emagmin = 0.1
# 
# COMPUTE MODEL OF THE OBSERVATIONS
# FROM THE SED
#
# PUT PHOTOMETRIC OBSERVATIONS IN LOG SCALE OF lambda * Flambda
#
# REMEMBER THAT LOG (l Fl) = -0.4 m_AB -19.436 - LOG l
# MA IO VOGLIO PER COMODITÀ UNA SED CHE SIA 0. PER MAG_AB=17 @ 10^4 A
#
 lw = np.log10(w)
 le = 0.4*err
#
 lf = [0.] * len(f)
 lfm = [0.] * len(lf)
 lfm = phot_from_sed_AB(w, stdev, blo, bup, wsed, flsedA, flsedB, zf1, ratio, norm, slope, IGM_tau)
 restot = 0.
 nphot = 0
 for ll in range (0, len(w)):
### ATTENZIONE AGLI ERRORI TROPPO PICCOLI (SI PENSI ALLA VARIABILITÀ)
### SE L'ERRORE IN MAG È < emagmin ALLORA LO PORTO A emagmin
         sigma = emagmin if emagmin > le[ll] else le[ll]
         if f[ll]>0. :
         #  VOGLIO PER COMODITÀ UNA SED CHE SIA 0. PER MAG_AB=17 @ 10^4 A
             lf[ll] = -0.4*f[ll] -lw[ll] + 10.8
             res = (lf[ll]-lfm[ll])/sigma
             nphot += 1
         elif f[ll]<0. :
         # SE f<0. SI TRATTA DI UN UPPER LIMIT
             lf[ll] = 0.4*f[ll] -lw[ll] + 10.8
             if lfm[ll] > lf[ll] :
                 res = (lf[ll]-lfm[ll])/sigma
                 nphot += 1
             else:
                 res = 0.
         restot += abs(res) 
#
# DIVIDO PER IL NUMERO DI PUNTI -1 SU CUI HO FATTO IL CONFRONTO
# MENO IL NUMERO DI PARAMETRI npar
 redux = nphot-1-npar if (nphot-1-npar) > 0 else 1
# AGGIUNGO LA CORREZIONE PER PASSARE DA MAD a SIGMA ??? 0.79788456
 restot = restot/redux/0.79788456
 return restot

###########################################################################################
def residualv(params, wsed, flsedA, flsedB, pixsed, w, stdev, blo, bup, f, err, IGM_tau):
#
### PARAMETRI INIZIALI
 zf1 = params['zf1'].value
 ratio = params['ratio'].value
 norm = params['norm'].value
 slope = params['slope'].value
 npar = 4
##  MINIMUM ERROR IN THE PHOTOMETRY - E.G. DUE TO VARIABILITY
 emagmin = 0.1
#  
# COMPUTE MODEL OF THE OBSERVATIONS
# FROM THE SED
#
# PUT PHOTOMETRIC OBSERVATIONS IN LOG SCALE OF lambda * Flambda
#
# REMEMBER THAT LOG (l Fl) = -0.4 m_AB -19.436 - LOG l
# MA IO VOGLIO PER COMODITÀ UNA SED CHE SIA 0. PER MAG_AB=17 @ 10^4 A
#
 lw = np.log10(w)
#  VOGLIO PER COMODITÀ UNA SED CHE SIA 0. PER MAG_AB=17 @ 10^4 A
 le = 0.4*err
#
 lf = [0.] * len(f)
 lfm = [0.] * len(lf)
 res = [0.] * len(lf)
 lfm = phot_from_sed_AB(w, stdev, blo, bup, wsed, flsedA, flsedB, zf1, ratio, norm, slope, IGM_tau)
 for ll in range (0, len(w)):
### ATTENZIONE AGLI ERRORI TROPPO PICCOLI (SI PENSI ALLA VARIABILITÀ)
### SE L'ERRORE IN MAG È < emagmin ALLORA LO PORTO A emagmin
         sigma = emagmin if emagmin > le[ll] else le[ll]
         if f[ll]>0. :
         #  VOGLIO PER COMODITÀ UNA SED CHE SIA 0. PER MAG_AB=17 @ 10^4 A
             lf[ll] = -0.4*f[ll] -lw[ll] + 10.8
             res[ll] = (lf[ll]-lfm[ll])/sigma
         elif f[ll]<0. :
         # SE f<0. SI TRATTA DI UN UPPER LIMIT
              lf[ll] = 0.4*f[ll] -lw[ll] + 10.8
              if lfm[ll] > lf[ll] :
                  res[ll] = (lf[ll]-lfm[ll])/sigma
              else:
                  res[ll] = 0.         
 return res
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
from astropy.table import QTable, Column, MaskedColumn
from numpy.polynomial import Polynomial
from lmfit import minimize, Parameters, fit_report
import time
##################   SET INITIAL OPTIONS   ###################
nbands=41
w = [None] * nbands
fl = [None] * nbands
err = [None] * nbands
errm = [None] * nbands
plotall = False
plotcum = True
# LA TABELLA CON LA FOTOMETRIA DEI QSOs  
# scegliere se usare tutti i quasar, la metà a z più basso o la metà a z più alto
z_selection='ALL'
if z_selection=='HIGHZ':
    phototab='./QSO_Bright_HighZ.fits'
elif z_selection=='LOWZ':
    phototab='./QSO_Bright_LowZ.fits'
elif z_selection=='ALL':
    phototab='./QSO_Bright2.fits'
elif z_selection=='TEST':
    phototab='./QSOreduv.fits'
else:
    raise(Exception('Please select a subset'))
    
# LA TABELLA CON LA SED APPROSSIMATA
SEDtab='QSOsedA.dat'
wsed =  []
flsed = []
efsed = []
pixsed= []
rmsed = []
npsed = []
p = Parameters()
#
sed = open(SEDtab, 'r')
for line in sed:
     line = line.strip()
     columns = line.split()
     wsed.append(float(columns[0]))
     flsed.append(float(columns[1]))
     efsed.append(float(columns[2]))
     pixsed.append(float(columns[3]))
     rmsed.append(float(columns[4]))
     npsed.append(float(columns[5]))
sed.close()
wsed = np.array(wsed)
flsedA = np.array(flsed)
# SECONDA SED DA LEGGERE
SEDtab='QSOsedB.dat'
wsed =  []
flsed = []
efsed = []
pixsed= []
rmsed = []
npsed = []
p = Parameters()
#
sed = open(SEDtab, 'r')
for line in sed:
     line = line.strip()
     columns = line.split()
     wsed.append(float(columns[0]))
     flsed.append(float(columns[1]))
     efsed.append(float(columns[2]))
     pixsed.append(float(columns[3]))
     rmsed.append(float(columns[4]))
     npsed.append(float(columns[5]))
sed.close()
wsed = np.array(wsed)
flsedB = np.array(flsed)
#
# Leggi la tabella fits con la fotometria 
hdu_list = fits.open(phototab, memmap=True)
#hdu_list.info()
#print(hdu_list[1].columns)
#
# OCCORRE LEGGERE RIGA PER RIGA?
#
t2 = Table(hdu_list[1].data)
hdu_list.close()

have_output=False
if have_output:
    #import the output table of an earlier run
    QSOfitfile='./PhotZ/Run9 - 2200A with unwise/QSO_PhotZ_out_e01_zfree.fits'
    hdu_list = fits.open(QSOfitfile, memmap=True)
    fittab = Table(hdu_list[1].data)
    hdu_list.close()
    
    unwise_w1=fittab['unwise_w1']
    norm_fit=fittab['norm_fit']
    
    #filter nan values
    W1=unwise_w1[~np.isnan(unwise_w1)]
    norm_clean=norm_fit[~np.isnan(unwise_w1)]
    
    #perform linear regression between W1 magnitude and norm_fit
    coef_norm=np.polyfit(W1,norm_clean,1)
    norm_func=np.poly1d(coef_norm)
else:
    #use default coefficients
    coef_norm=np.array([-0.29746476793722687, 4.47048318731566])
    norm_func=np.poly1d(coef_norm)
    
#remove objects without w1 magnitude
mask=np.isfinite(t2['unwise_w1'])
t2=t2[mask]

###
nrow = len(t2)
id = t2['qid']
z_spec = t2['z_spec']
AR = t2['RAd']
D = t2['DECd']
GA_b = t2['gaia_BP']
GA_g = t2['gaia_G']
GA_r = t2['gaia_RP']
J2M = t2['mag_J']
H2M = t2['mag_H']
K2M = t2['mag_K']
W1 =  t2['unwise_w1']
W2 =  t2['unwise_w2']
W3 =  t2['wise_w3']
W4 =  t2['wise_w4']
SM11_u = t2['Q1_SkyM11_mag_u']
SM11_v = t2['Q1_SkyM11_mag_v']
SM11_g = t2['Q1_SkyM11_mag_g']
SM11_r = t2['Q1_SkyM11_mag_r']
SM11_i = t2['Q1_SkyM11_mag_i']
SM11_z = t2['Q1_SkyM11_mag_z']
SM3_u = t2['Q1_SkyM3_mag_u']
SM3_v = t2['Q1_SkyM3_mag_v']
SM3_g = t2['Q1_SkyM3_mag_g']
SM3_r = t2['Q1_SkyM3_mag_r']
SM3_i = t2['Q1_SkyM3_mag_i']
SM3_z = t2['Q1_SkyM3_mag_z']
SDSS14_u = t2['SDSS_DR14Q_mag_u']
SDSS14_g = t2['SDSS_DR14Q_mag_g']
SDSS14_r = t2['SDSS_DR14Q_mag_r']
SDSS14_i = t2['SDSS_DR14Q_mag_i']
SDSS14_z = t2['SDSS_DR14Q_mag_z']
SDSS16_u = t2['SDSS_DR16Q_mag_u']
SDSS16_g = t2['SDSS_DR16Q_mag_g']
SDSS16_r = t2['SDSS_DR16Q_mag_r']
SDSS16_i = t2['SDSS_DR16Q_mag_i']
SDSS16_z = t2['SDSS_DR16Q_mag_z']
PS_g = t2['PanSTARRS1DR2_mag_g']
PS_r = t2['PanSTARRS1DR2_mag_r']
PS_i = t2['PanSTARRS1DR2_mag_i']
PS_z = t2['PanSTARRS1DR2_mag_z']
PS_y = t2['PanSTARRS1DR2_mag_Y']
#PsfMag or PetroMag? Y and H bands?
JVHS = t2['jPetroMag']
HVHS = t2['hPetroMag']
KVHS = t2['ksPetroMag']
# GALEX
NUV = t2['NUVmag']
# errors ATTENTO ALL'ORDINE!!
eGA_b = t2['gaia_BP_err']
eGA_g = t2['gaia_G_err']
eGA_r = t2['gaia_RP_err']
eJ2M = t2['mag_J_err']
eH2M = t2['mag_H_err']
eK2M = t2['mag_K_err']
eW1 =  t2['unwise_w1_err']
eW2 =  t2['unwise_w2_err']
eW3 =  t2['wise_w3_err']
eW4 =  t2['wise_w4_err']
eSM11_u = t2['Q1_SkyM11_mag_u_err']
eSM11_v = t2['Q1_SkyM11_mag_v_err']
eSM11_g = t2['Q1_SkyM11_mag_g_err']
eSM11_r = t2['Q1_SkyM11_mag_r_err']
eSM11_i = t2['Q1_SkyM11_mag_i_err']
eSM11_z = t2['Q1_SkyM11_mag_z_err']
eSM3_u = t2['Q1_SkyM3_mag_u_err']
eSM3_v = t2['Q1_SkyM3_mag_v_err']
eSM3_g = t2['Q1_SkyM3_mag_g_err']
eSM3_r = t2['Q1_SkyM3_mag_r_err']
eSM3_i = t2['Q1_SkyM3_mag_i_err']
eSM3_z = t2['Q1_SkyM3_mag_z_err']
eSDSS14_u = t2['SDSS_DR14Q_mag_u_err']
eSDSS14_g = t2['SDSS_DR14Q_mag_g_err']
eSDSS14_r = t2['SDSS_DR14Q_mag_r_err']
eSDSS14_i = t2['SDSS_DR14Q_mag_i_err']
eSDSS14_z = t2['SDSS_DR14Q_mag_z_err']
eSDSS16_u = t2['SDSS_DR16Q_mag_u_err']
eSDSS16_g = t2['SDSS_DR16Q_mag_g_err']
eSDSS16_r = t2['SDSS_DR16Q_mag_r_err']
eSDSS16_i = t2['SDSS_DR16Q_mag_i_err']
eSDSS16_z = t2['SDSS_DR16Q_mag_z_err']
ePS_g = t2['PanSTARRS1DR2_mag_g_err']
ePS_r = t2['PanSTARRS1DR2_mag_r_err']
ePS_i = t2['PanSTARRS1DR2_mag_i_err']
ePS_z = t2['PanSTARRS1DR2_mag_z_err']
ePS_y = t2['PanSTARRS1DR2_mag_Y_err']
eJVHS = t2['jPetroMagErr']
eHVHS = t2['hPetroMagErr']
eKVHS = t2['ksPetroMagErr']
eNUV = t2['e_NUVmag']
##
errm[0] = np.nanmedian(eGA_b)
errm[1] = np.nanmedian(eGA_g)
errm[2] = np.nanmedian(eGA_r)
errm[3] = np.nanmedian(eJ2M)
errm[4] = np.nanmedian(eH2M)
errm[5] = np.nanmedian(eK2M)
errm[6] = np.nanmedian(eW1)
errm[7] = np.nanmedian(eW2)
errm[8] = np.nanmedian(eW3)
errm[9] = np.nanmedian(eW4)
errm[10] = np.nanmedian(eSM11_u)
errm[11] = np.nanmedian(eSM11_v)
errm[12] = np.nanmedian(eSM11_g)
errm[13] = np.nanmedian(eSM11_r)
errm[14] = np.nanmedian(eSM11_i)
errm[15] = np.nanmedian(eSM11_z)
errm[16] = np.nanmedian(eSM3_u)
errm[17] = np.nanmedian(eSM3_v)
errm[18] = np.nanmedian(eSM3_g)
errm[19] = np.nanmedian(eSM3_r)
errm[20] = np.nanmedian(eSM3_i)
errm[21] = np.nanmedian(eSM3_z)
errm[22] = np.nanmedian(eSDSS14_u)
errm[23] = np.nanmedian(eSDSS14_g)
errm[24] = np.nanmedian(eSDSS14_r)
errm[25] = np.nanmedian(eSDSS14_i)
errm[26] = np.nanmedian(eSDSS14_z)
errm[27] = np.nanmedian(eSDSS16_u)
errm[28] = np.nanmedian(eSDSS16_g)
errm[29] = np.nanmedian(eSDSS16_r)
errm[30] = np.nanmedian(eSDSS16_i)
errm[31] = np.nanmedian(eSDSS16_z)
errm[32] = np.nanmedian(ePS_g)
errm[33] = np.nanmedian(ePS_r)
errm[34] = np.nanmedian(ePS_i)
errm[35] = np.nanmedian(ePS_z)
errm[36] = np.nanmedian(ePS_y)
errm[37] = np.nanmedian(eJVHS)
errm[38] = np.nanmedian(eHVHS)
errm[39] = np.nanmedian(eKVHS)
errm[40] = np.nanmedian(eNUV)
##
# Carica le lunghezze d'onda dei vari filtri
# usiamo il service
# http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=2MASS&asttype=
# NB - lambda mean!  not eff !!
# skymapper
# w = np.array([3496.61, 3837.71, 5099.45, 6157.30, 7778.37, 9161.71,
#  lambda eff
# w = np.array([3498.15, 3870.92, 4968.46, 6040.07, 7712.95, 9091.50,
# labda cen
w = np.array([
# GAIA
# 5319.87, 6719.55, 7939.10,
 5272.46, 6291.13, 7761.11,
# 2MASS
 12390.58, 16487.19, 21634.04,
# WISE
 # 33526.00, 46028.00, 115608.00, 220883.00,
 # Lambda Pivot
 33682.21, 46179.05, 120717.43, 221944.44,
 # 34655.16, 46443.00, 113081.34,
 # 132156.35,
 # 222228.81
 # Brown, Jarrett, Cluver 2014
 # 229562.36, 
# skymapper
 3503.42, 3828.56, 5033.64, 6174.52, 7806.40, 8873.64,
 3503.42, 3828.56, 5033.64, 6174.52, 7806.40, 8873.64,
# SDSS http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=SLOAN
 3565.05, 4700.33, 6174.48, 7533.63, 8781.69,
 3565.05, 4700.33, 6174.48, 7533.63, 8781.69,
# PAN-STARRS
 4900.12, 6241.28, 7563.76, 8690.10, 9644.63,
# VISTA VHS http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/filter-set
 #12520.00, 16450, 21470.00,
 #http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?&mode=search&search_text=vista
 12523.90, 16449.01, 21460.39,
# GALEX
 2297
                     ])
##
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
#BAND LOWER LIMIT
# PRENDIAMO GLI INTERVALLI IN CUI I FILTRI HANNO UNA TRASMISSIONE > 10% DEL PICCO
# ANCHE PER RIDURRE I TEMPI DI ESECUZIONE
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
#                279107.24,
#                288317.8
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

zfit_err = []
zfit_err =[-999. for x in range(nrow)]
rfit_err = []
rfit_err =[-999. for x in range(nrow)]
sfit_err = []
sfit_err =[-999. for x in range(nrow)]

corr_zr = []
corr_zr =[-999. for x in range(nrow)]
corr_zn = []
corr_zn =[-999. for x in range(nrow)]
corr_zs = []
corr_zs =[-999. for x in range(nrow)]
corr_rn = []
corr_rn =[-999. for x in range(nrow)]
corr_rs = []
corr_rs =[-999. for x in range(nrow)]
corr_sn = []
corr_sn =[-999. for x in range(nrow)]

z2fit = []
z2fit =[-999. for x in range(nrow)]
r2_fit = []
r2_fit =[-999. for x in range(nrow)]
s2_fit = []
s2_fit =[-999. for x in range(nrow)]
n2_fit = []
n2_fit =[-999. for x in range(nrow)]
chi2 = []
chi2 =[-999. for x in range(nrow)]
# CREATE IGM_TRANS E IGM_TAU DICTIONARIES
IGM_trans = {}
IGM_tau = {}
a_file = open("taueff_NzS70.txt")
for line in a_file:
    line = line.strip()
    value = line.split()
    iw = int(float(value[0]))
    for i in range(1, 71):
        key = (iw, i)
        IGM_trans[key] = float(value[i])
        IGM_tau[key] = -np.log10(IGM_trans[key])
a_file.close()
###
#   BIG LOOP <<<<  BIG LOOP ##
#   Grande loop su tutti gli oggetti nella tabella  #
"""correction = np.array(  [ [ np.nan for mm in range(len(wsed)) ] for nn in range(nrow) ] )
print(' Array = ', correction.ndim, correction.shape)
"""
start = time.time()
niter = nrow
# niter = 10
for i in range (nrow):
# RIEMPIO L'ARRAY FL CON LE MAG AB - ATTENTO ALL'ORDINE
 # GAIA Michael Weiler A&A 617, A138 (2018)
 fl[0] = GA_b[i]-25.362+25.3888
 fl[1] = GA_g[i]-25.6409+25.7455
 fl[2] = GA_r[i]-24.7600+25.1185
 # 2MASS http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
 fl[3] = J2M[i]+0.91
 fl[4] = H2M[i]+1.39
 fl[5] = K2M[i]+1.85
 # WISE https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html
 fl[6] = W1[i]
 fl[7] = W2[i]
 fl[8] = W3[i]
 fl[9] = W4[i]+0.1
 #fl[8] = W3[i]+5.174
# fl[20] = W4[i]+6.620
# https://ui.adsabs.harvard.edu/abs/2014PASA...31...49B/abstract
# Brown M.~J.~I., Jarrett T.~H., Cluver M.~E., 2014, PASA, 31, e049. doi:10.1017/pasa.2014.44
 #fl[9] = W4[i]+6.666
 # skymapper
 fl[10] = SM11_u[i]
 fl[11] = SM11_v[i]
 fl[12] = SM11_g[i]
 fl[13] = SM11_r[i]
 fl[14] = SM11_i[i]
 fl[15] = SM11_z[i]
 fl[16] = SM3_u[i]
 fl[17] = SM3_v[i]
 fl[18] = SM3_g[i]
 fl[19] = SM3_r[i]
 fl[20] = SM3_i[i]
 fl[21] = SM3_z[i]
 # SDSS
 fl[22] = SDSS14_u[i]
 fl[23] = SDSS14_g[i]
 fl[24] = SDSS14_r[i]
 fl[25] = SDSS14_i[i]
 fl[26] = SDSS14_z[i]
 fl[27] = SDSS16_u[i]
 fl[28] = SDSS16_g[i]
 fl[29] = SDSS16_r[i]
 fl[30] = SDSS16_i[i]
 fl[31] = SDSS16_z[i] 
 # PAN-STARRS
 fl[32] = PS_g[i]
 fl[33] = PS_r[i]
 fl[34] = PS_i[i]
 fl[35] = PS_z[i]
 fl[36] = PS_y[i]
 # VIS VHS http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/filter-set
 fl[37] = JVHS[i]+0.916
 fl[38] = HVHS[i]+1.366
 fl[39] = KVHS[i]+1.827
 # GALEX
 fl[40] = NUV[i]
#
 fl = np.array(fl)
# IN FL CI SONO MAGNITUDINI AB
# lfl = fl+np.log10(w)
# LE LASCIO ANCHE IN lfl, ma eventualmente posso cambiare
 lfl = np.array(fl)
# ERRORS - ATTENTO ALL'ORDINE !!
 # GAIA Michael Weiler A&A 617, A138 (2018)
 err[0] = eGA_b[i]
 err[1] = eGA_g[i]
 err[2] = eGA_r[i]
 # 2MASS http://www.astronomy.ohio-state.edu/~martini/usefuldataprint(err.html
 err[3] = eJ2M[i]
 err[4] = eH2M[i]
 err[5] = eK2M[i]
 # WISE https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html
 err[6] = eW1[i]
 err[7] = eW2[i]
 err[8] = eW3[i]
 err[9] = eW4[i]
 #SkyMapper
 err[10] = eSM11_u[i]
 err[11] = eSM11_v[i]
 err[12] = eSM11_g[i]
 err[13] = eSM11_r[i]
 err[14] = eSM11_i[i]
 err[15] = eSM11_z[i]
 err[16] = eSM3_u[i]
 err[17] = eSM3_v[i]
 err[18] = eSM3_g[i]
 err[19] = eSM3_r[i]
 err[20] = eSM3_i[i]
 err[21] = eSM3_z[i]
 #SDSS
 err[22] = eSDSS14_u[i]
 err[23] = eSDSS14_g[i]
 err[24] = eSDSS14_r[i]
 err[25] = eSDSS14_i[i]
 err[26] = eSDSS14_z[i]
 err[27] = eSDSS16_u[i]
 err[28] = eSDSS16_g[i]
 err[29] = eSDSS16_r[i]
 err[30] = eSDSS16_i[i]
 err[31] = eSDSS16_z[i]
 # PAN-STARRS
 err[32] = ePS_g[i]
 err[33] = ePS_r[i]
 err[34] = ePS_i[i]
 err[35] = ePS_z[i]
 err[36] = ePS_y[i]
 #VISTA VHS
 err[37] = eJVHS[i]
 err[38] = eHVHS[i]
 err[39] = eKVHS[i]
 # GALEX
 err[40] = eNUV[i]
##
# CHECK ERRORS - must be > emagmin.
 emagmin = 0.1
 for ee in range (0, len(err)):
     if not (err[ee] >emagmin):
         err[ee] = emagmin
 err = np.array(err)
#
# PRIMA DI PASSARE ALLA MINIMIZZAZIONE
# MEGLIO COSTRUIRE UN ARRAY  w, fl, err, stdo, blo, bup
# CON I SOLI DATI SIGNIFICATIVI
# ELIMINANDO I None
#
 wo = []
 lflo = []
 erro = []
 stdo  = []
 bloo = []
 bupo = []
 for m in range (0, len(w)):
     if np.isfinite(lfl[m]) and np.isfinite(err[m]) and abs(lfl[m])<50 and abs(lfl[m])>2 :
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
 
 #PERFORM FIT
 z_free=True
 if True:
    # PARAMETERS FOR THE MINIMIZATION  https://lmfit.github.io/lmfit-py/fitting.html
     params = Parameters()
     if z_free:
         params.add('zf1', value=2.54, min=1.4, max=6.1)
     else:
         params.add('zf1', value=(z_spec[i]+1.), vary=False)
     params.add('ratio', value=0.59, min=0.1, max=0.9)
     params.add('norm', value=(norm_func(W1[i])), vary=False)
     params.add('slope', value=0., min=-0.5, max=0.5)
     out = minimize(residual, params, args=(wsed, flsedA, flsedB, pixsed, wo, stdo, bloo, bupo, lflo, erro, IGM_tau),
                    nan_policy = 'omit', method='dual_annealing')

     z1fit[i] = out.params['zf1'].value
     r_fit[i] = out.params['ratio'].value
     n_fit[i] = out.params['norm'].value
     s_fit[i] = out.params['slope'].value
     chi[i] =  out.chisqr

### QUI INFILIAMO IL RE-FIT COL LEVENBERG ’leastsq’ out2 = minimize(residual, method='leastsq', params=out.params)
 ### https://lmfit.github.io/lmfit-py/confidence.html#label-confidence-advanced
### lmfit.report_fit(out2.params, min_correl=0.5)
 
 """togliere levenberg-marquadt e slope"""
 """provare a tenere slope e togliere LM e norm, fissandola in funzione di w1"""
    
    
 #REPEAT FIT WITH LEVENBERG-MARQUADT TO GET ERRORS (only if there are enough photometric points)
 leastsq=False
 if leastsq and len(params)<=len(lflo):
     params2 = Parameters()
     if z_free:
         params2.add('zf1', value=z1fit[i], min=1.4, max=6.1)   
     else:
         params2.add('zf1', value=(z_spec[i]+1.), vary=False)
     params2.add('ratio', value=r_fit[i], min=0.1, max=0.9)
     params2.add('norm', value=n_fit[i], min=-3., max=3.)
     params2.add('slope', value=s_fit[i], min=-0.5, max=0.5)
     out2 = minimize(residualv, params2, args=(wsed, flsedA, flsedB, pixsed, wo, stdo, bloo, bupo, lflo, erro, IGM_tau),
                     nan_policy = 'omit', method='leastsq')
     z1fit[i] = out2.params['zf1'].value
     r_fit[i] = out2.params['ratio'].value
     n_fit[i] = out2.params['norm'].value
     s_fit[i] = out2.params['slope'].value
     chi[i] =  out2.redchi
 
     zfit_err[i] = out2.params['zf1'].stderr
     rfit_err[i] = out2.params['ratio'].stderr
     nfit_err[i] = out2.params['norm'].stderr
     sfit_err[i] = out2.params['slope'].stderr
     
     if z_free and out2.params['zf1'].correl!=None:
         corr_zr[i] = out2.params['zf1'].correl['ratio']
         corr_zn[i] = out2.params['zf1'].correl['norm']
         corr_zs[i] = out2.params['zf1'].correl['slope']
     if out2.params['ratio'].correl !=None:
         corr_rn[i] = out2.params['ratio'].correl['norm'] 
         corr_rs[i] = out2.params['ratio'].correl['slope']
     if out2.params['slope'].correl !=None:
        corr_sn[i] = out2.params['slope'].correl['norm']
     
 if i % 10 == 0 :
    end = time.time()
    print((end-start), 'time elapsed')
    if i > 0:
        print('remaining time = ',((end-start)*float(niter-i)/float(i)/60.),' min')
    print('iter =',i,' z_spec=', z_spec[i],' z_fit=', (z1fit[i]-1.))
    print(fit_report(out))
    if leastsq:
        print(fit_report(out2))
#
 zfs = str(z1fit[i]-1.)
 rts = str(r_fit[i])
 sls = str(s_fit[i])
 p0s = str(n_fit[i])
 zsps = str(z_spec[i])
 schi = str(chi[i])
#
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
# wzsed sono le wavelengths osservate [redshiftate]
 wzsed = wsed*z1fit[i]
 lwzsed = np.log10(wzsed)
 
 #REPEAT FIT WITH FIXED Z
 repeat2=z_free
 if repeat2:
     params = Parameters()
     params.add('zf1', value=(z_spec[i]+1.), vary=False)         
     params.add('ratio', value=0.59, min=0.1, max=0.9)
     params.add('norm', value=(norm_func(W1[i])), vary=False)
     params.add('slope', value=0., min=-0.5, max=0.5)
     out = minimize(residual, params, args=(wsed, flsedA, flsedB, pixsed, wo, stdo, bloo, bupo, lflo, erro, IGM_tau),
                    nan_policy = 'omit', method='dual_annealing')
    #
     z2fit[i]  = out.params['zf1'].value
     r2_fit[i] = out.params['ratio'].value
     n2_fit[i] = out.params['norm'].value
     s2_fit[i] = out.params['slope'].value
     chi2[i] =  out.chisqr

     zfs2 = str(z2fit[i]-1.)
     sls2 = str(s2_fit[i])
     p0s2 = str(n2_fit[i])
     schi2 = str(chi2[i])
    # wzsed sono le wavelengths osservate [redshiftate]
     wzsed2 = wsed*z2fit[i]
     lwzsed2= np.log10(wzsed2)

#  PLOT RESULT

 if plotcum:
     print('iter =',i,'\n z_fit = ', z1fit[i]-1.,'\n ratio_fit = ', r_fit[i], '\n norm_fit = ', n_fit[i],
           '\n slope_fit = ', s_fit[i],'\n chi = ', chi[i])
     fig=plt.figure(figsize=(20,10))
     ax=fig.add_axes([0.1,0.1,0.85,0.85])
     ax.scatter(lw, lf, color='b')
     ax.errorbar(lw, lf, yerr=le, uplims=uplims, fmt='o', color='b')
     ax.set_xlabel('log_wave')
     ax.set_ylabel('log_flux')
     ax.scatter(lw, photfit, color='k', marker="+")
# COMPUTE MODEL SED
     lwsed = np.log10(wsed)
     lfsed = comp_mod_sed(wsed, flsedA, flsedB, r_fit[i], s_fit[i])
     lfp = np.array(lfsed) + n_fit[i] 
     lfbuf = lfp
     for pp in range (0, len(wsed)):
           if wsed[pp] < 1300.:
             iw = int(round(wsed[pp]))
             z10 = int(round(10.*(z1fit[i]-1.)))
             key = (iw, z10)
             labs =  IGM_tau[key]
             lfbuf[pp] -= labs
     lfp1 = lfbuf
     ax.plot(lwzsed, lfp1, color='r')

# COMPUTE MODEL SED AT Z_SPEC
     if repeat2:
         lfsed = comp_mod_sed(wsed, flsedA, flsedB, r2_fit[i], s2_fit[i])
         lfp = np.array(lfsed) + n2_fit[i]
         lfbuf = lfp
         for pp in range (0, len(wsed)):
               if wsed[pp] < 1300.:
                 iw = int(round(wsed[pp]))
                 z10 = int(round(10.*(z2fit[i]-1.)))
                 key = (iw, z10)
                 labs =  IGM_tau[key]
                 lfbuf[pp] -= labs
         lwsed = np.log10(wsed)
         lfp2 = lfbuf
         ax.plot(lwzsed2, lfp2, color='g')
     
#     title = '  z=' + str(z1fit[i]-1.)+ '  slope=' + str(s_fit[i]) +'  zpo =' + str(zpo)
     if repeat2:
         title = ' zf=' + zfs[0:5] + ' z_sp=' + zsps[0:5] + ' ratio=' +rts[0:5] + 'slope' + sls[0:5] + ' i=' + str(i) + ' qid=' + str(id[i]) + ' chi=', schi[0:7] + ' chi_sp=', schi2[0:7]
     else:
         title = ' zf=' + zfs[0:5] + ' z_sp=' + zsps[0:5] + ' ratio=' +rts[0:5] + 'slope' + sls[0:5] + ' i=' + str(i) + ' qid=' + str(id[i]) + ' chi=', schi[0:7]
     ax.set_title(title)
########
     plt.show()
     if z_selection=='HIGHZ':
         filename="./PhotZ/HIGHZe{}/QID_{}.png".format(str(emagmin).replace('.',''),id[i])
     elif z_selection=='LOWZ':
         filename="./PhotZ/LOWZe{}/QID_{}.png".format(str(emagmin).replace('.',''),id[i])
     elif z_selection=='ALL':
         if z_free and leastsq:
             filename="./PhotZ/ALLZe{}_zfree-norm/QID_{}.png".format(str(emagmin).replace('.',''),id[i])
         elif z_free and leastsq==False:
             filename="./PhotZ/ALLZe{}_zfree-norm_nols/QID_{}.png".format(str(emagmin).replace('.',''),id[i]) 
         else:
             filename="./PhotZ/ALLZe{}_zfixed-norm/QID_{}.png".format(str(emagmin).replace('.',''),id[i])
     elif z_selection=='TEST':
         filename="./PhotZ/TESTe{}/QID_{}.png".format(str(emagmin).replace('.',''),id[i])

     fig.savefig(filename)
     plt.clf()
     plt.close()
##
#  lfm = [0.] * len(lf)
# lfm = phot_from_sed(wo, stdev, blo, bup, wsed, lfsed, z1fit[i], slope, norm, IGM_trans):
#
# print('-------------------------------')
# print('Parameter    Value       Stderr')
# for name, param in out.params.items():
#    print('{:7s} {:11.5f} {:11.5f}'.format(name, param.value, param.stderr))
#
z_fit = np.array(z1fit)-1.
col_a = Column(name='z_fit', data=z_fit)
t2.add_column(col_a)
col_b = Column(name='ratio_fit', data=r_fit)
t2.add_column(col_b)
col_c = Column(name='norm_fit', data=n_fit)
t2.add_column(col_c)
col_d = Column(name='slope_fit', data=s_fit)
t2.add_column(col_d)
col_e = Column(name='chisq', data=chi)
t2.add_column(col_e)
if leastsq: #check whether the leastsq fit was performed
    col_a_err = Column(name='zfit_err', data=np.array(zfit_err,dtype='float'))
    t.add_column(col_a_err)
    col_b_err = Column(name='ratiofit_err', data=np.array(rfit_err,dtype='float'))
    t.add_column(col_b_err)
    col_c_err = Column(name='normfit_err', data=np.array(nfit_err,dtype='float'))
    t.add_column(col_c_err)
    col_d_err = Column(name='slopefit_err', data=np.array(sfit_err,dtype='float'))
    t.add_column(col_d_err)
    
    if z_free:
        t.add_column(Column(name='corr_zr', data=corr_zr))
        t.add_column(Column(name='corr_zn', data=corr_zn))
        t.add_column(Column(name='corr_zs', data=corr_zs))
    t.add_column(Column(name='corr_rn', data=corr_rn))
    t.add_column(Column(name='corr_rs', data=corr_rs))
    t.add_column(Column(name='corr_sn', data=corr_sn))

if repeat2:  #check whether the fit was repeated with fixed z
    col_f = Column(name='ratio2_fit', data=r2_fit)
    t2.add_column(col_f)
    col_g = Column(name='norm2_fit', data=n2_fit)
    t2.add_column(col_g)
    col_h = Column(name='slope2_fit', data=s2_fit)
    t2.add_column(col_h)
    col_i = Column(name='chisq2', data=chi2)
    t2.add_column(col_i)
if z_selection=='HIGHZ':
    t2.write('./PhotZ/QSO_PhotZ_HIGHZ.fits', overwrite=True)
elif z_selection=='LOWZ':
    t2.write('./PhotZ/QSO_PhotZ_LOWZ.fits', overwrite=True)
elif z_selection=='ALL':
    if z_free and leastsq:
        t2.write('./PhotZ/QSO_PhotZ_out_e{}_zfree-norm.fits'.format(str(emagmin).replace('.','')), overwrite=True)
    if z_free and leastsq==False:
        t2.write('./PhotZ/QSO_PhotZ_out_e{}_zfree-norm_nols.fits'.format(str(emagmin).replace('.','')), overwrite=True)
    else:
        t2.write('./PhotZ/QSO_PhotZ_out_e{}_zfixed-norm.fits'.format(str(emagmin).replace('.','')), overwrite=True)
elif z_selection=='TEST':
    t2.write('./PhotZ/QSOreduv_out.fits', overwrite=True)    
####   FINITO
