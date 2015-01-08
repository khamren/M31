import numpy as np

def sigclip_coadd_arr(lbin, spec, ivar, **kwargs):
    sigLim = kwargs.get('sigLim',3.5)
    verbose = kwargs.get('verbose',False)

    nlbin = len(lbin[0])
    nstars = len(lbin)

    coaddivar = np.nansum(ivar, axis = 0)
    coaddspec = np.nansum(spec*ivar, axis = 0)/coaddivar

    clipivar = []
    subspec = []

    for j in range(nlbin):
        if coaddivar[j] == 0:
            clipivar.append(0)
            subspec.append(np.nan)
        else:
            s = spec[:,j]
            iv = ivar[:,j]
            nsig = np.abs((s - coaddspec[j])*np.sqrt(iv))
            subspec.append(np.nansum(s[nsig < sigLim]*iv[nsig < sigLim]))
            clipivar.append(np.nansum(iv[nsig < sigLim]))

    clipspec = np.array(subspec)/np.array(clipivar)
    clipivar = np.array(clipivar)

    if verbose:
        return coaddspec, coaddivar, clipspec, clipivar
    else:
        return clipspec, clipivar

def sigclip_coadd_fits(fits,**kwargs):
    """
    Optional keywords are:
        calib - access the flux-calibrated tags, or not. Default is True.
        sigLim - sigma limit for clipping. Default is 3.5
        verbose - if False (default), out put just clipped spectrum and ivar. If True, output basic coadds as well

    """

    sigLim = kwargs.get('sigLim',3.5)
    calib = kwargs.get('calib',True)
    verbose = kwargs.get('verbose',False)

    if calib:
        ivartag = 'IVARCALIB'
        spectag = 'SPECCALIB'
    else:
        ivartag = 'IVAR'
        spectag = 'SPEC' 

    nlbin = len(fits.LBIN[0])
    nstars = len(fits)

    coaddivar = np.nansum(fits.field(ivartag), axis = 0)
    coaddspec = np.nansum(fits.field(spectag)*fits.field(ivartag), axis = 0)/coaddivar

    spec = fits.field(spectag)
    ivar = fits.field(ivartag)
    
    clipivar = []
    subspec = []

    for j in range(nlbin):

        if coaddivar[j] == 0:
            clipivar.append(0)
            subspec.append(np.nan)
        else:
            s = spec[:,j]
            iv = ivar[:,j]
            nsig = np.abs((s - coaddspec[j])*np.sqrt(iv))
            subspec.append(np.nansum(s[nsig < sigLim]*iv[nsig < sigLim]))
            clipivar.append(np.nansum(iv[nsig < sigLim]))

    clipspec = np.array(subspec)/np.array(clipivar)

    if verbose:
        return coaddspec, coaddivar, clipspec, clipivar
    else:
        return clipspec, clipivar
