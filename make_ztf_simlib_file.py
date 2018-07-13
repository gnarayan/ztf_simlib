#!/usr/bin/env python
import sys
import os
import numpy as np
import sqlite3
import pandas as pd
import astropy.table as at
from astroquery.irsa_dust import IrsaDust
import astropy.coordinates as coord
import astropy.units as u

def main():
    db = sqlite3.connect('test_schedule_v8_msip.db')
    table = pd.read_sql_query("SELECT * from SUMMARY", db)
    ind = table['subprogram'] == 'all_sky'
    msip = table[ind]

    release_date = '20180622'
    survey = 'ZTF_MSIP'
    filters = ''.join(np.unique(msip['filter']))
    user = 'gnarayan'
    host = 'grimnir.stsci.edu'
    comment = 'Based on ZTF observing log DB from Eric Bellm, Rahul Biswas on {}'.format(release_date)
    pixsize = 1.
    fields = np.unique(msip['fieldID'])
    nlibid = len(fields)

    outlines = []
    outlines.append('SURVEY: {}'.format(survey))
    outlines.append('FILTERS: {}'.format(filters))
    outlines.append('TELESCOPE: ZTF')
    outlines.append('USER: {}'.format(user))
    outlines.append('HOST: {}'.format(host))
    outlines.append('SKYSIG_UNIT: ADU_PER_SQARCSEC')
    outlines.append('PIXSIZE: {:0.1f}'.format(pixsize))
    outlines.append('NLIBID: {}'.format(nlibid))
    outlines.append('COMMENT: {}'.format(comment))
    outlines.append('BEGIN LIBGEN')

    for field in fields:
        outlines.append('# --------------------------------------------')
        # select from table, not MSIP in case some of the other programs
        # observe the same field this may not be useful since we don't have
        # access to non-MSIP data but in principle these observations have been
        # taken and could be used to classify the data

        outlines.append('LIBID: {}'.format(field))
        indf = (table['fieldID'] == field)

        # all the positions appear to be identical, so there's no way to
        # account for dithers or overlaps
        ra  = np.unique(table[indf]['fieldRA'])[0]
        dec = np.unique(table[indf]['fieldDec'])[0]

        coo = coord.SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
        dust = IrsaDust.get_query_table(coo, section='ebv')
        mwebv = dust['ext SandF mean'][0]


        nobs = len(table[indf])
        outlines.append('RA: {}    DEC: {}    NOBS: {}    PIXSIZE: {}    MWEBV: {}    FIELD: {}'.format(ra, dec, nobs, pixsize, mwebv, field))
        outlines.append('#                           CCD  CCD         PSF1 PSF2 PSF2/1')
        outlines.append('#     MJD      ID*NEXPOSE  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG')

        entries = at.Table.from_pandas(table[indf])

        for entry in entries:
            # get some quantities
            flt = entry['filter']
            skymag = entry['filtSkyBright']
            depth  = entry['fiveSigmaDepth']
            snr = 5.
            fwhm   = entry['FWHMeff']

            term1 = 2.0 * depth - skymag
            term2 = - (depth - skymag)

            # convert FWHM from arcsec to sigma_gaussian in pixels
            sigma_pixel = fwhm /2.35 /pixsize
            pixel_area = area = (1.51 * fwhm)**2
            arg = pixel_area * snr * snr

            # Background dominated limit assuming counts with system transmission only
            # is approximately equal to counts with total transmission
            zpt_approx = term1 + 2.5 * np.log10(arg)

            tmp = 10. **(-0.4 * term2)
            zpt_cor = 2.5 * np.log10(1. + 1. / (pixel_area * tmp))
            simlib_zptavg = zpt_approx + zpt_cor

            npix_asec = 1. / pixsize**2.
            skysig = np.sqrt((1.0 / npix_asec) * 10.**(-0.4 * (skymag - simlib_zptavg)))

            lst = ['S:',
                   "{0:5.4f}".format(entry['expMJD']),
                   "{0:10d}*2".format(entry['obsHistID']),
                   entry['filter'],
                   "{0:5.2f}".format(1.),                  # CCD Gain
                   "{0:5.2f}".format(0.25),                # CCD Noise
                   "{0:6.2f}".format(skysig),              # SKYSIG
                   "{0:4.2f}".format(sigma_pixel),         # PSF1
                   "{0:4.2f}".format(0.),                  # PSF2
                   "{0:4.3f}".format(0.),                  # PSFRatio
                   "{0:6.2f}".format(simlib_zptavg),       # ZPTAVG
                   "{0:6.3f}".format(0.005),               # ZPTNoise
                   "{0:+7.3f}".format(-99.)]               # MAG
            out = ' '.join(lst)
            outlines.append(out)
        outlines.append('END_LIBID: {}'.format(field))
    outlines = '\n'.join(outlines)

    with open('ztf_msip_simlib_{}.dat'.format(release_date), 'w') as f:
        f.write(outlines)






if __name__=='__main__':
    sys.exit(main())
