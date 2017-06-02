#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

import numpy as np
import scipy.interpolate as spi

class PMSTracks(object):
    """
    Docstring

    """
    def __init__(self, *args, **kwargs):
        """
        Docstring

        """
        self.tracks = kwargs.get('tracks', 'BHAC15')

        self.verbose = kwargs.get('verbose', False)

        if self.tracks == 'BHAC15':
            # self.reader = self.reader_BHAC15
            self.reader_BHAC15()
        else:
            raise ValueError("No valid reader method specified in pmstracks.")

    #
    # This method uses the scipy.interpolate.interp1d
    #   on ages for the closest mass tracks and then
    #   linearly interpolate the values between the two masses
    #   to get the result.
    #   This uses the mass tracks interpolation functions that
    #   I set up at the time of reading the table
    #   It has to be called specifying the mass and age for which we want the
    #   interpolation and using the correct dictionary label for the quantity
    #   that we want to get out.
    def interpolator_bilinear(self, mass, age, label, debug=False):
        #
        states = {0: 'No errors',
                  1: 'Problems with the tracks limits in mass',
                  2: 'Problems with the tracks limits in age',
                  3: 'Problems with the tracks limits in both age and mass'
                  }
        #
        # Find the two indices that bound the star
        im1, im2, m_status, ml_status = self._find_m1m2(mass)
        #
        # get the interpolated value
        int_value, i_status, il_status = self._get_intval(im1, im2, mass, age, label)
        #
        status = 0
        if m_status > 0:
            status += 1
        elif i_status > 0:
            status += 2
        #
        if debug:
            print('mass={0} mass[{1}]={2} mass[{3}]={4}'.format(mass,im1,self.mass[im1],im2,self.mass[im2]))
            print('    returned m_status={0} {1}'.format(m_status,ml_status))
            print('age={0} label={1} interpolated_value:{2}'.format(age,label,int_value))
            print('    returned i_status={0} {1}'.format(i_status, il_status))

        return int_value, status, states[status]

    #
    # This method returns the the interpolated value and a status
    #    for the requested age, given the two closest mass tracks.
    #    return the following states also:
    #       0: ok
    #       1: age is below the minimum age of one of the tracks
    #       2: age is above the maximum age of one of the tracks
    def _get_intval(self, im1, im2, mass, age, label):
        #
        states = { 0 : 'No errors',
                   1 : 'Age for lower mass below tracks limits', 2 : 'Age for lower mass above tracks limits',
                   3 : 'Age for higher mass below tracks limits', 6 : 'Age for higher mass above tracks limits',
                   4 : 'Age for both masses below tracks limits', 8 : 'Age for both masses above tracks limits',
                   7 : 'Unexpected age out of limits for both tracks', 5 : 'Unexpected age out of limits for both tracks'
                   }
        #
        int1_value, i1_status, l1_status = self._my_lint(im1, label, age)
        int2_value, i2_status, l2_status = self._my_lint(im2, label, age)
        status = i1_status+3*i2_status
        #
        if self.mass[im1] <= mass:
            int_value = int1_value
        elif self.mass[im2] >= mass:
            int_value = int2_value
        else:
            dm = mass - self.mass[im1]
            ddm = self.mass[im2] - self.mass[im1]
            ddy = int2_value - int1_value
            int_value = ddy/ddm*dm + int1_value
        #
        return int_value, status, states[status]

    #
    # used by _get_intval to extract the correct interpolation value and
    #    return an error or warning otherwise
    def _my_lint(self, im, label, age):
        #
        intlabel = label+'_int'
        #
        states = { 0 : 'No errors', 1 : 'Requested age for this mass below tracks limits', 2 : 'Requested age for this mass above tracks limits'}
        status = 0
        if age < ((self.tracks[im])['lage'])[0]:
            intval = ((self.tracks[im])[label])[0]
            status = 1
        elif age > ((self.tracks[im])['lage'])[-1]:
            intval = ((self.tracks[im])[label])[-1]
            status = 2
        else:
            intval = ((self.tracks[im])[intlabel])(age)
        #
        return intval, status, states[status]

    #
    # This method returns the two closest masses indices in the tracks
    #    for the requested mass. returns the status and a label
    #       0: ok
    #       1: mass is below the minimum mass of the tracks
    #       2: mass is above the maximum mass of the tracks
    def _find_m1m2(self, m):
        #
        states = { 0 : 'No errors', 1 : 'Requested mass below tracks limits', 2 : 'Requested mass above tracks limits'}
        status = 0
        if m < self.mass[0]:
            imin = 0
            imax = 0
            status = 1
        elif m > self.mass[-1]:
            imin = len(self.mass)-1
            imax = len(self.mass)-1
            status = 2
        else:
            imin = 0
            imax = len(self.mass)-1
            while imax - imin > 1:
                if m == self.mass[imin]:
                    imax = imin + 1
                elif m == self.mass[imax]:
                    imin = imax - 1
                else:
                    itry = imin+int((imax - imin)/2)
                    if m < self.mass[itry]:
                        imax = itry
                    else:
                        imin = itry
        #
        return imin, imax, status, states[status]

    #
    # This method uses the scipy.interpolate.interp2d
    #   to create the interpolation function, it has to be
    #   called using the correct dictionary label for the quantity
    #   that we want to get out.
    def interpolator2d(self, label):
        #
        return spi.interp2d(self.all_tracks['lmass'], self.all_tracks['lage'], self.all_tracks[label])

    #
    # This function reads the BHAC15 Evolutionary tracks
    def reader_BHAC15(self):
        self.infile_models = '/Users/ltesti/git_repository/diskpop/pmstracks/tracks/'+self.tracks+'/BHAC15_tracks+structure'
        #self.infile_models = './pmstracks/tracks/'+self.tracks+'/BHAC15_tracks+structure'
        if self.verbose:
            print("Reading file: {}".format(self.infile_models))

        doread = False
        mstar = []
        self.tracks = []
        age = []
        lum = []
        teff = []
        all_mass = []
        all_age = []
        all_lum = []
        all_teff = []
        newmass = True
        f = open(self.infile_models, 'r')
        for line in f.readlines():
            if (doread):
                if line[0] == '!':
                    pass
                elif line[0] == '\n':
                    if newmass != True:
                        self.tracks.append({'model_mass':mstar[-1], 'mass': mstar[-1]*np.ones(len(age)),
                                            'nage':len(age), 'lage':np.array(age),
                                            'llum':np.array(lum), 'llum_int': spi.interp1d(np.array(age),np.array(lum)),
                                            'teff':np.array(teff), 'teff_int': spi.interp1d(np.array(age),np.array(teff))})
                    newmass = True
                    age = []
                    lum = []
                    teff = []
                else:
                    columns = line.split()
                    if newmass:
                        mstar.append(float(columns[0]))
                        newmass = False
                    age.append(float(columns[1]))
                    lum.append(float(columns[3]))
                    teff.append(float(columns[2]))
                    all_mass.append(float(columns[0]))
                    all_age.append(float(columns[1]))
                    all_lum.append(float(columns[3]))
                    all_teff.append(float(columns[2]))

            else:
                if line[0] == '!':
                    doread = True
                    #nex += 1

        f.close()
        self.mass = np.array(mstar)
        #print('Mass = {}'.format(mstar))
        self.all_tracks = {'mass': np.array(all_mass), 'lage': np.array(all_age),
                           'llum': np.array(all_lum), 'teff': np.array(all_teff)}



