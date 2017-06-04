#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

import numpy as np
import scipy.interpolate as spi
import os
import glob
from . import TRACKS_DIR

class PMSTracks(object):
    """
    Docstring

    """
    def __init__(self, tracks='BHAC15', verbose=False):
        """
        Docstring

        """
        self.tracks_name = tracks

        self.verbose = verbose

        self.tracks_path = os.path.join(TRACKS_DIR, self.tracks_name)

        if self.tracks_name == 'BHAC15':
            self.infile_models = os.path.join(self.tracks_path, 'BHAC15_tracks+structure')
            self.reader = self.reader_BHAC15
        elif self.tracks_name == 'Siess00':
            self.infile_models = glob.glob(os.path.join(self.tracks_path, "*.hrd"))
            if self.verbose:
                print('{}'.format(self.infile_models))
            self.reader = self.reader_Siess00

        else:
            raise ValueError("No valid reader method specified in pmstracks.")
        #
        self.mass, self.tracks = self.reader()
        self.mass, self.tracks = self._sort_tracks()
        self.interp_age = self._tracks_age_interp()

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
        if self.mass[im1] >= mass:
            int_value = int1_value
        elif self.mass[im2] <= mass:
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
            intval = ((self.interp_age[im])[intlabel])(age)
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
    # This function is the reader for the Siess00 Evolutionary tracks
    def reader_Siess00(self):
        """
        Reader for Siess et al. (2000) track files

        The track files are downloaded from the Siess server and not edited
        The location of the tracks are defined by the attribute self.infile_models,
        which, in this case, is a list contaning all the files that need to be read.

        Parameters
        ----------
        None

        Returns
        -------
        mass, tracks

        mass: a numpy array containing the mass of each track

        tracks: a list of dictionaries containing all the tracks data
                'model_mass': mass for the track
                'mass': numpy array containing the same mass for each timestep (redundant, could be removed)
                'nage': number of time steps in the track (redundant, could be removed)
                'lage': numpy array of the time steps for this track
                'llum': numpy array of the log10(L/Lsun) for this track
                'teff': numpy array of the effective temperatures for this track

        Examples
        --------
        after defining self.infile_models so that is points to the file containing the tracks,
        running:

        self.mass, self.tracks = self.reader_BHAC15()

        will read the file and fill in the self.mass and self.tracks attributes

        """
        #
        # Define the file to read
        #
        mstar = []
        tracks = []
        for i_f in self.infile_models:
            age = []
            lum = []
            teff = []
            isfirst = True
            if self.verbose:
                print("Reading file: {}".format(i_f))
            f = open(i_f, 'r')
            for line in f.readlines():
                if line[0] != '#':
                    columns = line.split()
                    if isfirst:
                        mstar.append(float(columns[9]))
                        isfirst = False
                    age.append(np.log10(float(columns[10])))
                    lum.append(np.log10(float(columns[2])))
                    teff.append(float(columns[6]))
            f.close()
            tracks.append({'model_mass': mstar[-1], 'mass': mstar[-1] * np.ones(len(age)),
                           'nage': len(age), 'lage': np.array(age),
                           'llum': np.array(lum), 'teff': np.array(teff)})
            #
        return np.array(mstar), tracks

    #
    # This function is the reader for the BHAC15 Evolutionary tracks
    def reader_BHAC15(self):
        """
        Reader for Baraffe et al. (2015) track files

        The track files are downloaded from the Baraffe server and not edited
        The location of the tracks are defined by the attribute self.infile_models

        Parameters
        ----------
        None

        Returns
        -------
        mass, tracks

        mass: a numpy array containing the mass of each track

        tracks: a list of dictionaries containing all the tracks data
                'model_mass': mass for the track
                'mass': numpy array containing the same mass for each timestep (redundant, could be removed)
                'nage': number of time steps in the track (redundant, could be removed)
                'lage': numpy array of the time steps for this track
                'llum': numpy array of the log10(L/Lsun) for this track
                'teff': numpy array of the effective temperatures for this track

        Examples
        --------
        after defining self.infile_models so that is points to the file containing the tracks,
        running:

        self.mass, self.tracks = self.reader_BHAC15()

        will read the file and fill in the self.mass and self.tracks attributes

        """
        #
        # Define the file to read
        #
        if self.verbose:
            print("Reading file: {}".format(self.infile_models))
        #
        doread = False
        mstar = []
        tracks = []
        age = []
        lum = []
        teff = []
        newmass = True
        f = open(self.infile_models, 'r')
        dowrite = False
        for line in f.readlines():
            if (doread):
                if (line[0] == '!'):
                    if dowrite and (newmass != True):
                        tracks.append({'model_mass': mstar[-1], 'mass': mstar[-1] * np.ones(len(age)),
                                            'nage': len(age), 'lage': np.array(age),
                                            'llum': np.array(lum), 'teff': np.array(teff) })
                        dowrite = False
                        newmass = True
                        age = []
                        lum = []
                        teff = []
                    else:
                        pass
                elif line[0] == '\n':
                    pass
                else:
                    dowrite = True
                    columns = line.split()
                    if newmass:
                        mstar.append(float(columns[0]))
                        newmass = False
                    age.append(float(columns[1]))
                    lum.append(float(columns[3]))
                    teff.append(float(columns[2]))

            else:
                if line[0] == '!':
                    doread = True
        #
        f.close()
        #
        return np.array(mstar), tracks

    #
    # This method sorts the tracks in increasing mass and per age for each mass
    def _sort_tracks(self):
        """
        Used to sort the tracks by mass and then each track by age.

        Parameters
        ----------
        self.mass : numpy array
            Array with the values of the masses for each track.
        self.tracks : list of track dictionaries
            Each element contains a dictionary with the track data.

        Returns
        -------
        sort_mass : numpy array
            copy of self.mass sorted by increasing mass
        sort_tracks : list of track disctionaries
            copy of self.tracks sorted by mass and with the tracks resoted by increasing age

        Examples
        --------
        This method can be called after a track reader has filled self.mass and self.tracks:
        self.mass, self.tracks = self._sort_tracks()

        this will resort in place self.mass and self.tracks

        """
        msort = np.argsort(self.mass)
        sort_mass = self.mass[msort]
        sort_tracks = []
        for im in range(len(self.mass)):
            sort_tracks.append(self.tracks[msort[im]])
            isort = np.argsort((sort_tracks[im])['lage'])
            ((sort_tracks[im])['lage'])[:] = ((self.tracks[msort[im]])['lage'])[isort]
            ((sort_tracks[im])['llum'])[:] = ((self.tracks[msort[im]])['llum'])[isort]
            ((sort_tracks[im])['teff'])[:] = ((self.tracks[msort[im]])['teff'])[isort]

        return sort_mass, sort_tracks

    #
    # this method sets up the age interpolators
    def _tracks_age_interp(self):
        """
        Used to compute the interpolation functions for llum and teff as a function of age.
        Uses the scipy.interpolate.interp1d implementation of a linear interpolation, edges
        probelms need to be checked and cured separately.

        Parameters
        ----------
        self.mass : numpy array
            Array with the values of the masses for each track.
        self.tracks : list of track dictionaries
            Each element contains a dictionary with the track data.

        Returns
        -------
        interp_age : list of track l,t interpolation functions
            assumes that self.tracks have been sorted with _sort_tracks()

        Examples
        --------
        This method can be called after a track reader has filled self.mass and self.tracks, and
        _sort_tracks() has been used to resort them in place:

        self.interp_age = self._tracks_age_interp()

        self.interp_age is the list of dictionaries containing the llum and teff interpolators

        """
        #
        interp_age = []
        for im in range(len(self.mass)):
            interp_age.append({'llum_int': spi.interp1d((self.tracks[im])['lage'], (self.tracks[im])['llum']),
                               'teff_int': spi.interp1d((self.tracks[im])['lage'], (self.tracks[im])['teff'])})
        return interp_age


