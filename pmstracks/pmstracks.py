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
    This class defines the PMS tracks object and the methods to interpolate
    the tracks and obtain the values of the stellar parameters

    When the object is instantiated, it creates the tracks data by reading the
    specified set of tracks and creates the methods to interpolate the tracks

    the object can later be interrogated to return the interpolated value
    of the tracks.

    The track options are defined in __init__()

    The interpolation method is implemented in interpolator_bilinear()

    """
    def __init__(self, tracks='BHAC15', verbose=False):
        """
        When instantiated, the object creates a set of pms tracks reading the
        appropriate track files and the interpolators used to manipulate the
        tracks and extract the stellar properties for a given pair of mass and age.

        At the moment we have implemented readers for the following tracks:
         BHAC15  : Baraffe et al. 2015 tracks
         Siess00 : Siess et al. 2000 tracks

        Parameters
        ----------

        tracks : string

        string code for the tracks to be read, default is Baraffe et al. 2015 tracks

        verbose : boolean

        print out status messages, default is False

        """
        self.tracks_name = tracks

        self.verbose = verbose

        self.tracks_path = os.path.join(TRACKS_DIR, self.tracks_name)

        if self.tracks_name == 'BHAC15':
            self.infile_models = os.path.join(self.tracks_path, 'BHAC15_tracks.dat')
            self.reader = self.reader_bhac15
        elif self.tracks_name == 'Siess00':
            self.infile_models = glob.glob(os.path.join(self.tracks_path, "*.hrd"))
            if self.verbose:
                print('{}'.format(self.infile_models))
            self.reader = self.reader_siess00
        elif self.tracks_name == 'F16_std':
            self.infile_models = glob.glob(os.path.join(self.tracks_path, "*.trk"))
            if self.verbose:
                print('{}'.format(self.infile_models))
            self.reader = self.reader_feiden16_std
        elif self.tracks_name == 'F16_mag':
            self.infile_models = glob.glob(os.path.join(self.tracks_path, "*.ntrk"))
            if self.verbose:
                print('{}'.format(self.infile_models))
            self.reader = self.reader_feiden16_mag
        else:
            raise ValueError("No valid reader method specified in pmstracks.")
        #
        self.mass, self.tracks = self.reader()
        self.mass, self.tracks = self._sort_tracks()
        self.interp_age = self._tracks_age_interp()

    def __repr__(self):
        str = '{}'.format(self.tracks_name)
        return str

    def _pmsname(self):
        str = '{}'.format(self.tracks_name)
        return str
    
    #
    # This method determines the two closest tracks from a given (L, T) point
    #   the idea is to use this to identify the two closest tracks and then interpolate
    #   to get the stellar mass and age from tracks.
    def _iso_dist(lstar,tstar,mytrack):
        '''
        this utility function returns the minimum distance from a track, the corresponding signed
        luminosity and teff difference, and the indices of the two nearest point
        '''
        dl = (lstar-mytrack['llum'])
        dt = (np.log10(tstar)-np.log10(mytrack['teff']))
        dd = np.sqrt(dl**2.+dt**2.)
        ns = np.argsort(dd)
        return dd[ns[0]], dl[ns[0]], dt[ns[0]], ns[0], ns[1]
        
    def _get_trk_dist(self,lstar,tstar):
        '''
        this utility function returns the filled lists with minimum distance,
        luminosity and teff difference, and the indices of the two nearest point
        '''
        dd = []
        dl = []
        dt = []
        n0 = []
        n1 = []
        for i in range(len(self.mass)):
            ddt, dlt, dtt, n0t, n1t = self._iso_dist(lstar,tstar,self.tracks[i])
            dd.append(ddt)
            dl.append(dlt)
            dt.append(dtt)
            n0.append(n0t)
            n1.append(n1t)
        return dd, dl, dt, n0, n1

    def two_iso(self,lstar,tstar):
        '''
        this utility function returns the minimum distance from a track, the corresponding signed
        luminosity and teff difference, and the indices of the two nearest point
        '''
        ddt, dlt, dtt, n0t, n1t = self._get_trk_dist(lstar,tstar)
        

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
        """
        return the interpolated value for the specified label parameter

        This is the method that should be called from top level, once the
        class has been properly initialized, to obtain the interpolated value of the
        "label" parameter for a given mass and age.

        Parameters
        ----------

        mass : float

        value of the mass for which we want to have the interpolation
        (expected units: Msun)

        age : float

        value of the age for which we want to have the interpolation
        (expected units: Log10(age/Myr)

        label : string

        dictionary label for the quantity to be interpolated

        debug : boolean

        prints status messages for debug purposes

        Returns
        -------
        value, status, message

        value: result of the interpolation
        status: 0 if the input age is within the tracks, 1 to 3 out of tracks edges
        message: info message connected to the status

        Examples
        --------
        This is the top level method that should be called to get the interpolatd
        values:

        val, code_status, label_status = interpolator_bilinear(mass, age, label, debug=False)

        to find the interpolated value at the specified mass and age
        At the moment I have implemented 'llum' and 'teff' as labels, note that I am converting
        all tracks to have 'mass' in solar masses and 'age' in Log10(age/Myr) so implicitly this
        is what I am assuming. Furthermore, typically I have the 'llum' as Log10(L/Lsun) and
        'teff' in K.

        """
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
            print('mass={0} mass[{1}]={2} mass[{3}]={4}'.format(mass, im1, self.mass[im1], im2, self.mass[im2]))
            print('    returned m_status={0} {1}'.format(m_status, ml_status))
            print('age={0} label={1} interpolated_value:{2}'.format(age, label, int_value))
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
        """
        return the interpolated value between two tracks for the specified label parameter

        This method calls self._my_lint to interpolate each track in age and then interpolates
        to find the mass. Appropriate states are returned to flag cases where the requested mass
        and/or age are outside the track boundaries. The returned value correspond to the boundary.

        Parameters
        ----------
        im1, 1m2 : long

        indices for the tracks to be used

        mass : float

        value of the mass for which we want to have the interpolation

        age : float

        value of the age for which we want to have the interpolation

        label : string

        dictionary label for the quantity to be interpolated

        Returns
        -------
        value, status, message

        value: result of the interpolation
        status: 0 if the input age is within the tracks, 1 to 8 out of tracks edges
        message: info message connected to the status

        Examples
        --------
        after the self.mass and self.tracks are created and ordered, and after self.interp_age
        attributes are created, self._find_m1m2() can be used to find the two tracks to be used in
        the interpolation and then, finally, this method can be called with:

        val, code_status, label_status = self._get_intval(im1, im2, mass, age, label)

        to find the interpolated value at the specified age for the specified track index
        if masses/ages are out of bounds, the edges are returned and a proper message is returned

        """
        #
        states = {0: 'No errors',
                  1: 'Age for lower mass below tracks limits', 2: 'Age for lower mass above tracks limits',
                  3: 'Age for higher mass below tracks limits', 6: 'Age for higher mass above tracks limits',
                  4: 'Age for both masses below tracks limits', 8: 'Age for both masses above tracks limits',
                  7: 'Unexpected age out of limits for both tracks', 5: 'Unexpected age out of limits for both tracks'
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
        """
        return the interpolation value along a specified track for the appropriate label (parameter)

        The method checks whether the age is beyond the track limits and
        returns an appropriate message.

        This function assumes that self.tracks are appropriately ordered by increasing age and
        that self.interp_age contain the appropriate interpolation functions

        Parameters
        ----------
        im : long

        index for the track to be used

        label : string

        dictionary label for the quantity to be interpolated

        age : float

        value of the age for which we want to have the interpolation

        Returns
        -------
        value, status, message

        value: result of the interpolation
        status: 0 if the input age is within the tracks, 1 if below, 2 if above
        message: info message connected to the status

        Examples
        --------
        after the self.tracks and self.interp_age attributes are created , this method
        can be called with:

        val, code_status, label_status = self._my_lint(im, label, age)

        to find the interpolated value at the specified age for the specified track index
        if ages are out of bounds, the edges are returned and a proper message is returned

        """
        #
        intlabel = label+'_int'
        #
        states = {0: 'No errors', 1: 'Requested age for this mass below tracks limits',
                  2: 'Requested age for this mass above tracks limits'}
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
        """
        Returns the masses of the tracks that bracket the input mass

        The method checks whether the mass is beyond the track limits and
        returns an appropriate message.

        This function assumes that self.mass is ordered by increasing mass

        Parameters
        ----------
        m : float

        mass that we want to bracket

        Returns
        -------
        imin, imax, status, message

        imin: index within self.mass for the lower bracket
        imax: index within self.mass for the upper bracket
        status: 0 if the input mass is within the tracks, 1 if below, 2 if above
        message: info message connected to the status

        Examples
        --------
        after the self.mass and self.tracks attributes are created and ordered, this method
        can be called with:

        idx1, idx2, code_status, label_status = self._find_m1m2(mass)

        to find self.mass[idx1] <= mass <= self.mass[idx2]
        if mass is below self.mass[0], then idx1 = 0 and idx2 = 1, and the code_status = 1
        if mass is above self.mass[-1], then idx1 = len(self.mass)-2 and idx2 = len(self.mass)-1,
        and the code_status = 2

        """
        #
        states = {0: 'No errors', 1: 'Requested mass below tracks limits',
                  2: 'Requested mass above tracks limits'}
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

    def reader_feiden16_std(self):
        """
        Reader for Feiden et al. (2016) standard track files

        The track files are downloaded from the github 
        https://github.com/gfeiden/MagneticUpperSco 
        server and not edited
        The location of the tracks are defined by the attribute self.infile_models,
        which, in this case, is a list contaning all the files that need to be read.

        Parameters
        ----------

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

        Note: the reader produces tracks with:
        mass in Msun
        lage in Log10(age/Myr)
        llum in Log10(L/Lsun)
        teff in K

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
            if self.verbose:
                print("Reading file: {}".format(i_f))
            f = open(i_f, 'r')
            for line in f.readlines():
                if line[0] != '#':
                    columns = line.split()
                    age.append(np.log10(float(columns[0])))
                    lum.append(float(columns[3]))
                    teff.append(10**float(columns[1]))
                else:
                    if len(line) > 8:
                        if line[2] == 'M':
                            mstar.append(float(line[4:8]))
            f.close()
            tracks.append({'model_mass': mstar[-1], 'mass': mstar[-1] * np.ones(len(age)),
                           'nage': len(age), 'lage': np.array(age),
                           'llum': np.array(lum), 'teff': np.array(teff)})
            #
            if self.verbose:
                print("    Model Mass: {0}  nages={1}  age_start={2}  age_end={3}".format(mstar[-1],len(age),age[0],age[-1]))
        return np.array(mstar), tracks

    def reader_feiden16_mag(self):
        """
        Reader for Feiden et al. (2016) magnetic track files

        The track files are downloaded from the github 
        https://github.com/gfeiden/MagneticUpperSco 
        server and not edited
        The location of the tracks are defined by the attribute self.infile_models,
        which, in this case, is a list contaning all the files that need to be read.

        Parameters
        ----------

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

        Note: the reader produces tracks with:
        mass in Msun
        lage in Log10(age/Myr)
        llum in Log10(L/Lsun)
        teff in K

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
            if self.verbose:
                print("Reading file: {}".format(i_f))
            f = open(i_f, 'r')
            for line in f.readlines():
                if len(line) > 23:
                    if line[0] != '#':
                        columns = line.split()
                        age.append(np.log10(float(columns[2]))+9.)
                        lum.append(float(columns[3]))
                        teff.append(10**float(columns[6]))
                    else:
                        if line[2:12] == 'Total mass':
                            mstar.append(float(line[16:23]))
            if self.verbose:
                print("    Model Mass: {0}  nages={1}  age_start={2}  age_end={3}".format(mstar[-1],len(age),age[0],age[-1]))
            f.close()
            tracks.append({'model_mass': mstar[-1], 'mass': mstar[-1] * np.ones(len(age)),
                           'nage': len(age), 'lage': np.array(age),
                           'llum': np.array(lum), 'teff': np.array(teff)})
            #
        return np.array(mstar), tracks

    #
    # This function is the reader for the Siess00 Evolutionary tracks
    def reader_siess00(self):
        """
        Reader for Siess et al. (2000) track files

        The track files are downloaded from the Siess server and not edited
        The location of the tracks are defined by the attribute self.infile_models,
        which, in this case, is a list contaning all the files that need to be read.

        Parameters
        ----------

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

        Note: the reader produces tracks with:
        mass in Msun
        lage in Log10(age/Myr)
        llum in Log10(L/Lsun)
        teff in K

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
    def reader_bhac15(self):
        """
        Reader for Baraffe et al. (2015) track files

        The track files are downloaded from the Baraffe server and not edited
        The location of the tracks are defined by the attribute self.infile_models

        Parameters
        ----------

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

        Note: the reader produces tracks with:
        mass in Msun
        lage in Log10(age/Myr)
        llum in Log10(L/Lsun)
        teff in K

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
            if doread:
                if line[0] == '!':
                    if dowrite and (not newmass):
                        tracks.append({'model_mass': mstar[-1], 'mass': mstar[-1] * np.ones(len(age)),
                                       'nage': len(age), 'lage': np.array(age),
                                       'llum': np.array(lum), 'teff': np.array(teff)})
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
        sorted_mass : numpy array
            copy of self.mass sorted by increasing mass
        sorted_tracks : list of track disctionaries
            copy of self.tracks sorted by mass and with the tracks resoted by increasing age

        Examples
        --------
        This method can be called after a track reader has filled self.mass and self.tracks:
        self.mass, self.tracks = self._sort_tracks()

        this will resort in place self.mass and self.tracks

        """
        msort = np.argsort(self.mass)
        sorted_mass = self.mass[msort]
        sorted_tracks = []
        for im in range(len(self.mass)):
            sorted_tracks.append(self.tracks[msort[im]])
            isort = np.argsort((sorted_tracks[im])['lage'])
            ((sorted_tracks[im])['lage'])[:] = ((self.tracks[msort[im]])['lage'])[isort]
            ((sorted_tracks[im])['llum'])[:] = ((self.tracks[msort[im]])['llum'])[isort]
            ((sorted_tracks[im])['teff'])[:] = ((self.tracks[msort[im]])['teff'])[isort]

        return sorted_mass, sorted_tracks

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

    def plot_tracks(self, ax, ages=None, masses=None):
        """
        Plot evolutionary tracks
        :param ax:
        :param ages:
        :param masses:
        :return:
        """
        # TODO: ages grid are hardcoded. Should they change?
        if not ages:
            ages = np.linspace(5.2, 9.7, 50)
        if not masses:
            masses = np.array([0.011, 0.014, 0.083, 0.12, 0.17,
                               0.25, 0.3, 0.35, 0.77, 0.95, 1.05, 1.32])

        nages = len(ages)
        nmasses = len(masses)

        interp_ls00 = np.zeros((nmasses, nages))
        interp_ts00 = np.zeros((nmasses, nages))
        for i in range(nmasses):
            for j in range(nages):
                interp_ls00[i, j] = self.interpolator_bilinear(masses[i], ages[j], 'llum')[0]
                interp_ts00[i, j] = self.interpolator_bilinear(masses[i], ages[j], 'teff')[0]

        # plot evolutionary tracks
        for i in range(len(self.tracks)):
            ax.plot((self.tracks[i])['teff'], (self.tracks[i])['llum'],
                    color='k', linestyle='solid')

        # plot iso-mass tracks
        for i in range(nmasses):
            ax.plot(interp_ts00[i, :], interp_ls00[i, :],
                    color='red', linestyle='dotted')
