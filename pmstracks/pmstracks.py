#!/usr/bin/env python
# -*- coding: utf-8 -*-



class PMSTracks(object):
    """
    Docstring

    """
    def __init__(self, *args, **kwargs):
        """
        Docstring

        """
        self.tracks = kwargs.get('tracks', '')
        self.temp = kwargs.get('temp', 0.)
        self.evolve_method = kwargs.get('evolve_method', 0)

        if self.evolve_method == 'constant':
            self.evolve = self.constantStar
        else:
            raise ValueError("No valid evolve method specified in Star parameters.")

    def constantStar(self):
        pass


