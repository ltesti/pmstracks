from ._version import __version__

import os
import matplotlib

# To avoid Runtime Error
# RuntimeError: Python is not installed as a framework. The Mac OS X backend
# will not be able to function correctly if Python is not installed as a framework.
# See the Python documentation for more information on installing Python as a
# framework on Mac OS X. Please either reinstall Python as a framework, or try
# one of the other backends.
# If you are using (Ana)Conda please install python.app and replace the use of
# 'python' with 'pythonw'.
# See 'Working with Matplotlib on OSX' in the Matplotlib FAQ for more information.
# https://matplotlib.org/faq/osx_framework.html

if matplotlib.get_backend().lower() == 'macosx':
    matplotlib.use('TkAgg')

#print('{}'.format(__path__))
#TRACKS_DIR = os.path.join(os.path.split(__path__[0])[0],'tracks')
TRACKS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),'tracks')

print('tracks_dir: '+TRACKS_DIR)

from .pmstracks import PMSTracks

