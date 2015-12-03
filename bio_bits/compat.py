try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

try:
    from BytesIO import BytesIO
except ImportError:
    from io import BytesIO

from future.builtins import map, filter, zip

try:
    import unittest2 as unittest
except ImportError:
    import unittest


try:
    from functools import reduce
except:
    pass

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

try:
    from __builtin__ import open
except ImportError:
    from builtins import open

# Tests directory
from os.path import dirname
THIS = dirname(__file__)
