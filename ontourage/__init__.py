from . import _version
_versioneer_info = _version.get_versions()
if _versioneer_info['error']:
    raise RuntimeError(f'Could not create version string: {_versioneer_info["error"]}')
__version__ = _versioneer_info['version']

__author__ = "Peter Ebert"
__copyright__ = "Copyright 2021, Peter Ebert"
__license__ = "MIT"
__maintainer__ = "Peter Ebert"
__email__ = "peter.ebert@iscb.org"
__status__ = "Prototype"
