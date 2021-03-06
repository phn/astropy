===============================================
Building, Cython/C Extensions, and Distribution
===============================================

The build process currently uses the
`Distribute <http://packages.python.org/distribute/>`_ package to build and
install the astropy core (and any affiliated packages that use the template).
The user doesn't necessarily need to have `distribute` installed, as it will
automatically bootstrap itself using the ``distribute_setup.py`` file in the
source distribution if it isn't installed for the user.

Customizing setup/build for subpackages
---------------------------------------

As is typical, there is a single ``setup.py`` file that is used for the whole
`astropy` package.  To customize setup parameters for a given sub-package, a
``setup_package.py`` file can be defined inside a package, and if it is present,
the setup process will look for the following functions to customize the build
process:

* :func:`get_package_data`
    This function, if defined, should return a dictionary mapping the name of
    the subpackage(s) that need package data to a list of data file paths
    (possibly including wildcards) relative to the path of the package's source
    code.  e.g. if the source distribution has a needed data file
    ``astropy/wcs/tests/data/3d_cd.hdr``, this function should return
    ``{'astropy.wcs.tests:'['data/3d_cd.hdr']}``. See the ``package_data``
    option of the  :func:`distutils.core.setup` function.

    It is recommended that all such data be in a directory named "data" inside
    the package within which it is supposed to be used, and package data should
    be accessed via the `astropy.config.data.get_data_filename` and
    `astropy.config.data.get_data_fileobj` functions.

* :func:`get_extensions`
    This provides information for building C or Cython extensions. If defined,
    it should return a list of `distutils.core.Extension` objects controlling
    the Cython/C build process (see below for more detail).

* :func:`get_legacy_alias`
    This function allows for the creation of `shims` that allow a
    subpackage to be imported under another name.  For example,
    `astropy.io.fits` used to be available under the namespace
    `pyfits`.  For backward compatibility, it is helpful to have it
    still importable under the old name.  Under most circumstances,
    this function should call `astropy.setup_helpers.add_legacy_alias`
    to generate a legacy module and then return what it returns.

The `astropy.setup_helpers` modules includes a :func:`update_package_files`
function which automatically searches the given source path for
``setup_package.py`` modules and calls each of the above functions, if they
exist.  This makes it easy for affiliated packages to use this machinery in
their own ``setup.py``.

.. _building-c-or-cython-extensions:

C or Cython Extensions
----------------------

Astropy supports using C extensions for wrapping C libraries and Cython for
speeding up computationally-intensive calculations. Both Cython and C extension
building can be customized using the :func:`get_extensions` function of the
``setup_package.py`` file. If defined, this function must return a list of
`distutils.core.Extension` objects. The creation process is left to the
subpackage designer, and can be customized however is relevant for the
extensions in the subpackage.

While C extensions must always be defined through the :func:`get_extensions`
mechanism, Cython files (ending in ``.pyx``) are automatically located and
loaded in separate extensions if they are not in :func:`get_extensions`. For
Cython extensions located in this way, headers for numpy C functions are
included in the build, but no other external headers are included. ``.pyx``
files present in the extensions returned by :func:`get_extensions` are not
included in the list of extensions automatically generated extensions. Note
that this allows disabling a Cython file by providing an extension that
includes the Cython file, but giving it the special `name` 'cython_skip'. Any
extension with this package name will not be built by ``setup.py``.

.. note::

    If an :class:`~distutils.core.Extension` object is provided for Cython
    source files using the :func:`get_extensions` mechanism, it is very
    important that the ``.pyx`` files be given as the `source`, rather than the
    ``.c`` files generated by Cython.

Installing C header files
`````````````````````````

If your C extension needs to be linked from other third-party C code,
you probably want to install its header files along side the Python module.

    1) Create an `include` directory inside of your package for
       all of the header files.

    2) Use the :func:`get_package_data` hook in `setup_package.py` to
       install those header files.  For example, the `astropy.wcs`
       package has this::

           def get_package_data():
               return {'astropy.wcs': ['include/*.h']}

Preventing importing at build time
----------------------------------

In rare cases, some packages may need to be imported at build time.
Unfortunately, anything that requires a C or Cython extension or
processing through 2to3 will fail to import until the build phase has
completed.  In those cases, the `_ASTROPY_SETUP_` variable can be used
to determine if the package is being imported as part of the build and
choose to not import problematic modules.  `_ASTROPY_SETUP_` is
inserted into the builtins, and is `True` when inside of astropy's
`setup.py` script, and `False` otherwise.

For example, suppose there is a subpackage ``foo`` that needs to
import a module called ``version.py`` at build time in order to set
some version information, and also has a C extension, ``process``,
that will not be available in the source tree.  In this case,
``astropy/foo/__init__.py`` would probably want to check the value of
`_ASTROPY_SETUP_` before importing the C extension::

    if not _ASTROPY_SETUP_:
        from . import process

    from . import version

Release
-------

There have been no releases of the core package yet, so this isn't fully
defined. Some important items that are expected to be involved:

* The release process should be done by automated script.
* Documentation is currently hosted at http://readthedocs.org/docs/astropy,
  and tested at https://jenkins.shiningpanda.com/astropy/job/astropy-doc/,
  automatically updated with every commit to github master.
* C files generated by Cython should never be stored in source, but instead
  generated during the release process.  The setup script will always use these
  generated C files for release versions (instead of the .pyx Cython file).

There is a central `setup.py`.  It defines which Python packages to
install.  Each package does not have its own standalone `setup.py`.

Each package that needs to build C extensions has a module
`setup_package.py` that contains a function `get_extensions()` which
returns a list of `distutils.core.Extension` objects defining any
extensions to be built.

There are a set of helper functions for commonly occurring things when
building C extensions (e.g. finding the Numpy headers and library) in
`astropy.setup_helpers`.

Future directions
-----------------

We plan to switch to a newer packaging scheme when it's more stable, the
upcoming standard library `packaging` module, derived from the
`distutils2 <http://packages.python.org/Distutils2/library/distutils2.html>`_
project.  Until it's working right, however, we will be using `distribute` and
`distutils`.
