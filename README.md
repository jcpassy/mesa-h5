The HDF5 framework for MESA
===========================

This is the new implementation of [HDF5](https://support.hdfgroup.org/HDF5/) support for MESA
based on the former ``SE`` format. It provides ``mppnp`` compatibility.


Requirements
------------

You must have successfully installed a release of [MESA](http://mesa.sourceforge.net/).

Latest version of the MESA code on which this implementation has been tested: **10398**.

Installation
------------

Once MESA has been installed:
* copy the ``run_star_extras.f`` and ``mesa_hdf5_*.inc`` files located in the ``src`` folder
into the ``src`` folder in your work directory
* clean, compile, and enjoy!

Customization
-------------

The ``mesa_hdf5_params.inc`` file contains the different setup parameters of the HDF5 outputs.
You can modify them as you wish, but remember to clean and recompile for your changes to be applied.

Testing
-------

If you would like to run an example, you can use the inlist and columns files located in the ``test`` folder.

Authors
-------

[Jean-Claude Passy](https://github.com/jcpassy)

Falk Herwig

Samuel Jones

Copyright
---------

Copyright (c) 2014 - 2018, NuGrid Team. All rights reserved.

License
-------

BSD-3-Clause