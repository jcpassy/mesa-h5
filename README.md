MESA_H5
=======

This is the new implementation of [HDF5](https://support.hdfgroup.org/HDF5/) support for MESA
based on the former ``SE`` format. It provides ``mppnp`` compatibility. EK: test

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

Profiles in the HDF5 file
-------------------------

The user must specify which MESA profiles should be included in the HDF5. This is done in a file (default name ``hdf5_profile_columns.list`` similar to the ``profile_columns.list`` MESA file.

Note that the profiles in the HDF5 file **MUST** be a subset of the profiles written out in the MESA ASCII files. So whatever contains ``hdf5_profile_columns.list`` should also be in ``profile_columns.list``.

Customization
-------------

The ``mesa_hdf5_params.inc`` file contains the different setup parameters of the HDF5 outputs. You can modify them as you wish, but remember to clean and recompile for your changes to be applied. For example, for productions runs you will most likely set
```
       integer, parameter :: hdf5_num_mod_output
```
to a value of 100 or more depending on how long you will expect your run to be. Typically you would have 20-30 files per run. Putting 1000 cycles or time steps into one file is a practical number for longer runs.

Testing
-------

If you would like to run an example, you can use the inlist and columns files located in the ``test`` folder.

Known issues and tips
---------------------

You can use the `nugridpy.nugridse` module to work with the mesa_h5 hdf5 output files. Here are some notes:

* if the simulation has been stopped the last file may be corrupted and you may get some error when trying to extract or access data in that file. Just delete that last corrputed file.
* it is advisable to include `photo_interval = 100` in your inlist controls and add use in `src/mesa_hdf5_params.inc` the same value for `integer, parameter :: hdf5_num_mod_output = 100`

Authors
-------

[Jean-Claude Passy](https://github.com/jcpassy)

[Falk Herwig](https://github.com/fherwig)

Copyright
---------

Copyright (c) 2014 - 2018, NuGrid Team. All rights reserved.

License
-------

BSD-3-Clause "New" or "Revised" License (see LICENSE.md)