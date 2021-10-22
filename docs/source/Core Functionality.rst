Core Functionality
==================

These are the files in the L5NeuronSimulation folder that are used throughout the other folders.

.. toctree::
    :maxdepth: 1
    :hidden:
  
    clustering
    raster_maker

biophys_components
------------------

Notice that there are two different folders of both mechanisms and templates.
When running a simulation, you will need to remove the prefix (either L2-3 or L5)
from the two folders you want to be active, depending on which cell you want.
Remove the prefix of the cell you want to be using.

Files
-----

:doc:`clustering.py <clustering>`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Contains the classes and functions used to build functional groups and synaptic clustering.

run.py
^^^^^^

.. automodule:: L5NeuronSimulation.run
   :members:


:doc:`raster_maker.py <raster_maker>`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Contains the functions and class (SonataWriter) necessary for generating and saving the input spike rasters.

my_plotting.py
^^^^^^^^^^^^^^

.. automodule:: L5NeuronSimulation.my_plotting
   :members:
