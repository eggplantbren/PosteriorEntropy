TransitInfo
===========

(c) 2018 Brendon J. Brewer and others. LICENCE: MIT.

Usage
=====

First, clone the repository recursively:

`git clone --recursive https://github.com/eggplantbren/PlanetInfo`

Compile the C++:

`make`

Run the executable:

`./main`

View the results (can be done while `main` is still running):

`python postprocess.py`

The last line of output shows the measured conditional entropy
H(theta | data), i.e., the expected entropy of the posterior.
The example model is specified in include/Demo.h.