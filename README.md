PosteriorEntropy
================

This is an interface to [InfoNest](https://github.com/eggplantbren/InfoNest)
that simplifies the process of calculating H(theta | data).

(c) 2018 Brendon J. Brewer. LICENCE: MIT.

Usage
=====

First, clone the repository recursively:

`git clone --recursive https://github.com/eggplantbren/PosteriorEntropy`

Then compile the C++.

`make`

Run the executable:

`./main`

View the results (can be done while `main` is still running):

`python postprocess.py`

The last line of output shows the measured (differential) conditional entropy
H(theta | data), i.e., the expected entropy of the posterior.
The example model is specified in include/Demo.h. The parameter of
interest in the demo is the log width of the transit, and the
prior entropy is H(theta) = 1.419 nats.
