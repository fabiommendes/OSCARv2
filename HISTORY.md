# History

This file describes the history of editions and strategies used to evolve OSCAR 
code.

## Consolidate Python 2 code

[1b882b1720ab723cda2ca33113d23f700618bc96]

Create new non-stochastic and single-run tests. Default test run at rate of one
test every minute.

Fix bug at OSCAR_lite() function that mistakenly used `is` to check equality
of two strings.

[1b882b1720ab723cda2ca33113d23f700618bc96]
Consolidate all files into a single script to avoid problems with readfile() 
when converting to Python 3. For some reason tests are running at a noticiable
slower speed.

[8cc787be95de190c47e4fe98b4a3eebdf8f0b4dd]
Create a very simple modularization strategy. Basically we are replacing execfile()
with loading strings of text and calling exec()

[fede8d7f177e4ba6f1e5f25b6c16ed2d40039ca3]
Convert script to Python 3 and execute Black. Python 3 is now running 
approximately twice as slower as Python 2 and black formatted code is about
50% as big with 120 chars lines.

[529683bce6e9b9818359a69d5b6c57d7eea535b5]
Split data loading into several modules. Now, what used to be inside OSCAR-loadD
is moved to oscar.oscar_data.*. The names of modules are temporary e we still 
expect a lot of refactoring. It takes ~3.3 seconds to load all data modules.

[??]
Eliminate all exec's from code. Cleaned parameter, data and the OSCAR_lite 
function from all exec() statements. Now code us much more amenable to 
refactoring. Removing all exec's() produced a noticiable speed improvement to
about 45s runtime.

[??]
Remove all direct references from csv module. Abstracted the data loading 
procedures to call some OSCAR-specific functions that (for now) call the
CSV module, but can be later abstracted to use more efficient data storage
methods.

[??]
Add a caching scheme to load_data()-like functions to cache results from loading
the same dataset with the same parameters.

[??]
Remove all star imports and explicitly load variables from other modules.

[??]
Start refactoring architecture: group by role instead of element: historical 
data, scenarios, chemestry, etc.


## Extract data

## Instrumentation and tests


## Modularize

## Split functionalities into classes

## Create plugabble classes

## Optimize with Cython
