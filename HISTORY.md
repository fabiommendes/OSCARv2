# History

This file describes the history of editions and strategies used to evolve OSCAR 
code.

[1b882b1720ab723cda2ca33113d23f700618bc96]
Create new non-stochastic and single-run tests. Default test run at rate of one
test every minute.

Fix bug at OSCAR_lite() function that mistakenly used `is` to check equality
of two strings.

[fa1248439764ff5b18d7dd20099a4dd39d511069]
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

[41b9d4ee97ad107a3e61105ec0b554a7b789a4f3]
Eliminate all exec's from code. Cleaned parameter, data and the OSCAR_lite 
function from all exec() statements. Now code is much more amenable to 
refactoring. Removing all exec's() produced a noticeable speed improvement and 
reached an average of 45s runtime.

[46881d374a47dfcc50ba8d3002e15ef0ccf56049]
Remove all direct references from csv module. Abstracted the data loading 
procedures to call some OSCAR-specific functions that (for now) call the
CSV module, but can be later abstracted to use more efficient data storage
methods.

[70c8d02ee63e331e964ebe3140b730327b054667]
Add a caching scheme to load_data()-like functions to cache results from loading
the same dataset with the same parameters.

[f740923a5235fecf7b260dbdbd67f9fb8346e6fe]
Remove all star imports and explicitly load variables from other modules.

[6408acd62984e2dbf6c254fa2aba93b24ff2ced5]
Start refactoring architecture: group by role instead of element: historical 
data, scenarios, chemistry, etc.

[42b87ea9006823f653aa64c3c5f4a37ec34e61c6]
Change OSCAR_lite function and split implementation into several methods and
mixin classes.

[07a1e77e62cf954285e3003926b7740a18e69697]
Make module names simpler and modify them to the proposed desired final 
configuration.

[efb86ad7020800b585c48f8db2639e82f78b3d31]
Move configurations to a separate object. All modules now import the conf object
and use it under the "conf" namespace.

## Reorganize data around pandas dataframes and more agressive use of cache

## Instrumentation and tests

## Refactor simulation class to use composition instead of inheritance

## Create plugabble simulation classes

## Optimize with Cython
