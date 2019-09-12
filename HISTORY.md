# History

This file describes the history of editions and strategies used to evolve OSCAR 
code.

## Consolidate Python 2 code (doing)

[1b882b1720ab723cda2ca33113d23f700618bc96]

Create new non-stochastic and single-run tests. Default test run at rate of one
test every minute.

Fix bug at OSCAR_lite() function that mistakenly used `is` to check equality
of two strings.

[?]
Consolidate all files into a single script to avoid problems with readfile() 
when converting to Python 3. At this point, code cannot execute runs with 
different parameters anymore.  

## Convert to Python 3

## Extract data

## Instrumentation and tests

## Eliminate all exec's from code

## Modularize

## Split functionalities into classes

## Create plugabble classes

## Optimize with Cython
