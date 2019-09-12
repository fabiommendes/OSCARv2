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

[??]
Convert script to Python 3 and execute Black. Python 3 is now running 
approximately twice as slower as Python 2 and black formatted code is about
50% as big with 120 chars lines.

[??]
Eliminate all exec's from code


## Extract data

## Instrumentation and tests


## Modularize

## Split functionalities into classes

## Create plugabble classes

## Optimize with Cython
