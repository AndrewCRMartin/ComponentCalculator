ComponentCalculator
===================

## EARLY INITIAL CODE - VERY BASIC FUNCTIONALITY

### TODO:

- add command line options
- increase probability of altering component count if we are targeting
  low or high
- add sexual repro
- optimize probabilities
- optimize score function
- proper output display
- display values better
- default population size needs to change with E series
- adaptive mutation of values

-------------------------------------------------------------------------

(c) 2025, Andrew C.R. Martin

`ComponentCalculator` is a program to calculate combinations of
standard value resistors or capacitors that will make up a required
value. You can specify:

- Resistor or capacitor (default: capacitor)
- The maximum number of components to use (default: 5)
- The minimum number of components to use (default: 1)
- The required accuracy (default: 1%)
- The E-series of values (default: e24 for resistors, e6 for capacitors)
- Instead of an E-series, a file can be provided with allowed values

The method uses an evolutionary algorithm to find the best combination
of components by optimizing a function of the accuracy and the
similarity of values used. Having the values closer to one another
allows you to exploit the statistics of the distributions of values
for a given tolerance which improves the overall tolerance of a random
pair of components (see Douglas Self's book on Small Signal
Amplifiers).
