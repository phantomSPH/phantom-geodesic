# phantom-geodesic
geodesic integrator for relativistic orbits

Mainly used as a test code for the PHANTOM General Relativistic Smoothed Particle Hydrodynamics code, but functional as a standalone geodesic integrator

Please cite [Liptai & Price (2019), MNRAS 485, 819-842](https://ui.adsabs.harvard.edu/abs/2019MNRAS.485..819L) if you use this code

## Usage
```
make
./grtest
```

## Requirements

- gfortran compiler

## Options

Can change metric, setup or whether or not to use Newtonian gravity, e.g.

```
make NEWTONIAN=yes
```
And yes, compile-time options are a bit annoying. This would be nice to fix.

## Authors

(c) 2019- David Liptai & Daniel Price

Code was written mainly by David Liptai, with various modules imported from Phantom
