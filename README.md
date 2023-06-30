# Cherenkov detectors fast simulation using neural networks
The code used to generate the data for

> Derkach, D., Kazeev, N., Ratnikov, F., Ustyuzhanin, A., & Volokhova, A. (2020). Cherenkov detectors fast simulation using neural networks. Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, 952, 161804.

Forked from https://github.com/jmhardin/FastDIRC
# Building and running
Type

      make

to build.  Root should be the only dependency right now.

for a sample run try:
```bash
./dircfit -n <number of events> -E <Energy in GeV> -particle_theta <theta> -particle_phi <phi>
./dircfit -n 1000 -E 4.5 -particle_theta 4 -particle_phi 20
```  
This is "loop" mode and will run very quickly - the KDE is only simulated once.
This will produce "fitdirc.root".  This can be changed with the `-of` option.
There are also a lot of command line options to perturb the box geometry, they can be viewed in the source.

A signifigantly pared down driver file is also provided under `example_driver/` - type:
```bash
make example
```
to build it.  It accepts a subset of the commands of the regular file.  It also does not contain the experimental algorithms and input methods and so should be signifigantly easier to read.

To read from a csv file, use `-if foo.csv`

it expects the format:
```
Event_NUM PID_NUM BAR x_mm y_mm time_ns theta_deg phi_deg Energy_GeV
```
where event_num comes from the Montecarlo and pid is the geant PID. Bar is an indexing of the bars from the middle. 1 is the closest bar to the beam on the right, with 24 being the furthest right bar. -1 to -24 behave similarly on the left.  x and y are mm counts from the middle of the bar (NOT absolute coordinates).

Use can use the command flag `-force_kinematics` to overide the last 3 colums.  Note that force kinematics will override them even if no other options are passed, as they have default values.

The flag `-t <ns>` will cause events which are marked with a timestamp within that number of ns to produce hits on the PMT plane and "confound" the results.  Useful for testing windows and backgrounds.  Defaults to `-1`.

To fit a resolution, rung `graphicHistos.C(<file>,bool verbose, double Energy (GeV))` on the rootfile output from the simulation.
