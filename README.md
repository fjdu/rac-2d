## Install

Go to the ```src``` directory.  Inside it there is a ```makefile```.  You may
edit it for your own needs (but you don't have to).  Then run
```bash
    make
```
and an executable file with default name ```a.out``` will be generated in the
same directory.


## Run the code

There are a few input data files that are needed for the code to run.

The following files are compulsary:

    1. Configuration file.
    2. Chemical network.
    3. Initial chemical composition.
    4. Dust optical properties.

The following files are optional:

    1. Density structure.
    2. Enthalpy of formation of species.
    3. Molecular ransition data.
    4. Stellar spectrum.
    5. Points to output the intermediate steps of chemcial evolution.
    6. Species to output the intermediate steps of chemcial evolution.
    7. Species to check for grid refinement.

By default all these files are in the ```inp``` directory, though some of them
do not have to.  Go to this directory.  Edit the file ```configure.dat```.  It
has nearly 200 entries.  Some of them are for setting up the physics and
chemistry of the model, some are for setting up the running environment, while
others are switches telling the code whether or not it should execute some
specific tasks.  Details for editing the configure file are included below.

After you have get the configre file ready, and have all the needed files in
place, you can go to the directory on top of ```inp```, and then in a terminal
type in
```
    ./src/a.out ./inp/configure.dat
```
to start running the code.

## Contents of configure.dat

The configuration file is in the Fortran namelist format, so you may want to
set the language type for syntax highlighting of your editor to Fortran.
