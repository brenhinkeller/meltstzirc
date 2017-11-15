# meltstzirc
Source code for the alphaMELTS zircon saturation calculations used in the paper ["Temporal variation in relative zircon abundance throughout Earth history"](https://doi.org/10.7185/geochemlet.1721) in Geochemical Perspectives Letters by Keller, Boehnke, and Schoene

The code is written mostly in C, using MPI to distribute whole-rock compositions from an input file ([ignmajors.csv](ignmajors.csv)) to cores on different nodes, each of which then spins up an alphaMELTS simulation using a system() call, reads the output, and conducts Zr partitioning calculations using the mineral abundances reported by MELTS. 

Two versions of the main parallel code are provided: [meltsTzircParallelFull.c] which outputs zircon saturation state for each temperature step of the MELTS simulation, and [meltsTzircParallel.c] which outputs only the total mass saturated. 

Mineral/melt parition coefficients for Zr are obtained from the GERM database, as stored in [GERM.h], while a range of functions for interacting with alphaMELTS are provided in [runmelts.h]

