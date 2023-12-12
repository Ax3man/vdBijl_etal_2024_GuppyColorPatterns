# Selection decisions

In the finished folder are all the scripts used to select the fish used for the next generation of selection, and for sequencing. It also contains all the selection decisions in .csv format.

The pertinent information for data analysis (i.e. which selection direction the fish were in, and whether they were selected) can be loaded by sourcing `compile_decisions.R`.

The selection scripts were written to be used to update selection decisions during phenotyping, for logistical reasons. They also contain a lot of extra code to keep track of phenotyping performance etc.

Selection scripts are named in the form of `[replicate]_[generation].R` for males, and `fem_[replicate]_[generation].R` for females.

