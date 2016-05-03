# GenomicStore

This module is a place to put various Julia functions I need to utilise HDF5 as a store for genomic tracks. This (master) branch is composed of things that work for me for a specific application. A refactor into a more generally useful version of this is being done on the [restructure branch](https://github.com/nw11/GenomicStore.jl/tree/restructure).

GenomicStore aims to provides a convenient structure for accessing and storing genome annotation "tracks".

This is provided by interfacing with a [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) data store at the moment.
