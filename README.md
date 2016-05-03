# GenomicStore

GenomicStore is a Julia package for querying and storing genome annotation. These are typically organised as "tracks" that represent values associated with coordinates on a reference genome.

GenomicStore was initially built with support for a [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) based data store (see the [master branch here](https://github.com/nw11/GenomicStore.jl)).  The plan is to support other data store types based on other file formats, e.g. BigWig & BigBed. 

The supported HDF5 based data store uses the [Julia Data package](https://github.com/JuliaLang/JLD.jl) for convenience, so is denoted as a "JldGenomicStore".

Example usage is: 

```
genomic_store = JldGenomicStoreMetadata("/path/to/datastore.jld")
query_track(genomic_store,"chip-seq-track1","chr1", 1000000, 3000000) 
```

Support for other genomic stores would entail only defining 
 a new type of GenomicStore, for example:

```
genomic_store = BigWigGenomicStoreMetadata("/path/to/dir/of/bigwigs")
query_track(genomic_store,"chip-seq-track1","chr1", 1000000, 3000000)
```


## Intended basic API:

Create:
- save_track(genomic_store,filetype; kwargs... )

Retrieval:
- query_track(genomic_store,track_id,seq_id,start_pos,end_pos)
- query_track(genomic_store,track_ids,seq_ids::Array,starts_pos::Array,ends_pos::Array)

Update:
- save_track(genomic_store,track_id,seq_id,start_pos,end_pos,values)
 
(values must be of the same type of the type of track)

Delete:
- remove_track(genomic_store,track_id)
- remove_seq(genomic_store,track_id,seq_id)


## JLD GenomicStore
A JldGenomicStore can be bulk loaded by describing all the files to be loaded as tracks in a separate file. 
More details on this coming.

## Status
This branch is currently under significant development so this is pre-alpha at the moment.
More documentation coming. 
