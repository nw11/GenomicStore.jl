module GenomicStore

# package code goes here]
export get_track_interval
export get_chr_sizes_dict
export save_point_start_track
export save_interval_track
export save_bedgraph_track
export save_bed_track
export store_info_str
export store_info_dict
include(Pkg.dir( "GenomicStore","src","genomic-store.jl"))

export getdb
export delete_track!
include(Pkg.dir( "GenomicStore","src","crud.jl"))
end # module
