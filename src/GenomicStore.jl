module GenomicStore

export get_track_interval
export get_chr_sizes_dict
export save_point_start_track
export save_interval_track
export save_bedgraph_track
export save_bed_track
export store_info_str
export store_info_dict
include("genomic-store.jl")

export getdb
export delete_track!
include("crud.jl")

export store_samples
include("store.jl")
end # module
