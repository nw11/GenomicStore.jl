#=
  Notes: 
   might be better renamed as jld_genomic_store_bulk_load.jl
=#

include( Pkg.dir("GenomicStore","src","util.jl") )
include( Pkg.dir("GenomicStore","src","backends","jld_store.jl") )

function store_samples( genomic_store_path, sample_info_path, chr_sizes_path )

    track_df=read_sample_info(sample_info_path)
    for row=1:nrow(track_df)
        track_id     = track_df[row, :track_id]
        file_type    = track_df[row, :filetype]
        filepath    = track_df[row, :filename]
        if any(x->x == :coordshift,names(track_df))
            coord_shift  = track_df[row, :coordshift]
        else
            coord_shift=0
        end
        if any(x->x == :na_val,names(track_df))
             na_val = track_df[row, :na_val]
        else
            na_val=0.0
        end

        if any(x->x == :strandfilter,names(track_df))
            strand_filter = track_df[row, :strandfilter]
        else
            strand_filter = nothing
        end

        if !isabspath( filepath )
            file_path = abspath(filepath)
        end

        if track_exists(genomic_store_path,track_id)
            Lumberjack.info("==Found track $track_id in database: skip==")
            continue
        end

        Lumberjack.info("==Saving track $track_id to database==")
       
        # replace this with save_track()
        if file_type == "methpipe_bed_levels"
            store_methpipe_bed_points(filepath,
                                      genomic_store_path,
                                      track_id ,
                                      chr_sizes_path,
                                      na_val=na_val,
                                      coord_shift=coord_shift,
                                      measurement="levels")
        elseif file_type == "methpipe_bed_coverage"
            store_methpipe_bed_points(filepath ,
                                      genomic_store_path ,
                                      track_id ,
                                      chr_sizes_path,
                                      na_val=na_val,
                                      coord_shift=coord_shift,
                                      measurement="coverage")
        elseif file_type == "cpg_point"
            store_cpg_points(filepath ,
                                   genomic_store_path ,
                                   track_id ,
                                   chr_sizes_path,
                                   strand_filter_str="-",
                                   coord_shift= coord_shift)
        elseif file_type == "methpipe_bed_interval"
           store_methpipe_bed_intervals(filepath ,
                                   genomic_store_path ,
                                   track_id ,
                                   chr_sizes_path,
                                   coord_shift= coord_shift)
        elseif file_type == "bedgraph_cpg"
            store_bedgraph_cpg(filepath ,
                                   genomic_store_path ,
                                   track_id ,
                                   chr_sizes_path,
                                   coord_shift= coord_shift)
        else
            Lumberjack.info("Skip track - $track_id because filetype or storagetype not yet supported")
        end
    end
end


