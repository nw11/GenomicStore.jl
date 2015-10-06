#=
 store.jl

 A track can be in

 1. mask/interval format
 2. point format

 **mask/interval** format represents intervals, for example contained in a bed file  as consecutive 1's. The
 coordinates in between the intervals contain 0's.

 **point format** represents values associated with a particular base pair coordinate. Coordinates that
 do not contain values are given the same out of range number

 The aim of the io.jl is to parse various file formats and extract the right values and put these into the store
 in either of the above formats.

 Thus far the formats that are most relevant to deal with are:

 a) methpipe's bed format for cpgs
 chr1    3000826 +       CG      1       17
 chr1    3001006 +       CG      0.8     10
 chr1    3001017 +       CG      1       12
 chr1    3001276 +       CG      1       14
 chr1    3001628 +       CG      0.867   15

 b) methpipe's bed format for hypomethylated regions
 chr1    3670331 3672761 HYPO0   184     +
 chr1    3993551 3993737 HYPO1   5       +
 chr1    4490611 4490904 HYPO2   6       +
 chr1    4491370 4497651 HYPO3   275     +
 chr1    4571560 4572248 HYPO4   58      +

 c) tab delimited format describing CpG positions
 chr1    3000826 +
 chr1    3000827 -
 chr1    3001006 +
 chr1    3001007 -
 chr1    3001017 +

 Both inserting intervals and points

 Denote:
  a) as  methpipe_bed_point
  b) as  methpipe_bed_region
  c) as  cpg_point

 Note that methpoint_point has both floating point and integer values associated with it.

 REQUIRES API:
 store_methpipe_bed_points(value_type="levels") # can also be value_type="coverage"
 store_methpipe_bed_intervals()
 store_cpg_points()

 These functions are implemented in a store-backend.jl

 We also should aim for a general form :

 # Stores Points
 store_tab_file( coord_column_idxs=[1,2], value_column_idx=[5] )

 # Stores Intervals as mask
 #  - since value is empty and three column specification for the coordinate
 store_tab_file(coord_column_idx=[1,2,3],value_column_idx=[])

 # Stores Points as mask (1 or 0) because no value specified
 store_tab_file(coord_column_idx=[1,2],value_column_idx=[])

=#

using Compat
using HDF5,JLD
using GZip
using Lumberjack
using DataFrames
include("store-backend1.jl")

function read_sample_info(sample_info_path)
    try
        if !isfile(sample_info_path)
            Lumberjack.error("The sample_info path: $sample_info_path is not a valid file")
        end
        sample_info_df = readtable(sample_info_path, allowcomments=true)
        colnames=names(sample_info_df)
        Lumberjack.info("Column names in sample info: $colnames")
        required_colnames = [:filename,:filetype,:track_id]
        for required_colname in required_colnames
            colname_appears = any(x->x==required_colname,colnames)
            if !colname_appears
               Lumberjack.error("Required column $required_colname does not appear in the sample information")
            else
               Lumberjack.info("Found required column $required_colname in sample information")
            end
        end
        return sample_info_df
    catch e
        Lumberjack.error("$e\nError reading table: Check file exists and format correct.")
    end
end

function store_samples( genomic_store_path, sample_info_path, chr_sizes_path )
    # check for columns
    # store track with track id and type
    # GenomeStore(genome_store_path,sample_info_path, track_ids, AND OTHER THINGS)
    track_df=read_sample_info(sample_info_path)
    for row=1:nrow(track_df)
        track_id     = track_df[row, :track_id]
        file_type    = track_df[row, :filetype]
        #storage_type = track_df[row, :storagetype]
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

        Lumberjack.info("==Saving track $track_id to database==")

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
           # Lumberjack.info("Skip track - $track_id because filetype or storagetype not yet supported")
            store_cpg_points(filepath ,
                                   genomic_store_path ,
                                   track_id ,
                                   chr_sizes_path,
                                   strand_filter_char="-",
                                   coord_shift= coord_shift)
        elseif file_type == "methpipe_bed_interval"
           store_methpipe_bed_intervals(filepath ,
                                   genomic_store_path ,
                                   track_id ,
                                   chr_sizes_path,
                                   coord_shift= coord_shift)
        else
            Lumberjack.info("Skip track - $track_id because filetype or storagetype not yet supported")
        end
    end
end


