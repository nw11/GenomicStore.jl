#=
  jld_genomic_store.jl

  This module implements the user level query interface.

 =#
using Lumberjack

include( Pkg.dir("GenomicStore","src","backend","jld","jld_genomic_store_crud.jl") )
include( Pkg.dir("GenomicStore","src","util.jl") )
include( Pkg.dir("GenomicStore","src","readfile.jl"))


"""
    save_track(file,genomic_store,file::BedGraphMetadata)

    save_track() for bedGraph format

    This stores the value in the interval at the start position.
    It assumes that the start positions of the intervals do not overlap, since bedGraph
    is often used directly as a display format, this assumption is implied.

    The bedGraph files should be sorted.

    defaults:
        start_coord_shift=1 (bedgraph is zero based, and GenomicStore uses 1 based)
"""
function save_track(genomic_store::JldGenomicStore,file::BedGraphMetadata, track_id,chr_sizes_path;
                             start_coord_shift=1,
                             OUT_OF_RANGE_VAL=0.0,
                             gzip=false,
                             val_type="float32",)

   (seq_ids,starts,stops,scores) = read_file(file)
   bedgraph_file_length = length(seq_ids)
   Lumberjack.info("Read bedgraph file, total length: $bedgraph_file_length")
   _save_track(track_id,seq_ids,starts,stops,scores, chr_sizes_path,start_coord_shift=start_coord_shift, out_of_range_val=0.0,gzip,val_type=val_type)
end

function  _save_track(track_id,seq_ids,starts,stops,scores, chr_sizes_path;start_coord_shift=0, out_of_range_val=0.0,gzip=false,val_type=val_type)

    chr_sizes_dict=get_chr_sizes_dict(chr_sizes_path)
    db=getdb(genomic_store_path) # from crud.jl

    # initialise for the first sequence_id
    seq_id=seq_ids[1]
    seq_len  = chr_sizes_dict[seq_id]
    seq_vals = fill(OUT_OF_RANGE_VAL,seq_len)

    # initialise strucutres for checking sorted data
    seen_seq_ids = Set()
    last_start_pos = 0
    println("sequence id: $seq_id")
    for i=1:bedgraph_file_length
        # if there is a change in sequence_id
        if seq_id != seq_ids[i]
            println("change: $seq_id, $(seq_ids[i])")
            if( in(seq_id,seen_seq_ids ) )
                Lumberjack.error("Unordered sequence identifier found at $i ($seq_id)")
            end

            write_track(db,track_id,seq_id,seq_vals)

            # -- indicate we have processed the current seq_id
            push!(seen_seq_ids,seq_id)

            # -- reset last start position
            last_start_pos=0

            # -- update seq_id
            seq_id   = seq_ids[i]
            Lumberjack.info("Next sequence name: $seq_id")

            # -- check the id is valid --
            if !haskey(chr_sizes_dict,seq_id)
                Lumberjack.error("Invalid ID: ( $seq_id ), this file only saved up to line $(i-1)")
            end

            seq_len  = chr_sizes_dict[seq_id]
            seq_vals = fill(OUT_OF_RANGE_VAL,seq_len)
            Lumberjack.info("Next track to be $seq_len in length")
        end


        #--- check that the start positions are sorted
        if last_start_pos > starts[i]
            Lumberjack.error("Unsorted position found at $i ($last_start_pos > $(starts[i])")
        end
        # update last_start_pos
        last_start_pos = starts[i]

        #--- check that the stop position is not less than the start position
        if starts[i] > stops[i]
            Lumberjack.error("Start position greater than stop position line $i $(starts[i]) > $(stops[i])")
        end

        if val_type == "float32"
            try
                seq_vals[ starts[i] + start_coord_shift ] = Float32( scores[i] )
            catch y
                println("processing this line:\n")
                curr_row = join([ seq_ids[i], starts[i],stops[i], scores[i] ],"\t" )
                println(curr_row)
                println("shift: $start_coord_shift")
                Lumberjack.error(y)
            end
        end
        if val_type == "int32"
            seq_vals[ starts[i] + start_coord_shift ] = Int32( scores[i] )
        end
        if val_type == "int16"
            seq_vals[ starts[i] + start_coord_shift ] = Int16( scores[i] )
        end
        if val_type == "int8"
            seq_vals[ starts[i] + start_coord_shift ] = Int8( scores[i] )
        end
    end
    # write the final track
    write_track(db,track_id,seq_id,seq_vals)
    Lumberjack.info("Written track $seq_id")
end



function track_exists(genomic_store::JldGenomicStore, track_id)
   db=getdb(genomic_store.path)
   return _hastrack(db,track_id)
end

