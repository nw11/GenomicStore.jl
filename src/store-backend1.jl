  # This file provides functions for reading files into and hdf5a db
using Compat
using GZip
using Libz
using Lumberjack
using HDF5,JLD
include( Pkg.dir("GenomicStore","src","crud.jl") )

function read_gzip_file(filename)
   gzio=gzopen(filename,"r")
   line_num =0
   lines=ASCIIString[]
   for line in eachline(gzio)
        line=chomp(line)
        push!(lines,line)
        line_num +=1
        if line_num % 1000000 == 0
            Lumberjack.info("Procesed $line_num")
        end
    end
    return lines
end

function memory_read_file(filename)
    io = open(filename)
    Lumberjack.info("reading all")
    file=readall(io)
    Lumberjack.info("finished reading all")
    Lumberjack.info("split line")
    lines=split(file,'\n')
    Lumberjack.info("done split")
    close(io)
    return lines
end

function get_chr_sizes_dict(chrom_sizes_path)
    chrom_sizes=readdlm(chrom_sizes_path)
    chrom_sizes_dict = Dict{ASCIIString,Int64}()
    size(chrom_sizes,1)
    for i=1:size(chrom_sizes,1)
        chrom_sizes_dict[ chrom_sizes[i,1] ] = chrom_sizes[i,2]
    end
    return chrom_sizes_dict
end

function isgzip(filepath)
    (path,ext) = splitext(filepath)
    gzip = (ext==".gz" || ext==".gzip" ) ? true : false
    return gzip
end


#=
 memory_read_and_parse_point_file

 Reads into memory a file and parses a file in the format:

   chr1    3000826 +
   chr1    3000827 -

 returns tuple of 2 arrays
 1. seq_ids
 2. start positions
=#
function memory_read_and_parse_point_file(filename;strand_filter_char=nothing)
    io = open(filename)
    Lumberjack.info("reading all")
    file=readall(io)
    Lumberjack.info("finished reading all")
    Lumberjack.info("split line")
    lines=split(file,'\n')
    Lumberjack.info("done split")

    num_lines=length(lines)-1
    seq_ids = ASCIIString[] #fill("",num_lines)
    starts=   Int64[] #fill( 0, num_lines)
    num_kept = 0
    Lumberjack.info("parse and assign to array")
    #parseint takes for ever, so annoying about this.
    for idx=1:num_lines
        (seq_id,start,strand)= split(strip(lines[idx]),'\t')
        if strand == strand_filter_char
            continue
        end
        push!(seq_ids,seq_id)
        push!(starts,int64(start))
        if idx % 1000000 == 0
            println(idx)
        end
        num_kept+=1
    end
    Lumberjack.info("finished parsing integer and assigning to array ( number kept after filter if applied $num_kept). first seq_id is $(seq_ids[1]). Last seq_id $(seq_ids[end])")
    return (seq_ids,starts)
end



#=
 memory_read_and_parse_methpipe_cg_bed_coverage

 Reads into memory a file and parses a file in the format:

  chr1    3000826 +       CG      1       17
  chr1    3001006 +       CG      0.8     10

 extracting out the 6th column (the coverage column)

 returns tuple of 4 arrays

 1. seq_ids
 2, starts
 3. stops
 4. coverages

=#

function memory_read_and_parse_methpipe_cg_bed_coverage(filename;gzip=false)
    if gzip
        lines=read_gzip_file(filename)
    else
        lines = memory_read_file(filename)
    end
    num_lines=length(lines)-1
    seq_ids = fill("",num_lines)
    starts=fill( 0, num_lines)
    stops=fill( 0,  num_lines )
    coverages=fill(int32(0), num_lines )
    Lumberjack.info("parse and assign to array")
    #parseint takes for ever, and that is what is so annoying about this.
    for idx=1:num_lines
        (seq_id,start,strand,context,score,coverage)= split(lines[idx],'\t')
        try
            seq_ids[idx] = seq_id
            starts[idx]  = int64(start)
            stops[idx]   = int64(start) +1
            coverages[idx]  = int32(coverage)
            if idx % 1000000 == 0
                println(idx)
            end
       catch y
          println("processing this line($idx):\n")
          println(lines[idx])
          println("parsed to:")
          println("$seq_id $start $strand $context $score $coverage")
          error(y)
       end
    end
    Lumberjack.info("finished parsing integer and assigning to array")
    return (seq_ids,starts,stops,coverages)
end


#=
 memory_read_and_parse_methpipe_cg_bed_levels

 Reads into memory a file and parses a file in the format:

  chr1    3000826 +       CG      1       17
  chr1    3001006 +       CG      0.8     10

 extracting out the 6th column (the coverage column)

 returns tuple of 4 arrays

 1. seq_ids
 2, starts
 3. stops
 4. levels

=#


function memory_read_and_parse_methpipe_cg_bed_levels(filename;gzip=false)
    if gzip
        lines=read_gzip_file(filename)
    else
        lines = memory_read_file(filename)
    end
    num_lines=length(lines)-1
    seq_ids = fill("",num_lines)
    starts=fill( 0, num_lines)
    stops=fill( 0,  num_lines )
    scores=fill(float32(0.0), num_lines )
    Lumberjack.info("parse and assign to array")

    for idx=1:num_lines
        (seq_id,start,strand,context,score,coverage)= split(lines[idx],'\t')
        seq_ids[idx] = seq_id
        starts[idx]  = int64(start)
        stops[idx]   = int64(start) +1
        scores[idx]  = float32(score)
        if idx % 1000000 == 0
            println(idx)
        end
    end
    Lumberjack.info("finished parsing integer and assigning to array")
    return (seq_ids,starts,stops,scores)
end



#=
 memory_read_and_parse_interval_file

 This function expects a file to have the first three columns to contain
 a sequence id, start position, and stop/end position

 Initially used to parse files like below:

 chr1    3670331 3672761 HYPO0   184     +
 chr1    3993551 3993737 HYPO1   5       +

 Returns a tuple of vectors

 1. sequence_ids
 2. start_positions
 3. stop_positions

=#
function memory_read_and_parse_interval_file(filename; only_start=false,strand_filter_char=nothing)
    io = open(filename)
    Lumberjack.info("reading all")
    file=readall(io)
    Lumberjack.info("finished reading all")
    Lumberjack.info("split line")
    lines=split(file,'\n')
    Lumberjack.info("done split")

    num_lines=length(lines)-1 # lines -1 as it splits the final row
    seq_ids = fill("",num_lines)

    starts=fill( 0, num_lines)
    if !only_start
       stops=fill( 0,  num_lines )
    end
    Lumberjack.info("parse and assign to array")

    for idx=1:num_lines
        if only_start
            (seq_id,start)= split(lines[idx],'\t')
        else
            (seq_id,start,stop)= split(lines[idx],'\t')
        end
        seq_ids[idx]=seq_id
        starts[idx]=int64(start)
        if !only_start
            stops[idx]=int64(stop)
        end

        if idx % 1000000 == 0
            println(idx)
        end
    end
    Lumberjack.info("finished parsing integer and assigning to array")
    if only_start
        return (seq_ids,starts)
    end
    return  (seq_ids,starts,stops)
end


#=

  memory_read_and_parse_bedgraph

  This function expects to have the first three columns to contain
  seq id,start position, and stop/end postion.

  e.g.

  chr1 3670331 3670333 0.05

  Returns a tuple of vectors

  1. sequence_ids
  2. start_positions
  3. stop_positions
  4. score

=#

function memory_read_and_parse_bedgraph(filename;gzip=false)
    if gzip
        lines=read_gzip_file(filename)
    else
        lines = memory_read_file(filename)
    end
    #io = open(filename)
    #Lumberjack.info("reading all")
    #file=readall(io)
    #Lumberjack.info("finished reading all")
    #Lumberjack.info("split line")
    #lines=split(file,'\n')
    #Lumberjack.info("done split")
    # now we know what size thigns should be
    # get size and two columns and then
    # lines -1 as it splits the final row
    num_lines=length(lines)-1
    seq_ids = fill("",num_lines)
    starts=fill( 0, num_lines)
    stops=fill( 0,  num_lines )
    scores=fill(float32(0.0), num_lines )
    Lumberjack.info("parse and assign to array")
    #parseint takes for ever, and that is what is so annoying about this.
    for idx=1:num_lines
        (seq_id,start,stop,score)= split(lines[idx],'\t')
        seq_ids[idx] = seq_id
        starts[idx]  = int64(start)
        stops[idx]   = int64(stop)
        scores[idx]  = float32(score)
        if idx % 1000 == 0
            println(idx)
        end
    end
    Lumberjack.info("finished parsing and assigning to array")
    return (seq_ids,starts,stops,scores)
end

function read_and_parse_bedgraph(filename;gzip=false)
    if gzip
        stream=ZlibInflateInputStream(open(file),reset_on_end=true)
    else
        stream= open(file)
    end
    seq_ids=UTF8String[]
    starts=Int64[]
    stops=Int64[]
    scores=Float32[]
    count=0
    for line in eachline(stream)
           (seq_id,start,stop,score)=split(line,'\t')
           push!(starts,parse(Int64,start) )
           push!(scores,parse(Float32,score))
           push!(stops,parse(Int64,stop))
           push!(seq_ids,seq_id)
           count+=1
           if count % 1000000 == 0
             Lumberjack.info("Read $count")
           end
       end
    Lumberjack.info("finished parsing and assigning to array")
    return (seq_ids,starts,stops,scores)
end

#=

save_bed_track

 This function saves "bed" like files to hdf5 storage. The name should more accurately
 be "save_bedlike_track"

 It handles saving values specified in

 a) bedgraph format
 b) methpipe_cpg_bed
     i) extracting the methylation levels value
     ii) extracting the coverage levels

 Current requirement that input file is pre-sorted
=#

function save_bed_track(genomic_store_path,input_file,track_id,chr_sizes_path;
                             start_coord_shift=0,
                             OUT_OF_RANGE_VAL=0.0,
                             gzip=false,
                             val_type="float32",
                             bedtype="bedgraph")
    chr_sizes_dict=get_chr_sizes_dict(chr_sizes_path)
    #--- sort in temporary dir
    #--- load
    if bedtype == "bedgraph_cpg"
        (seq_ids,starts,stops,scores)=read_and_parse_bedgraph(input_file,gzip=gzip)
    elseif bedtype == "methpipe_cpg_bed_levels"
        (seq_ids,starts,stops,scores)=memory_read_and_parse_methpipe_cg_bed_levels(input_file,gzip=gzip)
    elseif bedtype == "methpipe_cpg_bed_coverage"
        (seq_ids,starts,stops,scores)=memory_read_and_parse_methpipe_cg_bed_coverage(input_file,gzip=gzip)
    else
         Lumberjack.error("Invalid bedtype specified $bedtype")
    end
    bedgraph_file_length = length(seq_ids)
    Lumberjack.info("Read bedgraph file, total length: $bedgraph_file_length")

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
                seq_vals[ starts[i] + start_coord_shift ] = float32( scores[i] )
            catch y
                println("processing this line:\n")
                curr_row = join([ seq_ids[i], starts[i],stops[i], scores[i] ],"\t" )
                println(curr_row)
                println("shift: $start_coord_shift")
                Lumberjack.error(y)
            end
        end
        if val_type == "int32"
            seq_vals[ starts[i] + start_coord_shift ] = int32( scores[i] )
        end
        if val_type == "int16"
            seq_vals[ starts[i] + start_coord_shift ] = int16( scores[i] )
        end
        if val_type == "int8"
            seq_vals[ starts[i] + start_coord_shift ] = int8( scores[i] )
        end
    end
    # write the final track
    write_track(db,track_id,seq_id,seq_vals)
    Lumberjack.info("Written track $seq_id")
end


#=
 save_point_start_track

 This function saves a "bed" like file containing just a start position and strand to hdf5 storage.
 The tab delimited format describing CpG positions looks like this:
 chr1    3000826 +
 chr1    3000827 -

 Current requirement that input file is pre-sorted

=#

function save_start_point_track(genomic_store_path,point_file,track_id,chr_sizes_path;start_coord_shift=0,strand_filter_char=nothing)
    chr_sizes_dict=get_chr_sizes_dict(chr_sizes_path)

    (seq_ids,starts)=memory_read_and_parse_point_file(point_file, strand_filter_char=strand_filter_char)
    bedgraph_file_length = length(seq_ids)
    Lumberjack.info("Read bedgraph file, total length: $bedgraph_file_length")

    db=getdb(genomic_store_path)

    # initialise for the first sequence_id
    seq_id=seq_ids[1]
    seq_len  = chr_sizes_dict[seq_id]
    seq_vals = fill(int8(0),seq_len)

    # initialise structures for checking sorted data
    seen_seq_ids = Set()
    last_start_pos = 0

    for i=1:bedgraph_file_length
        if seq_id != seq_ids[i]

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

            # -- check we have sorted data by chr
            seq_len  = chr_sizes_dict[seq_id]
            seq_vals = fill(int8(0),seq_len)
            Lumberjack.info("Next track to be $seq_len in length")
        end

        #--- check that the start positions are sorted
        if last_start_pos > starts[i]
            Lumberjack.error("Unsorted position found at $i ($last_start_pos > $(starts[i])")
        end
        # update last_start_pos
        last_start_pos = starts[i]

        seq_vals[ starts[i] + start_coord_shift ] = 1
    end
    # write the final track
    write_track(db,track_id,seq_id,seq_vals)
end


#=

 save_interval_track

# this assumes each interval needs to be recorded
# as existing or not existing - use 1 or 0
# this function just records the start position at the moment.
# interval_label_type may be either "binary" uses 1 or 0 to indicate an interval covers a basepair
# or multiclass_int8 which allows an 256 possible integer values (-128 .. 127)

 Current requirement that input file is pre-sorted
=#
function save_interval_track(genomic_store_path,interval_file,track_id,chr_sizes_path;
                             start_coord_shift=0,strand_filter_char=nothing, interval_label_type="binary")

    chr_sizes_dict=get_chr_sizes_dict(chr_sizes_path)

    if interval_label_type == "binary"
        (seq_ids,starts,stops)=memory_read_and_parse_interval_file(interval_file, strand_filter_char=strand_filter_char)
    elseif interval_label_type == "multiclass_int8"
        (seq_ids,starts,stops,scores)=memory_read_and_integer_parse_bedgraph(interval_file)
    else
        Lumberjack.error("Interval Label Type not supported")
    end

    bedgraph_file_length = length(seq_ids)
    Lumberjack.info("Read bedgraph file, total length: $bedgraph_file_length")

    db=getdb(genomic_store_path)

    # initialise for the first sequence_id
    seq_id=seq_ids[1]
    seq_len  = chr_sizes_dict[seq_id]
    seq_vals = fill(int8(0),seq_len)

    # initialise structures for checking sorted data
    seen_seq_ids = Set()
    last_start_pos = 0

    for i=1:bedgraph_file_length
        if seq_id != seq_ids[i]

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

            # -- check we have sorted data by chr
            seq_len  = chr_sizes_dict[seq_id]
            seq_vals = fill(int8(0),seq_len)
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

        starting_pos=starts[i] + start_coord_shift
        stopping_pos=stops[i]  + start_coord_shift
        if interval_label_type == "binary"
            score=1
        elseif interval_label_type == "multiclass_int8"
            score = scores[i]
        end
        for pos=starting_pos:stopping_pos
            seq_vals[ pos ] = score
        end
    end
    # write the final track
    write_track(db,track_id,seq_id,seq_vals)
end



# -- Exposed methods

#=
 store_methpipe_bed_points
   exported method
=#

function store_methpipe_bed_points(filepath::String,
                                   genomic_store_path::String,
                                   track_id::String,
                                   chr_sizes_path::String;
                                   na_val=0,
                                   coord_shift=0,
                                   measurement="levels")

    gzip = isgzip(filepath) ? true : false

    if measurement == "levels"
        val_type = "float32"
        bedtype = "methpipe_cpg_bed_levels"
        na_val=float32(na_val)
    elseif measurement == "coverage"
        val_type = "int32"
        bedtype = "methpipe_cpg_bed_coverage"
        na_val=int32(na_val)
    end

    save_bed_track(genomic_store_path,filepath,track_id,chr_sizes_path,
                                start_coord_shift=coord_shift, OUT_OF_RANGE_VAL=na_val, gzip=gzip,
                                val_type=val_type,
                                bedtype=bedtype
                           )
end


#=
 store_cpg_points
  default is to filter out lines on negative strand
=#
function store_cpg_points(filepath::String,
                                   genomic_store_path::String,
                                   track_id::String,
                                   chr_sizes_path::String;
                                   strand_filter_char="-",
                                   coord_shift=0)

  save_start_point_track(genomic_store_path, filepath,track_id, chr_sizes_path,
                         start_coord_shift=coord_shift
                         )

end

#=
 store_methpipe_bed_intervals
=#
function store_methpipe_bed_intervals(filepath::String,
                                   genomic_store_path::String,
                                   track_id::String,
                                   chr_sizes_path::String;
                                   coord_shift=0)

   save_interval_track(genomic_store_path,filepath,track_id,chr_sizes_path,
                       start_coord_shift=0, interval_label_type="binary")

end

#=
 store_bedgraph
=#

function store_bedgraph_cpg( filepath::String,
                         genomic_store_path::String,
                         track_id::String,
                         chr_sizes_path::String;
                         na_val=0,
                         coord_shift=0,
                         measurement="levels")

    gzip = isgzip(filepath) ? true : false
    val_type = "float32"
    na_val=float32(na_val)
    save_bed_track(genomic_store_path,filepath,track_id,chr_sizes_path,
                   start_coord_shift=coord_shift, OUT_OF_RANGE_VAL=na_val,
                   gzip=gzip, val_type=val_type, bedtype="bedgraph_cpg" )

end
