using HDF5,JLD
using GZip
using Lumberjack

import HDF5.dump

function dump(io::IO, x::Union(HDF5File, HDF5Group), n::Int, indent, display_num=10)
    println(typeof(x), " len ", length(x))
    if n > 0
        i = 1
        for k in names(x)
            print(io, indent, "  ", k, ": ")
            v = o_open(x, k)
            dump(io, v, n - 1, string(indent, "  "))
            close(v)
            if i > display_num
                println(io, indent, "  ...")
                break
            end
            i += 1
        end
    end
end

function store_info_dict(genomic_store_path)
  d=Dict{String, Array{Symbol} }()
  fid=h5open(genomic_store_path,"r")
  for seq in names(fid)
     track_array=names( fid[seq] )
     d[seq]=track_array
  end
  return d
end

function tracks(genomic_store_path)
  # query tracks
end

function store_info_str(genomic_store_path)
  buf=IOBuffer()
  fid=h5open(genomic_store_path,"r")
  dump(buf,fid,10,"\n",20)
  close(fid)
  return takebuf_string(buf)
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

function read_gzip_file(filename)
   gzio=gzopen(filename,"r")
   line_num =0
   lines=ASCIIString[]
   for line in eachline(gzio)
        line=chomp(line)
        push!(lines,line)
        line_num +=1
        if line_num % 1000000 == 0
            Lumberjack.info("full on speed like road runner $line_num already!")
        end
    end
    #pop!(lines) # clears out the final row which due to a bug in GZlib gives and empty string
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


# This file has first three columns chr start strand or something else
# but we ignore them.
function memory_read_and_integer_parse_point_file(filename;strand_filter_char=nothing)
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

# -
# - This file has first three columns chr start stop
#-  it could have others
# - but we ignore them.
# -
function memory_read_and_integer_parse_interval_file(filename; only_start=false,strand_filter_char=nothing)
    io = open(filename)
    Lumberjack.info("reading all")
    file=readall(io)
    Lumberjack.info("finished reading all")
    Lumberjack.info("split line")
    lines=split(file,'\n')
    Lumberjack.info("done split")
    # now we know what size thigns should be
    # get size and two columns and then
    # lines -1 as it splits the final row
    num_lines=length(lines)-1
    seq_ids = fill("",num_lines)
    starts=fill( 0, num_lines)
    if !only_start
       stops=fill( 0,  num_lines )
    end
    Lumberjack.info("parse and assign to array")
    #parseint takes for ever, so annoying about this.
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

# the integer_parse part of this method name needs to be renamed
function memory_read_and_integer_parse_methpipe_cg_bed(filename;gzip=false)
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
    #parseint takes for ever, and that is what is so annoying about this.
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

function memory_read_and_integer_parse_methpipe_cg_bed_coverage(filename;gzip=false)
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

function memory_read_and_integer_parse_bedgraph(filename;gzip=false)
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
        if idx % 1000000 == 0
            println(idx)
        end
    end
    Lumberjack.info("finished parsing integer and assigning to array")
    return (seq_ids,starts,stops,scores)
end

# writes a track - that must have a sequence id
function write_track(file_handle,seq_id,track_id,values)
    Lumberjack.info("Writing $seq_id")
    write(file_handle, bytestring("$seq_id/$track_id"), values)
    Lumberjack.info("Finished writing $seq_id")
end


#- We have different tracks for chromosomes so the files do not have
#- to be consistently sorted by chromosome

function save_bedgraph_track(genomic_store_path,input_file,track_id,chr_sizes_path;
                             start_coord_shift=0,
                             OUT_OF_RANGE_VAL=0.0,
                             gzip=false,
                             val_type="float32",
                             bedtype="bedgraph")

    return save_bed_track(genomic_store_path,input_file,track_id,chr_sizes_path,
                             start_coord_shift=start_coord_shift,
                             OUT_OF_RANGE_VAL=OUT_OF_RANGE_VAL,
                             gzip=gzip,
                             val_type=val_type,
                             bedtype="bedgraph")
end

function save_bed_track(genomic_store_path,input_file,track_id,chr_sizes_path;
                             start_coord_shift=0,
                             OUT_OF_RANGE_VAL=0.0,
                             gzip=false,
                             val_type="float32",
                             bedtype="bedgraph")
    chr_sizes_dict=get_chr_sizes_dict(chr_sizes_path)
    #--- sort in temporary dir
    #--- load
    if bedtype == "bedgraph"
        (seq_ids,starts,stops,scores)=memory_read_and_integer_parse_bedgraph(input_file,gzip=gzip)
    elseif bedtype == "methpipe_cpg_bed"
        (seq_ids,starts,stops,scores)=memory_read_and_integer_parse_methpipe_cg_bed(input_file,gzip=gzip)
    elseif bedtype == "methpipe_cpg_bed_coverage"
        (seq_ids,starts,stops,scores)=memory_read_and_integer_parse_methpipe_cg_bed_coverage(input_file,gzip=gzip)
    else
         Lumberjack.error("Invalid bedtype specified $bedtype")
    end
    bedgraph_file_length = length(seq_ids)
    Lumberjack.info("Read bedgraph file, total length: $bedgraph_file_length")
    if !isfile(genomic_store_path)
        fid=jldopen(genomic_store_path,"w",compress=true)
    else
        fid=jldopen(genomic_store_path,"r+",compress=true)
    end
    seq_id="FIRSTROW"
    if seq_ids[1] == "FIRSTROW"
       seq_id="FIRSTROW1"
    end
    seq_vals=nothing
    #seen_seq_id = Dict()
    #curr_seq_id=""
    for i=1:bedgraph_file_length
        if seq_id != seq_ids[i]
            if i != 1
                write_track(fid,seq_id,track_id,seq_vals)
            end
            seq_id   = seq_ids[i]
            Lumberjack.info("Next sequence name: $seq_id")

            # -- check the id is valid --

            if !haskey(chr_sizes_dict,seq_id)
                Lumberjack.error("Invalid ID: ( $seq_id ), this file only saved up to line $(i-1)")
            end

            # -- TO DO check we have sorted data by chr
            #if curr_seq_id != seq_id
            #    seen_seq_id[curr_seq_id]="DONE"
            #    curr_seq_id = seq_id
            #    seen_seq_id[curr_seq_id]="DOING"
            #end

            seq_len  = chr_sizes_dict[seq_id]
            seq_vals = fill(OUT_OF_RANGE_VAL,seq_len)
            Lumberjack.info("Next track to be $seq_len in length")
        end

        #--- will slow things down a little but check that the id is sorted, bare in mind a chr change.
        if val_type == "float32"
            try
                seq_vals[ starts[i] + start_coord_shift ] = float32( scores[i] )
            catch y
                println("processing this line:\n")
                curr_row = join([ seq_ids[i], starts[i],stops[i], scores[i] ],"\t" )
                println(curr_row)
                println("shift: $start_coord_shift")
                error(y)
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
    write_track(fid,seq_id,track_id,seq_vals)
    Lumberjack.info("Written track $seq_id")
    close(fid)
end

function get_track_interval( genomic_store_path, id, seq_id, start_pos, stop_pos )
    fid=h5open(genomic_store_path,"r")
    if has(fid,"$seq_id/$id")
        if start_pos == 0 && stop_pos == 0
            track_slice=h5read(genomic_store_path,bytestring("$seq_id/$id"))
        else
            track_slice=h5read(genomic_store_path,bytestring("$seq_id/$id"),(start_pos:stop_pos,))
        end
        close(fid)
        return track_slice
    end
    close(fid)
    return []
end

# this assumes each interval needs to be recorded
# as existing of not existing - use 1 or 0
# this function just records the start position at the moment.
function save_point_start_track(genomic_store_path,point_file,track_id,chr_sizes_path;start_coord_shift=0,strand_filter_char=nothing)
    chr_sizes_dict=get_chr_sizes_dict(chr_sizes_path)
    #--- sort in temporary dir
    #--- load
    (seq_ids,starts)=memory_read_and_integer_parse_point_file(point_file, strand_filter_char=strand_filter_char)
    bedgraph_file_length = length(seq_ids)
    Lumberjack.info("Read bedgraph file, total length: $bedgraph_file_length")

    if !isfile(genomic_store_path)
        fid=jldopen(genomic_store_path,"w",compress=true)
    else
        fid=jldopen(genomic_store_path,"r+",compress=true)
    end
    seq_id="FIRSTROW"
    if seq_ids[1] == "FIRSTROW"
       seq_id="FIRSTROW1"
    end
    seq_vals=nothing
    for i=1:bedgraph_file_length
        if seq_id != seq_ids[i]
            if i != 1
                write_track(fid,seq_id,track_id,seq_vals)
            end
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
        #---TODO will slow things down a little but check that the id is sorted baring in mind a chr change.
        seq_vals[ starts[i] + start_coord_shift ] = 1
    end
    # write the final track
    write_track(fid,seq_id,track_id,seq_vals)
    close(fid)
end

# this assumes each interval needs to be recorded
# as existing or not existing - use 1 or 0
# this function just records the start position at the moment.
# interval_label_type may be either "binary" uses 1 or 0 to indicate an interval covers a basepair
# or multiclass_int8 which allows an 256 possible integer values (-128 .. 127)
function save_interval_track(genomic_store_path,interval_file,track_id,chr_sizes_path;
                             start_coord_shift=0,strand_filter_char=nothing, interval_label_type="binary")

    chr_sizes_dict=get_chr_sizes_dict(chr_sizes_path)
    #--- sort in temporary dir
    #--- load
    if interval_label_type == "binary"
        (seq_ids,starts,stops)=memory_read_and_integer_parse_interval_file(interval_file, strand_filter_char=strand_filter_char)
    elseif interval_label_type == "multiclass_int8"
        (seq_ids,starts,stops,scores)=memory_read_and_integer_parse_bedgraph(interval_file)
    else
        Lumberjack.error("Interval Label Type not supported")
    end

    bedgraph_file_length = length(seq_ids)
    Lumberjack.info("Read bedgraph file, total length: $bedgraph_file_length")

    if !isfile(genomic_store_path)
        fid=jldopen(genomic_store_path,"w",compress=true)
    else
        fid=jldopen(genomic_store_path,"r+",compress=true)
    end
    seq_id="FIRSTROW"
    if seq_ids[1] == "FIRSTROW"
       seq_id="FIRSTROW1"
    end
    seq_vals=nothing
    for i=1:bedgraph_file_length
        if seq_id != seq_ids[i]
            if i != 1
                write_track(fid,seq_id,track_id,seq_vals)
            end
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
        #---TODO will slow things down a little but check that the id is sorted baring in mind a chr change.
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
    write_track(fid,seq_id,track_id,seq_vals)
    close(fid)
end
