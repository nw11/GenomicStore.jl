# This file provides functions for reading files into and hdf5a db
using Compat
using GZip
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

#=

save_bed_track

 This function saves "bed" like files to hdf5 storage. The name should more accurately
 be "save_bedlike_track"

 It handles saving values specified in

 a) bedgraph format
 b) methpipe_cpg_bed
     i) extracting the methylation levels value
     ii) extracting the coverage levels
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
    if bedtype == "bedgraph"
        (seq_ids,starts,stops,scores)=memory_read_and_parse_bedgraph(input_file,gzip=gzip)
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

    seq_id="FIRSTROW"
    #if seq_ids[1] == "FIRSTROW"
    #   seq_id="FIRSTROW1"
    #end
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
    write_track(db,track_id,seq_id,seq_vals)
    Lumberjack.info("Written track $seq_id")
end




function isgzip(filepath)
    (path,ext) = splitext(filepath)
    gzip = (ext==".gz" || ext==".gzip" ) ? true : false
    return gzip
end


#=
 store_methpipe_bed_points
   exported method
=#

function store_methpipe_bed_points(filepath::ASCIIString,
                                   genomic_store_path::ASCIIString,
                                   track_id::ASCIIString,
                                   chr_sizes_path::ASCIIString;
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

    save_bed_track(genomic_store_path,filepath,track_id,chr_sizes_path;
                                start_coord_shift=coord_shift, OUT_OF_RANGE_VAL=na_val, gzip=gzip,
                                val_type=val_type,
                                bedtype=bedtype
                           )
end


function store_cpg_points(filepath::ASCIIString,
                                   genomic_store_path::ASCIIString,
                                   track_id::ASCIIString,
                                   chr_sizes_path::ASCIIString;
                                   na_val=0,
                                   coord_shift=0)

end

#=
 store_methpipe_bed_intervals
=#
function store_methpipe_bed_intervals(filepath::ASCIIString,
                                   genomic_store_path::ASCIIString,
                                   track_id::ASCIIString,
                                   chr_sizes_path::ASCIIString;
                                   na_val=0,
                                   coord_shift=0)

end


