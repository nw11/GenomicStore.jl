#=
   readfile.jl

   Implements read_file and read_file_coordinates methods 
   for file types definedin in filetypes.jl
 
   read_file  returns the genomic coordinates
   and associated data (such as a corresponding data value)
 
   read_file_coordinates  returns just the genomic coordinates
   without any associated data value.

   These are not general methods for reading all the data in these formats,
   they return only data that is relevant for storing at genomic coordinates.
  
   Currently the files are read into memory, in the future this might be done
   more lazy like. 

   This might eventually utilise other packages that do these
   things more efficiently.
     
=#

using Lumberjack

include( Pkg.dir("GenomicStore","src","filetypes.jl"))



"""
  read_file_coordinates ( ChrStartStrandFileMetadatai, strand_filter="-" )
  This function expects to have the first three columns to contain
  seq id,start position, and strand.
  e.g.
  chr1    3000826 +
  chr1    3000827 -
  Returns a tuple of vectors
  1. sequence_ids
  2. start_positions
"""

function read_file_coordinates( file::ChrStartStrandFileMetadata, strand_filter="-"  )
    filename=file.path
    if file.gzip
        # reset on end must be done for bgzip files
        stream=ZlibInflateInputStream(open(filename),reset_on_end=true)
    else
        stream= open(filename)
    end
    seq_ids=UTF8String[]
    starts=Int64[]
    count=0
    pass_filter_count=-0
    for line in eachline(stream)
        (seq_id,start,strand)=split(line,'\t')
        strand=chomp(strand)
        count+=1
        if count % 1000000 == 0
           Lumberjack.info("Lines read: $count")
        end
        if strand == strand_filter_str
            continue
        end
        pass_filter_count +=1
        push!(starts,parse(Int64,start) )
        push!(seq_ids,seq_id)
    end
    Lumberjack.info("Finished parsing and assigning to array.\n
    Number of lines in file: $count\n
    Number of lines passing filter: $pass_filter_count\n
    First seq_id: $(seq_ids[1])\n
    Last seq_id: $(seq_ids[end])")
    return (seq_ids,starts)
end


"""
 read_file( file::MethpipeBEDMetadata; value_column=5 )
 Reads into memory a file and parses a file in the format:
  chr1    3000826 +       CG      1       17
  chr1    3001006 +       CG      0.8     10
 extracting out either the 5th or 6th column (the coverage column)
 returns tuple of 4 arrays
 1. seq_ids
 2. starts
 3. stops
 4. levels/coverage
"""
function read_file(file::MethpipeBEDMetadata; value_column=5)
    
    filename=file.path
    if file.gzip
        # reset on end must be done for bgzip files
        stream=ZlibInflateInputStream(open(filename),reset_on_end=true)
    else
        stream= open(filename)
    end
    seq_ids=UTF8String[]
    starts=Int64[]
    stops=Int64[]
    if  value_column == 5
        scores=Float32[]
        T=Float32
    elseif value_column == 6
        scores=Int32[]
        T=Int32
    else
        Lumberjack.error("invalid column number (valid: 5 or 6)")
    end
    count=0
    for line in eachline(stream)
        # This splits the line into the following:
        #(seq_id,start,strand,context,score,coverage) = split(line,'\t')
        split_line=split(line,'\t')
        count+=1
        if count % 1000000 == 0
           Lumberjack.info("Lines read: $count")
        end
        push!(seq_ids,split_line[1])
        start_int=parse(Int64,split_line[2])
        push!(starts,start_int)
        push!(stops,start_int+1)
        push!(scores,parse(T,split_line[value_column]))
    end
    Lumberjack.info("Finished parsing and assigning to array.\n
    Number of lines in file: $count\n
    First seq_id: $(seq_ids[1])\n
    Last seq_id: $(seq_ids[end])")
    return (seq_ids,starts,stops,scores)
end

"""
read_file_coordinates(file::BEDFileMetadatai, only_start=false,strand_filter_char=nothing)
 This function expects a file to have the first three columns to contain
 a sequence id, start position, and stop/end position
 Initially used to parse files like below:
 chr1    3670331 3672761 HYPO0   184     +
 chr1    3993551 3993737 HYPO1   5       +
 Returns a tuple of vectors
 1. sequence_ids
 2. start_positions
 3. stop_positions
"""

function read_file_coordinates(file::BEDMetadata; only_start=false,strand_filter_char=nothing)
    if file.gzip
        # reset on end must be done for bgzip files
        stream=ZlibInflateInputStream(open(filename),reset_on_end=true)
    else
        stream= open(filename)
    end
    seq_ids=UTF8String[]
    starts=Int64[]
    stops=Int64[]
    count=0
    for line in eachline(stream)
        # This splits the line into the following:
        #(seq_id,start,stop,id,cpg_count,strand) = split(line,'\t')
        split_line=split(line,'\t')
        count+=1
        if count % 1000000 == 0
           Lumberjack.info("Lines read: $count")
        end
        push!(seq_ids,split_line[1])
        push!(starts,parse(Int64,split_line[2]))
        push!(stops, parse(Int64,split_line[3]))
    end
    Lumberjack.info("Finished parsing and assigning to array.\n
    Number of lines in file: $count\n
    First seq_id: $(seq_ids[1])\n
    Last seq_id: $(seq_ids[end])")
    if only_start
        return (seq_ids,starts)
    end
    return  (seq_ids,starts,stops)
end


"""
  read_file(file::BedGraphMetadata)
  This function expects to have the first three columns to contain
  seq id,start position, and stop/end postion.
  e.g.
  chr1 3670331 3670333 0.05
  Returns a tuple of vectors
  1. sequence_ids
  2. start_positions
  3. stop_positions
  4. score
"""
function read_file(file::BedGraphMetadata)
    filename=file.path
    if file.gzip
        # reset on end must be done for bgzip files
        stream=ZlibInflateInputStream(open(filename),reset_on_end=true)
    else
        stream= open(filename)
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













