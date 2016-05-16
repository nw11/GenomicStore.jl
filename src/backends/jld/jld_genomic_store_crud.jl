#=
jld_genomic_store_crud.jl

Hold basic create, read, write, update functions using the JLD
module

 =#

using HDF5,JLD
using Lumberjack
import JLD.JldFile


"""
JldGenomicStore

Metadata for a GenomicStore that uses the JLD interface to HDF5
to store genomic tracks
"""
type JldGenomicStore
    path::AbstractString
    description::AbstractString
    compressed::Bool
    mmapped::Bool
end


# hdf5-blosc, hdf5-gzip hdf5-mmap, hdf5-uncompressed (uncompressed)
function getdb(dbpath::AbstractString; dbtype = "hdf5", description="" )
    # check for suffix
    if dbtype == "hdf5" # default with compression
        return JldGenomicStore(dbpath,description,true,false)
    elseif dbtype == "hdf5-blosc"
        return JldGenomicStore(dbpath,description,true,false)
    elseif dbtype == "hdf5-gzip"
        return JldGenomicStore(dbpath,description,true,false)
    elseif dbtype == "hdf5-mmap"
        return JldGenomicStore(dbpath,description,false,true)
    elseif dbtype == "hdf5-uncompressed"
       return JldGenomicStore(dbpath,description,false,false)
    end
end


"""
write_track

Sets up filehandle and writes values to a genomic track

# Arguments
   * `genomic_store::JldGenomicStore`: backend store type
   * `seq_id::AbstractString`: identifying name of genomic sequence
   * `track_id::AbstractString`: identifying name of track
   * `values::Number`
   * `overwrite::Boolean` : keyword argument indicating whether to overwrite the track
# Return
  * JldBackend
"""
function write_track{T <: Number}(genomic_store::JldGenomicStore,
                                  track_id::AbstractString, seq_id::AbstractString,
                                  values::Vector{T}; overwrite=false )
    genomic_store_path = genomic_store.path
    fid=nothing
    if !isfile(genomic_store_path)
        fid=jldopen(genomic_store_path,"w",compress=genomic_store.compressed)
    else
        fid=jldopen(genomic_store_path,"r+",compress=genomic_store.compressed)
    end

    if !has(fid.plain,"$seq_id/$track_id") & !overwrite
        Lumberjack.info("Writing $seq_id")
        _write_track(fid, track_id, seq_id,values)
        Lumberjack.info("Finished writing $seq_id")
    else
        Lumberjack.info("Track: $seq_id/$track_id exists ... skipping")
    end
    close(fid)
end



function delete_track!(genomic_store::JldGenomicStore,track_id::AbstractString, seq_id::AbstractString)
    genomic_store_path = genomic_store.path
    fid=jldopen(genomic_store_path,"r+",compress=genomic_store.compressed)
    Lumberjack.info("delete $seq_id")
    _delete_track!(fid, track_id, seq_id)
    Lumberjack.info("Finished removing $seq_id")
    close(fid)
end

"""
read_track

Sets up filehandle and writes values to a genomic track

# Arguments
   * `genomic_store::JldGenomicStore`: backend store type
   * `track_id::AbstractString`: identifying name of track
   * `seq_id::AbstractString`: identifying name of genomic sequence
   * `start_pos::Integer`
   * `stop_pos::Integer`
# Return
  * Array
"""


function read_track(genomic_store::JldGenomicStore,
                                 track_id::AbstractString,  seq_id::AbstractString,
                                 start_pos::Integer, stop_pos::Integer)
    genomic_store_path = genomic_store.path
    fid=nothing
    try
        fid=h5open(genomic_store_path,"r")
    catch e
        Lumberjack.error("$e\n Trying to open $genomic_store_path")
    end
    if has(fid,"$seq_id/$track_id")
        if start_pos == 0 && stop_pos == 0
            track_slice=h5read(genomic_store_path,bytestring("$seq_id/$track_id"))
        else
            track_slice=h5read(genomic_store_path,bytestring("$seq_id/$track_id"),(start_pos:stop_pos,))
        end
        close(fid)
        return track_slice
    end
    close(fid)
    return []
end


function update_track!{T <: Number}(genomic_store::JldGenomicStore,
                                  track_id::AbstractString,  seq_id::AbstractString,
                                  start_pos::AbstractString, stop_pos::AbstractString,
                                  values::Vector{T})
    Lumberjack.info("Updating $seq_id")
    genomic_store_path = genomic_store.path
    fid=jldopen(genomic_store_path,"r+",compress=genomic_store.compressed)
    _udpate_track!(fid,track_id , seq_id , start_pos,stop_pos,values)
    Lumberjack.info("Finished updating $seq_id")
    close(fid)
end

function _hastrack(genomic_store::JldGenomicStore,track_id::AbstractString)
    genomic_store_path = genomic_store.path
    fid=nothing
    try
        fid=h5open(genomic_store_path,"r")
    catch e
        Lumberjack.error("$e\n Trying to open $genomic_store_path")
    end
    for n in names(fid)
        if has(fid,"$n/$track_id")
            close(fid)
            return true
        end
    end
    close(fid)
    return false
end

function _write_track{T <: Number}(file_handle::JldFile, track_id::AbstractString, seq_id::AbstractString, values::Vector{T} )
    try
        write(file_handle, bytestring("$seq_id/$track_id"), values)
    catch err
        Lumberjack.error("Problem writing track: $err")
    end
end

function _delete_track!(file_handle::JldFile,track_id::AbstractString, seq_id::AbstractString)
    try
        delete!(file_handle,bytestring("$seq_id/$track_id"))
    catch err
        Lumberjack.error("Problem removing track: $err")
    end
end

function _update_track!{T <: Number}(file_handle::JldFile, track_id::AbstractString, seq_id::AbstractString,
                                      start_pos::Integer,stop_pos::Integer,values::Vector{T} )
    dset = file_handle["$seq_id/$track_id"]
    if stop_pos - start_pos != length(values)
        error("coordinates do not match length of values")
    end
    try
        dset[start_pos:stop_pos] = values
    catch err
        Lumberjack.error("Problem updating track: $err")
    end
end

function _list_tracks()
   # to be implemented
end

