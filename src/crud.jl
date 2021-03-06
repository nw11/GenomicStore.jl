using HDF5,JLD
using Docile
using Lumberjack

import JLD.JldFile
# if for whatever reason the backend needs to change
# then adding write/delete/update methods here with a
# different backend / file_handle type would be added

abstract StoreBackend

# possibly allows for filehandle in here and read/write privilidge?
# this is to avoid closing and opening handles all the time for
# many individual writes? bad idea? probably so don't put in for now.
# writes should generally be done more batch like - everything modified
# in memory as much as possible - then written to a track
# We want to try different forms of hdf5 e.g. compressed and mmapped, these can't be mixed

type JldStoreBackend <: StoreBackend
    path::AbstractString
    description::AbstractString
    compressed::Bool
    mmapped::Bool
end

# hdf5-blosc, hdf5-gzip hdf5-mmap, hdf5-uncompressed (uncompressed)
function getdb(dbpath::AbstractString; dbtype = "hdf5", description="" )
    # check for suffix
    if dbtype == "hdf5" # default with compression
        return JldStoreBackend(dbpath,description,true,false)
    elseif dbtype == "hdf5-blosc"
        return JldStoreBackend(dbpath,description,true,false)
    elseif dbtype == "hdf5-gzip"
        return JldStoreBackend(dbpath,description,true,false)
    elseif dbtype == "hdf5-mmap"
        return JldStoreBackend(dbpath,description,false,true)
    elseif dbtype == "hdf5-uncompressed"
       return JldStoreBackend(dbpath,description,false,false)
    end
end

"""
Sets up filehandle and writes values to a genomic track

# Arguments
   * `file_handle::JldBackend`: backend store type
   * `seq_id::String`: identifying name of genomic sequence
   * `track_id::String`: identifying name of track
   * `values::Number`
   * `overwrite::Boolean` : keyword argument indicating whether to overwrite the track

# Return
  * JldBackend
"""
function write_track{T <: Number}(backend_store::JldStoreBackend,
                                  track_id::AbstractString, seq_id::AbstractString,
                                  values::Vector{T}; overwrite=false )
    genomic_store_path = backend_store.path
    fid=nothing
    if !isfile(genomic_store_path)
        fid=jldopen(genomic_store_path,"w",compress=backend_store.compressed)
    else
        fid=jldopen(genomic_store_path,"r+",compress=backend_store.compressed)
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

function delete_track!(backend_store::JldStoreBackend,track_id::AbstractString, seq_id::AbstractString)
    genomic_store_path = backend_store.path
    fid=jldopen(genomic_store_path,"r+",compress=backend_store.compressed)
    Lumberjack.info("delete $seq_id")
    _delete_track!(fid, track_id, seq_id)
    Lumberjack.info("Finished removing $seq_id")
    close(fid)
end

function read_track(backend_store::JldStoreBackend,
                                 track_id::AbstractString,  seq_id::AbstractString,
                                 start_pos::Integer, stop_pos::Integer)
    genomic_store_path = backend_store.path
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


function update_track!{T <: Number}(backend_store::JldStoreBackend,
                                  track_id::AbstractString,  seq_id::AbstractString,
                                  start_pos::AbstractString, stop_pos::AbstractString,
                                  values::Vector{T})
    Lumberjack.info("Updating $seq_id")
    genomic_store_path = backend_store.path
    fid=jldopen(genomic_store_path,"r+",compress=backend_store.compressed)
    _udpate_track!(fid,track_id , seq_id , start_pos,stop_pos,values)
    Lumberjack.info("Finished updating $seq_id")
    close(fid)
end

function _istrack(backend_store::JldStoreBackend,track_id::AbstractString)
    genomic_store_path = backend_store.path
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

"""
  _hastrack( backend_store::JldStoreBackend,track_id::AbstractString )

  A track is a set of intervals or points associated with a sequence.
  So we look through all sequence names for evidence that a track exists.
  Should replace _istrack() _istrack will be deprecated in favour of this
  method.
"""

function _hastrack( backend_store::JldStoreBackend,track_id::AbstractString  )
    # _hastrack
    _istrack(backend_store,track_id)
end

"""
   _hasseqtrack(  backend_store::JldStoreBackend,track_id::AbstractString, seq_id::AbstractString )
   A track is a set of intervals or points associated with a sequence.
   This checks the existance of a seq_id/track_id path.
"""

function _hasseqtrack(  backend_store::JldStoreBackend,track_id::AbstractString, seq_id::AbstractString )
    genomic_store_path = backend_store.path
    fid=nothing
    try
         fid=h5open(genomic_store_path,"r")
    catch e
        Lumberjack.error("$e\n Trying to open $genomic_store_path")
    end
    if has(fid,"$seq_id/$track_id")
       close(fid)
       return true
    end
    close(fid)
    return false
end



# backend store - hdf5
# call these _write_track_hdf5  (crud-hdf5.jl?)
# catch errors here to allow logging of errors
# with specific messages relevant to these functions

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

"""
    _list_tracks

    Data organised as:

    "chr1"
     -> TRACKID1
     -> TRACKID2

    "chr2"
      -> TRACKID1
      -> TRACKID2

     Check whether a track_id exists under a sequence id.


"""
function _list_sequence_ids(backend_store::JldStoreBackend,track_id::AbstractString)
    genomic_store_path = backend_store.path
    fid=nothing
    try
        fid=h5open(genomic_store_path,"r")
    catch e
        Lumberjack.error("$e\n Trying to open $genomic_store_path")
    end
    seq_ids=[]
    for n in names(fid)
        if has(fid,"$n/$track_id")
            push!(seq_ids,n)
        end
    end
   close(fid)
   return seq_ids
end


function _list_tracks(backend_store::JldStoreBackend)
end


# a search function - e.g. indices or values  > than some value or matching a value .. prob not at this level - better done somewhere else
# specialise could be a missing value thing
# missing value tracks are very similar to storing intervals as 1 or 0
# this should all be done at a higher level than crud.
