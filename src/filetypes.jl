
abstract FileMetadata

abstract IntervalFileMetadata

abstract PointFileMetadata


"""
 MethpipeBEDMetadata

 Metadata associated with a Methpipe version of a "BED" file,
 which looks like this: 
 
 CHECK - is this methpipe or bismark - it's Bismark I think because
 it is for CG lvels.
 
  chr1    3000826 +       CG      1       17
  chr1    3001006 +       CG      0.8     10
"""

type MethpipeBEDMetadata
    path::UTF8String
    gzip::Bool
end

"""
   BEDMetadata

   Metadata associated with a BED file 
"""

type BEDMetadata
    path::UTF8String
    gzip::Bool
end

"""
   BedGraphMetadata

   Metadata associated with a Bedgraph file
"""

type BedGraphMetadata
    path::UTF8String 
    gzip::Bool
end

"""
   IntervalGeneralFileMetadata

   Metadata associated with any file that starts with the 
   fields  seq_id,start,end fields then followed by any 
   other number of tab delimited fields 
"""
type IntervalGeneralFileMetadata
    path::UTF8String
end


"""
   PointGeneralFileMetadata

   Metadata associated with any file that starts with the 
   fields  seq_id,start fields then followed by any other 
   number of tab delimited fields 
"""
type PointGeneralFileMetadata
     path::UTF8String
end


"""
  ChrStartStrandFileMetadata

  Metadata associated with a file which starts with
  the fields seq_id, start,strand.
 
  Yes. It has got to this. We just end up describing 
  bioinformatic files by what is in each column. I don't 
  know what else to call this, it's there, I use it because some
  program dumps it out. 

"""

type ChrStartStrandFileMetadata
    path::UTF8String
    gzip::Bool
end



#=
 Notes
  PointGeneralFile is subsumed by IntervalGeneralFile, but we will see how that goes.
  It looks like it might be useful to separate these two forms.
=#
