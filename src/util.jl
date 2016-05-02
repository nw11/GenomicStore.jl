
using Lumberjack
using DataFrames

function get_chr_sizes_dict(chrom_sizes_path)
    chrom_sizes=readdlm(chrom_sizes_path)
    chrom_sizes_dict = Dict{ASCIIString,Int64}()
    size(chrom_sizes,1)
    for i=1:size(chrom_sizes,1)
        chrom_sizes_dict[ chrom_sizes[i,1] ] = chrom_sizes[i,2]
    end
    return chrom_sizes_dict
end

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


function isgzip(filepath)
    (path,ext) = splitext(filepath)
    gzip = (ext==".gz" || ext==".gzip" ) ? true : false
    return gzip
end
