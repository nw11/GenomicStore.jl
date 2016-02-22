function get_chr_sizes_dict(chrom_sizes_path)
    chrom_sizes=readdlm(chrom_sizes_path)
    chrom_sizes_dict = Dict{ASCIIString,Int64}()
    size(chrom_sizes,1)
    for i=1:size(chrom_sizes,1)
        chrom_sizes_dict[ chrom_sizes[i,1] ] = chrom_sizes[i,2]
    end
    return chrom_sizes_dict
end
