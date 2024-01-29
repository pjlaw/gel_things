enhancers = read_tsv("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/perturbSeq/enhancers.bed", col_names=c("chr", "start", "end"))
mutspot = read_tsv("mutspot_indel.tsv")

mutspot_enh = apply(mutspot, 1, function(row) { row_chrom=row[[2]]; row_start=as.numeric(row[[3]]); row_end=as.numeric(row[[4]]); 
    overlaps = enhancers %>% filter(chr==row_chrom & start<=row_end & end >=row_start) %>% mutate(enh=glue("{chr}_{start}_{end}")) %>% pull(enh); 
    for (enhancer in overlaps) {return(paste(row_chrom, row_start, row_end, enhancer, sep="#")) } })

mutspot_enh = tibble(res=mutspot_enh)
mutspot_enh = mutspot_enh %>% separate(res, sep="#", into=c("chrom", "start", "end", "enhancer"), convert=T)
mutspot2 = mutspot %>% left_join(mutspot_enh, by=c("chrom", "start", "end"))
#remove duplicate matches in the same enhancer region, keep minimum
mutspot2 = mutspot2 %>% group_by(enhancer) %>% slice(which.min(pval))
write_tsv(mutspot2, "mutspot_indel_enh.tsv")

