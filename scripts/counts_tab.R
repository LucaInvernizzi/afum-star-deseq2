temp <- snakemake@input
myfiles <- list()

for (i in seq(from = 1, to = length(temp))){
        
        tmpfile <- read.table(file=temp[[i]],
                              header=T)
        myfiles[[i]] <- tmpfile
}

counts_col <- lapply(myfiles, function(x) x[7])
counts_tab <- do.call(cbind, counts_col) 

rownames(counts_tab) <- myfiles[[1]][,1]
colnames(counts_tab) <- sub("star_align.", "", colnames(counts_tab))
colnames(counts_tab) <- sub(".Aligned.sortedByCoord.out.bam", "", colnames(counts_tab))

write.table(counts_tab, file = snakemake@output[[1]])