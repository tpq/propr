RNA.seq           <- read.csv("./data-raw/data/RNA.seq.csv",        header=T, skip=1)
microarray        <- read.csv("./data-raw/data/microarray.csv",     header=T, skip=3)
names(microarray) <- sub("T", "timepoint", names(microarray))
go                <- read.csv("./data-raw/data/Complex_annotation", header=T, sep="\t")
names(go)[6]      <- "Systematic.name"

# Average the sums all the copies-per-cell ("cpc") counts for MM1 and MM2,
# treating any NAs as 0
tmp <- data.frame(Systematic.name=RNA.seq$Systematic.name,
                  RNA.seq=rowSums(
                    RNA.seq[,grep("MM[12].*cpc.*", names(RNA.seq))],
                    na.rm=TRUE
                    )/2
                  )

# Drop any mRNAs that have a zero count in the RNA-seq
tmp <- subset(tmp, RNA.seq > 0)

# Do an inner join of Abs and the microarray data based on the Systematic names
tmp <- merge(tmp, microarray, by="Systematic.name")

# Now use the relative abundances at each microarray timepoint to multiply
# the initial mRNA copies per cell. Remove any rows that contain NAs
multipliers   <- as.matrix(tmp[, grep("timepoint", names(tmp))])
Abs           <- data.frame(tmp$RNA.seq * multipliers)
rownames(Abs) <- tmp[,"Systematic.name"]
Abs           <- na.omit(Abs)
Abs.t         <- as.data.frame(t(Abs))

Rel   <- sweep(Abs,2,colSums(Abs, na.rm=TRUE),"/")
Rel.t <- as.data.frame(t(Rel))

marg.abs <- Abs.t

devtools::use_data(marg.abs)
