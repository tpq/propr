# Prepare count data
counts <- read.delim("data-raw/toad.counts", row.names = 1)
caneToad.counts <- t(counts)

# Prepare labels
groups <- read.delim("data-raw/toad.groups", row.names = 1)
caneToad.groups <- t(groups)
caneToad.groups[grepl("WA", caneToad.groups)] <- "WA"
caneToad.groups[grepl("QLD", caneToad.groups)] <- "QLD"
caneToad.groups <- caneToad.groups[1,]

# Build propr object
library(propr)
keep <- apply(caneToad.counts, 2, function(x) sum(x >= 10) >= 10)
caneToad.counts <- caneToad.counts[, keep]
rho <- perb(caneToad.counts)
best.995 <- rho[">", .995]
top <- simplify(best.995)

devtools::use_data(caneToad.counts, caneToad.groups, top)
