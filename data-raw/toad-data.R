# Prepare count data
counts <- read.delim("data-raw/toad.counts", row.names = 1)
caneToad.counts <- t(counts)

# Prepare labels
groups <- read.delim("data-raw/toad.groups", row.names = 1)
caneToad.groups <- t(groups)
caneToad.groups[grepl("WA", caneToad.groups)] <- "WA"
caneToad.groups[grepl("QLD", caneToad.groups)] <- "QLD"
caneToad.groups <- caneToad.groups[1,]

# Use data
devtools::use_data(caneToad.counts, caneToad.groups)
