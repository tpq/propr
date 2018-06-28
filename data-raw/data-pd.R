library(propr)
data(caneToad.counts)
data(caneToad.groups)

# Remove features with too many low counts or any zero counts
keep <- apply(caneToad.counts, 2, function(x) sum(x >= 40) >= 20 & all(x != 0))
caneToad.sub <- caneToad.counts[, keep]

pd <- propd(caneToad.sub, caneToad.groups) # propd object last built during propr merge
pd.dd <- setDisjointed(pd)
pd.ee <- setEmergent(pd)

# Subset top 1000 theta values only
pd.dd@results <- pd.dd@results[pd.dd@results$theta < sort(pd.dd@results$theta)[1001], ]
pd.ee@results <- pd.ee@results[pd.ee@results$theta < sort(pd.ee@results$theta)[1001], ]

# Rebuild pd.d with fewer counts
counts.d <- pd.dd@counts[, sort(union(pd.dd@results$Partner, pd.dd@results$Pair))]
pd.d <- setDisjointed(propd(counts.d, caneToad.groups))
pd.d@results <- pd.d@results[pd.d@results$theta < sort(pd.d@results$theta)[1001], ]#.
identical(round(pd.d@results$theta, 14), round(pd.dd@results$theta, 14)) # check subset

# Rebuild pd.e with fewer counts
counts.e <- pd.ee@counts[, sort(union(pd.ee@results$Partner, pd.ee@results$Pair))]
pd.e <- setEmergent(propd(counts.e, caneToad.groups))
pd.e@results <- pd.e@results[pd.e@results$theta < sort(pd.e@results$theta)[1001], ]#.
identical(round(pd.e@results$theta, 14), round(pd.ee@results$theta, 14)) # check subset

devtools::use_data(pd.d, pd.e)
