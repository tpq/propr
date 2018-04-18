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
pd.dd@theta <- pd.dd@theta[pd.dd@theta$theta < sort(pd.dd@theta$theta)[1001], ]
pd.ee@theta <- pd.ee@theta[pd.ee@theta$theta < sort(pd.ee@theta$theta)[1001], ]

# Rebuild pd.d with fewer counts
counts.d <- pd.dd@counts[, sort(union(pd.dd@theta$Partner, pd.dd@theta$Pair))]
pd.d <- setDisjointed(propd(counts.d, caneToad.groups))
pd.d@theta <- pd.d@theta[pd.d@theta$theta < sort(pd.d@theta$theta)[1001], ]#.
identical(round(pd.d@theta$theta, 14), round(pd.dd@theta$theta, 14)) # check subset

# Rebuild pd.e with fewer counts
counts.e <- pd.ee@counts[, sort(union(pd.ee@theta$Partner, pd.ee@theta$Pair))]
pd.e <- setEmergent(propd(counts.e, caneToad.groups))
pd.e@theta <- pd.e@theta[pd.e@theta$theta < sort(pd.e@theta$theta)[1001], ]#.
identical(round(pd.e@theta$theta, 14), round(pd.ee@theta$theta, 14)) # check subset

devtools::use_data(pd.d, pd.e)
