library(propriety)

data(iris)
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

permrange <- 2^(1:12)
cutrange <- seq(.95, 1, .01)
res <- sapply(permrange,
              function(perms){

                # Test original permutation method
                set.seed(1)
                theta <- propriety:::calculateTheta(counts, group)
                ptheta <- propriety:::permuteTheta(counts, group, p = perms)
                pt <- propriety:::calculateFDR(theta, ptheta, cutoff = cutrange)

                # Test prime permutation method
                set.seed(1)
                pd <- propd(counts, group, p = perms, cutoff = cutrange)

                # Calculate percent difference
                diff <- (pt$FDR - pd@fdr$FDR) / ((pt$FDR + pd@fdr$FDR) / 2)
                diff[is.na(diff)] <- 0
                return(diff)
              }
)

colnames(res) <- permrange
rownames(res) <- cutrange
df <- reshape2::melt(res)
colnames(df) <- c("cutoff", "perms", "diff")
df$cutoff <- as.factor(df$cutoff)
library(ggplot2)
png("data-raw/compare-permute.png")
ggplot(df, aes(perms, diff, color = cutoff)) + geom_point() + geom_line() + theme_bw() +
  xlab("Number of Permutations") + ylab("Percent Difference in FDR between Old and New Methods") +
  ggtitle("The Effect of Faster Permutations on FDR")
dev.off()
