M <- read.table("data-raw/gene_matrix.txt")
counts <- M
group = c(rep("A", 54), rep("B", 98-54))

time <- c(3, 10, 50, 100, 300, 500)
res <-
  lapply(time,
         function(k){

           ct <- counts[, 1:k]

           a <- microbenchmark::microbenchmark(
             propriety:::permuteTheta_old(ct, group, p = 20), times = 1L
           )
           b <- microbenchmark::microbenchmark(
             propriety:::permuteTheta_naive(ct, group, p = 20), times = 1L
           )
           c <- microbenchmark::microbenchmark(
             propriety:::permuteTheta(ct, group, p = 20), times = 1L
           )
           d <- microbenchmark::microbenchmark(
             propriety:::permuteTheta_prime(ct, group, p = 20), times = 1L
           )
           e <- microbenchmark::microbenchmark(
             propriety:::permuteTheta_false(ct, group, p = 20), times = 1L
           )

           data.frame(a$time/10^9, b$time/10^9, c$time/10^9,
                      d$time/10^9, e$time/10^9)
         })

df <- do.call(rbind, res)
df <- cbind(df, data.frame(time))
png("data-raw/bench-permute.png")
plot(y = df[,1] / 60, x = df$time, col = "red",
     ylab = "Minutes", xlab = "Total Features",
     main = "Theta Permutation Run-Time (20 Iter)")
points(y = df[,2] / 60, x = df$time, col = "orange")
# points(y = df[,3] / 60, x = df$time, col = "green")
points(y = df[,4] / 60, x = df$time, col = "blue")
# points(y = df[,5] / 60, x = df$time, col = "violet")
dev.off()
