M <- read.table("data-raw/gene_matrix.txt")
counts <- M
group = c(rep("A", 54), rep("B", 98-54))

time <- c(3, 10, 50, 100, 300, 500)
res <-
  lapply(time,
         function(k){

           ct <- counts[, 1:k]

           a <- microbenchmark::microbenchmark(
             propriety:::alphaTheta(ct, group, alpha = .1), times = 1L
           )

           b <- microbenchmark::microbenchmark(
             propriety:::alphaTheta_old(ct, group, alpha = .1), times = 1L
           )

           data.frame(a$time/10^9, b$time/10^9)
         })

df <- do.call(rbind, res)
df <- cbind(df, data.frame(time))
png("data-raw/bench-alphaTheta.png")
plot(y = df[,2] / 60, x = df$time, col = "red",
     ylab = "Minutes", xlab = "Total Features", type = "b",
     main = "Alpha Theta Calculation Run-Time (100 Samples)")
points(y = df[,1] / 60, x = df$time, col = "blue", type = "b")
dev.off()
