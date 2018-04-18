M <- read.table("data-raw/gene_matrix.txt")
counts <- M
group = c(rep("A", 54), rep("B", 98-54))

time <- c(3, 10, 50, 100, 300, 500, 1000, 2000, 4000, 8000, 10000)
res <-
  lapply(time,
         function(k){

           ct <- counts[, 1:k]

           a <- microbenchmark::microbenchmark(
             propriety:::calculateTheta(ct, group), times = 1L
           )

           data.frame(a$time/10^9)
         })

df <- do.call(rbind, res)
df <- cbind(df, data.frame(time))
png("data-raw/bench-calculate.png")
plot(y = df[,1] / 60, x = df$time, col = "red",
     ylab = "Minutes", xlab = "Total Features", type = "b",
     main = "Theta Calculation Run-Time (100 Samples)")
dev.off()
