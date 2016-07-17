# Set number of features in counts matrix
i <- 8000
X <- t(data.frame("a" = abs(rnorm(i)), "b" = abs(rnorm(i)), "c" = abs(rnorm(i)),
                  "d" = abs(rnorm(i)), "e" = abs(rnorm(i)), "f" = abs(rnorm(i))))

# Determine peak RAM needed for phi calculation
out <- lapply(list(X),
              function(df){

                usage <- data.frame()

                for(n in c(10, 100, 500, 1000, 2000, 4000, 6000, 8000)){

                  counts.subset <- X[, 1:n]
                  z <- function() {phit(counts.subset)}
                  max <- miSciTools::peakRAM(z)
                  usage <- rbind(usage, data.frame("n" = n, "max" = max))
                }

                return(usage)
              }
)

# Model additional peak RAM predictions
usage <- out[[1]]
plot(max ~ n, usage)
plot(sqrt(max) ~ n, usage)
fit <- lm(sqrt(max) ~ n, usage)
df <- data.frame("n" = c(1, 100, 500, 1000, 2000, 4000, 8000, 16000,
                         24000, 32000, 64000, 100000))
rownames(df) <- df$n
p <- predict(fit, df)^2 # predict
final <- data.frame(round(p))

# Format table
knitr::kable(final, format = "markdown")
