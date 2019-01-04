d <- read.table("gpa.txt",header=TRUE)

set.seed(536)                      # Set the random number seed so results
                                   # can be exactly replicated.

d <- d[sample(1:nrow(d),15),]      # Pretend we only have 15 observations

plot(d$hsgpa,d$cgpa)
test.statistic <- cor(d$hsgpa,d$cgpa)
test.statistic
cor.test(d$hsgpa,d$cgpa)           # Test the null hypothesis of no correlation

fm <- lm( d$cgpa ~ d$hsgpa )       # Equivalently, test if slope = 0 in SLR
summary(fm)
abline(fm)


n.permutations <- 10000            # Number of permutation samples
draws <- numeric(n.permutations)
e <- d
for ( i in 1:n.permutations ) {
  e[,"hsgpa"] <- sample(e[,"hsgpa"],nrow(e))   # Permute the high school GPA
  draws[i] <- cor(e$hsgpa,e$cgpa)              # Compute the correlation
}

hist(draws)                                    # Display results
abline(v=abs(test.statistic))
abline(v=-abs(test.statistic))

                                               # Compute the two-sided p.value
p.value = ( sum( draws > abs(test.statistic) ) + sum( draws < -abs(test.statistic) ) ) / length(draws)

p.value


