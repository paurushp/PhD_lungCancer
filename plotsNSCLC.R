ggplot(d3, aes(x = factor(Threshold), y = Score, fill=Measure)) + geom_bar(position = "dodge")+xlab("Threshold")+ ylab("Fraction in %")
