PRPC <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/infer_velocity_pseudotime/PRPC.csv")
from = c(59, 70, 76, 79, 87, 91, 100, 103, 116, 137, 141, 142, 162, 165)
to = c(8, 10, 10, 10, 13, 13, 15, 15,
       15, 19, 19, 19, 23, 23)
PRPC$Week <- plyr::mapvalues(PRPC$Days, from = from, to = to)

res0 <- cor.test(PRPC$Week,PRPC$velocity_pseudotime, method = "pearson")

res1 <- t.test(PRPC[(PRPC$Region == "Macula")&(PRPC$Week == 10),]$velocity_pseudotime, PRPC[(PRPC$Region == "Peripheral")&(PRPC$Week == 10),]$velocity_pseudotime)

res2 <- t.test(PRPC[(PRPC$Region == "Macula")&(PRPC$Week == 13),]$velocity_pseudotime, PRPC[(PRPC$Region == "Peripheral")&(PRPC$Week == 13),]$velocity_pseudotime)

res3 <- t.test(PRPC[(PRPC$Region == "Macula")&(PRPC$Week == 15),]$velocity_pseudotime, PRPC[(PRPC$Region == "Peripheral")&(PRPC$Week == 15),]$velocity_pseudotime)

res4 <- t.test(PRPC[(PRPC$Region == "Macula")&(PRPC$Week == 19),]$velocity_pseudotime, PRPC[(PRPC$Region == "Peripheral")&(PRPC$Week == 19),]$velocity_pseudotime)

res5 <- t.test(PRPC[(PRPC$Region == "Macula")&(PRPC$Week == 23),]$velocity_pseudotime, PRPC[(PRPC$Region == "Peripheral")&(PRPC$Week == 23),]$velocity_pseudotime)

p_values <- c(res0$p.value,res1$p.value, res2$p.value, res3$p.value, res4$p.value, res5$p.value)
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")