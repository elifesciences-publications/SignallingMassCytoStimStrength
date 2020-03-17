## script to write to file data underlying cydar plots

library(cydar)
load('data/cydar_results_singlet_activation.RData')

tab <- res.all$table
tab$qval.all <- qval.all

tab <- cbind(tab, intensities(cd2))

write.table(tab, file='data/DA_results.csv', quote=FALSE, sep=',')