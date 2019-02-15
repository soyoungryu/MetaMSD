
msstat.data = read.csv(file="/[DIRECTORYHERE]/LabelFree.csv", header=TRUE, sep=",")
names(msstat.data)
msstat.data.format = data.frame(Protein=msstat.data$Protein,
								Sign=msstat.data$Tvalue,
								Pvalue=msstat.data$pvalue)

write.table(msstat.data.format, file="/[DIRECTORYHERE]/Dataset1.txt", col.names=TRUE,row.names=FALSE, quote=FALSE, sep="\t")





