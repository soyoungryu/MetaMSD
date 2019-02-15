###################
#Specify Directory#
###################
directory="/Users/soyoungryu/Box\ Sync/BoxFolderOrg/Research/MetaAnalysis/Manuscript/SoftwareDeposit/Data/Simulation/SampleData/"

###########
#Parameter#
###########
##parameters from mouse data
mu=2.42; sigma=0.30
alpha=1.63; beta=-0.90; eta=0.50

#########################
#User Defined Parameters#
#########################
#varying parameters
##parameters from mouse data
mu=2.42; sigma=0.30
alpha=1.63; beta=-0.90; eta=0.50

#varying parameters
n=c(6,6) #two group comparison; n per group
nDataset=length(n) #combine 2 experiments
proportion=0.75 #percent of protein.n will be quantified for a particular dataset
wt=c(sqrt(2*n[1]), sqrt(2*n[2]))

n.sim=1000
q.value.cutoff=0.05
protein.n=2000
delta = array(c(0), dim=protein.n)
increment=ceiling(protein.n*0.05)
delta[1:increment]=0.58 #1.5 fold
delta[(increment+1):(2*increment)]=-0.58
delta[(2*increment+1):(3*increment)]=1 #2 fold
delta[(3*increment+1):(4*increment)]=-1
delta[(4*increment+1):(5*increment)]=2
delta[(5*increment+1):(6*increment)]=-2

gammaA=exp(rnorm(protein.n, mu, sigma))
gammaB=gammaA+delta
lambdaA = lambdaB = array(c(NA), dim=protein.n)
for (i in 1:protein.n){
lambdaA[i]=exp(rnorm(1, alpha+beta*log(gammaA[i]),eta))
lambdaB[i]=exp(rnorm(1, alpha+beta*log(gammaB[i]),eta))
}



############
#Simulation#
############

for (k in 1:n.sim){
	cat("Simulation",k, "\n")
DEprotein=NULL
p.list=stat.list=matrix(c(NA), nrow=protein.n, ncol=nDataset)
for (j in 1:nDataset){
Intensity.GroupA=matrix(c(NA), nrow=protein.n, ncol=n[j])
Intensity.GroupB=matrix(c(NA), nrow=protein.n, ncol=n[j])
p=stat=array(c(NA), dim=protein.n)
selected.protein.index=sample(1:protein.n)[1:(protein.n*proportion)]
for (i in 1:protein.n){
	if(i %in% selected.protein.index){
	Intensity.GroupA[i,] = rnorm(n[j], gammaA[i], lambdaA[i])
	Intensity.GroupB[i,] = rnorm(n[j], gammaB[i], lambdaB[i])
	test=t.test(Intensity.GroupA[i,], Intensity.GroupB[i,])
	p[i]=test$p.value; stat[i]=test$statistic
	}
}
p.list[,j] = p
stat.list[,j] = stat
}

DEprotein = list(p=p.list, stat=stat.list)

for (j in 1:ncol(DEprotein$p)){
	tmp=data.frame(Protein=paste("Protein", c(1:length(DEprotein$p[,j])), "Name", sep=""),
				   Sign=ifelse(DEprotein$stat[,j]>0,1,-1),
				   Pvalue=DEprotein$p[,j])
	write.table(tmp, file=paste(directory,"Simulation",k,"_Dataset", j, ".txt", sep=""), col.names=TRUE, row.names=FALSE, quote=FALSE)
}
}



