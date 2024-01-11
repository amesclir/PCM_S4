## read morphological data from file
anole.morphology<-read.csv("anole.data.csv",row.names=1)
## read ecological states from file
anole.ecomorph<-read.csv("ecomorph.csv",row.names=1,stringsAsFactors=TRUE)

## read in the phylogeny with ecomorph mapped
## using phytools::read.simmap as this is not a regular tree but a tree with 6 discrete traits mapped onto the phylogeny
ecomorph.tree<-read.simmap(file="anolis.mapped.nex",format="nexus",version=1.5)
ecomorph.tree

## check to see if data and phylogeny match using
## geiger::name.check
chk<-name.check(ecomorph.tree,anole.morphology)
summary(chk)

## trim our input data to match the tree using
## negative indexing
ecomorph.data<-anole.morphology[-which(rownames(anole.morphology)%in%chk$data_not_tree),]
name.check(ecomorph.tree,ecomorph.data)

library(OUwie)

## run phylogenetic PCA and print the results
pca<-phyl.pca(ecomorph.tree,ecomorph.data)
print(pca)

## create our OUwie data frame
ouwie.data<-data.frame(Genus_species=rownames(scores(pca)),Reg=anole.ecomorph[rownames(scores(pca)),],X=as.numeric(scores(pca)[,1]))
head(ouwie.data,n=10)

cols<-setNames(rainbow(n=6),levels(anole.ecomorph[,1]))
plot(ecomorph.tree,cols,lwd=2,ftype="i",fsize=0.4,ylim=c(-4,82),outline=TRUE)
add.simmap.legend(colors=cols,prompt=FALSE,x=0,y=-2,vertical=FALSE,fsize=0.9)

## fit standard, one-rate Brownian model\
fitBM<-OUwie(ecomorph.tree,ouwie.data,model="BM1",simmap.tree=TRUE)
fitBM

## fit multi-rate Brownian model
fitBMS<-OUwie(ecomorph.tree,ouwie.data,model="BMS",simmap.tree=TRUE,root.station=FALSE)
fitBMS

# fit multi-regime OU model
fitOUM<-OUwie(ecomorph.tree,ouwie.data,model="OUM",simmap.tree=TRUE,root.station=FALSE)
fitOUM

## extracting AIC scores
aic<-setNames(c(fitBM$AIC,fitBMS$AIC,fitOUM$AIC),c("BM1","BMS","OUM"))
aic
## compute Akaike weights
aic.w(aic)

## get the tip state for ecomorph for each species
tips<-getStates(ecomorph.tree,"tips")
## set these tip states to have the colors using the
## color scheme 
tip.cols<-cols[tips]
## plot tree with adjacent barplot using
## phytools::plotTree.barplot
plotTree.barplot(ecomorph.tree,scores(pca)[,1],args.plotTree=list(fsize=0.4),args.barplot=list(col=tip.cols,xlab=expression(paste("PC1 (",""%down%"","limbs, ",""%down%"","lamellae)",sep="")),cex.lab=0.8))
## add an informative legend
legend("topright",levels(anole.ecomorph[,1]),pch=22,pt.bg=cols,pt.cex=1.5,cex=0.9)
