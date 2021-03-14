### PCA for NAM 10/10/2019 ###

#######################
#####Load packages ####
#######################
library(edgeR)
library("DESeq2")
library(factoextra)     # Package visual for PCA
library("FactoMineR")
library(biomaRt)        # For BiomaRt annotation
library(openxlsx)
library(varhandle) # used to unfactor
library(relaimpo)

### Load base directory ###
base.dir = ".../analysis/NAM_group/"

### Annotation files ###
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
g2g <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name","start_position","end_position"), mart = ensembl)


###############################
###### Read in individuals ####
###############################
dds_counts=read.table(paste0(base.dir,"dds_counts_individual.txt"))


#############################################
########Color PCA for 'individuals' #########
#############################################
####group individuals by their 'reps' ####
group = group = c("Low", rep("High",3), rep("Low",3), "High", rep("Low",4),
                  rep("High",5), rep("Low",2), rep("High",2), "Low", "High") # Create group listings


####Start PCA plot ####
data=as.data.frame(t(log(dds_counts+1,2)))
pca_data<-PCA(data, graph=TRUE,scale.unit=FALSE)        # Run PCA
png(filename = paste0(base.dir, "/BERD 2019/PCA_all_samples2.png")) # Set png destination
fviz_pca_ind(pca_data, col.ind = group  ,mean.point = FALSE, # Plot PCA
             repel = TRUE)                                                                      
dev.off()  # png device off
label = "none" 
top=as.data.frame(head(pca_data$var$contrib,15))
top = cbind(top, g2g[match(rownames(top) , g2g$ensembl_gene_id),2]) # Annotate results with genes
top=top[,c(6,1:5)]
colnames(top)=c("Variable","PC1","PC2","PC3","PC4","PC5")
write.xlsx(top,paste0(base.dir,"/BERD 2019/Top_15variable_contributors_all_samples.xlsx"))



#########################
### Linear Regression ###
#########################
### Linear regression for predicting Base Coats levels ###
base.dir="/.../analysis/NAM_group/"

############################
### Predicting BaseCoats ###
############################
### Obtain model regression model ###
data=read.xlsx(paste0(base.dir,"Correlation/Top_cor_0.5_pearsons*.xlsx"))
data=data[1:23,-c(1,26,27)] # remove Ensemble ids, r and absoluate value of r columns and limiting variables to 30 (equal to samples)
rownames(data)=data[,1]
data=data[,-1]
data=as.data.frame(t(data)) # make variables the columns
# Checking variables to make sure they appear appropriate for linear regression
png(paste0(base.dir,"BERD 2019/Regression_variables_plots.png"))
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
plot(data$Baseline_COATS, data$`RP11-536I6.2`)
plot(data$Baseline_COATS,data$PAQR9)
plot(data$Baseline_COATS,data$SFRP2)
plot(data$Baseline_COATS,data$AC008132.12)
dev.off()
# all seem to have a linear relationship with BaseCoats
#Checking multicollinearity
x=data[,-1]
library(GGally)
png(paste0(base.dir,"BERD 2019/Regression_collinearity_plots.png"),height=2400,width=2400)
ggpairs(x)
dev.off()

# create null model
model.null = lm(Baseline_COATS ~ 1, 
                data=data)

# create full model
model.full = lm(Baseline_COATS ~ .,
                data=data)
# Use stepwise regression to determine final model
step(model.null,
     scope = list(upper=model.full),
     direction="both",
     test="Chisq",
     data=data)
# obtain final model
model.final = lm(formula = Baseline_COATS ~ `RP11-536I6.2` + NRG2 + `RP11-320M16.1` + 
                   PAQR9 + AC008132.12 + SFRP2 + `CTD-2199O4.6` + `RP11-377G16.2` + 
                   DCBLD1 + PCDHGA12 + NUDT10, data = data)

summary(model.final)

#Check diagnostics for final model
png(paste0(base.dir,"BERD 2019/Regression_diagnostic_plots.png"))
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
plot(model.final)
dev.off()


## Obtain relative importance of variables
t=calc.relimp(model.final,type=c("lmg","last","first","pratt"),rela=TRUE)
png(paste0(base.dir,"BERD 2019/Regression_Relative_Importance_plots.png"),width=1200)
barplot(t$lmg,ylab="Relative Importance",xlab="Variable",main="Relative Importance of Explanatory Variables")
dev.off()
## obtain numerical values of relative importance
t2=as.data.frame(t$lmg)
write.xlsx(t2,paste0(base.dir,"BERD 2019/Regression_Relative_Importance_values.xlsx"),row.names=TRUE,col.names=TRUE)

