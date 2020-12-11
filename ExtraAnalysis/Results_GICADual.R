setwd('~/Downloads/')
library(gtools)

files <- dir()

files <- mixedsort(files)


load(files[1])
output$ARIs
ARI <- output$ARIs

for(i in 2:length(files)){
  load(files[i])
  ARI <- rbind(ARI,output$ARIs)
}
ARI
apply(ARI, 2, mean)
apply(ARI, 2, sd)


V <- c(500,2000)                # voxels
Q <- c(2,5,20)                  # components
R <- c(2,4)                     # clusters
#D <- c('square', 'nonsquare')   # dimension mixing matrix
D <- c('square')   # dimension mixing matrix
E <- c(0.05, 0.2, 0.4)          # noise level

design <- expand.grid(V=V, Q = Q,
                      R = R, D = D, E=E)

res <- cbind(design, ARI)

library(doBy)
summary_by(data = res, formula = hcl+hcls+pam+km+pcakm~.)
summary_by(data = res, formula = hcl+hcls+pam+km+pcakm~V)
summary_by(data = res, formula = hcl+hcls+pam+km+pcakm~Q)
summary_by(data = res, formula = hcl+hcls+pam+km+pcakm~R)
summary_by(data = res, formula = hcl+hcls+pam+km+pcakm~E)

# G2
summary_by(data = res, formula = hclG2+hclsG1+pamG2+kmG2+pcakmG2~.)
summary_by(data = res, formula = hclG2+hclsG1+pamG2+kmG2+pcakmG2~V)
summary_by(data = res, formula = hclG2+hclsG1+pamG2+kmG2+pcakmG2~Q)
summary_by(data = res, formula = hclG2+hclsG1+pamG2+kmG2+pcakmG2~R)
summary_by(data = res, formula = hclG2+hclsG1+pamG2+kmG2+pcakmG2~E)
