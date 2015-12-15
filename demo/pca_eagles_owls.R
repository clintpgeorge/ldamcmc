library(ldamcmc)
library(lsa)


# Loads data

data("eagles-owls") # See help("eagles-owls")

K <- length(class.labels) # the number of topics 
V <- length(vocab) # the vocabulary size 

doc_tf <- calc_doc_tf(num.docs, V, ds$wid+1, ds$did+1); # tf-matrix 
colnames(doc_tf) <- vocab
rownames(doc_tf) <- 1:num.docs
doc_tf <- t(doc_tf)
doc_tfidf <- lw_logtf(doc_tf) * gw_idf(doc_tf)



# apply PCA
# This is based on SVD. So it's similar to PCA, if we center the data 
pca.res <- prcomp(t(doc_tfidf), retx=T, center=F)
summary(pca.res)
plot(pca.res)


# Sets cex properties for the line plots 
cex.axis <- 1.3;
cex.lab <- 1.5;
cex <- 1; 
cex.legend <- 1.3; 

# Set colors for different classes 
color.map <- rep("red", num.docs)
color.map[docs.metadata[,"category"] == "Category:Eagles"] <- "blue"


op <- par(bg = "white")
trellis.device(pdf, file="eagles-owls-pca-space.pdf", height=5.5, width=7, onefile=T)
plot(pca.res$x[,1], pca.res$x[,2], main="PCA", xlab="Component-1", ylab="Component-2", 
     col=color.map, pch=16, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis);
par(op)
dev.off() 



# 
# ####################### PCA from the scratch ###################################
# 
# 
# # Obtain data in a matrix
# Xoriginal <- t(doc_tfidf)
# 
# # Center the data so that the mean of each row is 0
# rm <- rowMeans(Xoriginal)
# X <- Xoriginal - matrix(rep(rm, dim(Xoriginal)[2]), nrow=dim(Xoriginal)[1])
# 
# # Calculate P
# A <- X %*% t(X)
# E <- eigen(A, TRUE)
# P <- t(E$vectors)
# 
# # Find the new data and standard deviations of the principal components
# newdata = P %*% X
# sdev = sqrt(diag((1/(dim(X)[2]-1)* P %*% A %*% t(P))))
# 

