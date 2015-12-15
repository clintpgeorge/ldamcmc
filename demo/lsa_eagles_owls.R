library(ldamcmc)
library(lsa)

# Loads data --------------------------------------------------------------
# Here, the dataset chosen is eagles-owls. We can load any available data for  
# this purpose, e.g., wt, felines, cats, canis, etc. 

data("eagles-owls") # See help("eagles-owls")

K <- length(class.labels) # the number of topics 
V <- length(vocab) # the vocabulary size 

doc_tf <- calc_doc_tf(num.docs, V, ds$wid+1, ds$did+1); # tf-matrix 
colnames(doc_tf) <- vocab
rownames(doc_tf) <- 1:num.docs
doc_tf <- t(doc_tf)
doc_tfidf <- lw_logtf(doc_tf) * gw_idf(doc_tf)

# Computes the LSA subspace 
lsa.space = lsa(doc_tfidf, dims=K)

# # document semantic space
# lsa.space$dk # D x K matrix 
# 
# # term semantic space 
# lsa.space$tk # V x K matrix 

# Display important terms 
lsa.beta <- lsa.space$tk
rownames(lsa.beta) <- vocab
colnames(lsa.beta) <- c("Component-1", "Component-2")
calc_top_topic_words(t(lsa.beta), vocab, num.words=20, num.digits=2)


x <- lsa.space$dk[,1]
y <- lsa.space$dk[,2]

color.map <- rep("red", num.docs)
color.map[docs.metadata[,"category"] == "Category:Eagles"] <- "blue"

# Sets cex properties for the line plots 
cex.axis <- 1.3;
cex.lab <- 1.5;
cex <- 1; 
cex.legend <- 1.3; 

op <- par(bg = "white")
trellis.device(pdf, file="eagles-owls-lsa-space.pdf", height=5.5, width=7, onefile=T)
plot(x,y, main="LSA Semantic Space", xlab="Component-1", ylab="Component-2", 
     col=color.map, pch=16, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis);
par(op)
dev.off() 





