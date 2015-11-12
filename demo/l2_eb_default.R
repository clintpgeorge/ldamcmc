eb.eta <- .46 
eb.alpha <- .09 
K <- 4

heb <- c(eb.eta, eb.alpha)
hdr <- c(1/K, 1/K) # Gensim package (2010)
hda <- c(.1, .1) # Asuncion et al. (2009) 
hdg <- c(.1, 50/K) # Griffiths and Steyvers (2004)

l2.heb.hdr <- sqrt(sum((heb - hdr)^2)) 
l2.heb.hda <- sqrt(sum((heb - hda)^2)) 
l2.heb.hdg <- sqrt(sum((heb - hdg)^2))

cat("(", heb[1], ", ", heb[2], ")[", 0.0000, "]", sep="")
cat("(", hdr[1], ", ", hdr[2], ")[", l2.heb.hdr, "]", sep="")
cat("(", hda[1], ", ", hda[2], ")[", l2.heb.hda, "]", sep="")
cat("(", hdg[1], ", ", hdg[2], ")[", l2.heb.hdg, "]", sep="")

