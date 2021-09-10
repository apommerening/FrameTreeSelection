# Agreement et al. for frame-tree selection. Ap - 07.05.2021, updated on 31.05.2021.

rm(list = ls())
options(digits = 4, width = 200)

# install.packages("irr", dep = T)
library(irr)
library(moments)
library(spatstat)

# Sys.setenv("LANGUAGE"="EN")

# Input file path
filePath <- "C:/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/FrameTrees/Data/"

# Output file path
outFilePath <- "C:/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/FrameTrees/Results/"

removeRowNames <- function(myVector) {
  xx <- NULL
  for (i in 1 : length(myVector))
    xx[i] <- myVector[[i]]
  return(xx)
}

# Loading data
# NDF
forests1 <- c("Ae2011", "Ardross2012", "Ardross2013", "Bin2010", "CannockChase2012", "CannockChase2013",
              "CannockChase2014-1", "CannockChase2014-2", "CannockChase2014-3",
              "Craigvinean2015-1", "Craigvinean2015-2", "Crychan2010", "Crychan2013",
              "Dalby2011", "Glentress2013", "Haldon2014", "LochArd2015", "PeckettStone2011", 
              "Dean2016-1", "Dean2016-2", "Dean2016-3", "Dean2017", "Dean2018", "Tummel2019")

# HDF, F
forests2 <- c("Ae2011", "Ardross2012", "Ardross2013", "Bin2010", "BlackIsle2013", "CannockChase2012", "CannockChase2013",
             "CannockChase2014-1", "CannockChase2014-2", "CannockChase2014-3", "Craigvinean2013", 
             "Craigvinean2015-1", "Craigvinean2015-2", "Crychan2010", "Crychan2013", "Dalby2011",
             "Glentress2013", "Haldon2014", "LochArd2015", "PeckettStone2011", "Dean2016-1", "Dean2016-2",  
             "Dean2016-3", "Dean2017", "Dean2018", "Tummel2019")


r.vector <- n.vector <- xforest <- exper <- kappa_F <- kappa_CHS <- rv.AV <- r.v <- P0 <- Pmax <- B.AV <- 
cv.dbh <- sk.dbh <-  cv.dbh <- cv.AV <- sk.AV <- cvhd.AV <- skhd.AV <- c.mean <- meanhd.AV <- H <- V <- Dg <- 
Dgk <- overlinen.i <- NULL

length(forests1) + 2 * length(forests2)

# i <- 52 # length(forests1)
for (i in 1 : (length(forests1) + 2 * length(forests2))) {
  if(i < (length(forests1) + 1)) {
      exper  <- "ThinningEx1" # NDF
      xforest <- forests1[i]
    } else if((i >= (length(forests1) + 1)) & (i < (length(forests1) + (length(forests2) + 1)))){
      exper <- "ThinningEx2" # HDF
      xforest <- forests2[i - length(forests1)]
    } else{
      exper <- "FrameEx2" # HDF
      xforest <- forests2[i - (length(forests1) + length(forests2))]
    }
    dFile <- paste(filePath, xforest, exper, ".txt", sep = "")
    tdata <- read.table(dFile, header = T) 
    # names(tdata)
    r <- ncol(tdata)
    dbh <- tdata$dbh
    ht <- tdata$h
    v <- tdata$V
    x <- tdata$X
    y <- tdata$Y
    nr <- tdata$Treenumber
    tdata <- tdata[,10 : r]

    area <- 0.1
    if((i == 1) | (i == length(forests1) + 1) | (i == (length(forests1) + (length(forests2) + 1))))
       area <- 0.133
    r <- ncol(tdata)
    r.vector[i] <- r
    n <- nrow(tdata)
    n.vector[i] <- n
    
    # Sort foresters according to s
    s <- colSums(tdata)
    tdata <- tdata[order(s, decreasing = T)]
    s <- colSums(tdata)
    overlinen.i[i] <- mean(s)
    
    # Sort trees according to marking frequencies
    rs <- removeRowNames(rowSums(tdata))
    tdata <- tdata[order(rs, decreasing = T), ]
    dbh <- dbh[order(rs, decreasing = T)]
    ht <- ht[order(rs, decreasing = T)]
    v <- v[order(rs, decreasing = T)]
    x <- x[order(rs, decreasing = T)]
    y <- y[order(rs, decreasing = T)]
    nr <- nr[order(rs, decreasing = T)]
    
    if(i >= (length(forests1) + (length(forests2) + 1))) {
      n.frame <- s / area
      distance <- sqrt(10000 / n.frame) # round(sqrt(10000 / n.frame))
      file2 <- paste(filePath, xforest, "ThinningEx2", ".txt", sep = "")
      data2 <- read.table(file2, header = T)
      range(x)
      range(y)
      (minx <- min(x))
      (miny <- min(y))
      xmax <- round(max(x))
      ymax <- round(max(y))
      xmax * ymax / 10000
      xwindow <- owin(c(0, xmax), c(0, ymax))
      degree <- NULL
      for (j in 1 : r) {
        # j <- 1
        dummy <- data.frame(Treenumber = nr, Ft = tdata[,j], x = x, y = y)
        index <- grep(names(tdata)[j], names(data2), value = FALSE)
        dummyy <- data2[, c(1, index)]
        dummy <- merge(dummy, dummyy, all.x = TRUE, sort = TRUE) # Data order is modified here!
        dummyx <- dummy[c(1 : 4)]
        dummyy <- dummy[c(1, 5, 3, 4)]
        rowV <- dummyx[c(2)]
        dummyx <- dummyx[rowV == 1,]
        myDataX <- ppp(dummyx$x, dummyx$y, window = xwindow)
        rowV <- dummyy[c(2)]
        dummyy <- dummyy[rowV == 1,]
        myDataY <- ppp(dummyy$x, dummyy$y, window = xwindow)
        myResult <- nncross(myDataX, myDataY, what = c("dist", "which"), k = 1 : (length(dummyy$Treenumber) - 1))
        competitors <- list()
        distances <- whichInd <- matrix(0, nrow = length(myResult$dist.1), ncol = length(dummyy$Treenumber) - 1)
        distances <- myResult[, 1 : (length(dummyy$Treenumber) - 1)]
        whichInd <- myResult[, length(dummyy$Treenumber) : (2 * (length(dummyy$Treenumber) - 1))]
        for (k in 1 : length(dummyx$Ft)) {
          myIndices <- whichInd[k, which(distances[k, ] < distance[[j]])]
          myIndices <- unlist(myIndices)
          test <- c()
          for (z in 1 : length(myIndices)) 
            test <- c(test, myIndices[[z]])
          competitors[[k]] <- dummyy$Treenumber[test]
          rm(test)
        }
        for (k in 1 : length(dummyx$Treenumber)) {
          if(length(competitors[[k]]) > 0) {
            fg <-  pi * (data2$dbh[which(dummyx$Treenumber[k] == data2$Treenumber)] / 200)^2
            sneighbours <- 0
            for (z in 1 : length(competitors[[k]]))
              sneighbours <- sneighbours + pi * (data2$dbh[which(competitors[[k]][z] == data2$Treenumber)] / 200)^2
            sneighbours <- sneighbours / fg
            # if((j == 1) & (k == 1))
            if(is.null(degree))
              degree <- data.frame(dbh = data2$dbh[which(dummyx$Treenumber[k] == data2$Treenumber)], deg = sneighbours) 
            else {
              dummy <- c(data2$dbh[which(dummyx$Treenumber[k] == data2$Treenumber)], sneighbours)
              degree <- rbind(degree, dummy) }
          }
        }
        
      }
      rm(data2, distance, file2, competitors, dummyx, dummyy, index, sneighbours, fg, dummy, myDataX, myDataY, myResult, distances, whichInd)
      # }
        limx <- 30
        if(max(degree$dbh) > 30)
          limx <- round(max(degree$dbh) + 1, digit = 1)
        limy <- 30
        if((max(degree$deg) > 10) & (max(degree$deg) < 15))
        limy <- 15
        else if((max(degree$deg) > 15) & (max(degree$deg) < 20))
          limy <- 20
        else if((max(degree$deg) > 20) & (max(degree$deg) < 25))
          limy <- 25
        trendl <- ksmooth(x = degree$dbh, y = degree$deg, kernel = "normal", bandwidth = 12)
        loss.L2 <- function(abdn, dbh, degg) {
          di <- abdn[1] * dbh^abdn[2]
          dev <- (degg - di)^2
          return(sum(dev, na.rm = TRUE))
        }  
        
        abdn0 <- c(0.32, 0.789)
        abdn.L2 <- optim(abdn0, loss.L2, dbh = degree$dbh, degg = degree$deg, control = list(maxit = 30000, temp = 2000, trace = F, REPORT = 500))
        pdf(file = paste(outFilePath, xforest, exper, "CompetitionDegree.pdf", sep = ""))
        par(mar = c(2, 3, 0.5, 0.5), mfrow = c(1, 1))
        plot(degree$dbh, degree$deg, xlab = "", ylab = "", xlim = c(0, limx), ylim = c(0, limy), main = "", pch = 16, axes = F) # xlim = c(0.4, 1), ylim = c(0.4, 1)
        # lines(trendl$x, trendl$y, lwd = 4, lty = 1, col = "red")
        curve(abdn.L2$par[1] * x^abdn.L2$par[2], from = min(degree$dbh), to = max(degree$dbh), lwd = 4, lty = 1, col = "red", add = T)
        axis(1, cex.axis = 1.8)
        axis(2, cex.axis = 1.8, las = 1) 
        abline(h = 1, lwd = 1, lty = 2)
        box(lwd = 2)
        dev.off()
    }
      
    # Mean conformity number
    c <- t(t(as.matrix(tdata))%*%rs) / s
    # C <- NULL
    # for(j in 1 : dim(tdata)[2])
    #   C[j] <- sum(rs[1 : s[j]], na.rm = T) / s[j]
    # c <- c / C
    c.mean[i] <- mean(c, na.rm = T)
    # if(mean(c, na.rm = T) == Inf) stop("Inf in conformity number.\n")
    rm(c)# , C)
    
    # Fleiss
    kappa_F[i] <- kappam.fleiss(tdata, exact = FALSE, detail = FALSE)$value # Fleiss
    kappa_CHS[i] <- kappam.fleiss(tdata, exact = TRUE, detail = FALSE)$value # CHS
    
    # Percentage marks per person
    numberMarks <- s / n 
    
    # Rater bar chart (active behaviour of raters)
    rlabels <- seq(1, r, 1)
    if(r > 9)
      pdf(file = paste(outFilePath, xforest, exper, "Rater.pdf", sep = ""), width = 12)
    else
      pdf(file = paste(outFilePath, xforest, exper, "Rater.pdf", sep = ""))
    par(mar = c(2, 3.2, 0.6, 0.5))
    mp <- barplot(numberMarks, beside = F, xlab = "", ylab = "", main = "", axes = F, ylim = c(0, 0.7), col = "white", names.arg ="") #las = 1, cex.axis = 1.7, 
    if(r < 16)
      axis(1, at = mp, labels = rlabels, cex.axis = 1.7) # at = mp
    else{
    for(j in 1 : length(mp))  # length(mp)
      axis(1, at = mp[j], labels = rlabels[j], cex.axis = 1.5) # at = mp[j]
    }
    axis(2, cex.axis = 1.7, las = 1) 
    box(lwd = 2)
    dev.off()
    rm(mp, rlabels, numberMarks)
    
    # Marking bar chart ("passive behaviour" of the trees) 
    table1 <- table(rs)
    table1 <- table1 / n
    # nr <- table1
    nr <- 0 # vector()
    nr <- c(nr, rep(0, r))
    for (j in 1 : length(nr)) 
      for (k in 1 : length(table1))
        if((j - 1) == as.numeric(names(table1[k])))
          nr[j] <- table1[[k]]

    if(is.na(xforest) | is.na(exper)) stop("NA in name.\n")
   
    rlabels <- seq(0, r, 1) # length(nr) - 1
    if(length(nr) > 9) # length(nr)
      pdf(file = paste(outFilePath, xforest, exper, "Marking.pdf", sep = ""), width = 12)
    else
      pdf(file = paste(outFilePath, xforest, exper, "Marking.pdf", sep = ""))
    par(mar = c(2, 3.9, 0.6, 0.5))
    mp <- barplot(nr, beside = T, xlab = "", ylab = "", main = "", col = "white", ylim = c(0, 0.9), names.arg = "", axes = F)
    if(length(nr) < 16) # length(nr)
      axis(1, at = mp, labels = rlabels, cex.axis = 1.7) # at = mp
    else{
      for(j in 1 : length(mp))  # length(mp)
        axis(1, at = mp[j], labels = rlabels[j], cex.axis = 1.5) # at = mp[j]
    }
    axis(2, cex.axis = 1.7, las = 1)
    box(lwd = 2)
    dev.off()
    rm(table1, nr, mp, rlabels)


    # Coefficient of variation of stem diameters
    cv.dbh[i] <- sd(dbh) / mean(dbh)
    sk.dbh[i] <- skewness(dbh)
    
    # Proportion of trees with 0's only
    P0[i] <- length(rs[rs == 0]) / n
    
    # Proportion of 20% largest classes
    t <- r - round(r * 0.2)
    Pmax[i] <- length(rs[rs >= t]) / n 

    # Coefficient of variation of rater behaviour
    ns <- colSums(tdata) / n
    r.v[i] <- sd(ns) / mean(ns)

    # SG ratio
    ba <-  pi * (dbh / 200)^2
    BA.total <- sum(ba)
    N.total <- length(dbh)
    BA.sel <- NA
    for (j in 1 : r)
      BA.sel[j] <- sum(ba[tdata[, j] == 1])
    N.sel <- NULL
    for (j in 1 : r)
      N.sel[j] <- length(ba[tdata[, j] == 1])
    rBA <- BA.sel / BA.total
    rN <- N.sel / N.total
    B.AV[i] <- mean(rN / rBA)
    rm(BA.total, BA.sel, N.sel, rN, rBA)
    
    # Entropy/Simpson index of rater behaviour
    nh1 <- colSums(tdata)
    nh1 <- nh1 / sum(n)
    nh2 <- 1 - nh1
    Hn <- nh1^2 + nh2^2 # -(nh1 * log(nh1)) + (nh2 * log(nh2))
    
    bh1 <- colSums(tdata * ba)
    bh1 <- bh1 / sum(ba)
    bh2 <- 1 - bh1
    Hb <- bh1^2 + bh2^2 # -(bh1 * log(bh1)) + (bh2 * log(bh2))
    H[i] <- mean(Hn / Hb)
    rm(nh1, nh2, Hn, bh1, bh2, Hb, ba)
    
    # Magin's k factor
    vh <- colSums(tdata * v) / colSums(tdata)
    V[i] <- mean(vh) / mean(v)
    dg.s <- sqrt(colSums(tdata * dbh^2)/ colSums(tdata))
    dg.r <- sqrt(colSums((1 - tdata) * dbh^2)/ colSums((1 - tdata)))
    Dg[i] <- mean(dg.s) / sqrt(sum(dbh^2)/length(dbh))
    Dgk[i] <- mean(dg.s) / mean(dg.r)
    rm(vh, dg.s, dg.r)
    
    # Diameter coefficient of variation and skewness of selected trees
    # AV
    dbh.res <- hd.res <- list()
    cv.d <- sk.d <- cv.hd <- sk.hd <- mean.hd <- NULL
    for (j in 1 : r) {
      dbh.res[[j]] <- dbh[tdata[, j] == 1] 
      cv.d[j] <- sd(dbh.res[[j]]) / mean(dbh.res[[j]]) 
      sk.d[j] <- skewness(dbh.res[[j]], na.rm = T)
      hd.res[[j]] <- ht[tdata[, j] == 1] / (dbh[tdata[, j] == 1] / 100)
      cv.hd[j] <- sd(hd.res[[j]]) / mean(hd.res[[j]]) 
      sk.hd[j] <- skewness(hd.res[[j]], na.rm = T)
      # mean.hd[j] <- mean(hd.res[[j]])
    }
    cv.AV[i] <- mean(cv.d) 
    sk.AV[i] <- mean(sk.d)
    cvhd.AV[i] <- mean(cv.hd) 
    skhd.AV[i] <- mean(sk.hd)
    meanhd.AV[i] <- mean(unlist(hd.res)) # mean(mean.hd) mean(unlist(List1))
    
    cat("Number: ", i, " Forest: ", xforest, " Experiment:", exper, " Number markers: ", r, " Number of trees: ", n, "\n") 
    # cat("Forest: ", xforest[i], " Experiment:", exper[i], " Number markers: ", r, " Number of trees: ", n, " Excess: ", excess[i], " Abs. diff. ", diff1[i], 
    #     " Trees not in AV/AP: ", diff2[i], " ChiHat1: ", ChiHat1[i], " ChiHat2: ", ChiHat2[i], " rv.AV: ", rv.AV[i], " rv.AP: ", rv.AP[i],
    #     " initBA: ", initBA[i], " AV_BA: ", AV_BA[i], " AP_BA: ", AP_BA[i], " AV trees selected: ", ii, " AP trees selected: ", length(refVotes), 
    #     " B.AV: ", B.AV[i], " B.AP: ", B.AP[i], " cv.AV: ", cv.AV[i], " cv.AP: ", cv.AP[i], " sk.AV: ", sk.AV[i], " sk.AP: ", sk.AP[i], "\n")
    rm(dbh, v, x, y, ht, dbh.res, hd.res, cv.d, sk.d, cv.hd, sk.hd, mean.hd, tdata, dFile, r, s, rs, t, ns)
}

plot(kappa_F, kappa_CHS)
abline(0, 1)
cor(kappa_F, kappa_CHS)

plot(B.AV, V)
points(B.AV, Dg, col = "red")
points(B.AV, Dgk, col = "blue")
points(B.AV, H, col = "green")
cor(B.AV, V)
cor(B.AV, Dg)
cor(B.AV, Dgk)
cor(B.AV, H)

pdf(file = paste(outFilePath, "P0.pdf", sep = ""))
my.labels <- c("Frame", "Low", "Crown")
par(lab = c(length(my.labels), 5, 7), mar = c(2.7, 3.2, 0.5, 0.5))#, mgp = c(3, 1.1, 0))
boxplot(P0[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))], P0[1 : length(forests1)], P0[(length(forests1) + 1) : (length(forests1) + length(forests2))],  
        cex = 0.6, las = 1, cex.axis = 1.5, lwd = 2, axes = FALSE, ylim = c(0, 1)) # , notch = TRUE
axis(side = 1, at = 1 : 3, labels = my.labels, lwd = 2, cex.axis = 1.8, padj = 0.5)
axis(side = 2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2) 
dev.off()

# Paired t-test
y1 <- P0[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))] # F trees
y2 <- P0[(length(forests1) + 1) : (length(forests1) + length(forests2))] # Crown thinning
length(y1)
length(y2)
t.test(y1, y2, paired = TRUE) 

diff <- y2 - y1
hist(diff)
qqnorm(diff)
qqline(diff)
shapiro.test(diff) # large p values indicate that null hypothesis of normality cannot be rejected

pdf(file = paste(outFilePath, "Pm.pdf", sep = ""))
my.labels <- c("Frame", "Low", "Crown")
par(lab = c(length(my.labels), 5, 7), mar = c(2.7, 3.2, 0.5, 0.5))
boxplot(Pmax[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))], Pmax[1 : length(forests1)], Pmax[(length(forests1) + 1) : (length(forests1) + length(forests2))],  
        cex = 0.6, las = 1, cex.axis = 1.5, lwd = 2, axes = FALSE, ylim = c(0, 0.4)) # , notch = TRUE
axis(side = 1, at = 1 : 3, labels = my.labels, lwd = 2, cex.axis = 1.8, padj = 0.5)
axis(side = 2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2) 
dev.off()

# Paired t-test
y1 <- Pmax[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))] # F trees
y2 <- Pmax[(length(forests1) + 1) : (length(forests1) + length(forests2))] # Crown thinning
length(y1)
length(y2)
t.test(y1, y2, paired = TRUE) 

diff <- y2 - y1
hist(diff)
qqnorm(diff)
qqline(diff)
shapiro.test(diff) # large p values indicate that null hypothesis of normality cannot be rejected

# PmP0 <- Pmax / P0
# my.labels <- c("Frame", "Low", "Crown")
# par(lab = c(length(my.labels), 5, 7), mar = c(2.7, 3.2, 0.5, 0.5))
# boxplot(PmP0[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))], PmP0[1 : length(forests1)], PmP0[(length(forests1) + 1) : (length(forests1) + length(forests2))],  
#         cex = 0.6, las = 1, cex.axis = 1.5, lwd = 2, axes = FALSE, ylim = c(0, 3)) # , notch = TRUE
# axis(side = 1, at = 1 : 3, labels = my.labels, lwd = 2, cex.axis = 1.8)
# axis(side = 2, las = 1, lwd = 2, cex.axis = 1.8)
# box(lwd = 2) 

pdf(file = paste(outFilePath, "Kappa.pdf", sep = ""))
my.labels <- c("Frame", "Low", "Crown")
par(lab = c(length(my.labels), 5, 7), mar = c(2.7, 3.2, 0.5, 0.5))
boxplot(kappa_F[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))], kappa_F[1 : length(forests1)], kappa_F[(length(forests1) + 1) : (length(forests1) + length(forests2))],  
        cex = 0.6, las = 1, cex.axis = 1.5, lwd = 2, axes = FALSE, ylim = c(0, 1)) # , notch = TRUE
axis(side = 1, at = 1 : 3, labels = my.labels, lwd = 2, cex.axis = 1.8, padj = 0.5)
axis(side = 2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2) 
dev.off()

# Paired t-test
y1 <- kappa_F[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))] # F trees
y2 <- kappa_F[(length(forests1) + 1) : (length(forests1) + length(forests2))] # Crown thinning
length(y1)
length(y2)
t.test(y1, y2, paired = TRUE) 

diff <- y2 - y1
hist(diff)
qqnorm(diff)
qqline(diff)
shapiro.test(diff) # large p values indicate that null hypothesis of normality cannot be rejected

# my.labels <- c("NDF", "HDF", "Z")
# par(lab = c(length(my.labels), 5, 7), mar = c(2.7, 3.2, 0.5, 0.5))
# boxplot(c.mean[1 : length(forests1)], c.mean[(length(forests1) + 1) : (length(forests1) + length(forests2))], c.mean[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))], 
#         cex = 0.6, las = 1, cex.axis = 1.5, lwd = 2, axes = FALSE, ylim = c(0, 10)) # , notch = TRUE
# axis(side = 1, at = 1 : 3, labels = my.labels, lwd = 2, cex.axis = 1.8)
# axis(side = 2, las = 1, lwd = 2, cex.axis = 1.8)
# box(lwd = 2) 

pdf(file = paste(outFilePath, "Rv.pdf", sep = ""))
my.labels <- c("Frame", "Low", "Crown")
par(lab = c(length(my.labels), 5, 7), mar = c(2.7, 3.2, 0.5, 0.5))
boxplot(r.v[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))], r.v[1 : length(forests1)], r.v[(length(forests1) + 1) : (length(forests1) + length(forests2))],  
        cex = 0.6, las = 1, cex.axis = 1.5, lwd = 2, axes = FALSE, ylim = c(0, 1)) # , notch = TRUE
axis(side = 1, at = 1 : 3, labels = my.labels, lwd = 2, cex.axis = 1.8, padj = 0.5)
axis(side = 2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2) 
dev.off()

# Paired t-test
y1 <- r.v[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))] # F trees
y2 <- r.v[(length(forests1) + 1) : (length(forests1) + length(forests2))] # Crown thinning
length(y1)
length(y2)
t.test(y1, y2, paired = TRUE) 

diff <- y2 - y1
hist(diff)
qqnorm(diff)
qqline(diff)
shapiro.test(diff) # large p values indicate that null hypothesis of normality cannot be rejected

pdf(file = paste(outFilePath, "B_AV.pdf", sep = ""))
my.labels <- c("Frame", "Low", "Crown")
par(lab = c(length(my.labels), 5, 7), mar = c(2.7, 3.2, 0.5, 0.5))
boxplot(B.AV[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))], B.AV[1 : length(forests1)], B.AV[(length(forests1) + 1) : (length(forests1) + length(forests2))],  
        cex = 0.6, las = 1, cex.axis = 1.5, lwd = 2, axes = FALSE, ylim = c(0, 2.5)) # , notch = TRUE
axis(side = 1, at = 1 : 3, labels = my.labels, lwd = 2, cex.axis = 1.8, padj = 0.5)
axis(side = 2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2) 
dev.off()

# Paired t-test
y1 <- B.AV[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))] # F trees
y2 <- B.AV[(length(forests1) + 1) : (length(forests1) + length(forests2))] # Crown thinning
length(y1)
length(y2)
t.test(y1, y2, paired = TRUE) 

diff <- y2 - y1
hist(diff)
qqnorm(diff)
qqline(diff)
shapiro.test(diff) # large p values indicate that null hypothesis of normality cannot be rejected

pdf(file = paste(outFilePath, "Cv_dbh.pdf", sep = ""))
my.labels <- c("Frame", "Low", "Crown")
par(lab = c(length(my.labels), 5, 7), mar = c(2.7, 3.2, 0.5, 0.5))
boxplot(cv.AV[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))], cv.AV[1 : length(forests1)], cv.AV[(length(forests1) + 1) : (length(forests1) + length(forests2))],  
        cex = 0.6, las = 1, cex.axis = 1.5, lwd = 2, axes = FALSE, ylim = c(0, 0.5)) # , notch = TRUE
axis(side = 1, at = 1 : 3, labels = my.labels, lwd = 2, cex.axis = 1.8, padj = 0.5)
axis(side = 2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2) 
dev.off()

# Paired t-test
y1 <- cv.AV[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))] # F trees
y2 <- cv.AV[(length(forests1) + 1) : (length(forests1) + length(forests2))] # Crown thinning
length(y1)
length(y2)
t.test(y1, y2, paired = TRUE) 

diff <- y2 - y1
hist(diff)
qqnorm(diff)
qqline(diff)
shapiro.test(diff) # large p values indicate that null hypothesis of normality cannot be rejected

# my.labels <- c("NDF", "HDF", "Z")
# par(lab = c(length(my.labels), 5, 7), mar = c(2.7, 3.8, 0.5, 0.5))
# boxplot(sk.AV[1 : length(forests1)], sk.AV[(length(forests1) + 1) : (length(forests1) + length(forests2))], sk.AV[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))], 
#         cex = 0.6, las = 1, cex.axis = 1.5, lwd = 2, axes = FALSE, ylim = c(-0.5, 1.5)) # , notch = TRUE
# axis(side = 1, at = 1 : 3, labels = my.labels, lwd = 2, cex.axis = 1.8)
# axis(side = 2, las = 1, lwd = 2, cex.axis = 1.8)
# box(lwd = 2) 

pdf(file = paste(outFilePath, "Hd.pdf", sep = ""))
my.labels <- c("Frame", "Low", "Crown")
par(lab = c(length(my.labels), 5, 7), mar = c(2.7, 3.5, 0.5, 0.5))
boxplot(meanhd.AV[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))], meanhd.AV[1 : length(forests1)], meanhd.AV[(length(forests1) + 1) : (length(forests1) + length(forests2))],  
        cex = 0.6, las = 1, cex.axis = 1.5, lwd = 2, axes = FALSE, ylim = c(0, 140)) # , notch = TRUE
axis(side = 1, at = 1 : 3, labels = my.labels, lwd = 2, cex.axis = 1.8, padj = 0.5)
axis(side = 2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2) 
dev.off()

# Paired t-test
y1 <- meanhd.AV[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))] # F trees
y2 <- meanhd.AV[(length(forests1) + 1) : (length(forests1) + length(forests2))] # Crown thinning
length(y1)
length(y2)
t.test(y1, y2, paired = TRUE) 

diff <- y2 - y1
hist(diff)
qqnorm(diff)
qqline(diff)
shapiro.test(diff) # large p values indicate that null hypothesis of normality cannot be rejected

# my.labels <- c("NDF", "HDF", "Z")
# par(lab = c(length(my.labels), 5, 7), mar = c(2.7, 3.8, 0.5, 0.5))
# boxplot(skhd.AV[1 : length(forests1)], sk.AV[(length(forests1) + 1) : (length(forests1) + length(forests2))], sk.AV[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))], 
#         cex = 0.6, las = 1, cex.axis = 1.5, lwd = 2, axes = FALSE, ylim = c(-0.5, 1.5)) # , notch = TRUE
# axis(side = 1, at = 1 : 3, labels = my.labels, lwd = 2, cex.axis = 1.8)
# axis(side = 2, las = 1, lwd = 2, cex.axis = 1.8)
# box(lwd = 2) 

kappaToClass <- function(kappa) { # Stoyan et al. (2017a, Table 4)
  kappa <- kappa + 1e-08
  xclass <- 5
  if(kappa <= 0.10) 
    xclass <- 1
  if((kappa > 0.10) & (kappa <= 0.33)) 
    xclass <- 2
  if((kappa > 0.33) & (kappa <= 0.50)) 
    xclass <- 3
  if((kappa > 0.50) & (kappa <= 0.67)) 
    xclass <- 4
  if((kappa > 0.67) & (kappa <= 0.90)) 
    xclass <- 5
  if(kappa > 0.90)
    xclass <- 6
  return(xclass)
}

kappaDistr <- matrix(0, nrow = 6, ncol = 3)
# for (i in 1 : length(kappa_F)) 
#   kappaDistr[kappaToClass(kappa_F[i]), 1] <- kappaDistr[kappaToClass(kappa_F[i]), 1] + 1
sum(kappaDistr)

kappaDistrL <-  matrix(0, nrow = 6, ncol = 1)
kappa_FL <- kappa_F[1 : length(forests1)]
for (i in 1 : length(kappa_FL)) 
  kappaDistrL[kappaToClass(kappa_FL[i])] <- kappaDistrL[kappaToClass(kappa_FL[i])] + 1
sum(kappaDistrL)

kappaDistrC <-  matrix(0, nrow = 6, ncol = 1)
kappa_FC <- kappa_F[(length(forests1) + 1) : (length(forests1) + length(forests2))]
for (i in 1 : length(kappa_FC)) 
  kappaDistrC[kappaToClass(kappa_FC[i])] <- kappaDistrC[kappaToClass(kappa_FC[i])] + 1
sum(kappaDistrC)

kappaDistrF <-  matrix(0, nrow = 6, ncol = 1)
kappa_FF <- kappa_F[(length(forests1) + length(forests2) + 1) : (length(forests1) + 2 * length(forests2))]
for (i in 1 : length(kappa_FF)) 
  kappaDistrF[kappaToClass(kappa_FF[i])] <- kappaDistrF[kappaToClass(kappa_FF[i])] + 1
sum(kappaDistrF)

kappaDistr[,1] <- kappaDistrL[,1]
kappaDistr[,2] <- kappaDistrC[,1]
kappaDistr[,3] <- kappaDistrF[,1]
names <- rep("", 6) # seq(1, 6, 1)

pdf(file = paste(outFilePath, "StackedAgreementScores.pdf", sep = ""))
par(mar = c(0.9, 2.7, 1, 0.5))
r <- rbind(kappaDistr[,1], kappaDistr[,2], kappaDistr[,3])
barplot(r, ylab = "", main = "", col = c("black", "red", "blue"), las = 1,  cex.axis = 1.7, lwd = 2, 
        cex.names = 1.7, names.arg = names, ylim = c(0, 45))
barplot(r, beside = FALSE, ylab = "", main = "", lwd = 1, col = c("black", "red", "blue"), density = 10,  add = TRUE, axes = FALSE) 
#legend("topright", c("Verbleibende Matrix", "Z", "Ausscheidende"),									# Add a legend 
#col = c("white", "black", "grey"), pt.bg = "black", pt.cex = 2, pt.lwd = 3, pch = c(20, 20, 20), bty = "n")
box(lwd = 2)
dev.off()

# Result file
# allResults <- data.frame(xforest, exper, r.vector, n.vector, kappa_F, kappa_CHS, cv.dbh, diff1, diff2, ChiHat1, ChiHat2, rv.AV, rv.AP, initBA, AV_BA, AP_BA, AVtrees, 
# APtrees, B.AV, B.AP, cv.AV, cv.AP, sk.AV, sk.AP)
# write.table(allResults, paste(outFilePath, "allResultsTable", round(h, digits = 2), ".txt", sep = ""))
# saveRDS(intersectionsAV, paste(outFilePath, "intersectionsAV", round(h, digits = 2), sep = ""))
# saveRDS(intersectionsAP, paste(outFilePath, "intersectionsAP", round(h, digits = 2), sep = ""))
# saveRDS(deservingness, paste(outFilePath, "deservingness", round(h, digits = 2), sep = ""))
