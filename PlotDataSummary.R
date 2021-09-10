# Approval voting and apportionment, calculation of basic forest data. Ap - 16.05.2020, updated on 26.05.2021.
rm(list = ls())
options(digits = 4, width = 200)

library(irr)
library(moments)

# Input file path
filePath <- "C:/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/FrameTrees/Data/"

# Output file path
outFilePath <- "C:/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/FrameTrees/Results/"

calcD100 <- function(myData, area) {
  numberD100 <- 100 * area
  inumberD100 <- trunc(numberD100)
  rest <- abs(numberD100 - inumberD100)
  dummy <- myData[order(myData, decreasing = TRUE)] [1:inumberD100]	
  d100 <- sum(dummy^2) 
  dummy <- myData[order(myData, decreasing = TRUE)] [1:(inumberD100 + 1)]	
  d100 <- d100 + (dummy * rest)^2
  d100 <- sqrt(sum(dummy^2)/numberD100)
  return(d100)
}

# Loading data
forests <- c("Ae2011", "Ardross2012", "Ardross2013", "Bin2010", "BlackIsle2013", "CannockChase2012", "CannockChase2013",
             "CannockChase2014-1", "CannockChase2014-2", "CannockChase2014-3", "Craigvinean2013", 
             "Craigvinean2015-1", "Craigvinean2015-2", "Crychan2010", "Crychan2013", "Dalby2011",
             "Glentress2013", "Haldon2014", "LochArd2015", "PeckettStone2011", "Dean2016-1", "Dean2016-2",  
             "Dean2016-3", "Dean2017", "Dean2018", "Tummel2019")

# Logarithmic height curve (from Jens' spreadsheets)
a <- c(6.2358, 3.7292, 3.7292, 6.6303, 3.2189, 3.1115, 3.1115, 4.4597, 4.4597, 4.4597, 5.0003, 4.5609, 4.5609, 
       1.7692, 3.7673, 4.6634, 6.9799, 4.0365, 7.6962, 8.6232, 3.854, 3.854, 3.854, 3.854, 5.369, 5.2645)
b <- c(-0.4107, 1.5411, 1.5411, -1.203, 1.338, 4.8628, 4.8628, 2.6998, 2.6998, 2.6998, -0.7003, 1.9505, 1.9505, 
       10.455, 5.3892, 3.2286, -0.9151, 5.044, -6.9394, -5.8544, 1.2584, 1.2584, 1.2584, 1.2584, -0.5729, -2.57)

mainSpecies <- NULL
N <- NULL
G <- NULL
dg <- NULL
h100 <- NULL
vd <- NULL
kd <- NULL
area <- 0.1

# i <- 8
for (i in 1 : length(forests)) {
    dFile <- paste(filePath, forests[i], "FrameEx2.txt", sep = "")
    tdata <- read.table(dFile, header = T) 

    if(i == 1)
      area <- 0.1332
    else 
      area <- 0.1
    
    # Calculate basal area,sph and dg
    ba <-  pi * (tdata$dbh / 200)^2
    G[i] <- sum(ba) / area
    dg[i] <- sqrt(sum(tdata$dbh^2)/length(tdata$dbh))
    N[i] <- length(tdata$dbh) / area
    rm(ba)

    # Top diameter and height
    d100 <- calcD100(tdata$dbh, area)
    h100[i] <- a[i] * log(d100) + b[i]
    rm(d100)
    
    # Main species
    stable <- table(tdata$Species)
    stable <- sort(stable, decreasing = TRUE)
    mainSpecies[i] <- names(stable[1])
    rm(stable)
    names(tdata)

    # Diameter coefficient of variation and skewness
    vd[i] <-  sd(tdata$dbh) / mean(tdata$dbh)
    kd[i] <- skewness(tdata$dbh)

    r <- ncol(tdata)
    tdata <- tdata[,10 : r]
    r <- ncol(tdata)
    n <- nrow(tdata)

    rm(dFile, tdata)
    cat("Forest: ", forests[i], " G: ", G[i], " dg: ", dg[i], " N: ", N[i], " h100: ", h100[i], " mainSpecies: ", mainSpecies[i], 
        " vd: ", vd[i], " kd: ", kd[i], " r: ", r, " n: ", n, "\n")
}

# Result file
basicResults <- data.frame(forests, G, dg, N, h100, mainSpecies, vd, kd)
write.table(basicResults, paste(outFilePath, "basicResultsTable.txt", sep = ""))
