library(ape)
library(doBy)
library(phytools)
library(quantreg)
library(diversitree)
library(apTreeshape)
library(RColorBrewer)

checkAndOpenTree <- function(filename){
  if (endsWith(filename, "nwk") || endsWith(filename, "newick")){
    write("Newick file detected", stderr())
    tree <- read.tree(filename)
    return(tree)
  } else {
    write("Tree must be in Newick format", stderr())
    quit(status = 1)
  }
}

prepTree <- function(tree, tip.date, random_resolve=TRUE){
  tree <- multi2di(tree, random = random_resolve)
  tree <- scaleTree(tree)
  tree <- rtt(t = tree, tip.dates = tip.date, objective = "correlation")
  return(tree)
}

boxlabel<-function(x, y, text, cex=1, bg="transparent", offset=0){
  w <- strwidth(text) * cex * 1.1
  h <- strheight(text) * cex * 1.4
  os <- offset * strwidth("W") * cex
  rect(x + os, y - 0.5 * h, x + w + os, y + 0.5 * h, col = bg, border = 0)
  text(x , y, text, pos =4, offset = offset, font = 3)
}

timeBinsToColors <- function(tips){
  x <- levels(tips)
  cols <- setNames(RColorBrewer::brewer.pal(length(x),"Set3"), x)
  return(cols)
}

getBinCount <- function(timeNum){
  n_bins <- round(sqrt(length(timeNum)))
  n_bins <- min(n_bins, 9)
  return(n_bins)
}

getLeafHeights <- function(tree){
  return(diag(vcv.phylo(tree)))
}

equalFreqBins <- function(timeNum, n_bins){
  n_leaves <- length(timeNum)
  n_repl <- floor(n_leaves/n_bins)
  n_plus <- sample(1:n_bins, n_leaves - n_repl*n_bins)
  n_rep <- rep(n_repl, n_bins)
  n_rep[n_plus] <- n_repl + 1
  timeNum[order(timeNum)] <- rep(seq.int(n_bins), n_rep)
  return(timeNum)
}

binTimeNum <- function(timeNum){
  timeNum <- equalFreqBins(timeNum, getBinCount(timeNum))
  states <- as.numeric(as.factor(timeNum))
  char <- as.factor(states)
  names(char) <- names(timeNum)
  names(states) <- names(timeNum)
  write(paste("Optimal number (max=9) of discrete time intervals is ", 
              length(levels(char)), sep=""), 
        stdout())
  level_names <- levels(char)
  binned_leaf_counts <- table(char)
  for (i in 1:length(levels(char))){
    write(paste("    Bin", level_names[i], 
                "contains", binned_leaf_counts[i], "leaves"), 
          stdout())
  }
  return(list("char"=char, "states"=states))
}

scaleTree <- function(tree){
  for (i in 1:length(tree$edge.length)){
    tree$edge.length[i] <-  max(0,tree$edge.length[i]) + 0.00001
  }
  return(tree)
}

ladderScores <- function(tree){
  tr <- ladderize(tree)
  b <- balance(tr)
  LT <- NULL
  for (i in 1:length(b[,1]))
  {
    left <- as.numeric(b[i,1])
    right <- as.numeric(b[i,2])
    ratio <- right/left
    LT <- c(LT,ratio)
  }
  LT <- sort(LT)
  return(LT)
}

subtips <- function(tree, node){
  if (node <- length(tree$tip.label)){
    ret <- 1
  }
  else{
    children <- which(tree$edge[,1] == node)
    ret <- (subtips(tree, 
                    tree$edge[children[1],2]) + subtips(tree,
                                                        tree$edge[children[2],2]))
  }
  return(ret)
}

# Modified from apTreeshape::aldous.test
aldous <- function (tree, xmin = 1, ...){
  if (identical(tree, NULL)) {
    stop("invalid tree", "\n")
  }
  if (xmin < 0.1) {
    xmin = 0.1
  }
  clade = smaller.clade.spectrum(tree)
  if (xmin > clade[1, 1]) {
    stop("'xmin' greater than the size of the tree")
  }
  xlab = "Size of parent clade (log scale)"
  ylab = "Size of smaller daughter clade (log scale)"
  plot(clade, pch = 20, xlab = xlab, log = "xy", ylab = ylab, 
       xlim = c(xmin, max(clade) + 10), ...)
  if (xmin >= 10) {
    x = xmin + 5
    text(x = x, y = 3/2 + 0.1, label = "PDA model", col = 4)
    text(x = x, y = x/4 - 1, label = "Yule model", col = 4)
  }
  else {
    x = 20
    text(x = x, y = 3/2 + 0.1, label = "PDA model", col = 4)
    text(x = x, y = x/4 + 2, label = "Yule model", col = 4)
  }
  mod = rq(clade[, 2] ~ clade[, 1])
  x = sort(clade[, 1])
  coef = coefficients(mod)
  abline(coef, untf = TRUE, col = 3)
  abline(0, 1/4, untf = TRUE)
  abline(3/2, 0, untf = TRUE)
  abline(v = 30, lty = 3)
}

formatData <- function(filename, tree){
  data <- read.table(filename, header = TRUE)
  data <- data[match(tree$tip.label, data$taxaId),]
  return(list("data"=data,"tree"=tree))
}

modelCharacterEvolution <- function(tree, timetable, treefile="TC.out", 
                                    boot=200, parsimony = F){
  fdata <- formatData(timetable, tree)
  data <- fdata$data
  n_bins <- getBinCount(data$timeNum)
  bins <- binTimeNum(data$timeNum)
  data$timeNum <- bins$states
  char <- bins$char
  states <- bins$states
  names(char) <- data$taxaId
  names(states) <- data$taxaId
  
  original_tree <- tree
  if (parsimony){
    tree <- compute.brlen(tree, 1)
  }
  
  phy <- chronos(tree, lambda=0.1)
  
  mkn <- make.mkn(phy, states, max(states))
  if (max(states) == 2)
    mkn <- constrain(mkn, q21 ~ 0)
  if (max(states) == 3)
    mkn <- constrain(mkn, q21 ~ 0, q31 ~ 0, q32 ~ 0)
  if (max(states) == 4)
    mkn <- constrain(mkn, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0)
  if (max(states) == 5)
    mkn <- constrain(mkn, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0, q51 ~ 0, q52 ~ 0, q53 ~ 0, q54 ~ 0)
  if (max(states) == 6)
    mkn <- constrain(mkn, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0, q51 ~ 0, q52 ~ 0, q53 ~ 0, q54 ~ 0, q61 ~ 0, q62 ~ 0, q63 ~ 0, q64 ~ 0, q65 ~ 0)
  if (max(states) == 7)
    mkn <- constrain(mkn, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0, q51 ~ 0, q52 ~ 0, q53 ~ 0, q54 ~ 0, q61 ~ 0, q62 ~ 0, q63 ~ 0, q64 ~ 0, q65 ~ 0, q71 ~ 0, q72 ~ 0, q73 ~ 0, q74 ~ 0, q75 ~ 0, q76 ~ 0)
  if (max(states) == 8)
    mkn <- constrain(mkn, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0, q51 ~ 0, q52 ~ 0, q53 ~ 0, q54 ~ 0, q61 ~ 0, q62 ~ 0, q63 ~ 0, q64 ~ 0, q65 ~ 0, q71 ~ 0, q72 ~ 0, q73 ~ 0, q74 ~ 0, q75 ~ 0, q76 ~ 0, q81 ~ 0, q82 ~ 0, q83 ~ 0, q84 ~ 0, q85 ~ 0, q86 ~ 0, q87 ~ 0)
  if (max(states) == 9)
    mkn <- constrain(mkn, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0, q51 ~ 0, q52 ~ 0, q53 ~ 0, q54 ~ 0, q61 ~ 0, q62 ~ 0, q63 ~ 0, q64 ~ 0, q65 ~ 0, q71 ~ 0, q72 ~ 0, q73 ~ 0, q74 ~ 0, q75 ~ 0, q76 ~ 0, q81 ~ 0, q82 ~ 0, q83 ~ 0, q84 ~ 0, q85 ~ 0, q86 ~ 0, q87 ~ 0, q91 ~ 0, q92 ~ 0, q93 ~ 0, q94 ~ 0, q95 ~ 0, q96 ~ 0, q97 ~ 0, q98 ~ 0)
  if (max(states) == 10)
    mkn <- constrain(mkn, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0, q51 ~ 0, q52 ~ 0, q53 ~ 0, q54 ~ 0, q61 ~ 0, q62 ~ 0, q63 ~ 0, q64 ~ 0, q65 ~ 0, q71 ~ 0, q72 ~ 0, q73 ~ 0, q74 ~ 0, q75 ~ 0, q76 ~ 0, q81 ~ 0, q82 ~ 0, q83 ~ 0, q84 ~ 0, q85 ~ 0, q86 ~ 0, q87 ~ 0, q91 ~ 0, q92 ~ 0, q93 ~ 0, q94 ~ 0, q95 ~ 0, q96 ~ 0, q97 ~ 0, q98 ~ 0, q1001 ~ 0, q1002 ~ 0, q1003 ~ 0, q1004 ~ 0, q1005 ~ 0, q1006 ~ 0, q1007 ~ 0, q1008 ~ 0, q1009 ~ 0)
  
  p <- runif(length(argnames(mkn)))
  fit.mkn <- find.mle(mkn, p)

  marginal_likelihood <- asr.marginal(mkn, coef(fit.mkn))
  t_marginal_likelihood <- t(marginal_likelihood)
  ancstates <- NULL
  for(i in 1:length(t_marginal_likelihood[,1])){
    a <- 1
    v <- t_marginal_likelihood[i,1]
    for(j in 2:length(t_marginal_likelihood[1,])){
      if (t_marginal_likelihood[i,j] > v)
      {
        a <- j
        v <- t_marginal_likelihood[i,j]
      }
    }
    ancstates <- c(ancstates,a)
  }
  
  internalnodeindices <- length(tree$tip.label) + (1:length(ancstates))
  internalnodestates <- as.data.frame(cbind(internalnodeindices,ancstates))
  names(internalnodestates) <- c("nodeId","state")
  leavesstates1 <- as.data.frame(cbind(1:length(tree$tip.label),tree$tip.label))
  temp <- as.data.frame(cbind(names(char),char))
  names(temp) <- c("taxaId","state")
  names(leavesstates1) <- c("nodeId","taxaId")
  leavesstates <- merge(temp, leavesstates1, by.x = "taxaId", by.y = "taxaId")
  leavesstates <- leavesstates[,c("nodeId","state")]
  nodeStates <- rbind(internalnodestates,leavesstates)
  temp1 <- as.data.frame(tree$edge)
  names(temp1) <- c("nodeId","child")
  changeStates1 <- merge(nodeStates, temp1, by.x = "nodeId", by.y = "nodeId")
  changeStates <- merge(changeStates1, nodeStates, by.x = "child", by.y = "nodeId")
  names(changeStates) <- c("child","nodeId","nodeState","childState")
  numChanges <- length(which(changeStates$nodeState!=changeStates$childState))
  actualNumChanges <- numChanges
  
  LT <- ladderScores(tree)
  randChanges <- NULL
  
  for (b in 1:boot){
    statesShuffled <- states
    names(statesShuffled) <- sample(names(statesShuffled))
    charShuffled <- as.factor(states)
    names(charShuffled) <- names(statesShuffled)
    ancstatesShuffled <- sample(ancstates)
    
    mknShuffled <- make.mkn(phy, statesShuffled, max(statesShuffled))
    if (max(statesShuffled) == 2)
      mknShuffled <- constrain(mknShuffled, q21 ~ 0)
    if (max(statesShuffled) == 3)
      mknShuffled <- constrain(mknShuffled, q21 ~ 0, q31 ~ 0, q32 ~ 0)
    if (max(statesShuffled) == 4)
      mknShuffled <- constrain(mknShuffled, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0)
    if (max(statesShuffled) == 5)
      mknShuffled <- constrain(mknShuffled, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0, q51 ~ 0, q52 ~ 0, q53 ~ 0, q54 ~ 0)
    if (max(statesShuffled) == 6)
      mknShuffled <- constrain(mknShuffled, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0, q51 ~ 0, q52 ~ 0, q53 ~ 0, q54 ~ 0, q61 ~ 0, q62 ~ 0, q63 ~ 0, q64 ~ 0, q65 ~ 0)
    if (max(statesShuffled) == 7)
      mknShuffled <- constrain(mknShuffled, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0, q51 ~ 0, q52 ~ 0, q53 ~ 0, q54 ~ 0, q61 ~ 0, q62 ~ 0, q63 ~ 0, q64 ~ 0, q65 ~ 0, q71 ~ 0, q72 ~ 0, q73 ~ 0, q74 ~ 0, q75 ~ 0, q76 ~ 0)
    if (max(statesShuffled) == 8)
      mknShuffled <- constrain(mknShuffled, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0, q51 ~ 0, q52 ~ 0, q53 ~ 0, q54 ~ 0, q61 ~ 0, q62 ~ 0, q63 ~ 0, q64 ~ 0, q65 ~ 0, q71 ~ 0, q72 ~ 0, q73 ~ 0, q74 ~ 0, q75 ~ 0, q76 ~ 0, q81 ~ 0, q82 ~ 0, q83 ~ 0, q84 ~ 0, q85 ~ 0, q86 ~ 0, q87 ~ 0)
    if (max(statesShuffled) == 9)
      mknShuffled <- constrain(mknShuffled, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0, q51 ~ 0, q52 ~ 0, q53 ~ 0, q54 ~ 0, q61 ~ 0, q62 ~ 0, q63 ~ 0, q64 ~ 0, q65 ~ 0, q71 ~ 0, q72 ~ 0, q73 ~ 0, q74 ~ 0, q75 ~ 0, q76 ~ 0, q81 ~ 0, q82 ~ 0, q83 ~ 0, q84 ~ 0, q85 ~ 0, q86 ~ 0, q87 ~ 0, q91 ~ 0, q92 ~ 0, q93 ~ 0, q94 ~ 0, q95 ~ 0, q96 ~ 0, q97 ~ 0, q98 ~ 0)
    if (max(statesShuffled) == 10)
      mknShuffled <- constrain(mknShuffled, q21 ~ 0, q31 ~ 0, q32 ~ 0, q41 ~ 0, q42 ~ 0, q43 ~ 0, q51 ~ 0, q52 ~ 0, q53 ~ 0, q54 ~ 0, q61 ~ 0, q62 ~ 0, q63 ~ 0, q64 ~ 0, q65 ~ 0, q71 ~ 0, q72 ~ 0, q73 ~ 0, q74 ~ 0, q75 ~ 0, q76 ~ 0, q81 ~ 0, q82 ~ 0, q83 ~ 0, q84 ~ 0, q85 ~ 0, q86 ~ 0, q87 ~ 0, q91 ~ 0, q92 ~ 0, q93 ~ 0, q94 ~ 0, q95 ~ 0, q96 ~ 0, q97 ~ 0, q98 ~ 0, q1001 ~ 0, q1002 ~ 0, q1003 ~ 0, q1004 ~ 0, q1005 ~ 0, q1006 ~ 0, q1007 ~ 0, q1008 ~ 0, q1009 ~ 0)

    aceShuffled <- asr.marginal(mknShuffled, coef(fit.mkn))
    t_marginal_likelihood <- t(aceShuffled)
    ancstatesShuffled <- NULL
    for(i in 1:length(t_marginal_likelihood[,1])){
      a <- 1
      v <- t_marginal_likelihood[i,1]
      for(j in 2:length(t_marginal_likelihood[1,]))
      {
        if (t_marginal_likelihood[i,j]>v)
        {
          a <- j
          v <- t_marginal_likelihood[i,j]
        }
      }
      ancstatesShuffled <- c(ancstatesShuffled,a)
    }
    
    internalnodeindices <- length(tree$tip.label) + (1:length(ancstatesShuffled))
    internalnodestates <- as.data.frame(cbind(internalnodeindices,
                                              ancstatesShuffled))
    names(internalnodestates) <- c("nodeId","state")
    leavesstates1 <- as.data.frame(cbind(1:length(tree$tip.label),
                                         tree$tip.label))
    temp <- as.data.frame(cbind(names(charShuffled),charShuffled))
    names(temp) <- c("taxaId","state")
    names(leavesstates1) <- c("nodeId","taxaId")
    leavesstates <- merge(temp, leavesstates1, 
                          by.x = "taxaId", by.y = "taxaId")
    leavesstates <- leavesstates[,c("nodeId","state")]
    nodeStates <- rbind(internalnodestates,leavesstates)
    temp1 <- as.data.frame(tree$edge)
    names(temp1) <- c("nodeId","child")
    changeStates1 <- merge(nodeStates, temp1, 
                           by.x = "nodeId", by.y = "nodeId")
    changeStates <- merge(changeStates1, nodeStates, 
                          by.x = "child", by.y = "nodeId")
    names(changeStates) <- c("child","nodeId","nodeState","childState")
    numChanges <- length(which(
      changeStates$nodeState != changeStates$childState))
    randChanges <- c(randChanges,numChanges)
    if (b%%(boot/10)==0){
      write(paste("shuffling ", b, "/", boot,
                  " (", (100*b/boot),"%)"), stderr())
    }
  }
  
  randTrees <- rmtree(N = boot, n = length(tree$tip.label), 
                      rooted = T, br = runif)
  LTrandTrees <- list()
  LTrand <- NULL
  for (i in 1:boot){
    randTrees[[i]] <- root(randTrees[[i]],1)
    randTrees[[i]] <- multi2di(randTrees[[i]], random = TRUE)

    LTrandTrees[[i]] <- ladderScores(randTrees[[i]])
    LTrand <- c(LTrand,length(which(LTrandTrees[[i]]<1))/length(LTrandTrees[[i]]))
    if (i%%(boot/10)==0){
      write(paste("shuffling ", i, "/", boot,
                  " (", (100*i/boot), "%)"), stderr())
    }
  }
  LTorig <- length(which(LT<1))/length(LT)
  
  
  ##########START PLOTTING AND OTHER STATS##########
  
  write("Plotting and Statistics", stderr())
  write("    -Plotting Tree",stderr())
  pdf("Ancestral_character_tree.pdf", width=15, height=15)
  plot(original_tree, show.tip.label=FALSE, cex=0.75)
  tiplabels(tree$tip.label, col=as.numeric(char[tree$tip.label]),
            frame="none", bg="white", cex=0.75)
  legend("bottomleft", lwd=3, paste("time point",1:max(states)),
         col=1:max(states))
  
  nodelabels(pie=t(marginal_likelihood),
             piecol=1:length(t(marginal_likelihood)[1,]), cex=0.5)
  title(main = "Ancestral Character Tree")
  dev.off()
 
  write("    -Calculating and Plotting TC", stderr())
  
  pdf("TC_plots.pdf", width=15, height=15)
  par(mfcol=c(2,4))
  
  hist(randChanges,xlim=c(min(randChanges,actualNumChanges)-1,
                          max(randChanges,actualNumChanges))+1,
       col="red",main="",xlab="# changes")
  abline(v=actualNumChanges,lwd=3,col="blue")
  legend("top",c("randomization","actual tree"),bg="white",
         col=c("red","blue"),pch=19,cex=1)
  
  TC <- (log(mean(randChanges))-log(actualNumChanges))/(log(max(randChanges))-log(length(levels(char))-1))
  zetaTC <- (actualNumChanges-mean(randChanges))/sqrt(var(randChanges)/length(randChanges))
  randChanges <- sort(randChanges)
  percTC <- round(100 * max(c(0, 
                              which(randChanges<actualNumChanges)
                              ))/length(randChanges),
                  2)
  
  write("    -Calculating and Plotting Staircase-ness", stderr())
  
  title(main = paste("TC score = ", round(TC,4),
              "( vs",boot,"random trees )","\n",
              " percentile = ", percTC))
  
  zetaLT <- (LTorig-mean(LTrand)) / sqrt(var(LTrand) / length(LTrand))
  LTrand <- sort(LTrand)
  percLT <- round(100 * min(c(length(LTrand), 
                            which(LTrand>LTorig))) / length(LTrand),2)
  
  rand1 <- stree(length(tree$tip.label), type = "left")
  theorMaxStairc <- length(which(ladderScores(rand1)<1))/length(ladderScores(rand1))
  hist(LTrand, 
       xlim=c(min(theorMaxStairc,LTorig,LTrand), 
              max(theorMaxStairc,LTorig,LTrand)),
       col="red", 
       xlab="scores",
       main=paste("staircase-ness = ",
                  round(LTorig,4),
                  "( vs",boot,"random trees )",
                  "\n"," percentile = ",
                  percLT)
       )
  abline(v=LTorig, lwd=3, col="blue")
  legend("right", lwd=2,c("tree","random distrib."), col=c("blue","red"))
  
  write("    -Calculating and Plotting root-to-tip distance vs time", stderr())
  
  alldist <- dist.nodes(tree)
  rttd <- NULL
  for (i in 1:length(tree$tip.label)){
    rttd=c(rttd, alldist[length(tree$tip.label)+1, i])
  }
  rttd <- as.data.frame(rttd)
  rttd <- cbind(tree$tip.label, rttd)
  names(rttd) <- c("taxaId","rttd")
  
  data1 <- merge(rttd, data, by.x="taxaId", by.y="taxaId")
  cortest <- cor.test(data1$timeNum,data1$rttd)
  pvalueCOR <- round(cortest$p.val,4)
  {if (pvalueCOR>0.0001) pvalueCOR=paste("=",round(pvalueCOR,4)) else pvalueCOR=paste("<=0.0001")}
  plot(data1$timeNum,data1$rttd,xlab="time", ylab="root-to-tip-dist")
  lines(lowess(data1$timeNum,data1$rttd), lwd=2,col="green")
  title(main = paste("root-to-tip distance vs. time","\n","Pearson's cor. = ",
              round(cortest$estim,4)," ( p-value ",pvalueCOR,")"))
  legend("topleft",c("lowess smoothing"),
         bg="white",col=c("green"),pch=19,cex=1)
  
  #write("After a few more functions", stderr())
  write("    -Calculating and Plotting root-to-tip distance vs time intervals", stderr())
  
  data2 <- as.data.frame(cbind(names(char),char))
  names(data2) <- c("taxaId","char")
  data3 <- merge(data1,data2,by.x="taxaId",by.y="taxaId")
  krustest <- kruskal.test(data3$rttd~data3$char)
  pvalueK <- round(krustest$p.val,4)
  {if (pvalueK>0.0001) pvalueK=paste("=",round(pvalueK,4)) else pvalueK=paste("<=0.0001")}
  
  boxplot(data3$rttd~data3$char, xlab="time points", ylab="root-to-tip-dist")
  title(main = paste("root-to-tip distance vs. time intervals","\n",
              "( p-value ",pvalueK,")"))
  
  write("    -Calculating and Plotting Aldous' Test", stderr())
  
  tree.shape <- as.treeshape(tree)
  aldous(tree.shape, xmin = min(tree.shape$merge))
  title(main = "Aldous' test")
  
  write("    -Calculating and Plotting Sackin's index for random Yule trees", stderr())
  
  statSY <- sackin(tree.shape,norm="yule")
  histSY <- sapply(rtreeshape(boot,tip.number=length(tree.shape$names),
                              model="yule"), FUN=sackin,norm="yule")
  histSY <- sort(histSY)
  percSY <- round(100*max(c(0,which(histSY<statSY)))/length(histSY),2)
  main <- paste("Sackin's index (vs. ", boot,
                " random Yule trees)", "\n percentile = ", percSY)
  xlab <- "Sackin's index"
  hist(histSY, freq=FALSE, main=main, 
       xlab=xlab, col="red", xlim=c(min(c(histSY,statSY)),
                                    max(c(histSY,statSY)))
       )
  abline(v=statSY,lwd=3,col="blue")
  
  write("    -Calculating and Plotting Sackin's index for random PDA trees", stderr())
  
  statSP <- sackin(tree.shape,norm="pda")
  histSP <- sapply(rtreeshape(boot,tip.number=length(tree.shape$names),
                              model="pda"), FUN=sackin, norm="pda")
  histSP <- sort(histSP)
  percSP <- round(100*max(c(0,which(histSP<statSP)))/length(histSP),2)
  main <- paste("Sackin's index (vs. ", boot, " random PDA trees)",
                "\n percentile = ", percSP)
  xlab <- "Sackin's index"
  hist(histSP, freq=FALSE, main=main, xlab=xlab, col="red")
  abline(v=statSP, lwd=3, col="blue")
  
  write("    -Calculating and Plotting Cherry count vs 200 random trees", stderr())
  
  spectre <- spectrum.treeshape(tree.shape)
  cherries <- spectre[length(spectre)]
  histCH <- NULL
  for (i in 1:boot){
    trq <- rtree(length(tree$tip.label))
    trq1 <- as.treeshape(trq)
    val <- spectrum.treeshape(trq1)
    histCH <- c(histCH,val[length(val)])
  }
  histCH <- sort(histCH)
  percCH <- round(100*max(c(0,which(histCH<cherries)))/length(histCH),2)
  hist(histCH,freq=FALSE,xlab="cherry counts",col="red",
       main=paste("Cherry count (vs. ",boot," random trees)",
                  "\n percentile = ",percCH), xlim=c(min(c(histCH,cherries)),
                                                       max(c(histCH,cherries))))
  abline(v=cherries,lwd=3,col="blue")
  dev.off()

  write(paste("TC score = ",round(TC,4)," ( percentile = ",percTC,")", sep=""),stdout())
  write(paste("Staircase-ness = ", round(LTorig,4), 
              " ( percentile = ", percLT,")", sep=""),stdout())
  write(paste("Pearson's cor. = ", round(cortest$estim,4), 
        " ( p-value ",pvalueCOR,")"),stdout())
  summary(tree.shape)
  ctestY <- colless.test(tree.shape, model = "yule")
  ctestP <- colless.test(tree.shape, model = "pda")
  stestY <- sackin.test(tree.shape, model = "yule")
  stestP <- sackin.test(tree.shape, model = "pda")
  likelihood.test(tree.shape, model = "yule", alternative="two.sided")
  likelihood.test(tree.shape, model = "pda", alternative="two.sided")
  paste("Sackin's index (vs. 300 random Yule trees) = ",
        statSY,", percentile = ",percSY)
  paste("Sackin's index (vs. 300 random PDA trees) =",
        statSP,", percentile = ",percSP)
  cherry(tree)
  paste("Cherry count (vs. 300 random trees) = ", 
        cherries, ", percentile = ", percCH)
  PybusGamma <- gammaStat(tree)
  twoTailPybusGamma <- 2*(1 - pnorm(abs(gammaStat(tree))))
  oneTailPybusGamma <- 1 - pnorm(abs(gammaStat(tree)))
  paste("Pybus' Gamma = ", PybusGamma, 
        " ( two-tail p-value = ",twoTailPybusGamma, 
        "; one-tail p-value = ",oneTailPybusGamma," )")

  x <- c("Set", "TimeRange", "TimeIntervals",
         "Tips", "PearsonRho", "TC", 
         "StaircaseNess", "Cherries", "PybusGamma",
         "Colless", "Sackin", "CollessYule",
         "CollessUnif", "SackinYule", "SackinUnif")
  y <- c(treefile, paste(min(data$timeNum), 
                         "-", max(data$timeNum),
                         "time_units"),
         length(levels(char)), 
         length(tree$tip.label), as.numeric(cortest$estim), 
         TC, LTorig,cherries, PybusGamma,
         colless(tree.shape), sackin(tree.shape),
         colless(tree.shape,norm="yule"), 
         colless(tree.shape,norm="pda"),
         sackin(tree.shape,norm="yule"),
         sackin(tree.shape,norm="pda"))
  names(y) <- x
  y <- as.data.frame(t(y))
  outputfilename <- paste(treefile,"phylotempo.txt",sep="_")
  write.table(y, file = outputfilename, sep="\t", row.names=F)
} # end modelCharEvo

randomTest <- function(){
  tree <- rtree(100)
  timeNum <- NULL
  heights <- getLeafHeights(tree)
  for (i in 1:length(tree$tip.label)){
    timeNum[i] <- rnorm(1, mean = heights[i])
  }
  timeNum <- as.integer(timeNum)
  tree <- prepTree(tree, tip.date = timeNum)
  df <- data.frame(timeNum)
  df$taxaId <- tree$tip.label
  write.table(df, file = "random_timetable.txt", row.names = FALSE)
  plot(tree)
  modelCharacterEvolution(tree, timetable = "random_timetable.txt", treefile = "TC.out", boot=200)
}

shankarappaTest <- function(){
  tree <- read.tree("set1_Shankarappa.tree")
  timetable <- "set1_Shankarappa.txt"
  fdata <- formatData(timetable, tree)
  data <- fdata$data
  tree <- prepTree(tree, tip.date = data$timeNum)
  modelCharacterEvolution(tree = tree, timetable = "set1_Shankarappa.txt", 
                        treefile = "set1", boot=200, parsimony = T)
}

temporalClustering <- function(tree_file, timetable_file, parsimony=F, bootstrap=200, output="TC", randomMulti2Di=T){
    tree <- read.tree(tree_file)
    fdata <- formatData(timetable_file, tree)
    data <- fdata$data
    tree <- prepTree(tree, tip.date = data$timeNum, random_resolve=randomMulti2Di)
    modelCharacterEvolution(tree = tree, timetable = timetable_file,
                        treefile = output, boot = bootstrap, parsimony = parsimony)
}
