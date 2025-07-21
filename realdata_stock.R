


source("realdata_func.R")

library(fastclime)
library(GGally)
library(igraph)
library(ggpubr)
library(network)

data('stockdata')


stockrate <- matrix(nrow=nrow(stockdata[["data"]])-1, ncol=ncol(stockdata[["data"]]))
for(i in 1:(nrow(stockdata[["data"]])-1)){
  stockrate[i,] <- log(stockdata[["data"]][i+1,]/stockdata[["data"]][i,])
}

n.vec <- c(252, nrow(stockrate)-252)
stockrate.all <- base::rbind(stockrate[(nrow(stockrate) - 251):nrow(stockrate),],
                             stockrate[1:(nrow(stockrate) - 252),])

stockinfo <- stockdata$info
ind.1 <- which(stockinfo[,2]=="Energy")
ind.2 <- which(stockinfo[,2]=="Consumer Staples") 
ind.3 <- which(stockinfo[,2]=="Utilities")
stockrate.p <- stockrate.all[,c(ind.1,ind.2,ind.3)]
stockinfo.p <- stockinfo[c(ind.1,ind.2,ind.3),]
#######
p <- ncol(stockrate.p)
#######
stock.type <- factor(stockinfo.p[,2])
stock.shape <- factor(stockinfo.p[,2])
indc <- levels(stock.type)

re <- te.trans.clime.rd.new(stockrate.p, 3.3, const = 0.5)
re.omega <- re$Omega
diag(re.omega) <- 0
for(i in 1:nrow(re.omega)){
  for(j in 1:ncol(re.omega)){
    if(re.omega[i,j]!=0){
      if(stockinfo.p[i,2] == stockinfo.p[j,2]){
        re.omega[i,j] <- 1
      }else{
        re.omega[i,j] <- 2
      }
    }
  }
}
nEn <- graph_from_adjacency_matrix(re.omega, mode="undirected", weight = TRUE)
igraph::V(nEn)$Type <- as.character(stock.type)
nEn.color <- E(nEn)$weight
nEn.color[nEn.color == 1] <- "grey"
nEn.color[nEn.color == 2] <- "red"
igraph::E(nEn)$color <- nEn.color
igraph::V(nEn)$shape <- stockinfo.p[,2]
library(randtoolbox)
pos1 <- halton(length(ind.1), 2)
pos2 <- halton(length(ind.2), 2)
pos2[,1] <- pos2[,1]+1
pos2[,2] <- pos2[,2]+1
pos3 <- halton(length(ind.3), 2)
pos3[,1] <- pos3[,1]+1
pos3[,2] <- pos3[,2]-1
pos <- rbind(pos1, pos2, pos3)


nEn.edge <- E(nEn)$weight
nEn.edge[nEn.edge == 1] <- 0.25
nEn.edge[nEn.edge == 2] <- 1

ggnet2(nEn, size=5, color = "Type", 
       edge.color = nEn.color, edge.size = nEn.edge,
       label = stockinfo.p[,1], label.size = 2, mode=pos, palette = "Set2")



plotnets <- function(data.all, const, n.vec, trans=TRUE){
  re <- te.trans.clime.rd.new(data.all,n.vec, const)
  if(trans){
    re.omega <- re$Omega
  }else{
    re.omega <- re$Omega.ori
  }
  diag(re.omega) <- 0
  for(i in 1:nrow(re.omega)){
    for(j in 1:ncol(re.omega)){
      if(re.omega[i,j]!=0){
        if(stockinfo.p[i,2] == stockinfo.p[j,2]){
          re.omega[i,j] <- 1
        }else{
          re.omega[i,j] <- 2
        }
      }
    }
  }
  nEn <- graph_from_adjacency_matrix(re.omega, mode="undirected", weight = TRUE)
  igraph::V(nEn)$Type <- as.character(stock.type)
  nEn.color <- igraph::E(nEn)$weight
  nEn.color[nEn.color == 1] <- "grey"
  nEn.color[nEn.color == 2] <- "red"
  igraph::E(nEn)$color <- nEn.color
  igraph::V(nEn)$shape <- stockinfo.p[,2]
  library(randtoolbox)
  pos1 <- halton(length(ind.1), 2)
  pos2 <- halton(length(ind.2), 2)
  pos2[,1] <- pos2[,1]+1
  pos2[,2] <- pos2[,2]+1
  pos3 <- halton(length(ind.3), 2)
  pos3[,1] <- pos3[,1]+1
  pos3[,2] <- pos3[,2]-1
  pos <- rbind(pos1, pos2,pos3)
  nEn.edge <- igraph::E(nEn)$weight
  nEn.edge[nEn.edge == 1] <- 0.25
  nEn.edge[nEn.edge == 2] <- 1
  
  ggnet2(nEn, size=5, color = "Type", 
         edge.color = nEn.color, edge.size = nEn.edge,
         label = stockinfo.p[,1], label.size = 2, mode=pos, palette = "Set2")
}


p1 <- plotnets(stockrate.p, 1.1, c(252, nrow(stockrate)-252))
p2 <- plotnets(stockrate.p, 1.1, c(252, nrow(stockrate)-252), trans=FALSE)
p1.1 <- plotnets(stockrate.p, 1.1*2, nrow(stockrate), trans=FALSE)
p3 <- plotnets(stockrate[,c(ind.1,ind.2,ind.3)], 
               0.95, nrow(stockrate))
p4 <- plotnets(stockrate[,c(ind.1,ind.2,ind.3)], 
               0.95, 1258)

p5 <- ggarrange(p1, p2, p1.1, common.legend=TRUE, ncol=3, legend = "bottom")
ggsave("p5.png", p5, width=100, height=50, unit="cm")


p6 <- ggarrange(p1, p2, common.legend=TRUE, ncol=2, legend = "bottom")
p6




