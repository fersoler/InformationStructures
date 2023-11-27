# This file contains some modifications to the file by S. Schreiber published in: 
#
# Schreiber, S. (2022). R Code for the article "Permanence via invasion graphs: 
# Incorporating community assembly into Modern Coexistence Theory" 
# by Josef Hofbauer and Sebastian J. Schreiber in the Journal of Mathematical 
# Biology. Zenodo. https://doi.org/10.5281/zenodo.7111753
#
# Modifications are indicated with ***FST***
# Author of the modifications: Fernando Soler Toscano
# Original code is in the file 'invasion_graph_main_functions.R'

###################
# MAIN DEFINITIONS
###################
  # As discussed in the paper, an invasion graph IG is based on an invasion schemer IS: a m x n matrix where m is the number of communities and n is the total number of species. As the invasion growth rate of a species supported by a community is zero, this code assumes that 0s in the IS correspond to species present in the community. The non-zeros correspond to the invasion growth rates of the missing species where <0 means the missing species can not invade and >0 means that the species can invade. As this code only makes use of the signs of the entries in IS, it suffices (but not required) for the entries of IS to be -1,0,1.

  #An invasion graph IG is a graph with m vertices with directed edges from community i to community j if IS[i,l]>0 for all the species l in community j but not community i, and IS[j,l]<0 for the species l in community i but not community j.

  # A -i community is a subset S of {1,2,..,,n} such that i is not in S and the invasion growth rates of of all other missing species are negative. In terms of the IG, it corresponds to a row whose ith-column is non-zero and all other non-zero column entries are negative.

#######################
# OVERVIEW OF FUNCTIONS
#######################

# LV.IS computes the invasion scheme for a Lotka Volterra system  dx/dt=x*(A%%x+b) with any number of species. The resulting IS can be passed to the IG.function.  This command assumes that the principle submatrices (i.e. submatrices by removing columns and the corresponding rows) are invertible.

# IG.function creates an invasion graph IG from an invasion scheme IS, checks whether the IG is acyclic, checks the permanence condition for all communities if it is acyclic, and identifies all saturated -i communities.

# IG.quick.function - like the IG.function but doesn't compute all of the -i communities and runs faster.

# IG.plot plots an invasion graph where the vertices are plotted at a height corresponding to the number of species in the community. The directed edges involving a single species invasion are drawn with a thicker edge and other directed edges with a lighter edge.

###################
# REQUIRED PACKAGES
###################
require('igraph')
require('gRbase')
require('dnet')

####################
# THE LV.IS FUNCTION
####################

# ***FST*** Modification of the original function to return not only the IS
# but also the communities.  

# Inputs: The interaction matrix A and intrinsic growth rate vector b for dx/dt=x*(A%%x+b). tolerance determines what numbers correspond to a numerical zero.
# Outputs: The interaction scheme IS with rows corresponding to communities and columns species. Entries are average per-capita growth rates.
LV.ISandComm=function(A,b,tolerance=1e-14){
  k=dim(A)[1]
  C=list() # communities (as equilibria)
  C[[1]]=rep(0,k)
  no.C=1
  # find all subcommunities
  for(i in 1:k){
    temp=combn(1:k,i)
    k2=dim(temp)[2]
    for(j in 1:k2){
      I=temp[,j]
      xtemp=solve(a = A[I,I],b=-b[I])
      if(min(xtemp)>0){
        no.C=no.C+1
        xtemp2=C[[1]]
        xtemp2[I]=xtemp
        C[[no.C]]=xtemp2
      }
    }
  }
  # create the invasion scheme matrix
  IS=matrix(NA,length(C),k)
  for(i in 1:length(C)){
    IS[i,]=A%*%C[[i]]+b
  }
  IS[which(abs(IS)<tolerance)]=0
  list(
    IS=IS,
    Comm=t(matrix(unlist(C),nrow=k))
  )
}

##########################
# THE IG.function FUNCTION
##########################

# ***FST*** Adapted to the new function LV.ISandComm

# Input: the invasion scheme matrix where rows correspond to communities, columns species, and entries average per-capita growth rates. Make sure that zeros correspond to species present in the community.
# Output: A list where
# IG is the invasion graph: a square matrix where rows and columns corresponds to communities, zeros correspond to no directed edge, and ones correspond to a directed edge from row entry to column entry.
# acyclic which is TRUE is the IG is acyclic
# permanent which lists all the permanent communities - only meaningful if IG is acyclic
# minus.i which lists all of the -i communities
# number.species which lists the species richness of all the communities
# composition which lists the species compositions of each community
# IS which returns the invasion scheme.
IG.functionComm=function(ISComm){
  IS = ISComm$IS
  Comm = ISComm$Comm
  n=dim(IS)[1] # number of communities
  composition=list() # a list to hold the composition of the communities
  # for(i in 1:n)composition[[i]]=which(IS[i,]==0) #  ***FST*** Original
  for(i in 1:n)composition[[i]]=which(Comm[i,]>0)  #  ***FST*** Modified
  # compute the invasion graph
  IG=matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      if(i!=j){
        b=composition[[j]]
        a=composition[[i]]
        c=setdiff(b,a)
        condition.1=1
        if(length(c)>0)condition.1=min(IS[i,c])
        c2=setdiff(a,b)
        condition.2=-1
        if(length(c2)>0)condition.2=max(IS[j,c2])
        if((condition.1>0)&(condition.2<0))IG[i,j]=1
      }
    }}
  # computing the graph to determine whether acyclic
  g=graph_from_adjacency_matrix(IG)
  # determine which communities are permanent
  # only sensible if the graph is acyclic
  permanent=matrix(TRUE,n,1)
  for(i in 1:n){
    for(j in 1:n){
      if((j!=i)&&all(composition[[j]] %in% composition[[i]])){if((max(IS[j,composition[[i]]])<=0)){permanent[i]=FALSE}}
    }
  }
  # identify the minus i communities
  minus.i=list()
  k=dim(IS)[2] # number of species
  for(i in 1:k){
    temp=c()
    for(j in 1:n){
      if((!is.element(i,composition[[j]]))&&(max(IS[j,-i])<=0))temp=c(temp,j)
    }
    minus.i[[i]]=temp
  }
  # count the number of species in each community
  # numbers of spp
  nums=numeric(n)
  for(i in 1:n)nums[i]=length(which(IS[i,]==0))
  
  return(list(IG=IG,acyclic=is_dag(g),permanent=permanent,minus.i=minus.i,number.species=nums,composition=composition,IS=IS))
}

######################
# THE plot.IG FUNCTION
######################
# Inputs: the output (out) from the IG.function command
# clear.out the list of communities that shouldn't be plotted
# bend.factor determines the convexity of the curve on which the vertices with the same species richess are plotted. Takes on values between 0 (no bending) to 1 (most bending)
# cols a vector of seven colors corresponding to
# cols[1] - edges S->T corresponding to set.diff(T,S)>1
# cols[2] - edges S->T corresponding to set.diff(T,S)<2
# cols[3] - background vertex color for an invasible -i community
# cols[4] - background vertex color for a non invasible -i community
# cols[5] - default color of vertex background
# cols[6] - default color of vertex labels and edge
# cols[7] - edge color for arrow leaving -i community
# impermanent.weight - what fraction of the nonstandard vertex background color to use for impermanent communities.
# multiple.invasion.edge.weight - weight on edges S->T where set.diff(T,S)>1

plot.IG=function(out,clear.out=c(),bend.factor=0.75,cols=c("lightgray","darkgray","darkolivegreen3","gold","white","black","darkolivegreen3"),impermanent.weight=0.172549,multiple.invasion.edge.weight=0.5){
  IG=out$IG
  IS=out$IS
  number.species=out$number.species
  minus.i=out$minus.i
  permanent=out$permanent
  composition=out$composition
  not.permanent=which(!out$permanent)
  # create graph
  g=graph_from_adjacency_matrix(IG)
  # get edge vertex info: column 1 = tail, 2 = tip
  edges_data_frame <- get.data.frame(g, what = "edges")
  # create edge widths based on one versus multiple invasion
  for(i in 1:length(edges_data_frame[,1])){
    E(g)$weight[i]=multiple.invasion.edge.weight
    E(g)$color[i]=cols[1]
    SS=composition[[edges_data_frame[i,1]]]
    TT=composition[[edges_data_frame[i,2]]]
    if(length(setdiff(TT,SS))<2){
      E(g)$weight[i]=2
      E(g)$color[i]=cols[2]}
    if(is.element(edges_data_frame[i,2],clear.out))E(g)$color[i]=NA
  }
  # create the reordered version of the graph (by edge weight that is used ultimately for the plotting
  g2 <- make_graph(as.vector(t(get.edgelist(g)[order(E(g)$weight),])))
  k=length(number.species) # number of communities
  # create vertex map
  xvals=numeric(k)
  yvals=numeric(k)
  lengths=numeric(max(number.species)+1)
  for(ii in 0:max(number.species)){
    lengths[ii+1]=length(which(number.species==ii))
  }
  counter=1
  for(i in 0:max(number.species)){
    for(j in 1:lengths[i+1]){
      xvals[counter]=j+(max(lengths)-lengths[i+1])/2
      x.temp=2*xvals[counter]/max(lengths)-1
      yvals[counter]=i-(1-(x.temp)^2)*bend.factor
      counter=counter+1}
  }

  # community names
  v.name=c(expression(symbol("\306")))
  for(i in 2:k){
    temp=composition[[i]]
    v.name=c(v.name,paste(unlist(temp),collapse=""))
  }
  for(i in clear.out)v.name[i]=""
  V(g)$name=v.name

  # minus.i community vertex colors and edge colors
  # for invadable -i, vertex is cols[3] and leaving edge is cols[7]
  # for saturated -i, vertex is cols[5]

  vcols=rep(cols[5],k)
  for(i in 1:dim(IS)[2]){
    for(j in minus.i[[i]]){
      if(IS[j,i]>0){
      vcols[j]=cols[3]
      for(jj in which(edges_data_frame$from==j)){
        SS=composition[j]
        TT=composition[edges_data_frame$to[jj]]
        if(length(setdiff(SS,TT))<2){
          E(g)$weight[jj]=2
          E(g)$color[jj]=cols[7]
          }}
      }else{vcols[j]=cols[4]}
    }
  }

  vertex.frame.cols=rep(cols[6],k)
  vertex.label.cols=rep(cols[6],k)
  vertex.shapes=rep("circle",k)
  for(i in clear.out){
    vcols[i]=NA
    vertex.frame.cols[i]=NA
    vertex.label.cols[i]=cols[5]
    for(j in which(edges_data_frame$to==i)){E(g)$color[j]=NA}
  }

  for(i in not.permanent){
    temp=(col2rgb(vcols[i])*impermanent.weight+col2rgb(cols[5])*(1-impermanent.weight))/255
    vcols[i]=rgb(temp[1],temp[2],temp[3])
    vertex.frame.cols[i]=cols[1]
    vertex.label.cols[i]=cols[2]
  }

  # assign properties to the reordered graph g2
  E(g2)$color <- E(g)$color[order(E(g)$weight)]
  E(g2)$weight<-E(g)$weight[order(E(g)$weight)]
  V(g2)$name <-v.name
  par(mar=c(0,0,0,0))
  plot(g2,vertex.color=vcols,edge.width=E(g2)$weight,
       vertex.frame.color=vertex.frame.cols,vertex.shape=vertex.shapes,vertex.label.color=vertex.label.cols,layout=cbind(xvals,yvals),vertex.label.cex=0.5,edge.arrow.size=0.4)
}

