#-----------------------------------------------------------------------------------
#        Function for Simulating Sample Data sets 
#       
#         Author: Eric J Walther
#         Last Modified: 9/7/2022
#
#   Simulation details:
#     1. Make the assumption that all species are everywhere
#     2. An individual is randomly drawn from the regional species 
#        pool and is "captured" based on a Bernoulli trail with defined 
#        capture probability p. This random draw can be equally balanced or based on
#        of group types represented in regional species pool. Default = equal
#     3. Currently this function only allows for two different collection methods.
#     4. The number of species caught (N) for each sampling event is independent
#        of collection method and is randomly drawn from Gamma(11,.1) 

#   Input objects:
#     SR - true regional species richness (integer)
#     n.type - number of functional types for grouping species (integer)
#     types_desc - description of function types. Must be same length as n.type (character)
#     p.cap_G1 - capture probabilities for each group type using collection method 1. Must be same length as n.type (numeric)
#     p.cap_G2 - capture probabilities for each group type using collection method 2. Must be same length as n.type (numeric)    
#     p.type - proportion of each functional type present in regional species pool SR (numeric)
#     n.samples - number of sampling events (integer)
#     G.prob - probability of gear type 1 being selected; default set to 0.65
#     equalsampling - determines if randomly draw species from pool is equally balances (==1) or based on proportional
#                     representation n regional species pool (==0). Default==1

#   Return objects in list:
#     1. SampleRecords - list of sampling event capture data matrices
#                        SpeciesID: unique identification code for captured species
#                        Count: number of individuals captured
#                        GT: Gear type used during sampling event
#                        Sample: Sampling event number
#     2. GearTypes - gear types used for each sampling event
#     3. Abundance - total number of individuals captured during sampling event
#     4. SpeciesID - Species ID table
#                    [,1] = unique species ID
#                    [,2] = species type
#                    [,3] = proportion of species group type in species pool
#     5. SpeciesObs - vector of all unique species detected across all sampling events
#     6. SpObservation - number of times species x was captured across sampling events
#     7. CaptureHistory - Capture history incidence matrix. Columns represent sampling events
#                        and rows represent species. Column 1 contains speciesID and last column
#                        contains species group type
#     8. Nsamples - number of sampling events in dataset
#     9. Ntypes - number of group types in dataset
#------------------------------------------------------------------------------------
sim.sample<-function(SR,n.type,types_desc,p.cap_G1,p.cap_G2,p.type,n.samples,G.prob=.65,equalsampling=2){
typeID<-1:n.type #ID code for each group

#assign functional group to species list
SP.ID <- matrix(c(1:SR,rep(NA,SR*3)),SR,4)
randvar <-  sample(typeID,size = SR,prob = p.type,replace = TRUE)
SP.ID[,2] <- randvar
for(i in 1:nrow(SP.ID)){
  SP.ID[i,3] <- p.type[SP.ID[i,2]]
}
for(t in 1:n.type){
  speciesgroup <- subset(SP.ID,SP.ID[,2]==t)
  sampledspp<-sample(speciesgroup[,1],nrow(speciesgroup),replace = FALSE)
  weights<-dlnorm(seq(1:nrow(speciesgroup)))
  spp.w<-cbind(sampledspp,weights)
  for(i in 1:nrow(spp.w)){
    SP.ID[which(SP.ID[,1]==speciesgroup[i,1]),4]<-spp.w[i,2]
  }
}

n.samples<-n.samples
surveycount<-1:n.samples
SampleHist<-NULL #empty vector to hold gear type used in each sampling event
samplelist<-vector(mode='list', length=n.samples) #empty list to store sampling data
N.caught<-NULL #empty vector to hold number of individuals 'caught' in each sample
#simulate sampling
for(i in 1:n.samples){
  if(rbinom(1,1,G.prob)==1){ #change this to alter proportion of gear types
    G<-p.cap_G1
    SampleHist[i]<-1
  }
  else{
    G<-p.cap_G2
    SampleHist[i]<-2
  }
  CH<-NULL
  N<-as.integer(rgamma(1,11,.1))
  N.caught[i]<-N
  if(N==0){  #currently 0 catch is not possible based on distribution of N. Could play around with dist if we want to include possiblity of 0s
    samplelist[[i]]<-data.frame(NA,0)
    next
  }
  else{
    N.c<-N+1
    c<-1
    while(c<N.c){
      equal<-equalsampling
      S.type<-sample(typeID,size = 1,prob = p.type)
      Speciestypetosample<-subset(SP.ID,SP.ID[,2]==S.type)
      s<-Speciestypetosample[ifelse(equal==1,which(Speciestypetosample[,1]==sample(Speciestypetosample[,1],1)),which(Speciestypetosample[,1]==sample(Speciestypetosample[,1],1,prob = Speciestypetosample[,4]))),]
      #s<-SP.ID[sample(SP.ID[,1],1),] #equally weighted
      #s<-SP.ID[sample(SP.ID[,1],1,prob =SP.ID[,3]),] #weighted based on abundance distribution
      caught<-s[1]*rbinom(1,1,G[s[2]])
      if(caught==0){
        next
      }
      else{
        CH[c]<-caught
      }
      c=length(CH)+1
    }
    catch<-as.data.frame(table(CH))
    names(catch)<-c("SpeciesID","count")
    samplelist[[i]]<-catch
  }
}
#compile sampling list
for(i in 1:n.samples){
  samplelist[[i]][,3]<-SampleHist[i]
  samplelist[[i]][,4]<-surveycount[i]
  names(samplelist[[i]])<-c("SpeciesID","Count","GT","Sample")
}

#create species lists
t_0<-unique(samplelist[[1]][,"SpeciesID"])
SL<-t_0
if(n.samples>1){
for(i in 2:n.samples){
  new.spp<-filter(samplelist[[i]],SpeciesID%nin%SL)[,"SpeciesID"]
  SL<-c(SL,new.spp)
}
}

#Create an incidence capture history matrix
CH.mat<-matrix(SL,ncol = n.samples+1,nrow = length(SL))
CH.mat[,2:(n.samples+1)]<-NA
for(i in 1:n.samples){
  for(r in 1:nrow(CH.mat)){
    CH.mat[r,(i+1)]<-ifelse(CH.mat[r,1]%in%samplelist[[i]][,"SpeciesID"],1,0)
  }
}
CH.mat.df<-as.data.frame(CH.mat)
for(i in 1:ncol(CH.mat.df)){
  CH.mat.df[,i]<-as.integer(CH.mat.df[,i])
}
CH.mat<-as.matrix(CH.mat.df)
obslist<-apply(CH.mat[,-1],1,sum) #output
SP.TyCH<-NULL
for(i in 1:nrow(CH.mat)){
  SP.TyCH[i]<-SP.ID[which(SP.ID[,1]==CH.mat[i,1]),2]
}
CH.matID<-cbind(CH.mat,SP.TyCH)

#output
return(list(SampleRecords = samplelist,
            GearTypes = SampleHist,
            Abundance = N.caught,
            SpeciesID = SP.ID,
            SpeciesObs = SL,
            SpObservation = obslist,
            CaptureHistory = CH.matID,
            Nsamples = n.samples,
            Ntypes = length(typeID)
            ))
}


