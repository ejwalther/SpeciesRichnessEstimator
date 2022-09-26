#-----------------------------------------------------------------------------------
#        Function for Simulating Sample Data sets from a multiple watersheds
#       (see single watershed file for argument details)
#         Author: Eric J Walther
#         Last Modified: 9/26/2022
#--------------------------------------------------------------------------------------
watershedsampling<-function(SP.ID,p.cap_G1,p.cap_G2,p.type,n.samples,G.prob=.65,equalsampling=2){
surveycount<-1:n.samples
SampleHist<-NULL #empty vector to hold gear type used in each sampling event
samplelist<-vector(mode='list', length=n.samples) #empty list to store sampling data
N.caught<-NULL 
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
  N<-N<-rnbinom(1,25,mu = 225)
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
      s<-Speciestypetosample[ifelse(equal==1,which(Speciestypetosample[,1]==sample(Speciestypetosample[,1],1)),
                                    which(Speciestypetosample[,1]==sample(Speciestypetosample[,1],1,prob = Speciestypetosample[,4]))),]
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




