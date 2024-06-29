
###############################################################################
# R code to produce results for manuscript 
# "Evidence for differences in patterns of temporal trends in meta-analyses of 
# diagnostic accuracy studies in the Cochrane Database of Systematic Reviews"
# (Manuscript revision 1, 28 June 2024)

# Authors/code written by: Jacqueline Murphy & Thomas Fanshawe
# Nuffield Department of Primary Care Health Sciences, University of Oxford
###############################################################################

# NOTE
# Code for bivariate meta-analysis using lme4/glmer (Method 2) is adapted from:
# Takwoingi Y, Dendukuri N, Schiller I, Rücker G, Jones HE, Partlett C, Macaskill P.
# Supplementary material 1 to Chapter 10: Code for undertaking meta-analysis. In: Deeks JJ,
# Bossuyt PM, Leeflang MM, Takwoingi Y (editors). Cochrane Handbook for Systematic
# Reviews of Diagnostic Test Accuracy. Version 2.0 (updated July 2023). Cochrane, 2023.
# Available from https://training.cochrane.org/handbook-diagnostic-test-accuracy/current



# DATA:  "CD..... .csv" (according to the CDSR study number)
# The analysis data files downloaded from the Cochrane Database of Systematic Reviews (CDSR) formatted with variable names:
# Code, Test, Study.ID, Year.of.study, TP, FP, TN, FN, Primary (1/0 to indicate primary meta-analysis for the review), Sens, Spec

# DATA:  "StudyInfo.csv"
# contains basic information about each review with variable names:
# Code	Author	Year	CochraneReviewGroup

###############################################################################

# Setup -----
library(ggplot2)
library(dplyr)
library(meta)
library(mada)
library(lme4)
library(optimx)
library(nlme)
library(astsa)
library(ncar)
library(gridExtra)
library(ggplotify)
library(cowplot)

rm(list=ls())
setwd("") 


###############################################################################
# Section 1. Extract basic study info -----
d<-dir()
d<-d[substr(d,1,2)=="CD"]
studynames<-substr(d,1,8)
ld<-length(studynames)

IS<-data.frame(matrix(NA,ld,13))
names(IS)<-c("Code","Test","n.studies","n.participants","date.first","date.last","yr.span",
             "n.deleted.studies", "n.studies.after.deletions","n.participants.after.deletions", 
             "date.first.after.deletions", "date.last.after.deletions",
             "yr.span.after.deletions")

i=1
for (i in 1:ld){
  z<-read.csv(d[i],header=TRUE,stringsAsFactors=FALSE)
  z<-z[,1:9] # Code	Test	Year of study	Primary TP FP TN FN 
  z$Sens<-z$TP/(z$TP+z$FN) # calculate sens and spec of each study within the selected MA
  z$Spec<-z$TN/(z$TN+z$FP)
  assign(substr(d[i],1,8),z[,1:11]) # name data set with review name
  z<-z[z$Primary==1,] #  keep only the primary meta analysis in this MA
  assign(paste0(substr(d[i],1,8),"p"),z[,1:11]) # name data set for the primary MA
  
  # output information about the primary MA
  IS$Code[i]<-z$Code[1] 
  IS$Test[i]<-z$Test[1] 
  IS$n.studies[i]<-length(z$Study.ID)
  IS$n.participants[i]<-cumsum(z$TP)[nrow(z)]+cumsum(z$FP)[nrow(z)]+cumsum(z$TN)[nrow(z)]+cumsum(z$FN)[nrow(z)]
  IS$date.first[i]<-min(z$Year.of.study)
  IS$date.last[i]<-max(z$Year.of.study)
  
  z2<-z[(z$TP+z$FN)>0 & (z$TN+z$FP)>0,]
  IS$n.deleted.studies[i]<-length(z$Study.ID) - length(z2$Study.ID)
  IS$n.studies.after.deletions[i]<-length(z2$Study.ID)
  IS$n.participants.after.deletions[i]<-cumsum(z2$TP)[nrow(z2)]+cumsum(z2$FP)[nrow(z2)]+cumsum(z2$TN)[nrow(z2)]+cumsum(z2$FN)[nrow(z2)]
  IS$date.first.after.deletions[i]<-min(z2$Year.of.study)
  IS$date.last.after.deletions[i]<-max(z2$Year.of.study)
  
}
IS$yr.span<-IS$date.last - IS$date.first
IS$yr.span.after.deletions<-IS$date.last.after.deletions - IS$date.first.after.deletions

rm(list = c("i","z","z2"))


## Summary of individual studies --------
studyinfo<-read.csv("StudyInfo.csv")[,1:4] %>% 
  relocate(c("CochraneReviewGroup","Author","Year"), .before="Code") %>%
  arrange(Code)

IS2<-data.frame(merge(studyinfo,IS,by.x = c("Code"))) %>%
  relocate(c("CochraneReviewGroup","Author","Year"), .before="Code") %>%
  rename("PrimaryTest"="Test") %>%
  arrange(CochraneReviewGroup, Author, Year)
colnames(IS2)

write.csv(IS2,"Results/IndividualReviewSummary.csv", row.names=FALSE)


## Aggregate summary of studies table --------
as1<-data.frame(table(IS2$Year))
colnames(as1)=c("Publication Year","Freq")
as2<-data.frame(table(IS2$CochraneReviewGroup))
colnames(as2)=c("Cochrane Review Group","Freq")
n.studies.after.deletions<-IS2 %>%
  summarize(med = Round(sum(n.studies.after.deletions),0))
as3a<-data.frame(" "="Number of studies in primary meta-analysis (N)",
                 "MedianIQR"=paste0(n.studies.after.deletions$med))
iqr.n.studies.after.deletions<-IS2 %>%
  summarize(med = Round(median(n.studies.after.deletions),0),
            Q1 = Round(quantile(n.studies.after.deletions, 0.25)),
            Q3 = Round(quantile(n.studies.after.deletions, 0.75)))
as3b<-data.frame(" "="Number of studies in primary meta-analysis (Median [IQR])",
                 "MedianIQR"=paste0(iqr.n.studies.after.deletions$med," [",iqr.n.studies.after.deletions$Q1,", ",iqr.n.studies.after.deletions$Q3,"]"))
n.participants.after.deletions<-IS2 %>%
  summarize(med = Round(sum(n.participants.after.deletions),0))
as4a<-data.frame(" "="Number of participants in primary meta-analysis (N)",
                 "MedianIQR"=paste0(n.participants.after.deletions$med))
iqr.n.participants.after.deletions <-IS2 %>%
  summarize(med = Round(median(n.participants.after.deletions)),
            Q1 = Round(quantile(n.participants.after.deletions, 0.25)),
            Q3 = Round(quantile(n.participants.after.deletions, 0.75)))
as4b<-data.frame(" "="Number of participants in primary meta-analysis (Median [IQR])",
                 "MedianIQR"=paste0(iqr.n.participants.after.deletions$med," [",iqr.n.participants.after.deletions$Q1,", ",iqr.n.participants.after.deletions$Q3,"]"))
iqr.yr.span.after.deletions <-IS2 %>%
  summarize(med = Round(median(yr.span.after.deletions)),
            Q1 = Round(quantile(yr.span.after.deletions, 0.25)),
            Q3 = Round(quantile(yr.span.after.deletions, 0.75)))
as5<-data.frame(" "="Time range of studies in primary meta-analysis (years) (Median [IQR])",
                "MedianIQR"=paste0(iqr.yr.span.after.deletions$med," [",iqr.yr.span.after.deletions$Q1,", ",iqr.yr.span.after.deletions$Q3,"]"))


write.table(data.frame(as1),"Results/AggregateReviewSummary.csv", sep=",", row.names=FALSE, append=FALSE)
cat("\n", file = "Results/AggregateReviewSummary.csv", append = TRUE)
write.table(data.frame(as2),"Results/AggregateReviewSummary.csv", sep=",", row.names=FALSE, append=TRUE)
cat("\n", file = "Results/AggregateReviewSummary.csv", append = TRUE)
write.table(data.frame(as3a),"Results/AggregateReviewSummary.csv", sep=",", row.names=FALSE, append=TRUE)
cat("\n", file = "Results/AggregateReviewSummary.csv", append = TRUE)
write.table(data.frame(as3b),"Results/AggregateReviewSummary.csv", sep=",", row.names=FALSE, append=TRUE)
cat("\n", file = "Results/AggregateReviewSummary.csv", append = TRUE)
write.table(data.frame(as4a),"Results/AggregateReviewSummary.csv", sep=",", row.names=FALSE, append=TRUE)
cat("\n", file = "Results/AggregateReviewSummary.csv", append = TRUE)
write.table(data.frame(as4b),"Results/AggregateReviewSummary.csv", sep=",", row.names=FALSE, append=TRUE)
cat("\n", file = "Results/AggregateReviewSummary.csv", append = TRUE)
write.table(data.frame(as5),"Results/AggregateReviewSummary.csv", sep=",", row.names=FALSE, append=TRUE)


rm(list = c("d","as1","as2","as3a","as3b","as4a","as4b","as5",
            "n.studies.after.deletions","iqr.n.studies.after.deletions",
            "iqr.n.participants.after.deletions",
            "iqr.yr.span.after.deletions","IS","IS2"))




###############################################################################
# Section 2. Sequential bivariate meta-analysis for each review -----
# Two methods are included: normal approximation (Reitsma 2005) and binomial likelihood (Chu 2006/Chu 2010)
# The option MAmethod allows user to select which is used for the analysis.

# If years are the same for more than one study, randomly order them (2nd argument to order function)
# poolfirst5 = 0 # first MA on only the first study
# poolfirst5 = 1 # fist MA pools the first 5 studies (but not ALL studies in the 5th pub year)
# poolfirst5 = 2 # DEFAULT first MA pools all studies up to the pub year of the 5th study

d<-ls()

seq.bma<-function(poolfirst5=2, MAmethod=2){

  d<-d[substr(d,1,2)=="CD" & substr(d,9,9)=="p"] # select only the data from primary MAs
  counter<-0
  for (i in 1:length(d)){
    counter<-counter+1
    x<-eval(str2expression(d[i]))
    
    cat("\n", x$Code[1]," \n") # , i," "
    
    x<-x[(x$TP+x$FN)>0 & (x$TN+x$FP)>0,] # remove studies where these sums are zero
    
    set.seed(240524) # set seed inside the seq.bma loop - reset seed for every new review so that the results do not depend on how many/which reviews are being processed using seq.bma
    x<-x[order(x$Year.of.study,runif(nrow(x))),] # random ordering of studies within a year
    
    # calculate cumulative sample size
    x$sampsize<-cumsum(x$TP) + cumsum(x$TN) + cumsum(x$FP) + cumsum(x$FN)
    
    # calculate quantities needed for Method 2 MA (lme4/glmer)
    x$n1 <- x$TP+x$FN
    x$n0 <- x$FP+x$TN
    x$true1 <- x$TP
    x$true0 <- x$TN
    x$recordid <- seq(1,nrow(x),1)
    
    # select the number of studies to be included in the first cMA
    if (poolfirst5==2) (x.start<-max((1:nrow(x))[x$Year.of.study==x$Year.of.study[5]])) 
    else if (poolfirst5==1) (x.start<-5) 
    else (x.start<-1)
    
    if (x$Code[1]=="CD013387" | x$Code[1]=="CD012768") (x.start <- x.start + 1) # start with the second cMA as the specificity results are very similar for the first cMA
    
    save.reitsma<-data.frame(matrix(NA,nrow(x)-x.start+1,17))
    names(save.reitsma)<-c("Code","n.studies",
                           "Sens","Logit.Sens","Logit.Sens.se","Sens.ci.l","Sens.ci.u",
                           "Spec","Logit.Spec","Logit.Spec.se","Spec.ci.l","Spec.ci.u",
                           "Order","Year","EndYear","YearRank","sampsize"
    )
    z<-0
    for (j in x.start:nrow(x)){
      
      # run the meta-analysis
      cat(j ," ")
      z<-z+1
      
      X = x[1:j,] # select the first j studies to be included in this MA
      
      save.reitsma[j-x.start+1,"Code"]<-x$Code[1] # Code
      save.reitsma[j-x.start+1,"n.studies"]<-j # n.studies
      
      # Calculate cumulative MA. The MAmethod parameter determines which results will be used for the final plots
      
      if (MAmethod==1) {
        colnames(save.reitsma)
        
        # Method 1: Reitsma 2005
        S<-summary(reitsma(X)) # run MA including studies up to the selected study j
        save.reitsma[j-x.start+1,c("Sens","Spec","Sens.ci.l","Spec.ci.l","Sens.ci.u","Spec.ci.u")]<-c(S$coefficients[3:4,c(1,5,6)])  
        save.reitsma[j-x.start+1,c("Logit.Sens", "Logit.Sens.se")]<-c(S$coefficients[1,1:2])  
        save.reitsma[j-x.start+1,c("Logit.Spec","Logit.Spec.se")]<-c(S$coefficients[2,1:2]) 
        
      }
      else {
        
        # Method 2: Chu 2006/Chu 2010 using lme4
        # prepare data (j studies)
        Y = reshape(X, direction="long", varying=list(c("n1", "n0"), c("true1",
                                                                       "true0")), timevar="sens", times=c(1,0), v.names=c("n","true"))
        Y = Y[order(Y$id),]
        Y$spec<- 1-Y$sens 
        
        # run the MA (j studies)
        
        if (Y$Code[1]=="CD013021"){ 
          # for this study specificity = 1 by definition (FP=0) so perform univariate RE logistic regression
          (MA_Y_sens = glmer(formula=cbind(true, n - true) ~ 0 + sens + (0+sens |Study.ID), data=Y, family=binomial, nAGQ=1,
                             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e4)))) #, verbose=1)
          S2_sens<-summary(MA_Y_sens)
          S2_sens
          
          # extracting relevant data
          Logit.Sens = S2_sens$coeff[1,1] 
          Logit.Sens.se = S2_sens$coeff[1,2] 
          Logit.Spec = NA 
          Logit.Spec.se = NA   
          
        }
        else{
          # choose a different optimisation option to solve convergence issues in these reviews
          if (Y$Code[1]=="CD012768" | Y$Code[1]=="CD012806" | 
              Y$Code[1]=="CD013186" | Y$Code[1]=="CD013346" | 
              Y$Code[1]=="CD013362"){
            (MA_Y = glmer(formula=cbind(true, n - true) ~ 0 + sens + spec + (0+sens + spec |Study.ID), data=Y, family=binomial, nAGQ=1,
                          control=glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))))
          }

          else{
            (MA_Y = glmer(formula=cbind(true, n - true) ~ 0 + sens + spec + (0+sens + spec |Study.ID), data=Y, family=binomial, nAGQ=1,
                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e4)))) 
          }
          
          S2<-summary(MA_Y)
          S2
          
          # extracting relevant data
          Logit.Sens = S2$coeff[1,1] 
          Logit.Sens.se = S2$coeff[1,2] 
          Logit.Spec = S2$coeff[2,1] 
          Logit.Spec.se= S2$coeff[2,2] 
        }
        
        Logit.Sens.ci.l = Logit.Sens-qnorm(0.975)*Logit.Sens.se
        Logit.Sens.ci.u = Logit.Sens+qnorm(0.975)*Logit.Sens.se
        Logit.Spec.ci.l = Logit.Spec-qnorm(0.975)*Logit.Spec.se
        Logit.Spec.ci.u = Logit.Spec+qnorm(0.975)*Logit.Spec.se
        
        Sens = exp(Logit.Sens)/(1+exp(Logit.Sens))
        Sens.ci.l = exp(Logit.Sens.ci.l)/(1+exp(Logit.Sens.ci.l))
        Sens.ci.u = exp(Logit.Sens.ci.u)/(1+exp(Logit.Sens.ci.u))
        Spec = exp(Logit.Spec)/(1+exp(Logit.Spec))
        Spec.ci.l = exp(Logit.Spec.ci.l)/(1+exp(Logit.Spec.ci.l))
        Spec.ci.u = exp(Logit.Spec.ci.u)/(1+exp(Logit.Spec.ci.u))
        
        if (Y$Code[1]=="CD013021"){ 
          Spec = 1
          Spec.ci.l = NA
          Spec.ci.u = NA
        }
        
        save.reitsma[j-x.start+1,"Sens"] = Sens 
        save.reitsma[j-x.start+1,"Logit.Sens"] = Logit.Sens
        save.reitsma[j-x.start+1,"Logit.Sens.se"] = Logit.Sens.se
        save.reitsma[j-x.start+1,"Sens.ci.l"] = Sens.ci.l
        save.reitsma[j-x.start+1,"Sens.ci.u"] = Sens.ci.u
        save.reitsma[j-x.start+1,"Spec"] = 1 - Spec # to match the output of method 1, need to modify these which are then re-converted later in the code
        save.reitsma[j-x.start+1,"Logit.Spec"] = -1*Logit.Spec # to match the output of method 1, need to modify these which are then re-converted later in the code
        save.reitsma[j-x.start+1,"Logit.Spec.se"] = -1*Logit.Spec.se # to match the output of method 1, need to modify these which are then re-converted later in the code
        save.reitsma[j-x.start+1,"Spec.ci.l"] = 1 - Spec.ci.l # to match the output of method 1, need to modify these which are then re-converted later in the code
        save.reitsma[j-x.start+1,"Spec.ci.u"] = 1 - Spec.ci.u # to match the output of method 1, need to modify these which are then re-converted later in the code
        
      }
      
      save.reitsma[j-x.start+1,"Order"]<-(z-1)
      save.reitsma[j-x.start+1,"Year"]<-x$Year.of.study[j]
      save.reitsma[j-x.start+1,"sampsize"]<-x$sampsize[j]
    }
    save.reitsma[,c(8,11,12)]<-1-save.reitsma[,c(8,11,12)]   
    save.reitsma[,9]<-(-1*save.reitsma[,9])   
    save.reitsma$EndYear<-1
    for (j in 1:(nrow(save.reitsma)-1)){
      if (save.reitsma$Year[j]==save.reitsma$Year[j+1]) (save.reitsma$EndYear[j]<-0)
    }
    
    # Reset the Order variable to be equal for all studies published in the same year
    # (and equal to the highest rank within those studies)
    years<-unique(save.reitsma$Year)
    for (j in 1:length(years)){
      save.reitsma$Order[save.reitsma$Year==years[j]]<-max(save.reitsma$Order[save.reitsma$Year==years[j]])
    }
    
    # calculate YearRank within each review
    save.reitsma$YearRank<-ifelse(save.reitsma$EndYear==1, cumsum(save.reitsma$EndYear),cumsum(save.reitsma$EndYear)+1)
    
    if (counter==1) (outp<-save.reitsma) else (outp<-rbind(outp,save.reitsma))
  }
  
  outp$remove<-ifelse((outp$Code=="CD009185" & outp$n.studies==15) | 
                        (outp$Code=="CD013483" & outp$n.studies==9) | 
                        (outp$Code=="CD014546" & (outp$n.studies==19 | outp$n.studies==20)),1,0)
  
  outp <- outp[outp$remove!=1,]
  
  outp
  
}


## Generate cMA results for each review ----
d<-ls()
save.reitsma<-seq.bma(poolfirst5=2, MAmethod=2) # for main analysis pool the first 5 studies
# see  sensitivity analysis for pooling all within a publication year (see end of code)

a<-unique(save.reitsma$Code)
n<-length(a)
first.last<-NULL
first.last<-data.frame(matrix(NA,n,24))
names(first.last)<-c("Code","n","n.EndYear","n.FirstPooled","Sens.FirstPooled","Sens.FirstPooled.ci.l","Sens.FirstPooled.ci.u",
                     "Logit.Sens.FirstPooled","Logit.Sens.SE.FirstPooled",
                     "Sens.Last","Sens.Last.ci.l","Sens.Last.ci.u",
                     "Logit.Sens.Last","Logit.Sens.SE.Last",
                     "Spec.FirstPooled","Spec.FirstPooled.ci.l","Spec.FirstPooled.ci.u",
                     "Logit.Spec.FirstPooled","Logit.Spec.SE.FirstPooled",
                     "Spec.Last","Spec.Last.ci.l","Spec.Last.ci.u",
                     "Logit.Spec.Last","Logit.Spec.SE.Last")
for (i in 1:n){
  first.last[i,]<-c(a[i],max(save.reitsma$n.studies[save.reitsma$Code==a[i]]),
                    nrow(save.reitsma[save.reitsma$Code==a[i] & save.reitsma$EndYear==1,]),
                    head(save.reitsma$n.studies[save.reitsma$Code==a[i]],1),
                    head(save.reitsma$Sens[save.reitsma$Code==a[i]],1),
                    head(save.reitsma$Sens.ci.l[save.reitsma$Code==a[i]],1),
                    head(save.reitsma$Sens.ci.u[save.reitsma$Code==a[i]],1),
                    head(save.reitsma$Logit.Sens[save.reitsma$Code==a[i]],1),
                    head(save.reitsma$Logit.Sens.se[save.reitsma$Code==a[i]],1),
                    tail(save.reitsma$Sens[save.reitsma$Code==a[i]],1),
                    tail(save.reitsma$Sens.ci.l[save.reitsma$Code==a[i]],1),
                    tail(save.reitsma$Sens.ci.u[save.reitsma$Code==a[i]],1),                
                    tail(save.reitsma$Logit.Sens[save.reitsma$Code==a[i]],1),
                    tail(save.reitsma$Logit.Sens.se[save.reitsma$Code==a[i]],1),
                    head(save.reitsma$Spec[save.reitsma$Code==a[i]],1),
                    head(save.reitsma$Spec.ci.l[save.reitsma$Code==a[i]],1),
                    head(save.reitsma$Spec.ci.u[save.reitsma$Code==a[i]],1),
                    head(save.reitsma$Logit.Spec[save.reitsma$Code==a[i]],1),
                    head(save.reitsma$Logit.Spec.se[save.reitsma$Code==a[i]],1),
                    tail(save.reitsma$Spec[save.reitsma$Code==a[i]],1),
                    tail(save.reitsma$Spec.ci.l[save.reitsma$Code==a[i]],1),
                    tail(save.reitsma$Spec.ci.u[save.reitsma$Code==a[i]],1),
                    tail(save.reitsma$Logit.Spec[save.reitsma$Code==a[i]],1),
                    tail(save.reitsma$Logit.Spec.se[save.reitsma$Code==a[i]],1))
}
colnames(first.last)
first.last[,2:24] <- sapply(first.last[,2:24],as.numeric)

# Calculate difference between first and last MA estimates
first.last<-first.last %>%
  mutate(Sens.Diff = Sens.Last - Sens.FirstPooled,
         Spec.Diff = Spec.Last - Spec.FirstPooled)

# format CIs
CIs<-first.last %>%
  mutate(Sens.FirstPooled.ci=paste0("(",sprintf("%.5s",Sens.FirstPooled.ci.l)," , ",sprintf("%.5s",Sens.FirstPooled.ci.u),")"),
         Sens.Last.ci=paste0("(",sprintf("%.5s",Sens.Last.ci.l)," , ",sprintf("%.5s",Sens.Last.ci.u),")"),
         Spec.FirstPooled.ci=paste0("(",sprintf("%.5s",Spec.FirstPooled.ci.l)," , ",sprintf("%.5s",Spec.FirstPooled.ci.u),")"),
         Spec.Last.ci=paste0("(",sprintf("%.5s",Spec.Last.ci.l)," , ",sprintf("%.5s",Spec.Last.ci.u),")")) %>%
  select("Code","Sens.FirstPooled","Sens.FirstPooled.ci","Sens.Last","Sens.Last.ci",
         "Spec.FirstPooled","Spec.FirstPooled.ci","Spec.Last","Spec.Last.ci","Sens.Diff","Spec.Diff") %>%
  mutate(Sens.FirstPooled=round(as.numeric(Sens.FirstPooled),3),
         Sens.Last=round(as.numeric(Sens.Last),3),
         Spec.FirstPooled=round(as.numeric(Spec.FirstPooled),3),
         Spec.Last=round(as.numeric(Spec.Last),3))

# merge into summary file
FL<-data.frame(merge(studyinfo,CIs,by = c("Code"))) %>%
  relocate(c("CochraneReviewGroup","Author","Year"), .before="Code") %>%
  arrange(Author, Year)




###############################################################################
# Section 3: Set up plotting functions -----

# Function: Plots evolution of diagnostic MA over time within each review in ROC space
# Black point=start (min 5 studies), red point=end (after last study added)
# add.cis add CIs for the first and last pooled estimates

par(mfcol=c(1,2),pty="s")
plot.biv<-function(use=1,print=FALSE,EndYear=TRUE,add.cis=FALSE, title=TRUE){
  n.use<-length(use)
  for (k in 1:n.use){
    a<-unique(save.reitsma$Code)
    sub<-save.reitsma[save.reitsma$Code==a[use[k]],]
    if (print) (print(sub))
    if (EndYear) (sub<-sub[sub$EndYear==1,])
    plot(1-sub$Spec[1],sub$Sens[1],xlim=c(0,1),ylim=c(0,1),pty="s",
         ylab="Sensitivity",xlab="Specificity",
         xaxt="n",yaxt="n",pch=19,cex=0.5)
    axis(1,seq(0,1,by=0.1),labels=seq(1,0,by=-0.1));axis(2,seq(0,1,by=0.1))
    for (j in 2:nrow(sub)){lines(1-sub$Spec[(j-1):j],sub$Sens[(j-1):j])}
    points(1-sub$Spec[nrow(sub)],sub$Sens[nrow(sub)],col=2,pch=19,cex=0.5)
    abline(0,1)
    if (add.cis){
      lines(x=rep(1-sub$Spec[1],2),y=c(sub$Sens.ci.l[1],sub$Sens.ci.u[1]),col=3)
      lines(x=1-c(sub$Spec.ci.l[1],sub$Spec.ci.u[1]),y=rep(sub$Sens[1],2),col=3)
      lines(x=rep(1-sub$Spec[nrow(sub)],2),y=c(sub$Sens.ci.l[nrow(sub)],sub$Sens.ci.u[nrow(sub)]),col=2)
      lines(x=1-c(sub$Spec.ci.l[nrow(sub)],sub$Spec.ci.u[nrow(sub)]),y=rep(sub$Sens[nrow(sub)],2),col=2)}
    if (title) (title(paste("Review",a[use[k]])))
    if (length(use)>1){
      cat(k,a[use[k]])
      readline()
    }
  }
}

#Setting EndYear=TRUE updates the sequential MA only at the end of each year
#Setting EndYear=FALSE updates the sequential MA after each new study (randomly ordering those published in the same year)

plot.sens<-function(use=1,which="both",EndYear=TRUE){
  n.use<-length(use)
  for (k in 1:n.use){
    a<-unique(save.reitsma$Code)
    sub<-save.reitsma[save.reitsma$Code==a[use[k]],]
    if (EndYear) (sub<-sub[sub$EndYear==1,])
    n<-nrow(sub)
    x.start<-min(sub$n.studies)
    x.end<-max(sub$n.studies)
    ymin<-0
    xmin=0
    xmax<-ceiling(x.end / 5) * 5
    a = seq(1,xmax,1)
    c = seq(xmin, xmax,5)
    if (which %in% c("both","sens")){
      plot(0,0,ylim=c(ymin,1),xlim=c(xmin,xmax),type="n",xaxt="n",ylab="Sensitivity",xlab="Number of studies")
      axis(1, at=a, labels=FALSE)
      axis(1, at=c, las=1)
      for (i in 1:n){
        points(sub$n.studies[i],sub$Sens[i],pch=19)
        lines(rep(sub$n.studies[i],2),c(sub$Sens.ci.l[i],sub$Sens.ci.u[i]))
      }}
    if (which %in% c("both","spec")){
      plot(0,0,ylim=c(ymin,1),xlim=c(xmin,xmax),type="n",xaxt="n",ylab="Specificity",xlab="Number of studies")
      axis(1, at=a, labels=FALSE)
      axis(1, at=c, las=1)
      for (i in 1:n){
        points(sub$n.studies[i],sub$Spec[i],pch=19)
        lines(rep(sub$n.studies[i],2),c(sub$Spec.ci.l[i],sub$Spec.ci.u[i]))
      }}
    cat(use[k]," ",a[use[k]])
    if (length(use)>1) (readline())
  }
}



###############################################################################
# Section 4. Trend regression analysis ----
# Function: Use Bagos 2009 method to fit a regression line to the cMAs and overlay onto plots

plot.bagos<-function(use=1,which="sens",data.sens=save.bagos.sens2,data.spec=save.bagos.spec2,
                     title=FALSE, bothtitles=FALSE, add.text=TRUE, add.ar=FALSE,EndYear=TRUE,
                     sens.linear=TRUE, spec.linear=TRUE){
  if (!(use %in% data.sens$Use)) (stop("No model fitted for this ID number"))
  if (which=="both") (which<-c("sens","spec"))
  for (k in 1:length(which)){
    plot.sens(use=use,which=which[k], EndYear)
    if (which[k]=="sens" & sens.linear==TRUE){
      curve(1/(1+exp(-(data.sens$REML.a[data.sens$Use==use]+
                         data.sens$REML.b[data.sens$Use==use]*x))),
            from=data.sens$n.first.MA[data.sens$Use==use],
            to=data.sens$n.studies[data.sens$Use==use],add=TRUE,col=2)
      if (add.text){
        if (add.ar==TRUE){
          p<-paste0("OR=",round(exp(data.sens$REML.b[data.sens$Use==use]),3),
                    " (",round(exp(data.sens$REML.b[data.sens$Use==use]-1.96*data.sens$REML.b.se[data.sens$Use==use]),3),
                    ",",round(exp(data.sens$REML.b[data.sens$Use==use]+1.96*data.sens$REML.b.se[data.sens$Use==use]),3),")")
          text(5,0.05, p,adj=0,cex=0.8)
        }
        if (add.ar==FALSE){
          p<-paste0("OR=",round(exp(data.sens$REML.b[data.sens$Use==use]),3),
                    " (",round(exp(data.sens$REML.b[data.sens$Use==use]-1.96*data.sens$REML.b.se[data.sens$Use==use]),3),
                    ",",round(exp(data.sens$REML.b[data.sens$Use==use]+1.96*data.sens$REML.b.se[data.sens$Use==use]),3),")")
          text(5,0.05, p,adj=0,cex=0.8)
        }
      }
    }
    if (title) (title(paste("Review",data.sens$Code[data.sens$Use==use])))
  }
  if (which[k]=="spec" & spec.linear==TRUE){
    curve(1/(1+exp(-(data.spec$REML.a[data.spec$Use==use]+
                       data.spec$REML.b[data.spec$Use==use]*x))),
          from=data.spec$n.first.MA[data.spec$Use==use],
          to=data.spec$n.studies[data.spec$Use==use],add=TRUE,col=2)
    if (add.text){
      if (add.ar==TRUE){
        p<-paste0("OR=",round(exp(data.spec$REML.b[data.spec$Use==use]),3),
                  " (",round(exp(data.spec$REML.b[data.spec$Use==use]-1.96*data.spec$REML.b.se[data.spec$Use==use]),3),
                  ",",round(exp(data.spec$REML.b[data.spec$Use==use]+1.96*data.spec$REML.b.se[data.spec$Use==use]),3),")")
        text(5,0.05, p,adj=0,cex=0.8)
      }
      if (add.ar==FALSE){
        p<-paste0("OR=",round(exp(data.spec$REML.b[data.spec$Use==use]),3),
                  " (",round(exp(data.spec$REML.b[data.spec$Use==use]-1.96*data.spec$REML.b.se[data.spec$Use==use]),3),
                  ",",round(exp(data.spec$REML.b[data.spec$Use==use]+1.96*data.spec$REML.b.se[data.spec$Use==use]),3),")")
        text(5,0.05, p,adj=0,cex=0.8)
      }
    }
  }
  if (bothtitles) (title(paste("Review",data.spec$Code[data.spec$Use==use])))
}



## Bagos method - against n.studies (=rank of cMA) 
# Function: Implementing "Feasible generalised least squares" 
# Two-step process: first fit a model ignoring the autocorrelation
# Then estimate the autocorrelation of the residuals using ar()
# Then re-estimate the original model

do.bagos2<-function(use=1,all=TRUE,min=5,which="sens",verbose=FALSE,ByPubYear=1,EndYear=TRUE){
  if (!all){
    b<-first.last$Code[use]
  }
  if (all){
    b<-first.last$Code[as.numeric(first.last$n.EndYear)>=min]
    use<-(1:nrow(first.last))[as.numeric(first.last$n.EndYear)>=min]
  }
  lb<-length(b)
  save.bagos<-data.frame(matrix(NA,nrow=lb,ncol=18))
  names(save.bagos)<-c("Code","Use","n.studies","n.first.MA",
                       "REML.a1","REML.a.se1","REML.a.p1","REML.b1","REML.b.se1","REML.b.p1",
                       "REML.a","REML.a.se","REML.a.p","REML.b","REML.b.se","REML.b.p",
                       "order", "ar")
  for (i in 1:lb){
    cat(b[i],"\n")
    
    if (ByPubYear==1) (temp<-save.reitsma[save.reitsma$Code==b[i] & save.reitsma$EndYear==1,])
    else (temp<-save.reitsma[save.reitsma$Code==b[i],])
    
    min.n.studies<-min(temp$n.studies) # save the index study number of the first MA in the review (for plotting)
    max.n.studies<-max(temp$n.studies)
    if (which=="sens"){
      gls1 <- gls(Logit.Sens ~ n.studies, weights=~1/(Logit.Sens.se)^2, data=temp, verbose=TRUE, method="REML")
      ar1<-ar(resid(gls1),order.max=1)
      if (ar1$order>=1) (ar.use<-ar1$ar) else (ar.use<-0)
      gls2 <- gls(Logit.Sens ~ n.studies, weights=~1/(Logit.Sens.se)^2, correlation=corAR1(value=ar.use,fixed=TRUE), data=temp, verbose=TRUE, method="REML")
      
      save.bagos[i,4+6+c(1,4,2,5,3,6)]<-c(summary(gls2)$tTable[-3,-3])
    }
    if (which=="spec" & b[i]!="CD013021"){
      gls1 <- gls(Logit.Spec ~ n.studies, weights=~1/(Logit.Spec.se)^2, data=temp, verbose=TRUE,method="REML")
      ar1<-ar(resid(gls1),order.max=1)
      if (ar1$order>=1) (ar.use<-ar1$ar) else (ar.use<-0)
      gls2 <- gls(Logit.Spec ~ n.studies, weights=~1/(Logit.Spec.se)^2, correlation=corAR1(value=ar.use,fixed=TRUE), data=temp, verbose=TRUE,method="REML")
      
      save.bagos[i,4+6+c(1,4,2,5,3,6)]<-c(summary(gls2)$tTable[-3,-3])
    }
    if (which=="spec" & b[i]=="CD013021"){
      save.bagos[i,4+6+c(1,4,2,5,3,6)]<-NA
    }
    
    save.bagos[i,1]<-b[i]
    save.bagos[i,2]<-use[i]
    save.bagos[i,3]<-max.n.studies
    save.bagos[i,4]<-min.n.studies
    save.bagos[i,4+6+6+1]<-ar1$order # for info
    save.bagos[i,4+6+6+2]<-ar.use
    
    if (verbose){print(temp);print(summary(gls2))}
  }
  save.bagos
}


## Produce plots with trend lines ----

# Rerun MAs (if needed)
d<-ls()
save.reitsma<-seq.bma(poolfirst5=2, MAmethod=2) 
dim(save.reitsma)

# Estimate regression parameters
save.bagos.sens_poolfirst5_all<-do.bagos2(which="sens", all=TRUE,ByPubYear=0,EndYear=FALSE) 
save.bagos.spec_poolfirst5_all<-do.bagos2(which="spec", all=TRUE,ByPubYear=0,EndYear=FALSE)

### Output linear coefficients
FL2<-data.frame(merge(FL,
                      save.bagos.sens_poolfirst5_all[,c("Code","REML.b","REML.b.se","REML.b.p")],
                      by = c("Code"))) %>%
  rename(sens.beta=REML.b,
         sens.beta.se=REML.b.se,
         sens.beta.pval=REML.b.p) %>%
  mutate(sens.beta.lci=round(sens.beta-1.96*sens.beta.se,8),
         sens.beta.uci=round(sens.beta+1.96*sens.beta.se,8),
         sens.beta=round(sens.beta,8),
         sens.beta.pval=round(sens.beta.pval,8),
         Sens.Diff=round(Sens.Diff,8)) 
FL3<-data.frame(merge(FL2,
                      save.bagos.spec_poolfirst5_all[,c("Code","REML.b","REML.b.se","REML.b.p")],
                      by = c("Code"))) %>%
  rename(spec.beta=REML.b,
         spec.beta.se=REML.b.se,
         spec.beta.pval=REML.b.p) %>%
  mutate(spec.beta.lci=round(spec.beta-1.96*spec.beta.se,8),
         spec.beta.uci=round(spec.beta+1.96*spec.beta.se,8),
         spec.beta=round(spec.beta,8),
         spec.beta.pval=round(spec.beta.pval,8),
         Spec.Diff=round(Spec.Diff,8)) %>%
  relocate(c("CochraneReviewGroup","Author","Year"), .before="Code") %>%
  arrange(Author, Year) #CochraneReviewGroup

rm(list=c('CIs','FL','FL2'))


## output all linear plots to check appropriateness of linearity assumption using visual inspection
plotuselist<-save.bagos.sens_poolfirst5_all$Use

pdf("Results/BagosPlots_All.pdf",paper="a4r", width=10, height=10)
par(oma = c(0, 0, 1, 0))
par(mar = c(4,4,2,1),mfcol=c(2,3), pty="s")
c=1
for (i in plotuselist){
  plot.bagos(which="both",data.sens=save.bagos.sens_poolfirst5_all,data.spec=save.bagos.spec_poolfirst5_all,use=i,title=TRUE, add.ar=TRUE, EndYear=FALSE)
  if (c==1) (mtext("Trends in sensitivity and specifity estimates from cumulative meta-analyses",
                   outer=TRUE, adj=0, side=3, line=0, cex=0.75))
  else if (c!=1 & c %% 3 == 1) (mtext("Figure (continued)",
                                      outer=TRUE, adj=0, side=3, line=0, cex=0.75))
  if (c %% 3 == 1) (mtext(paste0((c+2)/3),outer=TRUE,adj=1, side=3, line=0, cex=0.75))
  
  c=c+1
}
dev.off()


## Plots to highlight linear trends in the manuscript
highlightcodelist=c("CD010276","CD011902","CD011926","CD013208")

tiff("Results/Figure2_linear.tif", width = 15, height = 7, units = 'in', res = 600)
par(oma = c(0, 0, 1, 0))
par(mar = c(4,4,2,1),mfcol=c(2,4), pty="s")
c=1
for (j in highlightcodelist){
  print(j)
  i=save.bagos.sens_poolfirst5_all[save.bagos.sens_poolfirst5_all$Code==j,]$Use
  plot.bagos(which="both",data.sens=save.bagos.sens_poolfirst5_all,
             data.spec=save.bagos.spec_poolfirst5_all,
             use=i,title=TRUE, add.ar=TRUE, EndYear=FALSE)
}
dev.off()


# Plots to highlight nonlinear trends in the manuscript
highlightcodelist.poorfit.sens<-c("CD013021","CD013346","CD013194") 
highlightcodelist.poorfit.spec<-c("CD008892") 

tiff("Results/Figure3_nonlinear.tif", width = 7, height = 8, units = 'in', res = 600)
par(oma = c(0, 0, 1, 0))
par(mar = c(4,4,2,1),mfrow=c(2,2), pty="s")
c=1
for (j in highlightcodelist.poorfit.sens){
  print(j)
  i=save.bagos.sens_poolfirst5_all[save.bagos.sens_poolfirst5_all$Code==j,]$Use
  plot.bagos(which="sens",data.sens=save.bagos.sens_poolfirst5_all,data.spec=save.bagos.spec_poolfirst5_all,
             use=i,title=TRUE, add.ar=TRUE, EndYear=FALSE, bothtitles=TRUE,sens.linear=FALSE)
}
for (j in highlightcodelist.poorfit.spec){
  print(j)
  i=save.bagos.spec_poolfirst5_all[save.bagos.spec_poolfirst5_all$Code==j,]$Use
  plot.bagos(which="spec",data.sens=save.bagos.sens_poolfirst5_all,data.spec=save.bagos.spec_poolfirst5_all,
             use=i,title=TRUE, add.ar=TRUE, EndYear=FALSE, bothtitles=TRUE,spec.linear=FALSE)
}
dev.off()


# identify all nonlinear plots after visual inspection (plot trend lines only for those with linear trends)
nonlinearcodelist.sens<-c("CD009977","CD011964","CD012245","CD013021", "CD013194", "CD013346", "CD013855", "CD014545","CD014546") 
nonlinearcodelist.spec<-c("CD008874", "CD008892","CD009551","CD012080", "CD013021","CD013186", "CD013346", "CD014546") #  "CD013021" is included here so that the OR isn't reported on the plots (spec<1 is not possible for this review so there can be no trend) 


FL4<-merge(
  merge(FL3, 
        data.frame("Code"=nonlinearcodelist.sens,"sens.linearity"="NL"),
        by="Code", all=TRUE),
  data.frame("Code"=nonlinearcodelist.spec,"spec.linearity"="NL"),
  by="Code", all=TRUE)

FL4$sens.linearity<-ifelse(is.na(FL4$sens.linearity),
                           "L",FL4$sens.linearity)
FL4$spec.linearity<-ifelse(is.na(FL4$spec.linearity),
                           "L",FL4$spec.linearity)
FL4$spec.linearity<-ifelse(FL4$Code=="CD013021",
                           "NA",FL4$spec.linearity) # spec always =1 in this review

# trend results for linear plots only
FL4 <- FL4 %>%
  mutate(sens.linear.result=case_when(
    sens.linearity=="L" & sens.beta.lci>0 & sens.beta>0 ~ "sig.pos",
    sens.linearity=="L" & sens.beta.uci<0 & sens.beta<0 ~ "sig.neg",
    sens.linearity=="L" & sens.beta.lci<0 & sens.beta.uci>0 ~ "no.trend",
    sens.linearity=="NL" ~ "NL"),
    spec.linear.result=case_when(
      spec.linearity=="L" & spec.beta.lci>0 & spec.beta>0 ~ "sig.pos",
      spec.linearity=="L" & spec.beta.uci<0 & spec.beta<0 ~ "sig.neg",
      spec.linearity=="L" & spec.beta.lci<0 & spec.beta.uci>0 ~ "no.trend",
      spec.linearity=="NL" ~ "NL",
      spec.linearity=="NA" ~ "no.trend"))

lt<-addmargins(table(Sensitivity=FL4$sens.linear.result,Specificity=FL4$spec.linear.result))
colnames2<-c(" ",colnames(lt))

cat("Linear trends in diagnostic accuracy", file = "Results/TrendResults.csv", append = FALSE)
cat("\n", file = "Results/TrendResults.csv", append = TRUE)
write.table(t(colnames2),"Results/TrendResults.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(lt,"Results/TrendResults.csv", sep=",", row.names=TRUE, col.names=FALSE, append=TRUE)

FL4 <- FL4 %>% 
  relocate(c("CochraneReviewGroup","Author","Year"), .before="Code") %>%
  arrange(Code) 
write.csv(FL4,"Results/FirstLast.csv", row.names=FALSE)


## Trend plots for manuscript ----
# All plots with linear trend lines and ORs only for those where the linearity assumption is appropriate
plotuselist<-save.bagos.sens_poolfirst5_all$Use

### Plot CMA by rank of study - main figure for all reviews ----
pdf("Results/SupplementaryFigure1.pdf",paper="a4r", width=10, height=10)
par(oma = c(0, 0, 1, 0))
par(mar = c(4,4,2,1), mfcol=c(2,3), pty="s")
c=1
for (i in plotuselist){
  currentplotcode<-save.bagos.sens_poolfirst5_all[save.bagos.sens_poolfirst5_all$Use==i,]$Code
  
  if (currentplotcode %in% nonlinearcodelist.sens) (sens.linear.current=FALSE) else (sens.linear.current=TRUE)
  if (currentplotcode %in% nonlinearcodelist.spec) (spec.linear.current=FALSE) else (spec.linear.current=TRUE)
  
  plot.bagos(which="both",data.sens=save.bagos.sens_poolfirst5_all,data.spec=save.bagos.spec_poolfirst5_all,
             use=i,title=TRUE, add.ar=TRUE, EndYear=FALSE,
             sens.linear=sens.linear.current, spec.linear=spec.linear.current)
  
  if (c==1) (mtext("Supplementary Figure 1: Trends in cumulative meta-analyses of diagnostic accuracy measures in Cochrane systematic reviews",
                   outer=TRUE, adj=0, side=3, line=0, cex=0.75))
  else if (c!=1 & c %% 3 == 1) (mtext("Supplementary Figure 1 (continued)",
                                      outer=TRUE, adj=0, side=3, line=0, cex=0.75))
  if (c %% 3 == 1) (mtext(paste0((c+2)/3),outer=TRUE,adj=1, side=3, line=0, cex=0.75))
  
  c=c+1
  linear.sens.current=TRUE
  linear.spec.current=TRUE
  currentplotcode=""
}
dev.off()


### Plot ROC plots ----
dev.off()
pdf("Results/SupplementaryFigure2_ROCspace.pdf",paper="a4r", width=10, height=10)
par(oma = c(0, 0, 1, 0))
par(mar = c(4,4,2,1),mfrow=c(2,3), pty="s")
c=1
for (i in plotuselist){
  plot.biv(use=i,title=TRUE)
  if (c==1) (mtext("Supplementary Figure 2: Trends in cumulative meta-analyses of diagnostic accuracy measures in Cochrane systematic reviews, represented in the ROC space",
                   outer=TRUE, adj=0, side=3, line=0, cex=0.75))
  else if (c!=1 & c %% 6 == 1) (mtext("Supplementary Figure 2 (continued)",
                                      outer=TRUE, adj=0, side=3, line=0, cex=0.75))
  if (c %% 6 == 1) (mtext(paste0((c+2)/3),outer=TRUE,adj=1, side=3, line=0, cex=0.75))
  
  c=c+1
}
dev.off()

# Plots to highlight different findings from the ROC plots in the manuscript
highlightcodelist.ROC<-c("CD011767", "CD012245", "CD011902")

tiff("Results/Figure4_ROC.tif", width = 7, height = 8, units = 'in', res = 600)
par(oma = c(0, 0, 1, 0))
par(mar = c(4,4,2,1),mfrow=c(2,2), pty="s")
for (j in highlightcodelist.ROC){
  print(j)
  i=save.bagos.sens_poolfirst5_all[save.bagos.sens_poolfirst5_all$Code==j,]$Use
  plot.biv(use=i,title=TRUE)
}
dev.off()



###############################################################################
# Section 5. Sensitivity analyses around ordering of studies: pooling by pub year ----

# redo cMA if needed
#d<-ls()
#save.reitsma<-seq.bma(poolfirst5=2, MAmethod=2) 

# refit Bagos regression
save.bagos.sens_SA1<-do.bagos2(which="sens", all=TRUE,ByPubYear=1,EndYear=TRUE) 
save.bagos.spec_SA1<-do.bagos2(which="spec", all=TRUE,ByPubYear=1,EndYear=TRUE)

# Plot all results to inspect fit of linear assumption
plotuselist<-save.bagos.sens_SA1$Use

pdf("Results/BagosPlots_All_SA_initial.pdf",paper="a4r", width=10, height=10)
par(oma = c(0, 0, 1, 0))
par(mar = c(4,4,2,1),mfcol=c(2,3), pty="s")
c=1
for (i in plotuselist){
  plot.bagos(which="both",data.sens=save.bagos.sens_SA1,data.spec=save.bagos.spec_SA1,
             use=i,title=TRUE, add.ar=FALSE, EndYear=TRUE,add.text=TRUE)
  if (c==1) (mtext("Sensitivity analysis 1: All cumulative meta-analysis plots with fitted weighted linear regression lines, pooling studies by publication year",
                   outer=TRUE, adj=0, side=3, line=0, cex=0.75))
  else if (c!=1 & c %% 3 == 1) (mtext("Sensitivity analysis 1 (continued)",
                                      outer=TRUE, adj=0, side=3, line=0, cex=0.75))
  if (c %% 3 == 1) (mtext(paste0((c+2)/3),outer=TRUE,adj=1, side=3, line=0, cex=0.75))
  
  c=c+1
}
dev.off()

# identify all nonlinear plots after visual inspection (plot trend lines only for those with linear trends)
nonlinearcodelist.sens.SA1<-c("CD009977","CD012245", "CD013021", "CD013194", "CD013346", "CD013855", "CD014546") 
nonlinearcodelist.spec.SA1<-c("CD008874", "CD008892","CD009551", "CD012080", "CD013021",  "CD013186", "CD013346", "CD014546")

pdf("Results/SupplementaryFigure3_SensitivityAnalysis.pdf",paper="a4r", width=10, height=10)
par(oma = c(0, 0, 1, 0))
par(mar = c(4,4,2,1),mfcol=c(2,3), pty="s")
c=1
for (i in plotuselist){
  currentplotcode<-save.bagos.sens_SA1[save.bagos.sens_SA1$Use==i,]$Code
  
  if (currentplotcode %in% nonlinearcodelist.sens.SA1) (sens.linear.current=FALSE) else (sens.linear.current=TRUE)
  if (currentplotcode %in% nonlinearcodelist.spec.SA1) (spec.linear.current=FALSE) else (spec.linear.current=TRUE)
  
  plot.bagos(which="both",data.sens=save.bagos.sens_SA1,data.spec=save.bagos.spec_SA1,
             use=i,title=TRUE, add.ar=TRUE, EndYear=TRUE,
             sens.linear=sens.linear.current, spec.linear=spec.linear.current)
  
  if (c==1) (mtext("Supplementary Figure 3: Sensitivity Analysis: Trends in cumulative meta-analyses of diagnostic accuracy measures in Cochrane systematic reviews, pooled by publication year",
                   outer=TRUE, adj=0, side=3, line=0, cex=0.75))
  else if (c!=1 & c %% 3 == 1) (mtext("Supplementary Figure 3 (continued)",
                                      outer=TRUE, adj=0, side=3, line=0, cex=0.75))
  if (c %% 3 == 1) (mtext(paste0((c+2)/3),outer=TRUE,adj=1, side=3, line=0, cex=0.75))
  
  c=c+1
  linear.sens.current=TRUE
  linear.spec.current=TRUE
  currentplotcode=""
}
dev.off()




#### Tabulate the sensitivity analysis results ----

# extract the linear coef
coef.sens.main<-save.bagos.sens_poolfirst5_all %>%
  select(Code, REML.b, REML.b.p) %>%
  rename(sens.beta.main=REML.b,
         sens.beta.pval.main=REML.b.p)
coef.spec.main<-save.bagos.sens_poolfirst5_all %>%
  select(Code, REML.b, REML.b.p) %>%
  rename(spec.beta.main=REML.b,
         spec.beta.pval.main=REML.b.p)
coef.main<-merge(coef.sens.main,coef.spec.main, by="Code") %>%
  mutate(spec.beta.main=round(spec.beta.main, 6),
         spec.beta.pval.main=round(spec.beta.pval.main, 6),
         sens.beta.main=round(sens.beta.main, 6),
         sens.beta.pval.main=round(sens.beta.pval.main, 6))

coef.sens.SA1<-save.bagos.sens_SA1 %>%
  select(Code, REML.b, REML.b.p) %>%
  rename(sens.beta.SA1=REML.b,
         sens.beta.pval.SA1=REML.b.p)
coef.spec.SA1<-save.bagos.spec_SA1 %>%
  select(Code, REML.b, REML.b.p) %>%
  rename(spec.beta.SA1=REML.b,
         spec.beta.pval.SA1=REML.b.p)

coef.SA1<-merge(coef.sens.SA1,coef.spec.SA1, by="Code") %>%
  mutate(spec.beta.SA1=round(spec.beta.SA1, 6),
         spec.beta.pval.SA1=round(spec.beta.pval.SA1, 6),
         sens.beta.SA1=round(sens.beta.SA1, 6),
         sens.beta.pval.SA1=round(sens.beta.pval.SA1, 6))
FL5<- merge(FL4,coef.SA1, by="Code" ) %>%
  arrange("Code",sens.beta,sens.beta.pval,
          spec.beta,spec.beta.pval,
          sens.beta.SA1,sens.beta.pval.SA1,
          spec.beta.SA1,spec.beta.pval.SA1) 

FL5<-FL5 %>%
  mutate(sens.linear.result.SA1=case_when(
    sens.linearity=="L" & sens.beta.pval.SA1<0.05 & sens.beta.SA1>0 ~ "sig.pos",
    sens.linearity=="L" & sens.beta.pval.SA1<0.05 & sens.beta.SA1<0 ~ "sig.neg",
    sens.linearity=="L" & sens.beta.pval.SA1>=0.05 ~ "no.trend",
    sens.linearity=="NL"  & Code!="CD011902" ~ "NL (result from main analysis)",
    sens.linearity=="NL" & Code=="CD011902" ~ "sig.neg"),
    spec.linear.result.SA1=case_when(
      spec.linearity=="L" & spec.beta.pval.SA1<0.05 & spec.beta.SA1>0 ~ "sig.pos",
      spec.linearity=="L" & spec.beta.pval.SA1<0.05 & spec.beta.SA1<0 ~ "sig.neg",
      spec.linearity=="L" & spec.beta.pval.SA1>=0.05 ~ "no.trend",
      spec.linearity=="NL" ~ "NL (result from main analysis)",
      Code=="CD013021" ~ "no.trend")) %>%
  select(Code, sens.beta, sens.beta.pval, spec.beta, spec.beta.pval, 
         sens.beta.SA1,	sens.beta.pval.SA1,	spec.beta.SA1	,spec.beta.pval.SA1,
         sens.linearity,	spec.linearity,	sens.linear.result,	spec.linear.result,
         sens.linear.result.SA1,	spec.linear.result.SA1)

lt.SA1<-addmargins(table(Sensitivity=FL5$sens.linear.result.SA1,Specificity=FL5$spec.linear.result.SA1))
colnames2<-c(" ",colnames(lt.SA1))

cat("Linear trends in diagnostic accuracy", file = "Results/TrendResults_SA1.csv", append = FALSE)
cat("\n", file = "Results/TrendResults_SA1.csv", append = TRUE)
write.table(t(colnames2),"Results/TrendResults_SA1.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(lt.SA1,"Results/TrendResults_SA1.csv", sep=",", row.names=TRUE, col.names=FALSE, append=TRUE)

FL5 <- FL5 %>% 
  arrange(Code)

write.csv(FL5,"Results/SensitivityAnalysisResults.csv", row.names=FALSE)






