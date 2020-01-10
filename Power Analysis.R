# Power analysis 

# Parameters to set: 
# Attach - Correlation between attachment anxiety and avoidance 
# pAlone - Probability of being alone or with others
# AnxCoef - Coefficient for attachment anxiety 
# AvoidCoef - Coefficient for attachment avoidance 
# AloneCoef - Coefficient for solitude 
# n - Sample size
# k - Number of time-points
# e1sd - Standard deviation of level 1 error term
# e2sd - Standard deviation of level 2 error term
# y11 - Coefficient for interaction between anxiety and solitude
# y12 - Coefficient for interaction between avoidance and solitude
# y13 - Coefficient for interaction between anxiety and avoidance
# trials - Number of replications

# Load lmerTest
library(lmerTest)

# Function to generate one study
Study <- function(Attach, pAlone, AnxCoef, AvoidCoef, AloneCoef, 
                  n, k, e1sd, e2sd, y11, y12, y13){
  
  # Creating dataframe of simulated values
  # Attachment 
  # Note: Anxiety mean and SD values are based on existing attachment data 
  Anxiety = rnorm(n, mean=2, sd=1) 
  Error = rnorm(n, mean=0, sd=sqrt(1-Attach**2))
  Avoidance = mean(Anxiety) + Attach*Anxiety + Error
  AttachmentData <- data.frame(Anxiety, Avoidance)
  AttachmentData$ID <- c(1:n)
  
  # Sample Intercept
  # Note: Intercept mean is based on pilot data; Rating scale is 1-7
  Intercept <- rnorm(1, mean=3, sd=1)
  ID <-  c(1:n)
  InterceptData <- data.frame(ID, Intercept)
  
  # Alone or with others
  Alone <- rbinom(n*k, size=1, prob = pAlone)
  AloneData <- data.frame(Alone)
  
  # Level 1 error term 
  e1 <- rnorm(n*k, 0, e1sd) # Error terms have a mean of 0 and a SD that can be varied
  ID <-  c(rep((1:n),each=k))
  e1Data <- data.frame(ID, e1)
  
  # Level 2 error term
  e2 <- rnorm(n, 0, e2sd) 
  ID <-  c(1:n)
  e2Data <- data.frame(ID, e2)
  
  # Merge dataframes
  StudyData1 <- merge(AttachmentData, InterceptData, by="ID")
  ErrorData <- merge(e1Data, e2Data, by="ID")
  StudyData2 <- merge(StudyData1, ErrorData, by="ID")
  StudyData <- cbind.data.frame(StudyData2, AloneData)
  
  # Using "y" to symbol as "gamma"
  y00 = StudyData$Intercept[1]
  y01 = AnxCoef
  y02 = AvoidCoef
  y10 = AloneCoef
  y11 = y11
  y12 = y12
  y13 = y13
  
  # Mixed model for Valence
  StudyData$Valence = 
    y00 + 
    y01*StudyData$Anxiety + y02*StudyData$Avoidance +
    y10*StudyData$Alone   + 
    y11*StudyData$Anxiety*StudyData$Alone   +
    y12*StudyData$Avoidance*StudyData$Alone +
    y13*StudyData$Anxiety*StudyData$Avoidance +
    e1 + e2
  
  # Run MLM on data & extract p-values
  results = lmer(Valence ~ Alone*Anxiety + Alone*Avoidance + Anxiety*Avoidance + (1|ID), 
                 data = StudyData)
  coefs <- data.frame(coef(summary(results)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  
  # Save p-values
  Alone.p <- coefs$p.z[2]
  Anx.p <- coefs$p.z[3]
  Avoid.p <- coefs$p.z[4]
  AloneAnx.p <- coefs$p.z[5]
  AloneAvoid.p <- coefs$p.z[6]
  AnxAvoid.p <- coefs$p.z[7]
  
  # Print p-values
  cat(".")
  return(c(Alone.p, Anx.p, Avoid.p, AloneAnx.p, AloneAvoid.p, AnxAvoid.p))
  
}

# Check
Study(Attach=.3, pAlone=.5, AnxCoef=.2, AvoidCoef=.2, AloneCoef=.2, 
      n=100, k=10, e1sd=1, e2sd=.6, y11=.1, y12=.1, y13 = .1) 

# Function to replicate the study
Simulation <- function(Attach, pAlone, AnxCoef, AvoidCoef, AloneCoef, 
                       n, k, e1sd, e2sd, y11, y12, y13, trials){
  
  # Save output of all studies
  output <- replicate(trials,Study(Attach, pAlone, AnxCoef, AvoidCoef, 
                                   AloneCoef, n, k, e1sd, e2sd, y11, y12, y13)) 
  
  # Print output
  output <- t(output)
  
  Alone.p <- output[,1]
  Anx.p <- output[,2]
  Avoid.p <- output[,3]
  AloneAnx.p <- output[,4]
  AloneAvoid.p <- output[,5]
  AnxAvoid.p <- output[,6]
  
  # Optional graphing 
  graphing = 1
  if(graphing==1){
    # Graph p-value distribution of the main effects 
    Alone.histo <- hist(Alone.p, breaks=c(seq(0,1,.01)), plot=FALSE)
    Alone.histo$density <- Alone.histo$counts/(n/100) #converting to percentages rather than counts
    plot(Alone.histo,freq=FALSE,main="Solitude",xlab='P Values',ylab='Relative Frequency',col=c(rep('blue',5),rep('white',95)),abline(h=1,col='red'))
    
    Anx.histo <- hist(Anx.p, breaks=c(seq(0,1,.01)), plot=FALSE)
    Anx.histo$density <- Anx.histo$counts/(n/100) #converting to percentages rather than counts
    plot(Anx.histo,freq=FALSE,main="Attachment Anxiety",xlab='P Values',ylab='Relative Frequency',col=c(rep('blue',5),rep('white',95)),abline(h=1,col='red'))
    
    Avoid.histo <- hist(Avoid.p, breaks=c(seq(0,1,.01)), plot=FALSE)
    Avoid.histo$density <- Avoid.histo$counts/(n/100) #converting to percentages rather than counts
    plot(Avoid.histo,freq=FALSE,main="Attachment Avoidance",xlab='P Values',ylab='Relative Frequency',col=c(rep('blue',5),rep('white',95)),abline(h=1,col='red'))
    
    # And the interaction effects 
    AloneAnx.histo <- hist(AloneAnx.p, breaks=c(seq(0,1,.01)), plot=FALSE)
    AloneAnx.histo$density <- AloneAnx.histo$counts/(n/100) #converting to percentages rather than counts
    plot(AloneAnx.histo,freq=FALSE,main="Solitude x Anxiety",xlab='P Values',ylab='Relative Frequency',col=c(rep('blue',5),rep('white',95)),abline(h=1,col='red'))
    
    AloneAvoid.histo <- hist(AloneAvoid.p, breaks=c(seq(0,1,.01)), plot=FALSE)
    AloneAvoid.histo$density <- AloneAvoid.histo$counts/(n/100) #converting to percentages rather than counts
    plot(AloneAvoid.histo,freq=FALSE,main="Solitude x Avoidance",xlab='P Values',ylab='Relative Frequency',col=c(rep('blue',5),rep('white',95)),abline(h=1,col='red'))
    
    AnxAvoid.histo <- hist(AnxAvoid.p, breaks=c(seq(0,1,.01)), plot=FALSE)
    AnxAvoid.histo$density <- AnxAvoid.histo$counts/(n/100) #converting to percentages rather than counts
    plot(AnxAvoid.histo,freq=FALSE,main="Anxiety x Avoidance",xlab='P Values',ylab='Relative Frequency',col=c(rep('blue',5),rep('white',95)),abline(h=1,col='red'))
  }
  
  # Calculate the percentage of p-values below .05 for each distribution
  Alone.sig <- (sum(Alone.p<=.05)/trials)*100
  Anx.sig <- (sum(Anx.p<=.05)/trials)*100
  Avoid.sig <- (sum(Avoid.p<=.05)/trials)*100
  AloneAnx.sig <- (sum(AloneAnx.p<=.05)/trials)*100
  AloneAvoid.sig <- (sum(AloneAvoid.p<=.05)/trials)*100
  AnxAvoid.sig <- (sum(AnxAvoid.p<=.05)/trials)*100
  
  cat("Percent of significant results  ", "\n")
  
  cat(" Alone: ", Alone.sig, "\n")
  cat(" Anxiety: ", Anx.sig, "\n")
  cat(" Avoid: ", Avoid.sig, "\n")
  cat(" Alone*Anxiety: ", AloneAnx.sig, "\n")
  cat(" Alone*Avoid: ", AloneAvoid.sig, "\n")
  cat(" Anx*Avoid: ", AnxAvoid.sig, "\n")
  cat("\n")
  
}

# Vary parameters. Notes: 
# 1 - 'Attach' (correlation between attachment anxiety and avoidance) set to .3 
# based on existing data on attachment 
# 2 - 'pAlone' (probability of being alone or with others) is set to .5 because of no 
# prior data or reasoning to set it at another value
# 3 - 'AnxCoef', 'AvoidCoef', and 'AloneCoef' (main effects) set to .2 as a reasonable 
# estimate of effect sizes in psychology
# 4 - 'y11' and 'y12' and 'y13' (interaction effects) set to .1 as interaction effects 
# may be smaller than main effects 
# 5 - 'e1sd' and 'e2sd' (standard deviation of error terms) set to 1 and .6 respectively
# based on existing data where similar variables were analysized with MLM 
# (continuous outcome, binary level 1 predictor, continuous level 2 predictors) 

# Power with 200 participants and 32 repeated measurements
Simulation(Attach=.3, pAlone=.5,
           AnxCoef=.2, AvoidCoef=.2,
           AloneCoef=.2,
           n=200, k=32,
           e1sd=1, e2sd=.6, 
           y11=.1, y12=.1, y13=.1, trials=1000)

# Power with 200 participants and 25 repeated measurements
Simulation(Attach=.3, pAlone=.5,
           AnxCoef=.2, AvoidCoef=.2,
           AloneCoef=.2,
           n=200, k=25,
           e1sd=1, e2sd=.6, 
           y11=.1, y12=.1, y13=.1, trials=1000)

