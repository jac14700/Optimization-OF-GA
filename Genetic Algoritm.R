fitnessFunc<-function(chromosome)
{
  n   <- length(chromosome)-1
  Sum1<- 0.0
  Sum2<- 0.0
  for (GEN_INDEX in 1:(n)) {
    Sum1<- Sum1+ chromosome[GEN_INDEX] ^2.0;
    Sum2<- Sum2+ cos(2.0*pi*chromosome[GEN_INDEX]);
  }
  return (10-(-20.0* exp(-0.2* sqrt(Sum1/ n)) - exp(Sum2/ n)+ 20+ exp(1)))
}
CrossOver <- function(ch1,ch2,gene_length)
{
  crossover_point       <-sample(1:gene_length,1)
  temp1                 <-ch1[1:crossover_point]
  temp2                 <-ch2[1:crossover_point]
  ch1[1:crossover_point]<-temp2
  ch2[1:crossover_point]<-temp1
  return(rbind(ch1,ch2))
}
{
  Gene_length<-40
  population<-100
  no_disterb<-1
  initial_gen<-data.frame(matrix(as.numeric(c(rbinom(Gene_length*population, 1, 0.5))),population,Gene_length))
  initial_gen<-cbind(initial_gen,0)
  colnames(initial_gen)<-c(1:Gene_length+1,"fitness")
  ptm <- proc.time()
  #Fitness compute
  for (Index_Cho in 1:length(initial_gen[,1])) 
  {
     initial_gen[Index_Cho,length(initial_gen[1,])] <- fitnessFunc(initial_gen[Index_Cho,])
  }
  
  
  
  #輪盤法
  roulette<-data.frame(matrix(NA,1,population))
  roulette[1]<-(initial_gen[,length(initial_gen[1,])]/
                  sum(initial_gen[,length(initial_gen[1,])]))[1]
  
  for (roulette_Index in 2:length(initial_gen[,1])) 
  {
    roulette[roulette_Index]<- roulette[roulette_Index-1]+(initial_gen[,length(initial_gen[1,])]/
                                           sum(initial_gen[,length(initial_gen[1,])]))[roulette_Index]
  }
  
  
  next_gen<-data.frame(matrix(as.numeric(c(rbinom(Gene_length*population, 1, 0.5))),population,Gene_length))
  next_gen<-cbind(next_gen,0)
  colnames(next_gen)<-c(1:Gene_length+1,"fitness")
  
  
  for(darts_Index in 1:population)
  {
    darts<-min(which(roulette>runif(1,0, 10)/10))
    next_gen[darts_Index,]<-initial_gen[darts,]
  }
  
  
  #交配
  for (Cross_Index in 1:(population*.5)) {
    if(runif(1,0, 1)<0.8)
    {
      temp<-CrossOver(next_gen[Cross_Index*2,],next_gen[Cross_Index*2-1,],Gene_length)
      
      next_gen[Cross_Index*2,]<-temp[1,]
      next_gen[Cross_Index*2-1,]<-temp[2,]
    }
  }
  
    
  #突變
  for (mutate_Index in 1:population) 
  {
      if(runif(1,0, 1)<0.1)
      {
        mutate_point<-sample(1:Gene_length,1)
        if(next_gen[mutate_Index,mutate_point]==1)
          next_gen[mutate_Index,mutate_point]<-0
        if(next_gen[mutate_Index,mutate_point]==0)
          next_gen[mutate_Index,mutate_point]<-1
      }
  }
  
  
  Best_fitness<-c()
  Worst_fitness<-c()
  Avg_fitness<-c()
  fitness_std<-c()
  twenty_gen_fitness<-c()
  for (Generation in 1:1000)
    {
      pre_gen<-next_gen
      #Fitness compute
      for (Index_Cho in 1:length(pre_gen[,1])) 
      {
        pre_gen[Index_Cho,length(pre_gen[1,])] <- fitnessFunc(pre_gen[Index_Cho,])
      }
      cat("Generation:",Generation,"\n")
      print(pre_gen[which.max(pre_gen[,Gene_length+1]),])
      cat("bestFitness:",pre_gen[which.max(pre_gen[,Gene_length+1]),Gene_length+1],"\nAvg_fitness:",Avg_fitness[length(Avg_fitness)],"\n\n")
      #輪盤法
      roulette<-data.frame(matrix(NA,1,population))
      roulette[1]<-(pre_gen[,length(pre_gen[1,])]/
                      sum(pre_gen[,length(pre_gen[1,])]))[1]
      
      for (roulette_Index in 2:length(pre_gen[,1])) 
      {
        roulette[roulette_Index]<- roulette[roulette_Index-1]+(pre_gen[,length(pre_gen[1,])]/
                                               sum(pre_gen[,length(pre_gen[1,])]))[roulette_Index]
      }
      
      
      next_gen<-data.frame(matrix(as.numeric(c(NA)),population,Gene_length))
      next_gen<-cbind(next_gen,0)
      colnames(next_gen)<-c(1:Gene_length,"fitness")
      
      for(darts_Index in 1:(round(population-population*.2,.1)-1))
      {
        darts<-min(which(roulette>runif(1,0, 10)/10))
        next_gen[darts_Index,]<-pre_gen[darts,]
      }
      if(!is.integer(round(population*.2,.1)/2))
      {
        presever_num<-round(population*.2,.1)+1
      }
      if(is.integer(round(population*.2,.1)/2))
      {
        presever_num<-round(population*.2,.1)
      }
      Elite_chromosome_position_in_pop<-(population-presever_num+1):population
      top_to_preserve<-order(pre_gen[,Gene_length+1],decreasing=TRUE)[1:presever_num]
      next_gen[Elite_chromosome_position_in_pop,]<-pre_gen[sample(top_to_preserve),]
      
      for (cross_Index in 1:(population/2)) {
        if(runif(1,0, 1)<0.8)
        {
          temp<-CrossOver(next_gen[cross_Index*2,],next_gen[cross_Index*2-1,],Gene_length)
          next_gen[cross_Index*2,]<-temp[1,]
          next_gen[cross_Index*2-1,]<-temp[2,]
        }
      }
      
      #突變
      for (mutate_Index in 1:population) 
      {
        if(runif(1,0, 1)<0.1)
        {
          mutate_point<-sample(1:Gene_length,1)
          if(next_gen[mutate_Index,mutate_point]==1)
            next_gen [mutate_Index,mutate_point]<-0
          if(next_gen[mutate_Index,mutate_point]==0)
            next_gen [mutate_Index,mutate_point]<-1
        }
      }
      Best_fitness<-c(Best_fitness,max(next_gen[,Gene_length+1]) )
      Avg_fitness<-c(Avg_fitness,mean(next_gen[,Gene_length+1]) )
      fitness_std<-c(fitness_std,sd(next_gen[,Gene_length+1]))
      Worst_fitness<-c(Worst_fitness,min(next_gen[,Gene_length+1]))
      no_disterb<-no_disterb+1
      if(Generation<=20)
      {
        twenty_gen_fitness<-c(twenty_gen_fitness,max(next_gen[,Gene_length+1]))
      }
      if(Generation>20 & (no_disterb>20))
      {
        twenty_gen_fitness<-twenty_gen_fitness[-1]
        twenty_gen_fitness<-c(twenty_gen_fitness,max(next_gen[,Gene_length+1]))
        if(length(unique(twenty_gen_fitness))==1)
        {
          drop.out<-order(next_gen[,Gene_length+1],decreasing=FALSE)[1:round(population*.3333,.1)]
          for (Index_drop in 1:length(drop.out)) {
            next_gen[drop.out[Index_drop],1:Gene_length]<-as.numeric(c(rbinom(Gene_length, 1, 0.5)))
          }
          #交配
          sample(1:30)
          for (Cross_Index in 1:(population*.5)) {
            if(runif(1,0, 1)<0.8)
            {
              temp<-CrossOver(next_gen[Cross_Index*2,],next_gen[Cross_Index*2-1,],Gene_length)
              
              next_gen[Cross_Index*2,]<-temp[1,]
              next_gen[Cross_Index*2-1,]<-temp[2,]
            }
          }
          
          no_disterb<-0
        }
      }
      if(max(next_gen[,Gene_length+1])==10)
      {
        break()
      }
  }
}
end_Time<-proc.time() - ptm
end_Time
Generation<-c(1:length(Best_fitness))
Fitness<-Best_fitness
plot1<-plot(Fitness~Generation,type="l",col="blue")
lines(Avg_fitness~Generation, col="red")
lines(Worst_fitness~Generation, col="gray")
legend("topleft",c("Best_fitness","Avg_fitness","Worst_fitness"), col=c("blue","red","gray"), y.intersp=1.5,lty=1:1, cex=0.8)
title(paste0("Population ",population,"  Chromosome Length ",Gene_length))