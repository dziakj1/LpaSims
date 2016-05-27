rm(list=ls(all=TRUE));
#########################################################
# Set initial settings:
set.seed(700146);
N <- 1000;
nsim <- 1000;
# Loop over conditions;
for (measurement.quality in c("low","high")) {
  for (class.size.distribution in c("even","uneven")) {
    for (effect.size in c("small","medium","large")) {
      print(paste(measurement.quality,class.size.distribution,effect.size));
      output.path <- "C:\\Users\\jjd264\\Documents\\Sims-Lpa-Binary\\Datasets\\";
      #########################################################
      # Define parameters;
      true.item.means <- matrix(c(  1, 0, -1,
                                    1, 0, -1,
                                    1, 0, -1,
                                    1, 0, -1,
                                    1, 0, -1), ncol=3,byrow=TRUE);
      if (measurement.quality=="high") {
        true.item.sds <- matrix(sqrt(.5),nrow=5,ncol=3);
      }
      if (measurement.quality=="low") { 
        true.item.sds <- matrix(sqrt(.75),nrow=5,ncol=3);
      }
      if (class.size.distribution=="even") {
        true.class.sizes <- c(1,1,1)/3;
      }
      if (class.size.distribution=="uneven") {
        true.class.sizes <- c(.1,.3,.6);
      }
      if (effect.size == "large") {
        true.distal.prob.by.class <- c(.8,.5,.2);
      }
      if (effect.size == "medium") {
        true.distal.prob.by.class <- c(.7,.5,.3);
      }
      if (effect.size == "small") {
        true.distal.prob.by.class <- c(.55,.5,.45);
      }
      true.distal.prob <- true.distal.prob.by.class%*%true.class.sizes;
      setting.file.name <- paste(paste("settings",
                                       measurement.quality,
                                       class.size.distribution,
                                       effect.size,
                                       N, 
                                       sep="-"),".rdata",sep="");
      save(true.item.means,
           true.item.sds,
           true.class.sizes,
           true.distal.prob.by.class,
           true.distal.prob,
           file=paste(output.path,
                      setting.file.name,
                      sep=""));
      for (simulation in 1:nsim) {
        #########################################################    
        # Simulate datasets;
        data.file.name <- paste(paste("sim-data",
                                      measurement.quality,
                                      class.size.distribution,
                                      effect.size,
                                      N,
                                      simulation,
                                      sep="-"),".txt",sep="");
        trueclass <- sample(x=1:3,
                             size=N,
                             replace=TRUE,
                             prob=true.class.sizes);
        id <- 1:N;
        items <- matrix(NA,nrow=N,ncol=5);
        distal <- rep(NA,N);
        for (this.item in 1:5) {
          for (this.class in 1:3) {
            items[which(trueclass==this.class),this.item] <-
              round(rnorm(n=which(trueclass==this.class),
                          mean=true.item.means[this.item,this.class],
                          sd=true.item.sds[this.item,this.class]),6);
            distal[which(trueclass==this.class)] <-
              rbinom(n=which(trueclass==this.class),
                     size=1,
                     prob=true.distal.prob.by.class[this.class]);
          }
        } 
        #######################################################
        # Output the simulated dataset:
        simulated.dataset <- data.frame(cbind(id=id,
                                              y=distal,
                                              I1=items[,1],
                                              I2=items[,2],
                                              I3=items[,3],
                                              I4=items[,4],
                                              I5=items[,5],
                                              trueclass=trueclass));
        write.table(x=simulated.dataset,
                    quote=FALSE,
                    sep="\t",
                    col.names=TRUE,
                    row.names=FALSE, 
                    file=paste(output.path,data.file.name,sep=""));
      }
    }
  }
}

