rm(list=ls(all=TRUE));
#########################################################
# Set initial settings:
set.seed(200046);
memory.limit(200000);
input.path <- "C:\\Users\\jjd264\\Documents\\Sims-Lpa-Exponential\\Datasets\\";
working.path <- "C:\\Users\\jjd264\\Documents\\Sims-Lpa-Exponential\\";
LatentGOLD.path <- "C:\\Users\\jjd264\\DOCUME~1\\LATENT~1.1\\lg51.exe";
N <- 1000;
start.time <- Sys.time();
nsim <- 1000;
final.answers <- NULL;
measurement.quality.level.names <- c("low","high");
class.size.distribution.level.names <- c("even","uneven");
effect.size.level.names <- c("small","medium","large");
for (which.measurement.quality in 1:2) {
  for (which.class.size.distribution in 1:2) {
    for (which.effect.size in 1:3) {
      preliminary.answers <- NULL;
      measurement.quality <- measurement.quality.level.names[which.measurement.quality];
      class.size.distribution <- class.size.distribution.level.names[which.class.size.distribution];
      effect.size <- effect.size.level.names[which.effect.size];        
      setting.file.name <- paste(paste("settings",
                                       measurement.quality,
                                       class.size.distribution,
                                       effect.size,
                                       N, 
                                       sep="-"),".rdata",sep="");  
      print(paste(measurement.quality,class.size.distribution,effect.size));
      for (simulation in 1:nsim) {  
        print(simulation);
        # Summarize simulated data in order to check it;
        load(file=paste(input.path,setting.file.name,sep=""));
        data.file.name <- paste(paste("sim-data",
                                      measurement.quality,
                                      class.size.distribution,
                                      effect.size,
                                      N,
                                      simulation,
                                      sep="-"),".txt",sep="");  
        sim.data <- read.table(file=paste(input.path,data.file.name,sep=""),
                               header=TRUE,
                               sep="\t"); 
        oracle.item.means <- matrix(NA,5,3);
        oracle.distal.mean.by.class <- rep(NA,3);
        oracle.distal.mean.by.class.se <- rep(NA,3);
        items <- cbind(sim.data$I1,sim.data$I2,sim.data$I3,sim.data$I4,sim.data$I5);
        for (this.class in 1:3) {
          oracle.item.means[,this.class] <- apply(items[which(sim.data$trueclass==this.class),],2,mean);
          oracle.distal.mean.by.class[this.class] <- mean(sim.data$y[which(sim.data$trueclass==this.class)]);
          oracle.distal.mean.by.class.se[this.class] <- sd(sim.data$y[which(sim.data$trueclass==this.class)]) / 
            sqrt(length(which(sim.data$trueclass==this.class)));
        }
        oracle.class.sizes <- table(sim.data$trueclass)/N;
        oracle.distal.prob <- mean(sim.data$y);    
        # Get the working dataset ready for Latent GOLD;
        working.file.name <- "working.dat";
        if (file.exists(paste(working.path,working.file.name,sep=""))) {
          file.remove(paste(working.path,working.file.name,sep=""));
        }
        stopifnot(file.copy(from=paste(input.path,data.file.name,sep=""),
                            to=paste(working.path,working.file.name,sep="")));
        # Fit a non-inclusive LCA using Latent GOLD;
        if (file.exists(paste(working.path,"NoninclusivePostProbs.dat",sep=""))) {
          file.remove(paste(working.path,"NoninclusivePostProbs.dat",sep=""));
        }
        if (file.exists(paste(working.path,"FitNoninclusive.lst",sep=""))) {
          file.remove(paste(working.path,"FitNoninclusive.lst",sep=""));
        }
        shell(cmd=paste(LatentGOLD.path," ",working.path,"FitNoninclusive.lgs /b",sep=""));
        output <- readLines(con=paste(working.path,"FitNoninclusive.lst",sep=""));
        if (output[11]=="Estimation Warnings! See Iteration Detail\t") {
          these.preliminary.answers <- c(simulation = simulation,
                                         rep(NA,118));
        } else {
          # Get information on the entropy of the fitted model;
          stopifnot(strsplit(output[32],split="\t")[[1]][1]=="Entropy R-squared");
          noninclusive.entropy.Rsqd <- as.numeric(strsplit(output[32],split="\t")[[1]][2]);
          stopifnot(strsplit(output[35],split="\t")[[1]][1]=="Entropy");
          noninclusive.entropy.raw <- as.numeric(strsplit(output[35],split="\t")[[1]][2]);
          noninclusive.entropy.Ramaswamy <- 1 - (noninclusive.entropy.raw/(N*log(3)));
          # Read the estimated item means for the noninclusive LCA;
          stopifnot(output[190]=="Profile");
          temp <- strsplit(output[c(195,197,199,201,203)],2,split="\t");
          unsorted.item.means <- rbind(as.numeric(temp[[1]][c(2,4,6)]),
                                       as.numeric(temp[[2]][c(2,4,6)]),
                                       as.numeric(temp[[3]][c(2,4,6)]),
                                       as.numeric(temp[[4]][c(2,4,6)]),
                                       as.numeric(temp[[5]][c(2,4,6)]));
          # Determine the random class order for the noninclusive LCA;
          possible.permutations <- list(c(1,2,3),
                                        c(1,3,2),
                                        c(2,1,3),
                                        c(2,3,1),
                                        c(3,1,2),
                                        c(3,2,1));
          sum.sqd.error.by.permutation <- rep(NA,6);
          for (this.permutation in 1:6) {
            sum.sqd.error.by.permutation[this.permutation] <-
              sum((unsorted.item.means[,possible.permutations[[this.permutation]]]-true.item.means)^2);
          }
          best.permutation <- possible.permutations[[which.min(sum.sqd.error.by.permutation)]];
          sorted.item.means.noninclusive <- unsorted.item.means[,best.permutation];
          c1 <- best.permutation[1];
          c2 <- best.permutation[2];
          c3 <- best.permutation[3];
          # Latent GOLD analysis:  Modal None (Traditional three-step classify-analyze);
          shell(cmd=paste(LatentGOLD.path," ",working.path,"ModalNone.lgs /b",sep=""));
          output <- readLines(con=paste(working.path,"ModalNone.lst",sep=""));
          if(output[155]=="Cluster\tMean\ts.e.\t \t \t \t \t") {
            distal.mean.by.class.modal.none <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][2],
                                                            strsplit(output[156:158],split="\t")[[c2]][2],
                                                            strsplit(output[156:158],split="\t")[[c3]][2]));
            distal.mean.by.class.modal.none.se <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][3],
                                                               strsplit(output[156:158],split="\t")[[c2]][3],
                                                               strsplit(output[156:158],split="\t")[[c3]][3]));
          } else {
            distal.mean.by.class.modal.none <- rep(NA,3);
            distal.mean.by.class.modal.none.se <- rep(NA,3);
            print(paste("Modal none failed on ",simulation));
          }
          # Latent GOLD analysis:  Proportional None
          shell(cmd=paste(LatentGOLD.path," ",working.path,"ProportionalNone.lgs /b",sep=""));
          output <- readLines(con=paste(working.path,"ProportionalNone.lst",sep=""));
          if(output[155]=="Cluster\tMean\ts.e.\t \t \t \t \t") {
            distal.mean.by.class.proportional.none <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][2],
                                                                   strsplit(output[156:158],split="\t")[[c2]][2],
                                                                   strsplit(output[156:158],split="\t")[[c3]][2]));
            distal.mean.by.class.proportional.none.se <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][3],
                                                                      strsplit(output[156:158],split="\t")[[c2]][3],
                                                                      strsplit(output[156:158],split="\t")[[c3]][3]));
          } else {
            distal.mean.by.class.proportional.none <- rep(NA,3);
            distal.mean.by.class.proportional.none.se <- rep(NA,3);
            print(paste("Proportional none failed on ",simulation));
          }
          # Latent GOLD analysis:  Modal ML;
          shell(cmd=paste(LatentGOLD.path," ",working.path,"ModalML.lgs /b",sep=""));
          output <- readLines(con=paste(working.path,"ModalML.lst",sep=""));
          if(output[155]=="Cluster\tMean\ts.e.\t \t \t \t \t") {
            distal.mean.by.class.modal.ML <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][2],
                                                          strsplit(output[156:158],split="\t")[[c2]][2],
                                                          strsplit(output[156:158],split="\t")[[c3]][2]));
            distal.mean.by.class.modal.ML.se <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][3],
                                                             strsplit(output[156:158],split="\t")[[c2]][3],
                                                             strsplit(output[156:158],split="\t")[[c3]][3]));
          } else {
            distal.mean.by.class.modal.ML <- rep(NA,3);
            distal.mean.by.class.modal.ML.se <- rep(NA,3);
            print(paste("Modal ML failed on ",simulation));
          }
          # Latent GOLD analysis:  Proportional ML;
          shell(cmd=paste(LatentGOLD.path," ",working.path,"ProportionalML.lgs /b",sep=""));
          output <- readLines(con=paste(working.path,"ProportionalML.lst",sep=""));
          if(output[155]=="Cluster\tMean\ts.e.\t \t \t \t \t") {
            distal.mean.by.class.proportional.ML <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][2],
                                                                 strsplit(output[156:158],split="\t")[[c2]][2],
                                                                 strsplit(output[156:158],split="\t")[[c3]][2]));
            distal.mean.by.class.proportional.ML.se <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][3],
                                                                    strsplit(output[156:158],split="\t")[[c2]][3],
                                                                    strsplit(output[156:158],split="\t")[[c3]][3]));
          } else {
            distal.mean.by.class.proportional.ML <- rep(NA,3);
            distal.mean.by.class.proportional.ML.se <- rep(NA,3);
            print(paste("Proportional ML failed on ",simulation));
          }
          # Latent GOLD analysis:  Modal BCH;
          shell(cmd=paste(LatentGOLD.path," ",working.path,"ModalBCH.lgs /b",sep=""));
          output <- readLines(con=paste(working.path,"ModalBCH.lst",sep=""));
          if(output[155]=="Cluster\tMean\ts.e.\t \t \t \t \t") {
            distal.mean.by.class.modal.BCH <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][2],
                                                           strsplit(output[156:158],split="\t")[[c2]][2],
                                                           strsplit(output[156:158],split="\t")[[c3]][2]));
            distal.mean.by.class.modal.BCH.se <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][3],
                                                              strsplit(output[156:158],split="\t")[[c2]][3],
                                                              strsplit(output[156:158],split="\t")[[c3]][3]));
          } else {
            distal.mean.by.class.modal.BCH <- rep(NA,3);
            distal.mean.by.class.modal.BCH.se <- rep(NA,3);
            print(paste("Modal BCH failed on ",simulation));
          }
          # Latent GOLD analysis:  Proportional BCH;
          shell(cmd=paste(LatentGOLD.path," ",working.path,"ProportionalBCH.lgs /b",sep=""));
          if(output[155]=="Cluster\tMean\ts.e.\t \t \t \t \t") {
            distal.mean.by.class.proportional.BCH <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][2],
                                                                  strsplit(output[156:158],split="\t")[[c2]][2],
                                                                  strsplit(output[156:158],split="\t")[[c3]][2]));
            distal.mean.by.class.proportional.BCH.se <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][3],
                                                                     strsplit(output[156:158],split="\t")[[c2]][3],
                                                                     strsplit(output[156:158],split="\t")[[c3]][3]));
          } else {
            distal.mean.by.class.proportional.BCH <- rep(NA,3);
            distal.mean.by.class.proportional.BCH.se <- rep(NA,3);
            print(paste("Proportional BCH failed on ",simulation));
          }
          # Fit an inclusive LCA using Latent GOLD;
          if (file.exists(paste(working.path,"InclusivePostProbs.dat",sep=""))) {
            file.remove(paste(working.path,"InclusivePostProbs.dat",sep=""));
          }
          if (file.exists(paste(working.path,"FitInclusive.lst",sep=""))) {
            file.remove(paste(working.path,"FitInclusive.lst",sep=""));
          }
          shell(cmd=paste(LatentGOLD.path," ",working.path,"FitInclusive.lgs /b",sep=""));
          output <- readLines(con=paste(working.path,"FitInclusive.lst",sep=""));
          # Read the estimated item means for the inclusive LCA;
          if (output[222]=="Profile") {
            temp <- strsplit(output[c(227,229,231,233,235)],2,split="\t");
            unsorted.item.means <- rbind(as.numeric(temp[[1]][c(2,4,6)]),
                                         as.numeric(temp[[2]][c(2,4,6)]),
                                         as.numeric(temp[[3]][c(2,4,6)]),
                                         as.numeric(temp[[4]][c(2,4,6)]),
                                         as.numeric(temp[[5]][c(2,4,6)]));
            # Determine the random class order for the inclusive LCA;
            possible.permutations <- list(c(1,2,3),
                                          c(1,3,2),
                                          c(2,1,3),
                                          c(2,3,1),
                                          c(3,1,2),
                                          c(3,2,1));
            sum.sqd.error.by.permutation <- rep(NA,6);
            for (this.permutation in 1:6) {
              sum.sqd.error.by.permutation[this.permutation] <-
                sum((unsorted.item.means[,possible.permutations[[this.permutation]]]-true.item.means)^2);
            }
            best.permutation <- possible.permutations[[which.min(sum.sqd.error.by.permutation)]];
            c1 <- best.permutation[1];
            c2 <- best.permutation[2];
            c3 <- best.permutation[3];
            sorted.item.means.inclusive <- unsorted.item.means[,best.permutation];
            # Latent GOLD analysis:  Inclusive Modal None;
            shell(cmd=paste(LatentGOLD.path," ",working.path,"InclusiveModalNone.lgs /b",sep=""));
            output <- readLines(con=paste(working.path,"InclusiveModalNone.lst",sep=""));
            if (output[155]=="Cluster\tMean\ts.e.\t \t \t \t \t") {
              distal.mean.by.class.inclusive.modal.none <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][2],
                                                                        strsplit(output[156:158],split="\t")[[c2]][2],
                                                                        strsplit(output[156:158],split="\t")[[c3]][2]));
              distal.mean.by.class.inclusive.modal.none.se <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][3],
                                                                           strsplit(output[156:158],split="\t")[[c2]][3],
                                                                           strsplit(output[156:158],split="\t")[[c3]][3]));
            } else {
              distal.mean.by.class.inclusive.modal.none <- rep(NA,3);
              distal.mean.by.class.inclusive.modal.none.se <- rep(NA,3);
              print(paste("Inclusive modal failed on ",simulation));
            }
            # Latent GOLD analysis:  Inclusive Proportional None
            shell(cmd=paste(LatentGOLD.path," ",working.path,"InclusiveProportionalNone.lgs /b",sep=""));
            output <- readLines(con=paste(working.path,"InclusiveProportionalNone.lst",sep=""));
            if (output[155]=="Cluster\tMean\ts.e.\t \t \t \t \t") {
              distal.mean.by.class.inclusive.proportional.none <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][2],
                                                                               strsplit(output[156:158],split="\t")[[c2]][2],
                                                                               strsplit(output[156:158],split="\t")[[c3]][2]));
              distal.mean.by.class.inclusive.proportional.none.se <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][3],
                                                                                  strsplit(output[156:158],split="\t")[[c2]][3],
                                                                                  strsplit(output[156:158],split="\t")[[c3]][3]));
            } else {
              distal.mean.by.class.inclusive.proportional.none <- rep(NA,3);
              distal.mean.by.class.inclusive.proportional.none.se <- rep(NA,3);
              print(paste("Inclusive proportional failed on ",simulation));
            }
          } else {
              distal.mean.by.class.inclusive.modal.none <- rep(NA,3);
              distal.mean.by.class.inclusive.modal.none.se <- rep(NA,3);
              distal.mean.by.class.inclusive.proportional.none <- rep(NA,3);
              distal.mean.by.class.inclusive.proportional.none.se <- rep(NA,3);
              print(paste("Inclusive estimation failed on ",simulation));
          }
          # Fit a quadratic LCA using Latent GOLD;
          if (file.exists(paste(working.path,"QuadraticPostProbs.dat",sep=""))) {
            file.remove(paste(working.path,"QuadraticPostProbs.dat",sep=""));
          }
          if (file.exists(paste(working.path,"FitQuadratic.lst",sep=""))) {
            file.remove(paste(working.path,"FitQuadratic.lst",sep=""));
          }
          shell(cmd=paste(LatentGOLD.path," ",working.path,"FitQuadratic.lgs /b",sep=""));
          output <- readLines(con=paste(working.path,"FitQuadratic.lst",sep=""));
          if (output[246]=="Profile") {
            temp <- strsplit(output[c(251,253,255,257,259)],2,split="\t");
            unsorted.item.means <- rbind(as.numeric(temp[[1]][c(2,4,6)]),
                                         as.numeric(temp[[2]][c(2,4,6)]),
                                         as.numeric(temp[[3]][c(2,4,6)]),
                                         as.numeric(temp[[4]][c(2,4,6)]),
                                         as.numeric(temp[[5]][c(2,4,6)]));
            # Determine the random class order for the quadratic LCA;
            possible.permutations <- list(c(1,2,3),
                                          c(1,3,2),
                                          c(2,1,3),
                                          c(2,3,1),
                                          c(3,1,2),
                                          c(3,2,1));
            sum.sqd.error.by.permutation <- rep(NA,6);
            for (this.permutation in 1:6) {
              sum.sqd.error.by.permutation[this.permutation] <-
                sum((unsorted.item.means[,possible.permutations[[this.permutation]]]-true.item.means)^2);
            }
            best.permutation <- possible.permutations[[which.min(sum.sqd.error.by.permutation)]];
            c1 <- best.permutation[1];
            c2 <- best.permutation[2];
            c3 <- best.permutation[3];
            sorted.item.means.quadratic <- unsorted.item.means[,best.permutation];
            # Latent GOLD analysis:  Quadratic Modal None;
            shell(cmd=paste(LatentGOLD.path," ",working.path,"QuadraticModalNone.lgs /b",sep=""));
            output <- readLines(con=paste(working.path,"QuadraticModalNone.lst",sep=""));
            if (output[155]=="Cluster\tMean\ts.e.\t \t \t \t \t") {
              distal.mean.by.class.quadratic.modal.none <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][2],
                                                                        strsplit(output[156:158],split="\t")[[c2]][2],
                                                                        strsplit(output[156:158],split="\t")[[c3]][2]));
              distal.mean.by.class.quadratic.modal.none.se <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][3],
                                                                           strsplit(output[156:158],split="\t")[[c2]][3],
                                                                           strsplit(output[156:158],split="\t")[[c3]][3]));
            } else {
              distal.mean.by.class.quadratic.modal.none <- rep(NA,3);
              distal.mean.by.class.quadratic.modal.none.se <- rep(NA,3);
              print(paste("Quadratic modal failed on ",simulation));
            }
            # Latent GOLD analysis:  Quadratic Proportional None
            shell(cmd=paste(LatentGOLD.path," ",working.path,"QuadraticProportionalNone.lgs /b",sep=""));
            output <- readLines(con=paste(working.path,"QuadraticProportionalNone.lst",sep=""));
            if (output[155]=="Cluster\tMean\ts.e.\t \t \t \t \t") {
              distal.mean.by.class.quadratic.proportional.none <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][2],
                                                                               strsplit(output[156:158],split="\t")[[c2]][2],
                                                                               strsplit(output[156:158],split="\t")[[c3]][2]));
              distal.mean.by.class.quadratic.proportional.none.se <- as.numeric(c(strsplit(output[156:158],split="\t")[[c1]][3],
                                                                                  strsplit(output[156:158],split="\t")[[c2]][3],
                                                                                  strsplit(output[156:158],split="\t")[[c3]][3]));
            } else {
              distal.mean.by.class.quadratic.proportional.none <- rep(NA,3);
              distal.mean.by.class.quadratic.proportional.none.se <- rep(NA,3);
              print(paste("Quadratic proportional failed on ",simulation));
            }
          } else {
              distal.mean.by.class.quadratic.modal.none <- rep(NA,3);
              distal.mean.by.class.quadratic.modal.none.se <- rep(NA,3);
              distal.mean.by.class.quadratic.proportional.none <- rep(NA,3);
              distal.mean.by.class.quadratic.proportional.none.se <- rep(NA,3);
              print(paste("Quadratic estimation failed on ",simulation));
          }
          these.preliminary.answers <- c(simulation = simulation,
                                         Item1Oracle = oracle.item.means[1,],
                                         Item2Oracle = oracle.item.means[2,],
                                         Item3Oracle = oracle.item.means[3,],
                                         Item4Oracle = oracle.item.means[4,],
                                         Item5Oracle = oracle.item.means[5,],
                                         DistalByClassOracle = oracle.distal.mean.by.class,
                                         DistalByClassOracleSE = oracle.distal.mean.by.class.se,
                                         SizeByClassOracle = oracle.class.sizes,
                                         DistalOracle = oracle.distal.prob,
                                         NoninclusiveEntropyRamaswamy = noninclusive.entropy.Ramaswamy,
                                         NoninclusiveEntropyRaw = noninclusive.entropy.raw,
                                         NoninclusiveEntropyRsqd = noninclusive.entropy.Rsqd,
                                         Item1Noninclusive = sorted.item.means.noninclusive[1,],
                                         Item2Noninclusive = sorted.item.means.noninclusive[2,],
                                         Item3Noninclusive = sorted.item.means.noninclusive[3,],
                                         Item4Noninclusive = sorted.item.means.noninclusive[4,],
                                         Item5Noninclusive = sorted.item.means.noninclusive[5,],
                                         DistalByClassModalNone = distal.mean.by.class.modal.none,
                                         DistalByClassModalNoneSE = distal.mean.by.class.modal.none.se,
                                         DistalByClassProportionalNone = distal.mean.by.class.proportional.none,
                                         DistalByClassProportionalNoneSE = distal.mean.by.class.proportional.none.se,
                                         DistalByClassModalML = distal.mean.by.class.modal.ML,
                                         DistalByClassModalMLSE = distal.mean.by.class.modal.ML.se,
                                         DistalByClassProportionalML = distal.mean.by.class.proportional.ML,
                                         DistalByClassProportionalMLSE = distal.mean.by.class.proportional.ML.se,
                                         DistalByClassModalBCH = distal.mean.by.class.modal.BCH,
                                         DistalByClassModalBCHSE = distal.mean.by.class.modal.BCH.se,
                                         DistalByClassProportionalBCH = distal.mean.by.class.proportional.BCH,
                                         DistalByClassProportionalBCHSE = distal.mean.by.class.proportional.BCH.se,
                                         Item1Inclusive = sorted.item.means.inclusive[1,],
                                         Item2Inclusive = sorted.item.means.inclusive[2,],
                                         Item3Inclusive = sorted.item.means.inclusive[3,],
                                         Item4Inclusive = sorted.item.means.inclusive[4,],
                                         Item5Inclusive = sorted.item.means.inclusive[5,],
                                         DistalByClassModalInclusive = distal.mean.by.class.inclusive.modal.none,
                                         DistalByClassModalInclusiveSE = distal.mean.by.class.inclusive.modal.none.se,
                                         DistalByClassProportionalInclusive = distal.mean.by.class.inclusive.proportional.none,
                                         DistalByClassProportionalInclusiveSE = distal.mean.by.class.inclusive.proportional.none.se,
                                         DistalByClassModalQuadratic = distal.mean.by.class.quadratic.modal.none,
                                         DistalByClassModalQuadraticSE = distal.mean.by.class.quadratic.modal.none.se,
                                         DistalByClassProportionalQuadratic = distal.mean.by.class.quadratic.proportional.none,
                                         DistalByClassProportionalQuadraticSE = distal.mean.by.class.quadratic.proportional.none.se);
      }
      rownames(preliminary.answers) <- NULL;
      rownames(these.preliminary.answers) <- NULL;
      preliminary.answers <- rbind(preliminary.answers,these.preliminary.answers);
      save.image(paste("working",
                       which.measurement.quality,
                       which.class.size.distribution,
                       which.effect.size,
                       ".rdata",
                       sep=""));
      }
      these.final.answers <- c(  which.measurement.quality=which.measurement.quality,
                                 which.class.size.distribution=which.class.size.distribution,
                                 which.effect.size=which.effect.size,
                                 NoninclusiveEntropyRamaswamy = mean(preliminary.answers[,"NoninclusiveEntropyRamaswamy"]),
                                 NoninclusiveEntropyRaw = mean(preliminary.answers[,"NoninclusiveEntropyRaw"]),
                                 NoninclusiveEntropyRsqd = mean(preliminary.answers[,"NoninclusiveEntropyRsqd"]),
                                 WorkedOracle1 = mean(!is.na(preliminary.answers[,"DistalByClassOracle1"])),
                                 WorkedModalNone1 = mean(!is.na(preliminary.answers[,"DistalByClassModalNone1"])),
                                 WorkedModalML1 = mean(!is.na(preliminary.answers[,"DistalByClassModalML1"])),
                                 WorkedModalBCH1 = mean(!is.na(preliminary.answers[,"DistalByClassModalBCH1"])),
                                 WorkedModalInclusive1 = mean(!is.na(preliminary.answers[,"DistalByClassModalInclusive1"])),
                                 WorkedModalQuadratic1 = mean(!is.na(preliminary.answers[,"DistalByClassModalQuadratic1"])),
                                 WorkedProportionalNone1 = mean(!is.na(preliminary.answers[,"DistalByClassProportionalNone1"])),
                                 WorkedProportionalML1 = mean(!is.na(preliminary.answers[,"DistalByClassProportionalML1"])),
                                 WorkedProportionalBCH1 = mean(!is.na(preliminary.answers[,"DistalByClassProportionalBCH1"])),
                                 WorkedProportionalInclusive1 = mean(!is.na(preliminary.answers[,"DistalByClassProportionalInclusive1"])),
                                 WorkedProportionalQuadratic1 = mean(!is.na(preliminary.answers[,"DistalByClassProportionalQuadratic1"])),
                                 BiasOracle1 = mean(preliminary.answers[,"DistalByClassOracle1"] - true.distal.mean.by.class[1],na.rm=TRUE),
                                 BiasModalNone1 = mean(preliminary.answers[,"DistalByClassModalNone1"] - true.distal.mean.by.class[1],na.rm=TRUE),
                                 BiasModalML1 = mean(preliminary.answers[,"DistalByClassModalML1"] - true.distal.mean.by.class[1],na.rm=TRUE),
                                 BiasModalBCH1 = mean(preliminary.answers[,"DistalByClassModalBCH1"] - true.distal.mean.by.class[1],na.rm=TRUE),
                                 BiasModalInclusive1 = mean(preliminary.answers[,"DistalByClassModalInclusive1"] - true.distal.mean.by.class[1],na.rm=TRUE),
                                 BiasModalQuadratic1 = mean(preliminary.answers[,"DistalByClassModalQuadratic1"] - true.distal.mean.by.class[1],na.rm=TRUE),
                                 BiasProportionalNone1 = mean(preliminary.answers[,"DistalByClassProportionalNone1"] - true.distal.mean.by.class[1],na.rm=TRUE),
                                 BiasProportionalML1 = mean(preliminary.answers[,"DistalByClassProportionalML1"] - true.distal.mean.by.class[1],na.rm=TRUE),
                                 BiasProportionalBCH1 = mean(preliminary.answers[,"DistalByClassProportionalBCH1"] - true.distal.mean.by.class[1],na.rm=TRUE),
                                 BiasProportionalInclusive1 = mean(preliminary.answers[,"DistalByClassProportionalInclusive1"] - true.distal.mean.by.class[1],na.rm=TRUE),
                                 BiasProportionalQuadratic1 = mean(preliminary.answers[,"DistalByClassProportionalQuadratic1"] - true.distal.mean.by.class[1],na.rm=TRUE),
                                 BiasOracle2 = mean(preliminary.answers[,"DistalByClassOracle2"] - true.distal.mean.by.class[2],na.rm=TRUE),
                                 BiasModalNone2 = mean(preliminary.answers[,"DistalByClassModalNone2"] - true.distal.mean.by.class[2],na.rm=TRUE),
                                 BiasModalML2 = mean(preliminary.answers[,"DistalByClassModalML2"] - true.distal.mean.by.class[2],na.rm=TRUE),
                                 BiasModalBCH2 = mean(preliminary.answers[,"DistalByClassModalBCH2"] - true.distal.mean.by.class[2],na.rm=TRUE),
                                 BiasModalInclusive2 = mean(preliminary.answers[,"DistalByClassModalInclusive2"] - true.distal.mean.by.class[2],na.rm=TRUE),
                                 BiasModalQuadratic2 = mean(preliminary.answers[,"DistalByClassModalQuadratic2"] - true.distal.mean.by.class[2],na.rm=TRUE),
                                 BiasProportionalNone2 = mean(preliminary.answers[,"DistalByClassProportionalNone2"] - true.distal.mean.by.class[2],na.rm=TRUE),
                                 BiasProportionalML2 = mean(preliminary.answers[,"DistalByClassProportionalML2"] - true.distal.mean.by.class[2],na.rm=TRUE),
                                 BiasProportionalBCH2 = mean(preliminary.answers[,"DistalByClassProportionalBCH2"] - true.distal.mean.by.class[2],na.rm=TRUE),
                                 BiasProportionalInclusive2 = mean(preliminary.answers[,"DistalByClassProportionalInclusive2"] - true.distal.mean.by.class[2],na.rm=TRUE),
                                 BiasProportionalQuadratic2 = mean(preliminary.answers[,"DistalByClassProportionalQuadratic2"] - true.distal.mean.by.class[2],na.rm=TRUE),
                                 BiasOracle3 = mean(preliminary.answers[,"DistalByClassOracle3"] - true.distal.mean.by.class[3],na.rm=TRUE),
                                 BiasModalNone3 = mean(preliminary.answers[,"DistalByClassModalNone3"] - true.distal.mean.by.class[3],na.rm=TRUE),
                                 BiasModalML3 = mean(preliminary.answers[,"DistalByClassModalML3"] - true.distal.mean.by.class[3],na.rm=TRUE),
                                 BiasModalBCH3 = mean(preliminary.answers[,"DistalByClassModalBCH3"] - true.distal.mean.by.class[3],na.rm=TRUE),
                                 BiasModalInclusive3 = mean(preliminary.answers[,"DistalByClassModalInclusive3"] - true.distal.mean.by.class[3],na.rm=TRUE),
                                 BiasModalQuadratic3 = mean(preliminary.answers[,"DistalByClassModalQuadratic3"] - true.distal.mean.by.class[3],na.rm=TRUE),
                                 BiasProportionalNone3 = mean(preliminary.answers[,"DistalByClassProportionalNone3"] - true.distal.mean.by.class[3],na.rm=TRUE),
                                 BiasProportionalML3 = mean(preliminary.answers[,"DistalByClassProportionalML3"] - true.distal.mean.by.class[3],na.rm=TRUE),
                                 BiasProportionalBCH3 = mean(preliminary.answers[,"DistalByClassProportionalBCH3"] - true.distal.mean.by.class[3],na.rm=TRUE),
                                 BiasProportionalInclusive3 = mean(preliminary.answers[,"DistalByClassProportionalInclusive3"] - true.distal.mean.by.class[3],na.rm=TRUE),
                                 BiasProportionalQuadratic3 = mean(preliminary.answers[,"DistalByClassProportionalQuadratic3"] - true.distal.mean.by.class[3],na.rm=TRUE),
                                 RootMSEOracle1 = sqrt(mean((preliminary.answers[,"DistalByClassOracle1"] - true.distal.mean.by.class[1])^2,na.rm=TRUE)),
                                 RootMSEModalNone1 = sqrt(mean((preliminary.answers[,"DistalByClassModalNone1"] - true.distal.mean.by.class[1])^2,na.rm=TRUE)),
                                 RootMSEModalML1 = sqrt(mean((preliminary.answers[,"DistalByClassModalML1"] - true.distal.mean.by.class[1])^2,na.rm=TRUE)),
                                 RootMSEModalBCH1 = sqrt(mean((preliminary.answers[,"DistalByClassModalBCH1"] - true.distal.mean.by.class[1])^2,na.rm=TRUE)),
                                 RootMSEModalInclusive1 = sqrt(mean((preliminary.answers[,"DistalByClassModalInclusive1"] - true.distal.mean.by.class[1])^2,na.rm=TRUE)),
                                 RootMSEModalQuadratic1 = sqrt(mean((preliminary.answers[,"DistalByClassModalQuadratic1"] - true.distal.mean.by.class[1])^2,na.rm=TRUE)),
                                 RootMSEProportionalNone1 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalNone1"] - true.distal.mean.by.class[1])^2,na.rm=TRUE)),
                                 RootMSEProportionalML1 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalML1"] - true.distal.mean.by.class[1])^2,na.rm=TRUE)),
                                 RootMSEProportionalBCH1 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalBCH1"] - true.distal.mean.by.class[1])^2,na.rm=TRUE)),
                                 RootMSEProportionalInclusive1 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalInclusive1"] - true.distal.mean.by.class[1])^2,na.rm=TRUE)),
                                 RootMSEProportionalQuadratic1 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalQuadratic1"] - true.distal.mean.by.class[1])^2,na.rm=TRUE)),
                                 RootMSEOracle2 = sqrt(mean((preliminary.answers[,"DistalByClassOracle2"] - true.distal.mean.by.class[2])^2,na.rm=TRUE)),
                                 RootMSEModalNone2 = sqrt(mean((preliminary.answers[,"DistalByClassModalNone2"] - true.distal.mean.by.class[2])^2,na.rm=TRUE)),
                                 RootMSEModalML2 = sqrt(mean((preliminary.answers[,"DistalByClassModalML2"] - true.distal.mean.by.class[2])^2,na.rm=TRUE)),
                                 RootMSEModalBCH2 = sqrt(mean((preliminary.answers[,"DistalByClassModalBCH2"] - true.distal.mean.by.class[2])^2,na.rm=TRUE)),
                                 RootMSEModalInclusive2 = sqrt(mean((preliminary.answers[,"DistalByClassModalInclusive2"] - true.distal.mean.by.class[2])^2,na.rm=TRUE)),
                                 RootMSEModalQuadratic2 = sqrt(mean((preliminary.answers[,"DistalByClassModalQuadratic2"] - true.distal.mean.by.class[2])^2,na.rm=TRUE)),
                                 RootMSEProportionalNone2 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalNone2"] - true.distal.mean.by.class[2])^2,na.rm=TRUE)),
                                 RootMSEProportionalML2 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalML2"] - true.distal.mean.by.class[2])^2,na.rm=TRUE)),
                                 RootMSEProportionalBCH2 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalBCH2"] - true.distal.mean.by.class[2])^2,na.rm=TRUE)),
                                 RootMSEProportionalInclusive2 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalInclusive2"] - true.distal.mean.by.class[2])^2,na.rm=TRUE)),
                                 RootMSEProportionalQuadratic2 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalQuadratic2"] - true.distal.mean.by.class[2])^2,na.rm=TRUE)),
                                 RootMSEOracle3 = sqrt(mean((preliminary.answers[,"DistalByClassOracle3"] - true.distal.mean.by.class[3])^2,na.rm=TRUE)),
                                 RootMSEModalNone3 = sqrt(mean((preliminary.answers[,"DistalByClassModalNone3"] - true.distal.mean.by.class[3])^2,na.rm=TRUE)),
                                 RootMSEModalML3 = sqrt(mean((preliminary.answers[,"DistalByClassModalML3"] - true.distal.mean.by.class[3])^2,na.rm=TRUE)),
                                 RootMSEModalBCH3 = sqrt(mean((preliminary.answers[,"DistalByClassModalBCH3"] - true.distal.mean.by.class[3])^2,na.rm=TRUE)),
                                 RootMSEModalInclusive3 = sqrt(mean((preliminary.answers[,"DistalByClassModalInclusive3"] - true.distal.mean.by.class[3])^2,na.rm=TRUE)),
                                 RootMSEModalQuadratic3 = sqrt(mean((preliminary.answers[,"DistalByClassModalQuadratic3"] - true.distal.mean.by.class[3])^2,na.rm=TRUE)),
                                 RootMSEProportionalNone3 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalNone3"] - true.distal.mean.by.class[3])^2,na.rm=TRUE)),
                                 RootMSEProportionalML3 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalML3"] - true.distal.mean.by.class[3])^2,na.rm=TRUE)),
                                 RootMSEProportionalBCH3 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalBCH3"] - true.distal.mean.by.class[3])^2,na.rm=TRUE)),
                                 RootMSEProportionalInclusive3 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalInclusive3"] - true.distal.mean.by.class[3])^2,na.rm=TRUE)),
                                 RootMSEProportionalQuadratic3 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalQuadratic3"] - true.distal.mean.by.class[3])^2,na.rm=TRUE)),
                                 CoverageOracle1 = mean(abs(preliminary.answers[,"DistalByClassOracle1"]-true.distal.mean.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassOracleSE1"],na.rm=TRUE),
                                 CoverageModalNone1 = mean(abs(preliminary.answers[,"DistalByClassModalNone1"]-true.distal.mean.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassModalNoneSE1"],na.rm=TRUE),
                                 CoverageModalML1 = mean(abs(preliminary.answers[,"DistalByClassModalML1"]-true.distal.mean.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassModalMLSE1"],na.rm=TRUE),
                                 CoverageModalBCH1 = mean(abs(preliminary.answers[,"DistalByClassModalBCH1"]-true.distal.mean.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassModalBCHSE1"],na.rm=TRUE),
                                 CoverageModalInclusive1 = mean(abs(preliminary.answers[,"DistalByClassModalInclusive1"]-true.distal.mean.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassModalInclusiveSE1"],na.rm=TRUE),
                                 CoverageModalQuadratic1 = mean(abs(preliminary.answers[,"DistalByClassModalQuadratic1"]-true.distal.mean.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassModalQuadraticSE1"],na.rm=TRUE),
                                 CoverageProportionalNone1 = mean(abs(preliminary.answers[,"DistalByClassProportionalNone1"]-true.distal.mean.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassProportionalNoneSE1"],na.rm=TRUE),
                                 CoverageProportionalML1 = mean(abs(preliminary.answers[,"DistalByClassProportionalML1"]-true.distal.mean.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassProportionalMLSE1"],na.rm=TRUE),
                                 CoverageProportionalBCH1 = mean(abs(preliminary.answers[,"DistalByClassProportionalBCH1"]-true.distal.mean.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassProportionalBCHSE1"],na.rm=TRUE),
                                 CoverageProportionalInclusive1 = mean(abs(preliminary.answers[,"DistalByClassProportionalInclusive1"]-true.distal.mean.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassProportionalInclusiveSE1"],na.rm=TRUE),
                                 CoverageProportionalQuadratic1 = mean(abs(preliminary.answers[,"DistalByClassProportionalQuadratic1"]-true.distal.mean.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassProportionalQuadraticSE1"],na.rm=TRUE),
                                 CoverageOracle2 = mean(abs(preliminary.answers[,"DistalByClassOracle2"]-true.distal.mean.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassOracleSE2"],na.rm=TRUE),
                                 CoverageModalNone2 = mean(abs(preliminary.answers[,"DistalByClassModalNone2"]-true.distal.mean.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassModalNoneSE2"],na.rm=TRUE),
                                 CoverageModalML2 = mean(abs(preliminary.answers[,"DistalByClassModalML2"]-true.distal.mean.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassModalMLSE2"],na.rm=TRUE),
                                 CoverageModalBCH2 = mean(abs(preliminary.answers[,"DistalByClassModalBCH2"]-true.distal.mean.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassModalBCHSE2"],na.rm=TRUE),
                                 CoverageModalInclusive2 = mean(abs(preliminary.answers[,"DistalByClassModalInclusive2"]-true.distal.mean.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassModalInclusiveSE2"],na.rm=TRUE),
                                 CoverageModalQuadratic2 = mean(abs(preliminary.answers[,"DistalByClassModalQuadratic2"]-true.distal.mean.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassModalQuadraticSE2"],na.rm=TRUE),
                                 CoverageProportionalNone2 = mean(abs(preliminary.answers[,"DistalByClassProportionalNone2"]-true.distal.mean.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassProportionalNoneSE2"],na.rm=TRUE),
                                 CoverageProportionalML2 = mean(abs(preliminary.answers[,"DistalByClassProportionalML2"]-true.distal.mean.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassProportionalMLSE2"],na.rm=TRUE),
                                 CoverageProportionalBCH2 = mean(abs(preliminary.answers[,"DistalByClassProportionalBCH2"]-true.distal.mean.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassProportionalBCHSE2"],na.rm=TRUE),
                                 CoverageProportionalInclusive2 = mean(abs(preliminary.answers[,"DistalByClassProportionalInclusive2"]-true.distal.mean.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassProportionalInclusiveSE2"],na.rm=TRUE),
                                 CoverageProportionalQuadratic2 = mean(abs(preliminary.answers[,"DistalByClassProportionalQuadratic2"]-true.distal.mean.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassProportionalQuadraticSE2"],na.rm=TRUE),
                                 CoverageOracle3 = mean(abs(preliminary.answers[,"DistalByClassOracle3"]-true.distal.mean.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassOracleSE3"],na.rm=TRUE),
                                 CoverageModalNone3 = mean(abs(preliminary.answers[,"DistalByClassModalNone3"]-true.distal.mean.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassModalNoneSE3"],na.rm=TRUE),
                                 CoverageModalML3 = mean(abs(preliminary.answers[,"DistalByClassModalML3"]-true.distal.mean.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassModalMLSE3"],na.rm=TRUE),
                                 CoverageModalBCH3 = mean(abs(preliminary.answers[,"DistalByClassModalBCH3"]-true.distal.mean.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassModalBCHSE3"],na.rm=TRUE),
                                 CoverageModalInclusive3 = mean(abs(preliminary.answers[,"DistalByClassModalInclusive3"]-true.distal.mean.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassModalInclusiveSE3"],na.rm=TRUE),
                                 CoverageModalQuadratic3 = mean(abs(preliminary.answers[,"DistalByClassModalQuadratic3"]-true.distal.mean.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassModalQuadraticSE3"],na.rm=TRUE),
                                 CoverageProportionalNone3 = mean(abs(preliminary.answers[,"DistalByClassProportionalNone3"]-true.distal.mean.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassProportionalNoneSE3"],na.rm=TRUE),
                                 CoverageProportionalML3 = mean(abs(preliminary.answers[,"DistalByClassProportionalML3"]-true.distal.mean.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassProportionalMLSE3"],na.rm=TRUE),
                                 CoverageProportionalBCH3 = mean(abs(preliminary.answers[,"DistalByClassProportionalBCH3"]-true.distal.mean.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassProportionalBCHSE3"],na.rm=TRUE),
                                 CoverageProportionalInclusive3 = mean(abs(preliminary.answers[,"DistalByClassProportionalInclusive3"]-true.distal.mean.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassProportionalInclusiveSE3"],na.rm=TRUE),
                                 CoverageProportionalQuadratic3 = mean(abs(preliminary.answers[,"DistalByClassProportionalQuadratic3"]-true.distal.mean.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassProportionalQuadraticSE3"],na.rm=TRUE));
      save.image(file=paste("Preliminary",which.measurement.quality,which.class.size.distribution,which.effect.size,".rdata",sep=""));
      rm(preliminary.answers,these.preliminary.answers);
      final.answers <- rbind(final.answers,these.final.answers);
    }
  }
}
finish.time <- Sys.time();
final.answers <- data.frame(final.answers);
print(apply(final.answers,2,mean));
write.table(x=final.answers,file=paste("simanswers",nsim,".txt",sep=""),
            quote=FALSE,row.names=TRUE,col.names=TRUE);
save.image(paste("FinishedSimulations-",nsim,".rdata",sep=""));