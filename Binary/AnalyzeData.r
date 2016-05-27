rm(list=ls(all=TRUE));
#########################################################
# Set initial settings:
set.seed(700046);
input.path <- "C:\\Users\\jjd264\\Documents\\Sims-Lpa-Binary\\Datasets\\";
working.path <- "C:\\Users\\jjd264\\Documents\\Sims-Lpa-Binary\\";
LatentGOLD.path <- "C:\\Users\\jjd264\\DOCUME~1\\LATENT~1.1\\lg51.exe";
N <- 1000;
start.time <- Sys.time();
nsim <- 1000;
final.answers <- NULL;
all.preliminary.answers <- NULL;
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
        oracle.distal.prob.by.class <- rep(NA,3);
        oracle.distal.prob.by.class.se <- rep(NA,3);
        items <- cbind(sim.data$I1,sim.data$I2,sim.data$I3,sim.data$I4,sim.data$I5);
        for (this.class in 1:3) {
          oracle.item.means[,this.class] <- apply(items[which(sim.data$trueclass==this.class),],2,mean);
          oracle.distal.prob.by.class[this.class] <- mean(sim.data$y[which(sim.data$trueclass==this.class)]);
          oracle.distal.prob.by.class.se[this.class] <- sd(sim.data$y[which(sim.data$trueclass==this.class)]) / 
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
        stopifnot(output[152]=="Cluster\t0\ts.e.\t1\ts.e.\t \t \t");
        distal.prob.by.class.modal.none <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][4],
                                                        strsplit(output[153:155],split="\t")[[c2]][4],
                                                        strsplit(output[153:155],split="\t")[[c3]][4]));
        distal.prob.by.class.modal.none.se <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][5],
                                                           strsplit(output[153:155],split="\t")[[c2]][5],
                                                           strsplit(output[153:155],split="\t")[[c3]][5])); 
        # Latent GOLD analysis:  Proportional None
        shell(cmd=paste(LatentGOLD.path," ",working.path,"ProportionalNone.lgs /b",sep=""));
        output <- readLines(con=paste(working.path,"ProportionalNone.lst",sep=""));
        stopifnot(output[152]=="Cluster\t0\ts.e.\t1\ts.e.\t \t \t");
        distal.prob.by.class.proportional.none <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][4],
                                                               strsplit(output[153:155],split="\t")[[c2]][4],
                                                               strsplit(output[153:155],split="\t")[[c3]][4]));
        distal.prob.by.class.proportional.none.se <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][5],
                                                                  strsplit(output[153:155],split="\t")[[c2]][5],
                                                                  strsplit(output[153:155],split="\t")[[c3]][5]));
        # Latent GOLD analysis:  Modal ML;
        shell(cmd=paste(LatentGOLD.path," ",working.path,"ModalML.lgs /b",sep=""));
        output <- readLines(con=paste(working.path,"ModalML.lst",sep=""));
        stopifnot(output[152]=="Cluster\t0\ts.e.\t1\ts.e.\t \t \t");
        distal.prob.by.class.modal.ML <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][4],
                                                      strsplit(output[153:155],split="\t")[[c2]][4],
                                                      strsplit(output[153:155],split="\t")[[c3]][4]));
        distal.prob.by.class.modal.ML.se <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][5],
                                                         strsplit(output[153:155],split="\t")[[c2]][5],
                                                         strsplit(output[153:155],split="\t")[[c3]][5])); 
        # Latent GOLD analysis:  Proportional ML;
        shell(cmd=paste(LatentGOLD.path," ",working.path,"ProportionalML.lgs /b",sep=""));
        output <- readLines(con=paste(working.path,"ProportionalML.lst",sep=""));
        stopifnot(output[152]=="Cluster\t0\ts.e.\t1\ts.e.\t \t \t");
        distal.prob.by.class.proportional.ML <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][4],
                                                             strsplit(output[153:155],split="\t")[[c2]][4],
                                                             strsplit(output[153:155],split="\t")[[c3]][4]));
        distal.prob.by.class.proportional.ML.se <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][5],
                                                                strsplit(output[153:155],split="\t")[[c2]][5],
                                                                strsplit(output[153:155],split="\t")[[c3]][5]));
        # Latent GOLD analysis:  Modal BCH;
        shell(cmd=paste(LatentGOLD.path," ",working.path,"ModalBCH.lgs /b",sep=""));
        output <- readLines(con=paste(working.path,"ModalBCH.lst",sep=""));
        stopifnot(output[152]=="Cluster\t0\ts.e.\t1\ts.e.\t \t \t");
        distal.prob.by.class.modal.BCH <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][4],
                                                       strsplit(output[153:155],split="\t")[[c2]][4],
                                                       strsplit(output[153:155],split="\t")[[c3]][4]));
        distal.prob.by.class.modal.BCH.se <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][5],
                                                          strsplit(output[153:155],split="\t")[[c2]][5],
                                                          strsplit(output[153:155],split="\t")[[c3]][5])); 
        # Latent GOLD analysis:  Proportional BCH;
        shell(cmd=paste(LatentGOLD.path," ",working.path,"ProportionalBCH.lgs /b",sep=""));
        output <- readLines(con=paste(working.path,"ProportionalBCH.lst",sep=""));
        stopifnot(output[152]=="Cluster\t0\ts.e.\t1\ts.e.\t \t \t");
        distal.prob.by.class.proportional.BCH <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][4],
                                                              strsplit(output[153:155],split="\t")[[c2]][4],
                                                              strsplit(output[153:155],split="\t")[[c3]][4]));
        distal.prob.by.class.proportional.BCH.se <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][5],
                                                                 strsplit(output[153:155],split="\t")[[c2]][5],
                                                                 strsplit(output[153:155],split="\t")[[c3]][5]));
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
        stopifnot(output[204]=="Profile");
        temp <- strsplit(output[c(209,211,213,215,217)],2,split="\t");
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
        stopifnot(output[152]=="Cluster\t0\ts.e.\t1\ts.e.\t \t \t");
        distal.prob.by.class.inclusive.modal.none <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][4],
                                                                  strsplit(output[153:155],split="\t")[[c2]][4],
                                                                  strsplit(output[153:155],split="\t")[[c3]][4]));
        distal.prob.by.class.inclusive.modal.none.se <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][5],
                                                                     strsplit(output[153:155],split="\t")[[c2]][5],
                                                                     strsplit(output[153:155],split="\t")[[c3]][5])); 
        # Latent GOLD analysis:  Inclusive Proportional None
        shell(cmd=paste(LatentGOLD.path," ",working.path,"InclusiveProportionalNone.lgs /b",sep=""));
        output <- readLines(con=paste(working.path,"InclusiveProportionalNone.lst",sep=""));
        stopifnot(output[152]=="Cluster\t0\ts.e.\t1\ts.e.\t \t \t");
        distal.prob.by.class.inclusive.proportional.none <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][4],
                                                                         strsplit(output[153:155],split="\t")[[c2]][4],
                                                                         strsplit(output[153:155],split="\t")[[c3]][4]));
        distal.prob.by.class.inclusive.proportional.none.se <- as.numeric(c(strsplit(output[153:155],split="\t")[[c1]][5],
                                                                            strsplit(output[153:155],split="\t")[[c2]][5],
                                                                            strsplit(output[153:155],split="\t")[[c3]][5]));
        these.preliminary.answers <- c(simulation = simulation,
                                       Item1Oracle = oracle.item.means[1,],
                                       Item2Oracle = oracle.item.means[2,],
                                       Item3Oracle = oracle.item.means[3,],
                                       Item4Oracle = oracle.item.means[4,],
                                       Item5Oracle = oracle.item.means[5,],
                                       DistalByClassOracle = oracle.distal.prob.by.class, 
                                       DistalByClassOracleSE = oracle.distal.prob.by.class.se, 
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
                                       DistalByClassModalNone = distal.prob.by.class.modal.none, 
                                       DistalByClassModalNoneSE = distal.prob.by.class.modal.none.se, 
                                       DistalByClassProportionalNone = distal.prob.by.class.proportional.none, 
                                       DistalByClassProportionalNoneSE = distal.prob.by.class.proportional.none.se, 
                                       DistalByClassModalML = distal.prob.by.class.modal.ML, 
                                       DistalByClassModalMLSE = distal.prob.by.class.modal.ML.se, 
                                       DistalByClassProportionalML = distal.prob.by.class.proportional.ML, 
                                       DistalByClassProportionalMLSE = distal.prob.by.class.proportional.ML.se, 
                                       DistalByClassModalBCH = distal.prob.by.class.modal.BCH,
                                       DistalByClassModalBCHSE = distal.prob.by.class.modal.BCH.se, 
                                       DistalByClassProportionalBCH = distal.prob.by.class.proportional.BCH, 
                                       DistalByClassProportionalBCHSE = distal.prob.by.class.proportional.BCH.se,
                                       Item1Inclusive = sorted.item.means.inclusive[1,],
                                       Item2Inclusive = sorted.item.means.inclusive[2,],
                                       Item3Inclusive = sorted.item.means.inclusive[3,],
                                       Item4Inclusive = sorted.item.means.inclusive[4,],
                                       Item5Inclusive = sorted.item.means.inclusive[5,],
                                       DistalByClassModalInclusive = distal.prob.by.class.inclusive.modal.none, 
                                       DistalByClassModalInclusiveSE = distal.prob.by.class.inclusive.modal.none.se, 
                                       DistalByClassProportionalInclusive = distal.prob.by.class.inclusive.proportional.none, 
                                       DistalByClassProportionalInclusiveSE = distal.prob.by.class.inclusive.proportional.none.se );
        preliminary.answers <- rbind(preliminary.answers,these.preliminary.answers);
        all.preliminary.answers <- rbind(all.preliminary.answers,preliminary.answers);
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
                                 BiasOracle1 = mean(preliminary.answers[,"DistalByClassOracle1"] - true.distal.prob.by.class[1]),
                                 BiasModalNone1 = mean(preliminary.answers[,"DistalByClassModalNone1"] - true.distal.prob.by.class[1]),
                                 BiasModalML1 = mean(preliminary.answers[,"DistalByClassModalML1"] - true.distal.prob.by.class[1]),
                                 BiasModalBCH1 = mean(preliminary.answers[,"DistalByClassModalBCH1"] - true.distal.prob.by.class[1]),
                                 BiasModalInclusive1 = mean(preliminary.answers[,"DistalByClassModalInclusive1"] - true.distal.prob.by.class[1]), 
                                 BiasProportionalNone1 = mean(preliminary.answers[,"DistalByClassProportionalNone1"] - true.distal.prob.by.class[1]),
                                 BiasProportionalML1 = mean(preliminary.answers[,"DistalByClassProportionalML1"] - true.distal.prob.by.class[1]),
                                 BiasProportionalBCH1 = mean(preliminary.answers[,"DistalByClassProportionalBCH1"] - true.distal.prob.by.class[1]),
                                 BiasProportionalInclusive1 = mean(preliminary.answers[,"DistalByClassProportionalInclusive1"] - true.distal.prob.by.class[1]), 
                                 BiasOracle2 = mean(preliminary.answers[,"DistalByClassOracle2"] - true.distal.prob.by.class[2]),
                                 BiasModalNone2 = mean(preliminary.answers[,"DistalByClassModalNone2"] - true.distal.prob.by.class[2]),
                                 BiasModalML2 = mean(preliminary.answers[,"DistalByClassModalML2"] - true.distal.prob.by.class[2]),
                                 BiasModalBCH2 = mean(preliminary.answers[,"DistalByClassModalBCH2"] - true.distal.prob.by.class[2]),
                                 BiasModalInclusive2 = mean(preliminary.answers[,"DistalByClassModalInclusive2"] - true.distal.prob.by.class[2]), 
                                 BiasProportionalNone2 = mean(preliminary.answers[,"DistalByClassProportionalNone2"] - true.distal.prob.by.class[2]),
                                 BiasProportionalML2 = mean(preliminary.answers[,"DistalByClassProportionalML2"] - true.distal.prob.by.class[2]),
                                 BiasProportionalBCH2 = mean(preliminary.answers[,"DistalByClassProportionalBCH2"] - true.distal.prob.by.class[2]),
                                 BiasProportionalInclusive2 = mean(preliminary.answers[,"DistalByClassProportionalInclusive2"] - true.distal.prob.by.class[2]), 
                                 BiasOracle3 = mean(preliminary.answers[,"DistalByClassOracle3"] - true.distal.prob.by.class[3]),
                                 BiasModalNone3 = mean(preliminary.answers[,"DistalByClassModalNone3"] - true.distal.prob.by.class[3]),
                                 BiasModalML3 = mean(preliminary.answers[,"DistalByClassModalML3"] - true.distal.prob.by.class[3]),
                                 BiasModalBCH3 = mean(preliminary.answers[,"DistalByClassModalBCH3"] - true.distal.prob.by.class[3]),
                                 BiasModalInclusive3 = mean(preliminary.answers[,"DistalByClassModalInclusive3"] - true.distal.prob.by.class[3]), 
                                 BiasProportionalNone3 = mean(preliminary.answers[,"DistalByClassProportionalNone3"] - true.distal.prob.by.class[3]),
                                 BiasProportionalML3 = mean(preliminary.answers[,"DistalByClassProportionalML3"] - true.distal.prob.by.class[3]),
                                 BiasProportionalBCH3 = mean(preliminary.answers[,"DistalByClassProportionalBCH3"] - true.distal.prob.by.class[3]),
                                 BiasProportionalInclusive3 = mean(preliminary.answers[,"DistalByClassProportionalInclusive3"] - true.distal.prob.by.class[3]), 
                                 RootMSEOracle1 = sqrt(mean((preliminary.answers[,"DistalByClassOracle1"] - true.distal.prob.by.class[1])^2)),
                                 RootMSEModalNone1 = sqrt(mean((preliminary.answers[,"DistalByClassModalNone1"] - true.distal.prob.by.class[1])^2)),
                                 RootMSEModalML1 = sqrt(mean((preliminary.answers[,"DistalByClassModalML1"] - true.distal.prob.by.class[1])^2)),
                                 RootMSEModalBCH1 = sqrt(mean((preliminary.answers[,"DistalByClassModalBCH1"] - true.distal.prob.by.class[1])^2)),
                                 RootMSEModalInclusive1 = sqrt(mean((preliminary.answers[,"DistalByClassModalInclusive1"] - true.distal.prob.by.class[1])^2)), 
                                 RootMSEProportionalNone1 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalNone1"] - true.distal.prob.by.class[1])^2)),
                                 RootMSEProportionalML1 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalML1"] - true.distal.prob.by.class[1])^2)),
                                 RootMSEProportionalBCH1 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalBCH1"] - true.distal.prob.by.class[1])^2)),
                                 RootMSEProportionalInclusive1 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalInclusive1"] - true.distal.prob.by.class[1])^2)), 
                                 RootMSEOracle2 = sqrt(mean((preliminary.answers[,"DistalByClassOracle2"] - true.distal.prob.by.class[2])^2)),
                                 RootMSEModalNone2 = sqrt(mean((preliminary.answers[,"DistalByClassModalNone2"] - true.distal.prob.by.class[2])^2)),
                                 RootMSEModalML2 = sqrt(mean((preliminary.answers[,"DistalByClassModalML2"] - true.distal.prob.by.class[2])^2)),
                                 RootMSEModalBCH2 = sqrt(mean((preliminary.answers[,"DistalByClassModalBCH2"] - true.distal.prob.by.class[2])^2)),
                                 RootMSEModalInclusive2 = sqrt(mean((preliminary.answers[,"DistalByClassModalInclusive2"] - true.distal.prob.by.class[2])^2)), 
                                 RootMSEProportionalNone2 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalNone2"] - true.distal.prob.by.class[2])^2)),
                                 RootMSEProportionalML2 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalML2"] - true.distal.prob.by.class[2])^2)),
                                 RootMSEProportionalBCH2 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalBCH2"] - true.distal.prob.by.class[2])^2)),
                                 RootMSEProportionalInclusive2 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalInclusive2"] - true.distal.prob.by.class[2])^2)), 
                                 RootMSEOracle3 = sqrt(mean((preliminary.answers[,"DistalByClassOracle3"] - true.distal.prob.by.class[3])^2)),
                                 RootMSEModalNone3 = sqrt(mean((preliminary.answers[,"DistalByClassModalNone3"] - true.distal.prob.by.class[3])^2)),
                                 RootMSEModalML3 = sqrt(mean((preliminary.answers[,"DistalByClassModalML3"] - true.distal.prob.by.class[3])^2)),
                                 RootMSEModalBCH3 = sqrt(mean((preliminary.answers[,"DistalByClassModalBCH3"] - true.distal.prob.by.class[3])^2)),
                                 RootMSEModalInclusive3 = sqrt(mean((preliminary.answers[,"DistalByClassModalInclusive3"] - true.distal.prob.by.class[3])^2)), 
                                 RootMSEProportionalNone3 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalNone3"] - true.distal.prob.by.class[3])^2)),
                                 RootMSEProportionalML3 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalML3"] - true.distal.prob.by.class[3])^2)),
                                 RootMSEProportionalBCH3 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalBCH3"] - true.distal.prob.by.class[3])^2)),
                                 RootMSEProportionalInclusive3 = sqrt(mean((preliminary.answers[,"DistalByClassProportionalInclusive3"] - true.distal.prob.by.class[3])^2)), 
                                 CoverageOracle1 = mean(abs(preliminary.answers[,"DistalByClassOracle1"]-true.distal.prob.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassOracleSE1"]),
                                 CoverageModalNone1 = mean(abs(preliminary.answers[,"DistalByClassModalNone1"]-true.distal.prob.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassModalNoneSE1"]),
                                 CoverageModalML1 = mean(abs(preliminary.answers[,"DistalByClassModalML1"]-true.distal.prob.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassModalMLSE1"]),
                                 CoverageModalBCH1 = mean(abs(preliminary.answers[,"DistalByClassModalBCH1"]-true.distal.prob.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassModalBCHSE1"]),
                                 CoverageModalInclusive1 = mean(abs(preliminary.answers[,"DistalByClassModalInclusive1"]-true.distal.prob.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassModalInclusiveSE1"]),
                                 CoverageProportionalNone1 = mean(abs(preliminary.answers[,"DistalByClassProportionalNone1"]-true.distal.prob.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassProportionalNoneSE1"]),
                                 CoverageProportionalML1 = mean(abs(preliminary.answers[,"DistalByClassProportionalML1"]-true.distal.prob.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassProportionalMLSE1"]),
                                 CoverageProportionalBCH1 = mean(abs(preliminary.answers[,"DistalByClassProportionalBCH1"]-true.distal.prob.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassProportionalBCHSE1"]),
                                 CoverageProportionalInclusive1 = mean(abs(preliminary.answers[,"DistalByClassProportionalInclusive1"]-true.distal.prob.by.class[1]) < 1.96*preliminary.answers[,"DistalByClassProportionalInclusiveSE1"]),
                                 CoverageOracle2 = mean(abs(preliminary.answers[,"DistalByClassOracle2"]-true.distal.prob.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassOracleSE2"]),
                                 CoverageModalNone2 = mean(abs(preliminary.answers[,"DistalByClassModalNone2"]-true.distal.prob.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassModalNoneSE2"]),
                                 CoverageModalML2 = mean(abs(preliminary.answers[,"DistalByClassModalML2"]-true.distal.prob.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassModalMLSE2"]),
                                 CoverageModalBCH2 = mean(abs(preliminary.answers[,"DistalByClassModalBCH2"]-true.distal.prob.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassModalBCHSE2"]),
                                 CoverageModalInclusive2 = mean(abs(preliminary.answers[,"DistalByClassModalInclusive2"]-true.distal.prob.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassModalInclusiveSE2"]),
                                 CoverageProportionalNone2 = mean(abs(preliminary.answers[,"DistalByClassProportionalNone2"]-true.distal.prob.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassProportionalNoneSE2"]),
                                 CoverageProportionalML2 = mean(abs(preliminary.answers[,"DistalByClassProportionalML2"]-true.distal.prob.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassProportionalMLSE2"]),
                                 CoverageProportionalBCH2 = mean(abs(preliminary.answers[,"DistalByClassProportionalBCH2"]-true.distal.prob.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassProportionalBCHSE2"]),
                                 CoverageProportionalInclusive2 = mean(abs(preliminary.answers[,"DistalByClassProportionalInclusive2"]-true.distal.prob.by.class[2]) < 1.96*preliminary.answers[,"DistalByClassProportionalInclusiveSE2"]),
                                 CoverageOracle3 = mean(abs(preliminary.answers[,"DistalByClassOracle3"]-true.distal.prob.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassOracleSE3"]),
                                 CoverageModalNone3 = mean(abs(preliminary.answers[,"DistalByClassModalNone3"]-true.distal.prob.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassModalNoneSE3"]),
                                 CoverageModalML3 = mean(abs(preliminary.answers[,"DistalByClassModalML3"]-true.distal.prob.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassModalMLSE3"]),
                                 CoverageModalBCH3 = mean(abs(preliminary.answers[,"DistalByClassModalBCH3"]-true.distal.prob.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassModalBCHSE3"]),
                                 CoverageModalInclusive3 = mean(abs(preliminary.answers[,"DistalByClassModalInclusive3"]-true.distal.prob.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassModalInclusiveSE3"]),
                                 CoverageProportionalNone3 = mean(abs(preliminary.answers[,"DistalByClassProportionalNone3"]-true.distal.prob.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassProportionalNoneSE3"]),
                                 CoverageProportionalML3 = mean(abs(preliminary.answers[,"DistalByClassProportionalML3"]-true.distal.prob.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassProportionalMLSE3"]),
                                 CoverageProportionalBCH3 = mean(abs(preliminary.answers[,"DistalByClassProportionalBCH3"]-true.distal.prob.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassProportionalBCHSE3"]),
                                 CoverageProportionalInclusive3 = mean(abs(preliminary.answers[,"DistalByClassProportionalInclusive3"]-true.distal.prob.by.class[3]) < 1.96*preliminary.answers[,"DistalByClassProportionalInclusiveSE3"]));
      rm(preliminary.answers);
      final.answers <- rbind(final.answers,these.final.answers);
    }
  }
}
finish.time <- Sys.time();
distal.prob.estimates <- all.preliminary.answers[,c("DistalByClassOracle1",
                                                    "DistalByClassOracle2",
                                                    "DistalByClassOracle3",
                                                    "DistalByClassOracleSE1",
                                                    "DistalByClassOracleSE2",
                                                    "DistalByClassOracleSE3",
                                                    "DistalByClassModalNone2",
                                                    "DistalByClassModalNone3",
                                                    "DistalByClassModalNoneSE1",
                                                    "DistalByClassModalNoneSE2",
                                                    "DistalByClassModalNoneSE3",
                                                    "DistalByClassProportionalNone1",
                                                    "DistalByClassProportionalNone2",
                                                    "DistalByClassProportionalNone3",
                                                    "DistalByClassProportionalNoneSE1",
                                                    "DistalByClassProportionalNoneSE2",
                                                    "DistalByClassProportionalNoneSE3",
                                                    "DistalByClassModalML1",
                                                    "DistalByClassModalML2",
                                                    "DistalByClassModalML3",
                                                    "DistalByClassModalMLSE1",
                                                    "DistalByClassModalMLSE2",
                                                    "DistalByClassModalMLSE3",
                                                    "DistalByClassProportionalML1",
                                                    "DistalByClassProportionalML2",
                                                    "DistalByClassProportionalML3",
                                                    "DistalByClassProportionalMLSE1",
                                                    "DistalByClassProportionalMLSE2",
                                                    "DistalByClassProportionalMLSE3",
                                                    "DistalByClassModalBCH1",
                                                    "DistalByClassModalBCH2",
                                                    "DistalByClassModalBCH3",
                                                    "DistalByClassModalBCHSE1",
                                                    "DistalByClassModalBCHSE2",
                                                    "DistalByClassModalBCHSE3",
                                                    "DistalByClassProportionalBCH1",
                                                    "DistalByClassProportionalBCH2",
                                                    "DistalByClassProportionalBCH3",
                                                    "DistalByClassProportionalBCHSE1",
                                                    "DistalByClassProportionalBCHSE2",
                                                    "DistalByClassProportionalBCHSE3")];
final.answers <- data.frame(final.answers);  
print(apply((distal.prob.estimates<=0)|(distal.prob.estimates>=1),2,sum)); 
print(final.answers);
print(apply(final.answers,2,mean));
write.table(x=final.answers,file=paste("simanswers",nsim,".txt",sep=""),
                                       quote=FALSE,row.names=TRUE,col.names=TRUE);
save.image(paste("FinishedSimulations-",nsim,".rdata",sep=""));