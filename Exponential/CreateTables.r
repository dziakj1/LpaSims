rm(list=ls(all=TRUE));
load("Preliminary123.rdata");
temp.from.123 <- these.final.answers;
load("FinishedSimulations-1000.rdata");
final.answers <- rbind(final.answers[1:5,],temp.from.123,final.answers[6:11,]); # one condition was mistakenly not appended earlier;

########################################################
# What is the average estimate of entropy for each condition?
print(round(cbind(final.answers$which.measurement.quality,
              final.answers$which.class.size.distribution,
              final.answers$which.effect.size,
              final.answers$NoninclusiveEntropyRamaswamy),2));
#########################################################
# Summarize across classes

final.answers$CoverageOracle <- ((final.answers$CoverageOracle1)+(final.answers$CoverageOracle2)+(final.answers$CoverageOracle3))/3;
final.answers$CoverageModalNone <- ((final.answers$CoverageModalNone1)+(final.answers$CoverageModalNone2)+(final.answers$CoverageModalNone3))/3;
final.answers$CoverageProportionalNone <- ((final.answers$CoverageProportionalNone1)+(final.answers$CoverageProportionalNone2)+(final.answers$CoverageProportionalNone3))/3;
final.answers$CoverageModalML <- ((final.answers$CoverageModalML1)+(final.answers$CoverageModalML2)+(final.answers$CoverageModalML3))/3;
final.answers$CoverageProportionalML <- ((final.answers$CoverageProportionalML1)+(final.answers$CoverageProportionalML2)+(final.answers$CoverageProportionalML3))/3;
final.answers$CoverageModalBCH <- ((final.answers$CoverageModalBCH1)+(final.answers$CoverageModalBCH2)+(final.answers$CoverageModalBCH3))/3;
final.answers$CoverageProportionalBCH <- ((final.answers$CoverageProportionalBCH1)+(final.answers$CoverageProportionalBCH2)+(final.answers$CoverageProportionalBCH3))/3;
final.answers$CoverageModalInclusive <- ((final.answers$CoverageModalInclusive1)+(final.answers$CoverageModalInclusive2)+(final.answers$CoverageModalInclusive3))/3;
final.answers$CoverageProportionalInclusive <- ((final.answers$CoverageProportionalInclusive1)+(final.answers$CoverageProportionalInclusive2)+(final.answers$CoverageProportionalInclusive3))/3;
final.answers$CoverageModalQuadratic <- ((final.answers$CoverageModalQuadratic1)+(final.answers$CoverageModalQuadratic2)+(final.answers$CoverageModalQuadratic3))/3;
final.answers$CoverageProportionalQuadratic <- ((final.answers$CoverageProportionalQuadratic1)+(final.answers$CoverageProportionalQuadratic2)+(final.answers$CoverageProportionalQuadratic3))/3;

final.answers$AbsBiasOracle <- (abs(final.answers$BiasOracle1)+abs(final.answers$BiasOracle2)+abs(final.answers$BiasOracle3))/3;
final.answers$AbsBiasModalNone <- (abs(final.answers$BiasModalNone1)+abs(final.answers$BiasModalNone2)+abs(final.answers$BiasModalNone3))/3;
final.answers$AbsBiasProportionalNone <- (abs(final.answers$BiasProportionalNone1)+abs(final.answers$BiasProportionalNone2)+abs(final.answers$BiasProportionalNone3))/3;
final.answers$AbsBiasModalML <- (abs(final.answers$BiasModalML1)+abs(final.answers$BiasModalML2)+abs(final.answers$BiasModalML3))/3;
final.answers$AbsBiasProportionalML <- (abs(final.answers$BiasProportionalML1)+abs(final.answers$BiasProportionalML2)+abs(final.answers$BiasProportionalML3))/3;
final.answers$AbsBiasModalBCH <- (abs(final.answers$BiasModalBCH1)+abs(final.answers$BiasModalBCH2)+abs(final.answers$BiasModalBCH3))/3;
final.answers$AbsBiasProportionalBCH <- (abs(final.answers$BiasProportionalBCH1)+abs(final.answers$BiasProportionalBCH2)+abs(final.answers$BiasProportionalBCH3))/3;
final.answers$AbsBiasModalInclusive <- (abs(final.answers$BiasModalInclusive1)+abs(final.answers$BiasModalInclusive2)+abs(final.answers$BiasModalInclusive3))/3;
final.answers$AbsBiasProportionalInclusive <- (abs(final.answers$BiasProportionalInclusive1)+abs(final.answers$BiasProportionalInclusive2)+abs(final.answers$BiasProportionalInclusive3))/3;
final.answers$AbsBiasModalQuadratic <- (abs(final.answers$BiasModalQuadratic1)+abs(final.answers$BiasModalQuadratic2)+abs(final.answers$BiasModalQuadratic3))/3;
final.answers$AbsBiasProportionalQuadratic <- (abs(final.answers$BiasProportionalQuadratic1)+abs(final.answers$BiasProportionalQuadratic2)+abs(final.answers$BiasProportionalQuadratic3))/3;

final.answers$RootMSEOracle <- sqrt(((final.answers$RootMSEOracle1)^2+(final.answers$RootMSEOracle2)^2+(final.answers$RootMSEOracle3)^2)/3);
final.answers$RootMSEModalNone <- sqrt(((final.answers$RootMSEModalNone1)^2+(final.answers$RootMSEModalNone2)^2+(final.answers$RootMSEModalNone3)^2)/3);
final.answers$RootMSEProportionalNone <- sqrt(((final.answers$RootMSEProportionalNone1)^2+(final.answers$RootMSEProportionalNone2)^2+(final.answers$RootMSEProportionalNone3)^2)/3);
final.answers$RootMSEModalML <- sqrt(((final.answers$RootMSEModalML1)^2+(final.answers$RootMSEModalML2)^2+(final.answers$RootMSEModalML3)^2)/3);
final.answers$RootMSEProportionalML <- sqrt(((final.answers$RootMSEProportionalML1)^2+(final.answers$RootMSEProportionalML2)^2+(final.answers$RootMSEProportionalML3)^2)/3);
final.answers$RootMSEModalBCH <- sqrt(((final.answers$RootMSEModalBCH1^2)+(final.answers$RootMSEModalBCH2)^2+(final.answers$RootMSEModalBCH3)^2)/3);
final.answers$RootMSEProportionalBCH <- sqrt(((final.answers$RootMSEProportionalBCH1^2)+(final.answers$RootMSEProportionalBCH2)^2+(final.answers$RootMSEProportionalBCH3)^2)/3);
final.answers$RootMSEModalInclusive <- sqrt(((final.answers$RootMSEModalInclusive1)^2+(final.answers$RootMSEModalInclusive2)^2+(final.answers$RootMSEModalInclusive3)^2)/3);
final.answers$RootMSEProportionalInclusive <- sqrt(((final.answers$RootMSEProportionalInclusive1)^2+(final.answers$RootMSEProportionalInclusive2)^2+(final.answers$RootMSEProportionalInclusive3)^2)/3);
final.answers$RootMSEModalQuadratic <- sqrt(((final.answers$RootMSEModalQuadratic1)^2+(final.answers$RootMSEModalQuadratic2)^2+(final.answers$RootMSEModalQuadratic3)^2)/3);
final.answers$RootMSEProportionalQuadratic <- sqrt(((final.answers$RootMSEProportionalQuadratic1)^2+(final.answers$RootMSEProportionalQuadratic2)^2+(final.answers$RootMSEProportionalQuadratic3)^2)/3);



#########################################################
# Bias table:
temp <- rbind(final.answers$which.measurement.quality,
            final.answers$which.class.size.distribution,
            final.answers$which.effect.size,
            round(final.answers$AbsBiasModalNone,3),
            round(final.answers$AbsBiasProportionalNone,3),
            round(final.answers$AbsBiasModalML,3),
            round(final.answers$AbsBiasProportionalML,3),
            round(final.answers$AbsBiasModalBCH,3),
            round(final.answers$AbsBiasProportionalBCH,3),
            round(final.answers$AbsBiasModalInclusive,3),
            round(final.answers$AbsBiasProportionalInclusive,3),
            round(final.answers$AbsBiasModalQuadratic,3),
            round(final.answers$AbsBiasProportionalQuadratic,3),
            round(final.answers$AbsBiasOracle,3));
print(rbind(temp[,1:6],temp[,7:12]));

########################################################
# Zoom-in bias table;

subset.answers <- final.answers[which(final.answers$which.effect.size==3 & final.answers$which.class.size.distribution==1),];  # only consider large effects and even classes;


temp <- rbind(true.distal.mean.by.class,
              c(subset.answers$which.measurement.quality[1],subset.answers$which.class.size.distribution[1],subset.answers$which.effect.size[1]),
              c(subset.answers$BiasModalNone1[1],subset.answers$BiasModalNone2[1],subset.answers$BiasModalNone3[1]),
              c(subset.answers$BiasModalML1[1],subset.answers$BiasModalML2[1],subset.answers$BiasModalML3[1]),
              c(subset.answers$BiasModalBCH1[1],subset.answers$BiasModalBCH2[1],subset.answers$BiasModalBCH3[1]),
              c(subset.answers$BiasModalInclusive1[1],subset.answers$BiasModalInclusive2[1],subset.answers$BiasModalInclusive3[1]),
              c(subset.answers$BiasModalQuadratic1[1],subset.answers$BiasModalQuadratic2[1],subset.answers$BiasModalQuadratic3[1]),
              c(subset.answers$which.measurement.quality[2],subset.answers$which.class.size.distribution[2],subset.answers$which.effect.size[2]),
              c(subset.answers$BiasModalNone1[2],subset.answers$BiasModalNone2[2],subset.answers$BiasModalNone3[2]),
              c(subset.answers$BiasModalML1[2],subset.answers$BiasModalML2[2],subset.answers$BiasModalML3[2]),
              c(subset.answers$BiasModalBCH1[2],subset.answers$BiasModalBCH2[2],subset.answers$BiasModalBCH3[2]),
              c(subset.answers$BiasModalInclusive1[2],subset.answers$BiasModalInclusive2[2],subset.answers$BiasModalInclusive3[2]),
              c(subset.answers$BiasModalQuadratic1[2],subset.answers$BiasModalQuadratic2[2],subset.answers$BiasModalQuadratic3[2]) );

print(temp);



#########################################################
# Root MSE table:
temp <- rbind(final.answers$which.measurement.quality,
            final.answers$which.class.size.distribution,
            final.answers$which.effect.size,
            round(final.answers$RootMSEModalNone,3),
            round(final.answers$RootMSEProportionalNone,3),
            round(final.answers$RootMSEModalML,3),
            round(final.answers$RootMSEProportionalML,3),
            round(final.answers$RootMSEModalBCH,3),
            round(final.answers$RootMSEProportionalBCH,3),
            round(final.answers$RootMSEModalInclusive,3),
            round(final.answers$RootMSEProportionalInclusive,3),
            round(final.answers$RootMSEModalQuadratic,3),
            round(final.answers$RootMSEProportionalQuadratic,3),
            round(final.answers$RootMSEOracle,3));
print(rbind(temp[,1:6],temp[,7:12]));

#########################################################
# Coverage table:  
temp <- rbind(final.answers$which.measurement.quality,
              final.answers$which.class.size.distribution,
              final.answers$which.effect.size,
              round(final.answers$CoverageModalNone,3),
              round(final.answers$CoverageProportionalNone,3),
              round(final.answers$CoverageModalML,3),
              round(final.answers$CoverageProportionalML,3),
              round(final.answers$CoverageModalBCH,3),
              round(final.answers$CoverageProportionalBCH,3),
              round(final.answers$CoverageModalInclusive,3),
              round(final.answers$CoverageProportionalInclusive,3),
              round(final.answers$CoverageModalQuadratic,3),
              round(final.answers$CoverageProportionalQuadratic,3),
              round(final.answers$CoverageOracle,3));
print(rbind(temp[,1:6],temp[,7:12]));
