rm(list=ls(all=TRUE));
#########################################################
# Set initial settings:
set.seed(46000); 
N <- 1000;
nsim <- 500;
# Loop over conditions; 
for (class.size.distribution in c("even","uneven")) {
  for (effect.size in c("large","medium","small")) {
    print(paste(class.size.distribution,effect.size));
    output.path <- "C:\\Users\\jjd264\\Documents\\Sims-Lpa\\Datasets\\";
    #########################################################
    # Define parameters; 
    if (class.size.distribution=="even") {
      true.class.sizes <- c(1/3,1/3,1/3);
    }
    if (class.size.distribution=="uneven") {
      true.class.sizes <- c(.2,.3,.5);
    }
    if (effect.size == "large") {
      true.distal.prob.by.class <- c(.9,.3,.5);
    }
    if (effect.size == "medium") {
      true.distal.prob.by.class <- c(.7,.3,.5);
    }
    if (effect.size == "small") {
      true.distal.prob.by.class <- c(.5,.3,.4);
    }
    true.distal.prob <- true.distal.prob.by.class%*%true.class.sizes;
    contingency.table <- rbind(1-true.distal.prob.by.class,
                               true.distal.prob.by.class)
    
    temp <- chisq.test(contingency.table);
    o <- temp$observed;
    e <- temp$expected;
    w <- sqrt(sum(((o-e)^2)/e));
    print(c(true.distal.prob.by.class,true.distal.prob,w));
  } 
}
