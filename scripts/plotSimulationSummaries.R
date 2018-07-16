rm(list = ls())
library(rjson)
library(bdskytools)
library(yaml)

colmap <- list('I'=cblue, "R"=corange, 'E'=cpurple, 'S'=cgreen, "Y"=cred)
legmap <- list('I'='Infected', 'R'='Removed', 'E'='Exposed', 'S'='Susceptible', 'Y'="Sampled")

models <- list("constantRe.summary"="Constant", 
               "decreaseRe.summary"="Piecewise\nconstant",
               "seir.summary"      ="SEIR  \nmodel",
               "sirs.summary"      ="SIRS  \nmodel")


par(mfrow=c(length(models),1))
for (model in names(models)) {
  
  #####################################################################################
  # Read config file and simulation results
  pars  <- yaml.load_file(paste0("../config/",model,".cfg"))
  sims  <- fromJSON(file=paste0("../results/summaries/",model,".json"))
  
  
  
  #####################################################################################
  # Plot trajectories
  par(mar=c(5,4,4,4)+0.1)
  plot(1, type='n', xlim=range(sims$t), ylim=range(pretty(c(0,5000))), xlab='Time', ylab='Population', bty='o', las=1, cex.lab=1.5, cex.axis=1.2)
  grid(col=pal.dark(cgray))
  
  for (type in c("S", "I","E","Y")) {
      if (type %in% names(sims)) {
      
          mean  <- sims[[type]]$mean[[1]]
          upper <- sims[[type]]$mean[[1]]+sims[[type]]$std[[1]]
          lower <- sapply(sims[[type]]$mean[[1]]-sims[[type]]$std[[1]], max, 0)
        
          lines(sims$t, mean, lwd=1, col=pal.dark(colmap[[type]]))
          polygon(c(sims$t, rev(sims$t)), c(upper, rev(lower)), col=pal.dark(colmap[[type]],0.3), 
                  border=pal.dark(colmap[[type]]), lty=2, lwd=1)
      }
  }
  
  
  #####################################################################################
  # Plot Re
  par(new=TRUE)
  plot(1, type='n', xlim=range(sims$t), ylim=c(0,5), bty='n', axes=FALSE, xlab="", ylab="")
  axis(4, las=1)
  mtext(expression("R"[e]),4, line=2, cex.lab=1.5)

  if ('popSize' %in% names(pars)) {
      susMean  <- sims$S$mean[[1]]
      susUpper <- sims$S$mean[[1]] + sims$S$std[[1]]
      susLower <- sapply(sims$S$mean[[1]]-sims$S$std[[1]], max, 0)
      lines(sims$t, (pars$reproductiveNumber*susMean)/(pars$popSize), col=pal.dark(cred,0.3), lwd=2, lty=1)
      lines(sims$t, (pars$reproductiveNumber*susUpper)/(pars$popSize), col=pal.dark(cred,0.3), lwd=2, lty=2)
      lines(sims$t, (pars$reproductiveNumber*susLower)/(pars$popSize), col=pal.dark(cred,0.3), lwd=2, lty=2)
  } else
    if ('shiftTime' %in% names(pars)) {
      #abline(v=pars$shiftTime, lty=2, col=pal.dark(cred))
      lines(c(0, pars$shiftTime, max(sims$t)), c(pars$reproductiveNumber1, pars$reproductiveNumber2, pars$reproductiveNumber2), col=pal.dark(cred), lty=2, lwd=2, type='s')
    } else {
      lines(range(sims$t), rep(pars$reproductiveNumber, 2), col=pal.dark(cred), lty=2, lwd=2, type='s')
    }


  #####################################################################################
  # Set up legend
  pops <- c()
  cols <- c()
  for (i in c("E","I","Y")) {
    if (i %in% names(sims)) {
      pops <- c(pops, legmap[[i]])
      cols <- c(cols, pal.dark(colmap[[i]]))
    }
  }
  legend('topright', legend=c(pops, expression("R"[e])), col=c(cols, pal.dark(cred)), lty=c(rep(1,length(pops)),2))
}

