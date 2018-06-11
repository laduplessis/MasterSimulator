rm(list = ls())
library(rjson)
library(bdskytools)
library(yaml)

colmap <- list('I'=cblue, "R"=corange, 'E'=cpurple, 'S'=cgreen)
legmap <- list('I'='Infected', 'R'='Removed', 'E'='Exposed', 'S'='Susceptible')

models <- list("constantRe"="Constant", 
               "decreaseRe"="Piecewise\nconstant",
               "seir"      ="SEIR  \nmodel")#model <- "decreaseRe")
numtraces <- 25

model <- names(models)[[1]]

for (model in names(models)) {
    
    #####################################################################################
    # Read config file and simulation results
    pars  <- yaml.load_file(paste0("../config/",model,".cfg"))
    sims  <- fromJSON(file=paste0("../results/simulations/",model,".json"))
    
    # Subsample traces from results
    traces <- sample(1:length(sims$trajectories), numtraces, replace=FALSE)
    
    # Get x and y limits for plots
    maxtimes <- c()
    maxpops  <- c()
    for (i in traces) {cd 
        maxtimes <- c(maxtimes, max(sims$trajectories[[i]]$t))
        
        popmax <- c()
        types  <- c()
        for (type in names(sims$trajectories[[i]])) {
            if (type != 't' && type != 'R' && type != 'S') {
                types <- c(types, type)
                popmax <- c(popmax, max(sims$trajectories[[i]][[type]]))
            }
        }
        maxpops <- rbind(maxpops, popmax)
        rownames(maxpops)[nrow(maxpops)] <- i
        colnames(maxpops) <- types
    }
    maxT  <- max(maxtimes)
    maxP  <- max(maxpops)
    
    
    #####################################################################################
    # Plot trajectories
    par(mar=c(5,4,4,4)+0.1)
    plot(1, type='n', xlim=c(0,maxT), ylim=range(pretty(c(0,maxP))), xlab='Time', ylab='Population', bty='o', las=1, cex.lab=1.5, cex.axis=1.2)
    grid(col=pal.dark(cgray))
    
    for (i in traces) {
      traj  <- sims$trajectories[[i]]
      for (j in 1:length(names(traj))) {
          type <- names(traj)[j]
          if (type != 't' && type != 'Y' && type != 'R' && type != 'S')
             lines(sims$trajectories[[i]]$t, sims$trajectories[[i]][[type]], col=pal.dark(colmap[[type]],0.3), type='s', lwd=1)    
      }
      
      # Plot sampling events on trajectory of infenctious individuals
      samplingidxs  <- which(diff(traj$Y) == 1)+1
      samplingtimes <- traj$t[samplingidxs]
      points(samplingtimes, traj$I[samplingidxs], pch=20, col=pal.dark(cred,0.2), cex=0.8)
      # print(length(samplingidxs))
    }
    
    
    #####################################################################################
    # Plot Re
    par(new=TRUE)
    plot(1, type='n', xlim=c(0,maxT), ylim=c(0,5), bty='n', axes=FALSE, xlab="", ylab="")
    axis(4, las=1)
    mtext(expression("R"[e]),4, line=2, cex.lab=1.5)
    
    if ('popSize' %in% names(pars)) {
        for (i in traces) {
          traj          <- sims$trajectories[[i]]
          lines(traj$t, (pars$reproductiveNumber*traj$S)/(pars$popSize), col=pal.dark(cred,0.3), lwd=2, lty=2)
        }
    } else 
    if ('shiftTime' %in% names(pars)) {
        #abline(v=pars$shiftTime, lty=2, col=pal.dark(cred))
        lines(c(0, pars$shiftTime, maxT), c(pars$reproductiveNumber1, pars$reproductiveNumber2, pars$reproductiveNumber2), col=pal.dark(cred), lty=2, lwd=2, type='s')
    } else {
        lines(c(0, maxT), rep(pars$reproductiveNumber, 2), col=pal.dark(cred), lty=2, lwd=2, type='s')
    }
    
    
    #####################################################################################
    # Set up legend
    pops <- c()
    cols <- c()
    for (i in c("E","I")) {
        if (i %in% names(sims$trajectories[[1]])) {
            pops <- c(pops, legmap[[i]])
            cols <- c(cols, pal.dark(colmap[[i]]))
        }
    }
    legend('topright', legend=c(pops, "Samples", expression("R"[e])), col=c(cols, rep(pal.dark(cred),2)), lty=c(rep(1,length(pops)),0,2), pch=c(rep(NA,length(pops)),20,NA))
    
}

