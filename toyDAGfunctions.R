## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Sample DAGs with BiDAG
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
sampleDAGs <- function(inData, nDigraphs = 50, seed=101, dname="", ...){
  n <- ncol(inData) # number of variables
  stepsave <- 10*round(n*n*log(n)) ## thinning interval, 
  ## to improve the properties of the chain (sample independence)
  iterations <- nDigraphs*stepsave ## nDAGs is the desired size of the DAGs ensemble
  
  if(!dir.exists("./saveout")) {dir.create("./saveout")} ## create output folder if it doesn't exist
  if(!file.exists(paste0("./saveout/dagdraw", n, "seed", seed, dname, ".RData"))){
    # initialise the score object needed for later functions
    scoreObject <- scoreparameters(data=inData, ..., 
                                   nodeslabels = colnames(inData))
    set.seed(seed) 
    ## set the seed for the generation of random numbers (for reproducibility)
    
    # find the search space with iterative search
    itFit <- iterativeMCMC(scoreObject, 
                           startorder = sample(c(1:n)), 
                           scoreout = TRUE) ## find iterative search space
    searchSpace <- itFit$endspace
    
    # sample a starting DAG with order MCMC
    # the default length of the chain is 6*n^2/log(n)
    # orderSample$score is the score for the highest scoring DAG found
    orderSample <- orderMCMC(scoreObject, 
                             MAP = FALSE, 
                             startspace = searchSpace, 
                             startorder = sample(c(1:n)), 
                             chainout = TRUE)
    startDAG <- last(orderSample$traceadd$incidence) # extract selected (last) DAG
    startDAGscore <- last(orderSample$trace) # extract score for the selected DAG
    
    # sample an ensemble of nDAGs DAGs from the posterior by using partition MCMC
    partitionSample <- partitionMCMC(scoreObject, 
                                     startspace = searchSpace, 
                                     startDAG = startDAG, 
                                     iterations = iterations, 
                                     stepsave = stepsave)
    
    ## extract the collection of sampled DAGs
    sampledDAGs <- partitionSample$traceadd$incidence 
    DAGscores <-partitionSample$trace # extract their scores
    
    ## save the the collection of sampled DAGs and their scores to an .RData file
    save(sampledDAGs, DAGscores, scoreObject,
         file=paste0("./saveout/dagdraw", n, "seed", seed, dname, ".RData"))
  }
}  


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Compute intervention effects with Bestie
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
computeEffects <- function(n, seed=101, dname=""){
  
  if (!dir.exists("./saveout")) {dir.create("./saveout")} ## create output folder if it doesn't exist
  if(file.exists(paste0("./saveout/dagdraw", n, "seed", seed, dname, ".RData"))){
    # load the collection of sampled DAGs
    load(file = paste0("./saveout/dagdraw", n, "seed", seed, dname, ".RData"))
  
    # sample parameters and derive a matrix of intervention effects for each DAG, 
    # saving them all in a list;
    # First check if effects have already been estimated and saved to file
    if(!file.exists(paste0("./saveout/effects", n, "seed", seed, dname, ".RData"))){
      causalMats <- DAGintervention(sampledDAGs, scoreObject, sample=TRUE)

      # save estimated intervention effects to an .RData file
      save(causalMats,
           file=paste0("./saveout/effects", n, "seed", seed, dname, ".RData"))
    }
  } else(warning("Looking for causal effects without a DAG? Try runif()!"))
}
## There appear to be no sampled DAGs stored, you muppet! :)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Load collection of sampled DAGs and effects
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
loadsamples <- function(seeds, nn, dname = "") {
  nSeeds <- length(seeds)
  alldigraphs <- vector("list", nDAGs * nSeeds) # to store the graphs
  alleffs <- vector("list", nDAGs * nSeeds) # to store the matrices of effects
  for (nlevel in 1:nSeeds) {
    seednumber <- seeds[nlevel]
    ## Retrieve sampled DAGs - DAG chain
    load(file = paste0("./saveout/dagdraw", nn, "seed", seednumber, dname, ".RData"))
    alldigraphs[1:nDAGs + (nlevel - 1) * nDAGs] <- sampledDAGs[-1] # remove the starting point
    ## Retrieve estimated effects
    load(file = paste0("./saveout/effects", nn, "seed", seednumber, dname, ".RData"))
    alleffs[1:nDAGs + (nlevel - 1) * nDAGs] <- causalMats[-1] # remove the starting point
  }
  list(alldigraphs=alldigraphs, alleffs=alleffs)
}

# Define a function for summarising estimates of DAGs and effects from the structure learning procedure.
# No longer needed, integrated within the respective plotting functions
# combineEsts <- function(allmatrices, eff=FALSE) round(simplify2array(allmatrices), 8)
# sdcausalmat <- apply(allarray, c(1,2), sd)  ## no longer used
# dags4plot <- combineEsts(allmatrices = allDAGs)
# effects4plot <- combineEsts(allmatrices = allEffects, eff = TRUE)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Internal function to set-up plot style
## Define some default plotting settings for the effects
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
setPlot <- function(dt = NULL, ...){
  plot(dt, 
       axes=FALSE, 
       frame.plot=FALSE, 
       xlab="", 
       ylab="", 
       main="",
       ...)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Plot a summary DAG from a collection of samples
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
dagviz <- function(dags4plot, # roundprobs,
                   daglabs = colnames(dags4plot[[1]]), # character vector of names
                   nn = length(daglabs),
                   grouped_nodes = NULL, # list(sample(order(daglabs), 2)),
                   # each element in the above list 
                   # will have its nodes placed inside a group
                   ## It won't work with fewer than 2 nodes
                   node_colour = "#1e90ff", # node colour
                   font_colour = "#fffaf0", # font colour in nodes
                   edge_colour = "#68228b", # edge colour
                   edge_threshold = 0.1, # only show edges above this threshold
                   group_colour = "#add8e6", # colour for group of nodes
                   edge_width = 1.5, # edge width
                   title_text = "Directed Acyclic Graph\n",
                   # title text (use \n for line breaks)
                   style_mat = matrix(c(1,2), nrow=2*nn + 1, ncol=2*nn)[1:nn, 1:nn],
                   # checker-board for illustration
                   edge_styles = c("solid", "dashed")
                   ){
  ## calculate plotting parameters from the DAG data
  dagsarray <- round(simplify2array(lapply(dags4plot, as.matrix)), 8)
  edgeprobs <- apply(dagsarray, c(1,2), mean)
  roundprobs <- round(255*edgeprobs, 0) # rounded to 0-255 for line intensity
  
  ## Build graphviz code
  graphcommand <- paste0("digraph G { \n node [color=\"", 
                         node_colour, 
                         "\", style=filled, fontcolor=\"", 
                         font_colour, "\"]; \n ")
  if (length(grouped_nodes) > 0) {
    for (kk in 1:length(grouped_nodes)) {
      node_text1 <- paste0(daglabs[grouped_nodes[[kk]]], 
                           ";", collapse = " ")
      node_text2 <- paste0(daglabs[grouped_nodes[[kk]]], 
                           collapse = " ")
      graphcommand <- paste0(graphcommand,
                           "subgraph cluster", 
                           kk, 
                           " { \n bgcolor=\"", 
                           group_colour, 
                           "\" \n penwidth=0 \n", 
                           node_text1,
                           " \n {rank=same ", 
                           node_text2,"} \n } \n ") 
      }
    }
  for(ii in 1:nn){
    for(jj in 1:nn){
    if(edgeprobs[ii,jj] > edge_threshold){
      graphcommand<-paste0(graphcommand, 
                          daglabs[ii],
                          " -> ", 
                          daglabs[jj], 
                          " [color=\"", 
                          edge_colour, 
                          as.hexmode(roundprobs[ii,jj]),
                          "\", penwidth=", 
                          edge_width, 
                          ", style=", 
                          edge_styles[style_mat[ii, jj]],"] \n")
      }
    }
  }
  grViz(paste0(graphcommand,"labelloc=\"t\";\n label= \"", title_text, "\";}"))
}

## filename for gv file if desired, not needed with DiagrammeR
## gv_output <- "DAG.gv"
## write(graphcommand, file=gv_output)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## grViz() by itself produces a html output not visible in pdfs (even by setting always_allow_html: true in the YAML)
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Save summary DAG to a file and display on pdf
## ---------------------------------------------------------------------------------------------------------------------------------------------------------

displayDAG <- function(g2plot, figname = "dataDAG.png"){
  g2plot %>%
    export_svg %>%
    charToRaw %>%
    rsvg_png(figname, width = 12000, height = 12000)
  ## Increasing the size also increases the quality of the graph
  knitr::include_graphics(figname)
}

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## Plot a matrix of posterior distributions of effects
## ---------------------------------------------------------------------------------------------------------------------------------------------------------
plotEffects <- function(effects4plot, 
                        sortlabs, ## an order vector of indiced
                        ## would it be better if names are saved in the list of effect matrices?
                        ## then use match if given a vector of order labels to find the indices in the Effect matrices
                        sigcutoff = 0.025, 
                        title_text = "Distributions of Downstream Causal Effects\n"){
  orderingy = c(0, sortlabs) 
  ## Provide labels in sortlabs in a customised order if desired
  ### Pick an ordering which matches to some extent the node hierarchy on the DAG
  
  nn <- length(sortlabs)
  effsarray <- round(simplify2array(effects4plot), 8)
  aveffs <- apply(effsarray, c(1,2), mean) 
  # estimate of the posterior mean of the intervention effects
  roundeffs <- round(aveffs, 3)
  
  par(mfrow = c(nn+1, nn+1))
  par(oma = rep(0.2,4) + c(0,0,5,0))
  par(mar = rep(0,4))
  
  for(i in 0:nn+1){
    ii <- orderingy[i+1]
    for(j in 0:nn){
      jj <- orderingy[j+1]
      if(i==(nn+1) || j==0){
        setPlot(xlim=c(-1.1,1.1), ylim=c(0,1), col="darkorchid4")
        if(j==0 && i<(nn+1)){
          text(0,0.5, sortlabs[ii], cex=.6)
          }
        if(i==(nn+1) && j>0){
          text(0,0.5, sortlabs[jj], srt=90, cex=.6)
          }
        } else{
          if(!(roundeffs[ii,jj]==0 || roundeffs[ii,jj]==1) ){
            d <- density(effsarray[ii,jj,])
            setPlot(d, xlim=c(-0.1,0.6), col="dodgerblue")
            testquantiles <- quantile(effsarray[ii,jj,], c(sigcutoff, 1-sigcutoff))
            if(testquantiles[1]>0 || testquantiles[2]<0){
              u <- par("usr") # The coordinates of the plot area
              rect(u[1], u[3], u[2], u[4], col="#ffa07a44", border=NA)
              }
            polygon(d, col="dodgerblue", border="dodgerblue")
            abline(v=0, col="firebrick3", lty=2)
            u <- par("usr") # The coordinates of the plot area
            text( (u[1]+u[2])/2,(u[3]+u[4])/2, 
                  format(round(roundeffs[ii,jj], 2), nsmall = 2) )
            } else {
              setPlot(xlim=c(-0.1,0.6), ylim=c(0,1), col="darkorchid4")
              text(0, 0.5, roundeffs[ii,jj])
            }
        }
    }
  }
  mtext(title_text, outer = TRUE, cex = 0.9)
}
