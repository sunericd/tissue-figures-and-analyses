# generates simulated spatial transcriptomics using SRTsim GUI

library(SRTsim)


setwd("C:/Users/Eric Sun/Desktop/RESEARCH/GeneImputationTISSUE/SRTsim")


savename <- "main_two_conditions_same" # main_two_conditions_same, main_two_conditions_different


# run shiny app to specify simulation: https://xzhoulab.github.io/SRTsim/03_Reference_Free_Example/

shinySRT1 <- SRTsim_shiny()


# after exiting, run this to save:
simSRT1 <- Shiny2SRT(shinySRT1)
saveRDS(simSRT1, paste0(savename,".rds"))

# save counts matrix
write.csv(as.data.frame(as.matrix(simSRT1@simCounts)), paste0(savename,"_counts.txt"))

# save metadata
write.csv(simSRT1@simcolData, paste0(savename,"_meta.txt"))


# to regenerate from object
simSRT1 <- readRDS("main_two_conditions_different.rds")
newCT1  <- reGenCountshiny(simSRT1)
