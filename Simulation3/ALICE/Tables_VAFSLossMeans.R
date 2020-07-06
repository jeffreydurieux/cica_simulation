# Mon Jul  6 13:18:56 2020
# Author: Jeffrey Durieux, MSc


# table of VAFS/ loss and screes



vafs <- sapply(1:length(ListVAFDATA), function(i) ListVAFDATA[[i]]$VAF)

rep(1:20, )
vafs30 <- vafs[,seq(from = 1, to = 20, by = 2)]
vafs70 <- vafs[,seq(from = 2, to = 20, by = 2)]

m3 <- apply(vafs30, MARGIN = 1, mean)
sd3 <- apply(vafs30, MARGIN = 1, sd)

m7 <- apply(vafs70, MARGIN = 1, mean)
sd7 <- apply(vafs70, MARGIN = 1, sd)

cbind(comp = ListVAFDATA$rep_1CICA_simdata_0.30$n.comp,
      clus = ListVAFDATA$rep_1CICA_simdata_0.30$n.clus,
      m3, sd3)

cbind(comp = ListVAFDATA$rep_1CICA_simdata_0.30$n.comp,
      clus = ListVAFDATA$rep_1CICA_simdata_0.30$n.clus,
      m7, sd7)

