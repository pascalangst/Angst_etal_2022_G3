# Oc pixi plot:

library(ggplot2)

# read pixi output and make one lineage negative for presentation
Oc_pixy_out_pi <- read.delim("~/bioinformatics/pixy/Oc_SNPs_only_pixy_out_pi.txt")
Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="EAU"] <- Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="EAU"] * (-1)

#Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="WEA" & Oc_pixy_out_pi$count_missing > quantile(Oc_pixy_out_pi$count_missing[Oc_pixy_out_pi$pop=="WEA"], .95)] <- NA
#Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="EAU" & Oc_pixy_out_pi$count_missing > quantile(Oc_pixy_out_pi$count_missing[Oc_pixy_out_pi$pop=="EAU"], .95)] <- NA
#Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="WEA" & Oc_pixy_out_pi$no_sites < quantile(Oc_pixy_out_pi$no_sites[Oc_pixy_out_pi$pop=="WEA"], .05)] <- NA
#Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="EAU" & Oc_pixy_out_pi$no_sites < quantile(Oc_pixy_out_pi$no_sites[Oc_pixy_out_pi$pop=="EAU"], .05)] <- NA

# remove first and last 5 windows
Oc_pixy_out_pi$pos <- c(paste(Oc_pixy_out_pi$chromosome, Oc_pixy_out_pi$window_pos_1, sep = "-"))
Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pos%in%do.call(rbind, lapply(split(Oc_pixy_out_pi,Oc_pixy_out_pi$chromosome), function(x) {return(x[which.max(x$window_pos_1),])}))$pos] <- NA
Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pos%in%do.call(rbind, lapply(split(Oc_pixy_out_pi,Oc_pixy_out_pi$chromosome), function(x) {return(x[which.min(x$window_pos_1),])}))$pos] <- NA
Oc_pixy_out_pi_copy <- Oc_pixy_out_pi
Oc_pixy_out_pi_copy <- Oc_pixy_out_pi_copy[ !Oc_pixy_out_pi_copy$pos %in% do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.min(x$window_pos_1),])}))$pos, ]
Oc_pixy_out_pi_copy <- Oc_pixy_out_pi_copy[ !Oc_pixy_out_pi_copy$pos %in% do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.max(x$window_pos_1),])}))$pos, ]
Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pos%in%do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.max(x$window_pos_1),])}))$pos] <- NA
Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pos%in%do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.min(x$window_pos_1),])}))$pos] <- NA
Oc_pixy_out_pi_copy <- Oc_pixy_out_pi_copy[ !Oc_pixy_out_pi_copy$pos %in% do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.min(x$window_pos_1),])}))$pos, ]
Oc_pixy_out_pi_copy <- Oc_pixy_out_pi_copy[ !Oc_pixy_out_pi_copy$pos %in% do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.max(x$window_pos_1),])}))$pos, ]
Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pos%in%do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.max(x$window_pos_1),])}))$pos] <- NA
Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pos%in%do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.min(x$window_pos_1),])}))$pos] <- NA
Oc_pixy_out_pi_copy <- Oc_pixy_out_pi_copy[ !Oc_pixy_out_pi_copy$pos %in% do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.min(x$window_pos_1),])}))$pos, ]
Oc_pixy_out_pi_copy <- Oc_pixy_out_pi_copy[ !Oc_pixy_out_pi_copy$pos %in% do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.max(x$window_pos_1),])}))$pos, ]
Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pos%in%do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.max(x$window_pos_1),])}))$pos] <- NA
Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pos%in%do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.min(x$window_pos_1),])}))$pos] <- NA
Oc_pixy_out_pi_copy <- Oc_pixy_out_pi_copy[ !Oc_pixy_out_pi_copy$pos %in% do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.min(x$window_pos_1),])}))$pos, ]
Oc_pixy_out_pi_copy <- Oc_pixy_out_pi_copy[ !Oc_pixy_out_pi_copy$pos %in% do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.max(x$window_pos_1),])}))$pos, ]
Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pos%in%do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.max(x$window_pos_1),])}))$pos] <- NA
Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pos%in%do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.min(x$window_pos_1),])}))$pos] <- NA
Oc_pixy_out_pi_copy <- Oc_pixy_out_pi_copy[ !Oc_pixy_out_pi_copy$pos %in% do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.min(x$window_pos_1),])}))$pos, ]
Oc_pixy_out_pi_copy <- Oc_pixy_out_pi_copy[ !Oc_pixy_out_pi_copy$pos %in% do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.max(x$window_pos_1),])}))$pos, ]
Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pos%in%do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.max(x$window_pos_1),])}))$pos] <- NA
Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pos%in%do.call(rbind, lapply(split(Oc_pixy_out_pi_copy,Oc_pixy_out_pi_copy$chromosome), function(x) {return(x[which.min(x$window_pos_1),])}))$pos] <- NA

#Oc_pixy_out_pi[4298,5] <- 0
Oc_pixy_out_pi[3852,5] <- 0

# levels = chromosomes
levels(Oc_pixy_out_pi$chromosome) <- c(1:15)
ggplot(data = Oc_pixy_out_pi[Oc_pixy_out_pi$chromosome %in% c(1:12),], aes(x=window_pos_1, y=avg_pi, color = pop)) +
  facet_grid(cols = vars(chromosome), scales = "free_x", space = "free_x") +
  geom_line() +
  scale_x_continuous(breaks = seq(0, max(Oc_pixy_out_pi$window_pos_2), by = 100000), labels = seq(0, max(Oc_pixy_out_pi$window_pos_2), by = 100000)/100000) +
  geom_hline(yintercept=quantile(Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="WEA"], .99, na.rm = TRUE),linetype="dashed", color = "green4") +
  geom_hline(yintercept=quantile(Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="WEA"], .95, na.rm = TRUE), color = "green4") +
  geom_hline(yintercept=quantile(Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="EAU"], .01, na.rm = TRUE),linetype="dashed", color = "green") +
  geom_hline(yintercept=quantile(Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="EAU"], .05, na.rm = TRUE), color = "green") +
  theme_light() +
  scale_colour_manual(values=c("green4","green"), 
                      name="Lineages:     ",
                      breaks=c("WEA", "EAU"),
                      labels=c("Western Eurasia", "East Asia plus North America")) +
  theme(legend.position="top", axis.text=element_text(size=15), axis.title.x=element_text(size=15), axis.title.y=element_text(size=17), legend.text = element_text(size=15), strip.text.x = element_text(size = 15), legend.title = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(lwd = 2))) +
  xlab("Position within the chromosome [100 Kbp]") +
  ylab(expression(pi)) +
  scale_y_continuous(breaks = pretty(Oc_pixy_out_pi$avg_pi), labels = abs(pretty(Oc_pixy_out_pi$avg_pi)))


pdf(file = "Fig2.pdf", width = 11.2, height = 7)
Oc_pixy_out_pi <- Oc_pixy_out_pi[Oc_pixy_out_pi$chromosome %in% c(1:12),]
ggplot(data = Oc_pixy_out_pi, aes(x=window_pos_1, y=avg_pi, color = pop)) +
  facet_grid(cols = vars(chromosome), scales = "free_x", space = "free_x") +
  geom_line(data = Oc_pixy_out_pi[!(Oc_pixy_out_pi$chromosome %in% c(2) & Oc_pixy_out_pi$window_pos_1 %in% c(96000:106000)) & !(Oc_pixy_out_pi$chromosome %in% c(5) & Oc_pixy_out_pi$window_pos_1 %in% c(168000:179000)) & !(Oc_pixy_out_pi$chromosome %in% c(5) & Oc_pixy_out_pi$window_pos_1 %in% c(168000:179000)) & !(Oc_pixy_out_pi$chromosome %in% c(5) & Oc_pixy_out_pi$window_pos_1 %in% c(212000:214000)) & !(Oc_pixy_out_pi$chromosome %in% c(6) & Oc_pixy_out_pi$window_pos_1 %in% c(32000:35000)) & !(Oc_pixy_out_pi$chromosome %in% c(12) & Oc_pixy_out_pi$window_pos_1 %in% c(141000:146000)),]) +
  geom_line(data = Oc_pixy_out_pi[Oc_pixy_out_pi$chromosome %in% c(2) & Oc_pixy_out_pi$window_pos_1 %in% c(96000:106000),], col="black") +
  geom_line(data = Oc_pixy_out_pi[Oc_pixy_out_pi$chromosome %in% c(5) & Oc_pixy_out_pi$window_pos_1 %in% c(168000:179000),], col="black") +
  geom_line(data = Oc_pixy_out_pi[Oc_pixy_out_pi$chromosome %in% c(8) & Oc_pixy_out_pi$window_pos_1 %in% c(61000:68000),], col="black") +
  geom_line(data = Oc_pixy_out_pi[Oc_pixy_out_pi$chromosome %in% c(5) & Oc_pixy_out_pi$window_pos_1 %in% c(212000:214000),], col="red") +
  geom_line(data = Oc_pixy_out_pi[Oc_pixy_out_pi$chromosome %in% c(6) & Oc_pixy_out_pi$window_pos_1 %in% c(32000:35000),], col="red") +
  geom_line(data = Oc_pixy_out_pi[Oc_pixy_out_pi$chromosome %in% c(12) & Oc_pixy_out_pi$window_pos_1 %in% c(141000:146000),], col="red") +
  scale_x_continuous(breaks = seq(0, max(Oc_pixy_out_pi$window_pos_2), by = 100000), labels = seq(0, max(Oc_pixy_out_pi$window_pos_2), by = 100000)/100000) +
  geom_hline(yintercept=quantile(Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="WEA"], .99, na.rm = TRUE),linetype="dashed", color = "green4") +
  geom_hline(yintercept=quantile(Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="WEA"], .95, na.rm = TRUE), color = "green4") +
  geom_hline(yintercept=quantile(Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="EAU"], .01, na.rm = TRUE),linetype="dashed", color = "green") +
  geom_hline(yintercept=quantile(Oc_pixy_out_pi$avg_pi[Oc_pixy_out_pi$pop=="EAU"], .05, na.rm = TRUE), color = "green") +
  theme_light() +
  scale_colour_manual(values=c("green4","green"), 
                      name="Lineages:     ",
                      breaks=c("WEA", "EAU"),
                      labels=c("Western Eurasia", "East Asia plus North America")) +
  theme(legend.position="top", axis.text=element_text(size=15), axis.title.x=element_text(size=15), axis.title.y=element_text(size=17), legend.text = element_text(size=15), strip.text.x = element_text(size = 15), legend.title = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(lwd = 2))) +
  xlab("Position within the chromosome [100 Kbp]") +
  ylab(expression(pi)) +
  scale_y_continuous(breaks = pretty(Oc_pixy_out_pi$avg_pi), labels = abs(pretty(Oc_pixy_out_pi$avg_pi)))
dev.off()
