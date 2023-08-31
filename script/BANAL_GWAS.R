#!/usr/bin/env R

library(tidyverse)
library(broom)
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)
library(scales)
library("ggmsa")
library("ggrepel")

setwd("current_directory")


species.interest.v <- c("R.cornutus","R.ferrumequinum","R.shamel","R.pusillus","R.macrotis","R.affinis","R.sinicus","Pangolin","Hamster")

data.align.name <- 'input/ACE2_aligned.txt'
data.align <- read.table(data.align.name,header=T,check.names=F,colClasses= "character")
colnames(data.align)[1] <- "species"

data.name <- 'input/220922_ACE2_pseudovirus.txt'
data <- read.table(data.name,header=T)
data[,2:3] <- log(data[,2:3],10)
data <- data %>% mutate(diff = BANAL.236 - BANAL.52)
data <- data %>% inner_join(data.align,by="species") %>% arrange(desc(diff))

plot.df <- data %>% select(BANAL.52,BANAL.236,species)
plot.df <- plot.df %>% mutate(group = ifelse(species %in% species.interest.v,"include","exclude"))
plot.df.filtered <- plot.df %>% filter(group == "include")


g <- ggplot(plot.df,aes(y=BANAL.52,x=BANAL.236,label=species, color = group))
g <- g + geom_point(shape=20)
g <- g + geom_smooth(data = plot.df.filtered, method="lm")
#g <- g + geom_text_repel()
g <- g + theme_set(theme_classic(base_size = 10, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
)
g <- g + theme_set(theme_classic(base_size = 10, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
)
g <- g + xlab("BANAL-52 (log10)")+ ylab("BANAL-236 (log10)")
g <- g + scale_x_continuous(breaks=c(3.5,4,4.5,5,5.5)) + scale_y_continuous(breaks=c(2.5,3.5,4.5,5.5))
g <- g + scale_color_manual(breaks=c("include","exclude"),values=c("black","gray70"))


pdf.name <- 'output/BANAL52_236_infectivity.pdf'
pdf(pdf.name,width=4,height=3)
plot(g)
dev.off()


g <- ggplot(plot.df,aes(y=BANAL.52,x=BANAL.236,label=species, color = group))
g <- g + geom_point(shape=20)
g <- g + geom_smooth(data = plot.df.filtered, method="lm")
g <- g + geom_text_repel()
g <- g + theme_set(theme_classic(base_size = 10, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank()
)
g <- g + theme_set(theme_classic(base_size = 10, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank()
)
g <- g + xlab("BANAL-52 (log10)")+ ylab("BANAL-236 (log10)")
g <- g + scale_x_continuous(breaks=c(3.5,4,4.5,5,5.5),limits=c(4,5)) + scale_y_continuous(breaks=c(2.5,3.5,4.5,5.5),limits=c(3.5,5.5))
g <- g + scale_color_manual(breaks=c("include","exclude"),values=c("black","gray70"))


pdf.name <- 'output/BANAL52_236_infectivity.pos.pdf'
pdf(pdf.name,width=4,height=3)
plot(g)
dev.off()

g




data <- data %>% filter(species %in% species.interest.v)

res.df <- data.frame()
for(i in 5:(ncol(data)-1)){

pos <- colnames(data)[i]

temp.df <- data %>% select(BANAL.52,BANAL.236,diff) %>% mutate(aa = as.factor(data[,i]))

if(length(unique(temp.df$aa)) > 1){

fit.diff <- aov(diff ~ aa,data=temp.df) %>% tidy()
fit.BANAL.52 <- aov(BANAL.52 ~ aa,data=temp.df) %>% tidy()
fit.BANAL.236 <- aov(BANAL.236 ~ aa,data=temp.df) %>% tidy()

pval.diff <- fit.diff$p.value[1]
pval.BANAL.52 <- fit.BANAL.52$p.value[1]
pval.BANAL.236 <- fit.BANAL.236$p.value[1]

tmp.df.res <- data.frame(pos = pos, pval.BANAL.52 = pval.BANAL.52, pval.BANAL.236 = pval.BANAL.236, pval.diff = pval.diff)

res.df <- rbind(res.df,tmp.df.res)
}
}

res.df.filtered <- res.df %>% filter(pval.diff < 0.1)
data.selected <- data %>% select(all_of(c("species","BANAL.52","BANAL.236","diff",res.df.filtered$pos)))

data.selected <- data.selected %>% mutate(species = factor(species,levels = rev(data.selected$species)))
data.selected.scaled <- data.selected %>% select(BANAL.52,BANAL.236,diff) %>% scale() %>% as.data.frame() %>% mutate(species = data.selected$species)
data.selected.scaled <- data.selected.scaled %>% mutate(species = factor(species,levels = rev(data.selected$species)))

data.selected.scaled.long <- data.selected.scaled %>% gather(key = type,value=value,-species)

g1 <- ggplot(data.selected.scaled.long,aes(x=type,y=species,fill=value))
g1 <- g1 + geom_tile()
g1 <- g1 + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g1 <- g1 + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 strip.text = element_text(size=8)
)
g1 <- g1 + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 strip.text = element_text(size=8)
)
g1 <- g1 + scale_fill_gradientn(colors = brewer.pal(9, "BuPu"), limits = c(-1.5,1.5), oob = squish)
g1 <- g1 + xlab("") + ylab("")
g1

pdf.name <- 'output/BANAL52_236_heatmap.pdf'
pdf(pdf.name,width=4,height=4)
plot(g1)
dev.off()


data.seq <- data.selected %>% unite(col = seq,5:ncol(data.selected),sep = "")
seq.v <- data.seq$seq
names(seq.v) <- data.seq$species

AA_sequences <- AAStringSet(seq.v)

g1 <- ggmsa(AA_sequences,
           1, 8,
           color = "Chemistry_AA",
           font = "DroidSansMono",
           char_width = 0.5,
           seq_name = T)

res.df.filtered.long <- res.df.filtered %>% gather(key=type,value=value,-pos)
res.df.filtered.long <- res.df.filtered.long %>% mutate(pos = factor(pos,levels=res.df.filtered$pos))

g2 <- ggplot(res.df.filtered.long,aes(x=pos,y=type,fill=-log(value,10)))
g2 <- g2 + geom_tile()
g2 <- g2 + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g2 <- g2 + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 strip.text = element_text(size=8)
)
g2 <- g2 + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 strip.text = element_text(size=8)
)
g2 <- g2 + scale_fill_gradientn(colors = brewer.pal(9, "BuGn"), limits= c(0,2), oob = squish)
g2 <- g2 + xlab("") + ylab("")
g2

g3 <- g1 / g2 #+ plot_layout(ncol = 1, heights = c(3, 1))

pdf.name <- 'output/BANAL52_236_msa.pdf'
pdf(pdf.name,width=5,height=5)
plot(g3)
dev.off()




g <- ggplot(res.df,aes(x=as.numeric(pos),y=-log(pval.diff,10)))
g <- g + geom_point()
g <- g + geom_text_repel(data = res.df %>% filter(pval.diff < 0.1), 
                aes(x=as.numeric(pos),y=-log(pval.diff,10),label=pos))
g <- g + geom_hline(yintercept=1, linetype="dashed", color = "gray70", size=0.5)
g <- g + xlab("Position") + ylab("-log10(P)")

pdf.name <- 'output/BANAL52_236_manhattan.pdf'
pdf(pdf.name,width=4.5,height=3)
plot(g)
dev.off()



