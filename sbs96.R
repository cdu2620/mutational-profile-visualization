setwd('"/Users/cdu2620"')

library(tidyverse)
library(hrbrthemes)
library(scales)
library(ggpubr)
library(ggsci)
library(janitor)
library(factoextra)
library(ggplot2)
library(dygraphs)
library(rlist)

signatureFunction <- function(dataset, columnInput, log.show, reg.show, colorInput) {
   input <- read_csv(dataset,col_names = T)
  ## merge small number of signature to "Signature Subs-others" using default creteria minNumMutation 100
  
  
  studydata <- input %>% select(Samples,Study)
  similaritydata <- input %>% select(Samples,Similarity)
  data <- input %>% select(-Study,-Similarity) %>% gather(Signatures,Contribution,contains('Signature'))
  reg.show <- TRUE
  if (reg.show == TRUE) {
    minNumMutation <- 1000
    existsigs <- data %>% group_by(Signatures) %>% summarise(Total=(sum(Contribution))) %>% filter(Total>=minNumMutation) %>% pull(Signatures)
    # show the signatures merged to others
    data %>% group_by(Signatures) %>% summarise(Total=(sum(Contribution))) %>% filter(Total<minNumMutation)
    
    data <- data %>% mutate(Signatures=if_else(Signatures %in% existsigs,Signatures,'Signature Subs-others')) %>% group_by(Signatures,Samples)%>% summarise(Contribution=(sum(Contribution))) %>%  spread(Signatures,Contribution) 
  }
  
  if (log.show == TRUE) {
    minNumMutation <- 10
    existsigs <- data %>% group_by(Signatures) %>% summarise(Total=log(sum(Contribution))) %>% filter(Total>=minNumMutation) %>% pull(Signatures)
    
    # show the signatures merged to others
    data %>% group_by(Signatures) %>% summarise(Total=log(sum(Contribution))) %>% filter(Total<minNumMutation)
    
    data <- data %>% mutate(Signatures=if_else(Signatures %in% existsigs,Signatures,'Signature Subs-others')) %>% group_by(Signatures,Samples)%>% summarise(Contribution=log(sum(Contribution))) %>%  spread(Signatures,Contribution) 
  }
  
  sigcolor <- c(
    'Signature Subs-01'='#4a9855',
    'Signature Subs-02'='#e2a8ab',
    'Signature Subs-03'='#40004b',
    'Signature Subs-04'='#5aa1ca',
    'Signature Subs-05'='#305d39',
    'Signature Subs-06'='#785940',
    "Signature Subs-07a"='#6e70b7',
    "Signature Subs-07b"='#ff7f00',
    "Signature Subs-07c"='#fec44f',
    "Signature Subs-07d"='#846a2a',
    "Signature Subs-08"='#cab2d6',
    "Signature Subs-10a"='#8dd3c7',
    "Signature Subs-10b"='#5e4fa2',
    'Signature Subs-12'='#ffed6f',
    'Signature Subs-13'='#e41a1c',
    'Signature Subs-14'='#ffffbf',
    'Signature Subs-15'='#4d4d4d',
    #'Signature Subs-17'='#513276',
    'Signature Subs-17a'='#df4c7d',
    'Signature Subs-17b'='#2c7fb8',
    'Signature Subs-18'='#b3de69',
    'Signature Subs-20'='#b2182b',
    'Signature Subs-21'='#e6f5d0',
    'Signature Subs-25'='#35978f',
    'Signature Subs-28'='#de77ae',
    'Signature Subs-31'='#f781bf',
    'Signature Subs-33'='#b25d7e',
    'Signature Subs-36'='yellow',
    'Signature Subs-40'='#b15928',
    'Signature Subs-46'='#e6f598',
    'Signature Subs-42'='#ae017e',
    'Signature Subs-54'='#fcc5c0',
    'Signature Subs-others'='#cececa'
  )
  
  studyname <- unique(studydata$Study)
  studycolor <- pal_npg()(length(studyname))
  names(studycolor) <- studyname
  
  sigName <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
  colorGrad <- colorRampPalette(c("yellow", "purple"))
  colorGrad <- colorGrad(length(sigName))
  
  sigColors = c()
  for (i in 1:length(sigName)) {
    val <- sigName[[i]]
    sigColors[val] <- colorGrad[[i]]
  }
  
  names(sigColors) <-sigName
  sigcolorindex <- sigcolor[colnames(data)[-1]]
  sigcolor <- c(sigcolorindex,"white",studycolor)
  names(sigcolor) <- c(names(sigcolorindex),"        ",names(studycolor))
  
  sigcolor2 <- c(sigcolorindex, "white", sigColors)
  names(sigcolor2) <- c(names(sigcolorindex),"        ",names(sigColors))
  
  tmp <- data %>% adorn_percentages('row') 
  tmp2 <- tmp
  mdata <- as.matrix(tmp[,-1])
  rownames(mdata) <- tmp$Samples
  
  fviz_nbclust(mdata, kmeans, method = "gap_stat")
  
  color12 <- rev(c(
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#b15928"
  ))
  
  # propColors <- c('A'='#a6cee3',
  #                      'B'='#1f78b4',
  #                      'C'='#b2df8a',
  #                      'D'='#33a02c',
  #                      'E'='#fb9a99',
  #                      'F'='#e31a1c',
  #                      "G"='#fdbf6f',
  #                      "H"='#ff7f00',
  #                      "I"='#fec44f',
  #                      'J'='#cab2d6',
  #                      'K'='#6a3d9a')
  
  #user inputs a list of their own colors
  if (colorInput == TRUE) {
    my.colors <- readline(prompt = "Enter your colors:")
    
    x <- as.character(unlist(strsplit(my.colors, ",")))
    
    color12 = rev(c(x))
    
    #more user color inputs
    my.sigcolors <- readline(prompt = "Enter your signature colors as rgbs:")
    
    listSigColors <- as.character(unlist(strsplit(my.sigcolors, ",")))
    
    sigcolor <- c()
    matches <- append(sigcolor, colnames(tmp2))
    
    for (i in 2:length(matches)) {
      val <- matches[[i-1]]
      sigcolor[val] <- listSigColors[[i-1]]
    }
  }
  res <- hcut(mdata,k = 8, stand = TRUE,hc_func = 'hclust',hc_metric = 'euclidean',hc_method = 'ward.D2')
  
  tempAllGraphs <- list()
  #p1
  assign(paste0("p", length(tempAllGraphs) + 1), fviz_dend(res, rect = TRUE, cex = 0.5,k_colors = color12,lwd = 0.5,show_labels = FALSE)+scale_x_discrete(expand = c(0,0))+ theme(plot.margin=margin(b=-0.2,unit="cm"),title = element_blank()))
  tempAllGraphs <- list.append(tempAllGraphs, p1)
  
  #p2
  assign(paste0("p", length(tempAllGraphs) +1), data %>% gather(Signature,Weight,contains("Signature")) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,Weight,fill=factor(Signature,levels = names(sigcolor))))+geom_bar(stat="identity",position="fill",col="gray95",width = 1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12)+theme(legend.title = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1),panel.grid.major=element_line(),legend.position = "bottom",legend.box.background = element_blank(),legend.box.spacing = unit(-0.5,"cm"),legend.key = element_rect(size = 0),axis.ticks.y = element_line(colour = "black"),legend.key.size = unit(0.25, "cm"),legend.key.width =unit(1, "cm"))+scale_y_continuous(expand = c(0, 0),breaks = pretty_breaks())+xlab("")+scale_fill_manual(values = sigcolor,drop=FALSE)+guides(fill=guide_legend(nrow=2,byrow=TRUE,label.position = "bottom"))+ylab("Signature contribution\n"))
  tempAllGraphs <- list.append(tempAllGraphs, p2)
  
  #p3
  assign(paste0("p", length(tempAllGraphs) + 1), data %>% gather(Signature,Weight,contains("Signature")) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,Weight,fill=factor(Signature,levels = names(sigcolor2))))+geom_bar(stat="identity",position="fill",col="gray95",width = 1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12)+theme(legend.title = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1),panel.grid.major=element_line(),legend.position = "bottom",legend.box.background = element_blank(),legend.box.spacing = unit(-0.5,"cm"),legend.key = element_rect(size = 0),axis.ticks.y = element_line(colour = "black"),legend.key.size = unit(0.25, "cm"),legend.key.width =unit(1, "cm"))+scale_y_continuous(expand = c(0, 0),breaks = pretty_breaks())+xlab("")+scale_fill_manual(values = c(sigcolor, sigcolor2),drop=FALSE)+guides(fill=guide_legend(nrow=2,byrow=TRUE,label.position = "bottom"))+ylab("Signature contribution\n"))
  tempAllGraphs <- list.append(tempAllGraphs, p3)
  
  # get_proportion <- function(d) {
  #   require(tidyverse)
  #   
  #   d_gather <- d %>% 
  #     gather(key = "signature", value = "sig_count", -Samples)
  #   
  #   d_sum <- d_gather %>% 
  #     group_by(Samples) %>% 
  #     summarise(sample_sum = sum(sig_count))
  #   
  #   d_prop <- d_gather %>% 
  #     left_join(d_sum, by = "Samples") %>% 
  #     mutate(prop = sig_count / sample_sum) %>% 
  #     select(-sig_count, -sample_sum) %>% 
  #     mutate(signature = paste0(signature,"_prop")) %>% 
  #     spread(key = signature, value = prop)
  #   
  #   d %>% left_join(d_prop, by = "Samples")
  #   
  # }
  # 
  # pset3 <- get_proportion(data)
  # 
  # pset3 <- pset3[, c(1, 15:27)]
  
  #separates the proportion charts out by signature
  listNames <- list()
  listSigs <- list()
  for (i in 2:ncol(tmp2)) {
    listNames[[i]] <- colnames(tmp2[,i])
    listSigs[[i]] <- assign(paste0("x", i-1), tmp2[, c(1, i)])
  }
  listNames[[1]] <- NULL
  
  # listSigsY <- list()
  # for (i in 2:nrow(tmp2)) {
  #   listSigsY[[i]] <- assign(paste0("y", i), tmp[c(i),])
  # }
  
  for (i in 2:ncol(tmp2)) {
    listSigs[[i]][, 2][listSigs[[i]][, 2] >= 0 & listSigs[[i]][, 2] < .1 ] <- "A"
    listSigs[[i]][, 2][listSigs[[i]][, 2] >= .1 & listSigs[[i]][, 2] < .2 ] <- "B" 
    listSigs[[i]][, 2][listSigs[[i]][, 2] >= .2 & listSigs[[i]][, 2] < .3 ] <- "C"
    listSigs[[i]][, 2][listSigs[[i]][, 2] >= .3 & listSigs[[i]][, 2] < .4 ] <- "D"
    listSigs[[i]][, 2][listSigs[[i]][, 2] >= .4 & listSigs[[i]][, 2] < .5 ] <- "E"
    listSigs[[i]][, 2][listSigs[[i]][, 2] >= .5 & listSigs[[i]][, 2] < .6 ] <- "F"
    listSigs[[i]][, 2][listSigs[[i]][, 2] >= .6 & listSigs[[i]][, 2] < .7 ] <- "G"
    listSigs[[i]][, 2][listSigs[[i]][, 2] >= .7 & listSigs[[i]][, 2] < .8 ] <- "H"
    listSigs[[i]][, 2][listSigs[[i]][, 2] >= .8 & listSigs[[i]][, 2] < .9 ] <- "I"
    listSigs[[i]][, 2][listSigs[[i]][, 2] >= .9 & listSigs[[i]][, 2] < 1 ] <- "J"
    assign(paste0("x", i), listSigs[[i]])
  }
  
  # for (i in 2:nrow(tmp2)) {
  #   listSigsY[[i]][1,][listSigsY[[i]][1,] >= 0 & listSigsY[[i]][1,] < .1 ] <- "A"
  #   listSigsY[[i]][1,][listSigsY[[i]][1,] >= .1 & listSigsY[[i]][1,] < .2 ] <- "B"
  #   listSigsY[[i]][1,][listSigsY[[i]][1,] >= .2 & listSigsY[[i]][1,] < .3 ] <- "C"
  #   listSigsY[[i]][1,][listSigsY[[i]][1,] >= .3 & listSigsY[[i]][1,] < .4 ] <- "D"
  #   listSigsY[[i]][1,][listSigsY[[i]][1,] >= .4 & listSigsY[[i]][1,] < .5 ] <- "E"
  #   listSigsY[[i]][1,][listSigsY[[i]][1,] >= .5 & listSigsY[[i]][1,] < .6 ] <- "G"
  #   listSigsY[[i]][1,][listSigsY[[i]][1,] >= .6 & listSigsY[[i]][1,] < .7 ] <- "H"
  #   listSigsY[[i]][1,][listSigsY[[i]][1,] >= .7 & listSigsY[[i]][1,] < .8 ] <- "I"
  #   listSigsY[[i]][1,][listSigsY[[i]][1,] >= .8 & listSigsY[[i]][1,] < .9 ] <- "J"
  #   listSigsY[[i]][1,][listSigsY[[i]][1,] >= .9 & listSigsY[[i]][1,] < 1 ] <- "K"
  #   assign(paste0("y", i), listSigsY[[i]])
  # }
  listSigs[[1]] <- NULL
  testDf <- data.frame("A" = "0-.1", "B" = ".1-.2", "C" = ".2-.3",
                       "D" = ".3-.4", "E" = ".4-.5", "F" = ".5-.6", "G" = ".6-.7",
                       "H" = ".7-.8", "I" = ".8-.9", "J" = ".9-1")
  #p4
  assign(paste0("p", length(tempAllGraphs) +1), studydata %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,1,fill=factor(Study,levels = names(studycolor))))+geom_tile(col="black")+scale_fill_manual(values =studycolor)+theme_minimal()+theme(legend.position = "none",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-0.8,unit="cm"),title = element_blank())+ylim(c(0,2)))
  tempAllGraphs <- list.append(tempAllGraphs, p4)
  
  for (i in 1:length(listSigs)) {
    if (listNames[[i]] == columnInput) {
      mm1 <- as.matrix(listSigs[[i]])
      mm2 <- matrix(mm1, ncol = ncol(listSigs[[i]]), dimnames = NULL)
      result <- listSigs[[i]] %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,1,fill=factor(mm2[,2],levels = names(sigColors))))+geom_tile(col="black")+theme_minimal()+theme(legend.position = "none",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-0.8,unit="cm"),title = element_blank())+ylim(c(0,2))
      p7 <- result+geom_bar(stat="identity",position="fill",col="gray95",width = 1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12)+theme(legend.title = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1),panel.grid.major=element_line(),legend.position = "bottom",legend.box.background = element_blank(),legend.box.spacing = unit(-0.5,"cm"),legend.key = element_rect(size = 0),legend.key.size = unit(0.25, "cm"),legend.key.width =unit(1, "cm"), axis.ticks.y = element_blank(), axis.text.y = element_blank())+xlab("")+scale_fill_manual(values = c(sigcolor, sigcolor2),drop=FALSE)+guides(fill=guide_legend(nrow=2,byrow=TRUE,label.position = "bottom"))
      p3.5 <- listSigs[[i]] %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,1,fill=factor(mm2[,2],levels = names(sigColors))))+geom_tile(col="black")+scale_fill_manual(values =sigColors)+theme_minimal()+theme(legend.position = "none",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-0.8,unit="cm"),title = element_blank())+ylim(c(0,2))
      }
  }
  
  #p5
  assign(paste0("p", length(tempAllGraphs)+1), data %>% gather(Signature,Weight,contains("Signature")) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,Weight,fill=factor(Signature,levels = names(sigcolor))))+geom_bar(stat="identity",col="gray95",width=1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = "xy")+theme(legend.title = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_blank(),panel.grid.major.x=element_blank(),legend.position = "none",legend.box.spacing = unit(0,"cm"),plot.margin=margin(b=-1.3,t = 1,unit="cm"))+scale_y_continuous(expand = c(0, 0),breaks = pretty_breaks(),labels = comma)+scale_fill_manual(values = sigcolor,drop=FALSE)+guides(fill=guide_legend(nrow=2,byrow=TRUE))+xlab("")+ylab("Number of mutations \n"))
  tempAllGraphs <- list.append(tempAllGraphs, p5)
  
  listOfShown <- list()
  allShows <- list()
  allGraphs <- list()
  
   p1.show <- readline(prompt = "Show clusters? (TRUE/FALSE):")
  allShows <- list.append(allShows, p1.show)
  allGraphs <- list.append(allGraphs, p1)
   p4.show <- readline(prompt = "Show study colors? (TRUE/FALSE):")
  allShows <- list.append(allShows, p4.show)
  allGraphs <- list.append(allGraphs, p4)
   if (p4.show == TRUE) {
     p2.show <- TRUE
   }
  else {
    p2.show <- FALSE
  }
  allShows <- list.append(allShows, p2.show)
  allGraphs <- list.append(allGraphs, p2)
   p5.show <- readline(prompt = "Show proportion colors? (TRUE/FALSE):")
  allShows <- list.append(allShows, p5.show)
  allGraphs <- list.append(allGraphs, p5)
   p3.show <- readline(prompt = "Show prop legend? (TRUE/FALSE):")
  allShows <- list.append(allShows, p3.show)
  allGraphs <- list.append(allGraphs, p3)
  
  val = 1
  for (i in 1:length(allShows)) {
    if (allShows[[i]] == TRUE) {
        if (identical(allGraphs[[i]], p3)) {
          saveIndex <- i
        }
        listOfShown[[val]] <- allGraphs[[i]]
        val = val + 1
    }
  }
  if (p3.show == TRUE){
    listOfShown <- list.insert(listOfShown, saveIndex, p3.5)
  }
  ggarrange(plotlist=listOfShown,ncol=1,nrow = 6,align = 'v',heights = c(2,0.1,3,7,2,2))
  ggsave("SBS96_ludmil.pdf",height = 10,width = 20,device = cairo_pdf)
  
  # testList <- list()
  # for (i in 1:length(listNames)) {
  #   val <- paste0("outline", i)
  #   testList[[val]] <- checkbox(FALSE, listNames[[i]])
  # }
  # library(gWidgets2)
  # library(manipulate)
  # manipulate({p4
  #   if ("outline1" == TRUE) {
  #     p3.5
  #   }},
  #   testList
  # )
  # 
  # manipulate(
  #   boxplot(Freq ~ Class, data = Titanic, outline = outline),
  #   outline = checkbox(FALSE, "Show outliers"),)
  
  #notes for future:
   #further customization: match up the size of each graph
  #put in the data table as an image
}


#user puts in the dataset, the column they want looked at, if they want the graph to have a log transformation, and whether they're importing colors or using the default colors
signatureFunction('SBS_signatures_in_samples.csv', "Signature Subs-05", FALSE, TRUE, FALSE)

