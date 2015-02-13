#' Experimental function to plot phylogenies of the class 'phylo' alongside multiple dataframes which can be an alignments
#'
#' @param An object of the class "phylo" and a list of dataframe
#' @param palette_list a list of color pallettes to use, for discrete traits, the order of the colours must be according to factors 
#' @export
#' @author Joseph Hughes, Anton Camacho
#' @examples see plot_ggphydata_test.R

# assume that if it is strings, it will be discrete traits, otherwise treat it as continuous
plot_ggphydata<-function(phylo,dataframelist,plotsize=NULL,titles=NULL,Palette=NULL,tip_labels=F,tip_attribute=NULL,var_tip_labels=NULL,var_tip_colour=NULL,tip.dates=NULL,branch.unit=NULL){
# input tree and dataframe (an alignment can be a dataframe)
   ggphy<-phylo2ggphy(phylo, tip.dates = NULL, branch.unit = NULL,verbose = TRUE)
   rownames(ggphy[[1]])<-ggphy[[1]]$label
   
   # processing the matrices
   # keep all tree tips and add empty rows for missing info in dataframes, remove rows not found in tree tips
   # type corresponds to whether it is DNA, discrete or continuous as the color scale will be different
   alldata <- data.frame(Label=character(),
                 ColNb=character(), 
                 value=character(), 
                 Dataset=character(),
                 Type=character()) 
   p<-list()  
   size<-c(0)
   nbdatasets<-length(dataframelist)
   for (i in 1:nbdatasets){
     # need to melt all the different data frames together and then present them all as a facet_grid 
     if (inherits(dataframelist[[i]],"DNAbin")){
       # process DNAbin in a particular way
       m_dna<-as.character(dataframelist[[i]])
       size<-c(size,ncol(m_dna))
       # sort the DNAbin according to order of tree labels and add empty rows
       # also need to remove a row if labels are not present in phylo
       m_dna_empty<-merge(ggphy[[1]],m_dna,by="row.names",all.x=TRUE)[,-(1:5)]
       rownames(m_dna_empty)<-ggphy[[1]]$label
       colnames(m_dna_empty)<-1:ncol(m_dna_empty)
       longData<-melt(as.matrix(m_dna_empty))
       #write.table(longData,file="test2",row.names=TRUE,sep="\t")
      # zp1 <- ggplot(longData,aes(x = Var2, y = Var1, fill = value))
      # zp1 <- zp1 + geom_tile() + theme_bw() + xlab("") + ylab("")   
      # print(zp1)
      longData$Var2<-as.character(longData$Var2)
       longData$Dataset<-paste("Set",i,sep="")
       longData$Type<-"DNA"
       alldata<-rbind(alldata,longData)
       #print(head(alldata))
       #print(titles[i])
       longData$Var2<-as.numeric(as.character(longData$Var2))
       breaks<-seq(10,max(longData$Var2),10)
       
       p[[i]]<- ggplot(longData,aes(x = Var2, y = Var1, fill = value)) + geom_tile() + xlab(titles[i]) + scale_x_continuous(breaks = breaks) 
       if (!is.null(Palette[[i]])){
         p[[i]]<- p[[i]] + scale_fill_manual(values=Palette[[i]])
       }
       p[[i]]<- p[[i]] + theme(legend.position="right",axis.title.y=element_blank(),
           panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
           panel.grid.minor=element_blank(),plot.background=element_blank(),plot.margin = unit(x = c(0, 0, 0, 0), units = "mm"),legend.margin=unit(-0.6,"cm"))
     }else{
     # check if it is discrete/integer and check if it is continuous
       #print(dataframelist[[i]])
       m_df_empty<-merge(ggphy[[1]],dataframelist[[i]],by="row.names",all.x=TRUE)[,-(1:5)]
       m_df_empty<-as.data.frame(m_df_empty)
       
       #print(m_df_empty)
       numeric_check<-sapply(m_df_empty, is.numeric)
       #print(sum(numeric_check))
       if (sum(numeric_check)){
         #print("numeric")
         rownames(m_df_empty)<-ggphy[[1]]$label
         longNumeric<-melt(as.matrix(m_df_empty))
         size<-c(size,ncol(m_df_empty))
         #print(longNumeric)
         longNumeric$Dataset<-paste("Set",i,sep="")
         longNumeric$Type<-"numeric"
         alldata<-rbind(alldata,longNumeric)
         p[[i]]<- ggplot(longNumeric,aes(x = Var2, y = Var1, fill = value)) + geom_tile() + xlab(titles[i]) 
         if (!is.null(Palette[[i]])){
           if(length(Palette[[i]])==2){
           p[[i]]<- p[[i]] + scale_fill_gradient2(low=Palette[[i]][1],high=Palette[[i]][2])
           }else if(length(Palette[[i]])==3){
           p[[i]]<- p[[i]] + scale_fill_gradient2(low=Palette[[i]][1],mid=Palette[[i]][2],high=Palette[[i]][2])
           }else{
           print("Colour gradient for continuous dataframe must have 2 or 3 colours. Using default")
           }
         }
         
         p[[i]]<- p[[i]] + theme(legend.position="right",axis.title.y=element_blank(),
           panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
           panel.grid.minor=element_blank(),plot.background=element_blank(),plot.margin = unit(x = c(0, 0, 0, 0), units = "mm"),legend.margin=unit(-0.6,"cm"))
       }else if (!sum(numeric_check)){
         #print("discrete")
         rownames(m_df_empty)<-ggphy[[1]]$label
         longDiscrete<-melt(as.matrix(m_df_empty))
         size<-c(size,ncol(m_df_empty))
         #print(longDiscrete)
         longDiscrete$Dataset<-paste("Set",i,sep="")
         longDiscrete$Type<-"discrete"
         alldata<-rbind(alldata,longDiscrete)
         p[[i]]<- ggplot(longDiscrete,aes(x = Var2, y = Var1, fill = value)) + geom_tile()  + xlab(titles[i]) 
         if (!is.null(Palette[[i]])){
           p[[i]]<- p[[i]] + scale_fill_manual(values=Palette[[i]])
         }
         p[[i]]<- p[[i]] + theme(legend.position="right",axis.title.y=element_blank(),
           panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
           panel.grid.minor=element_blank(),plot.background=element_blank(),plot.margin = unit(x = c(0, 0, 0, 0), units = "mm"),legend.margin=unit(-0.6,"cm"))
       }
       
     }
     
   }
   
   # plotting the tree - code taken from Anton Camacho's plotggphy
   # need to modify y axis positions to align with matrix figures
   df_tip<-ggphy[[1]]
   df_node<-ggphy[[2]]
   df_edge<-ggphy[[3]]
   #print(df_tip)
   #print(df_node)
   #print(df_edge)
   is_x_date<-inherits(df_edge$x_beg,"Date")
   is_x_date<-inherits(df_edge$x_beg,"Date")
   
   if(!is.null(tip_attribute) & !is.null(var_tip_labels)){
    #merge df_tip with tip attributes
        tmp<-merge(df_tip,tip_attribute,by.x="label",by.y=var_tip_labels)
        df_tip<-tmp
   }
    
   #theme_set(theme_grey())
   theme_old<-theme_update(
                            axis.ticks.y = element_blank(),plot.margin = unit(x = c(0, 0, 0, 0), units = "mm"),legend.margin=unit(-0.6,"cm"),
                            axis.title.y = element_blank(),panel.grid.major.y = element_blank(),
                            panel.grid.minor.y = element_blank(),panel.background=element_blank())
    
    pp<-ggplot(df_edge)
    pp<-pp+geom_segment(data=df_edge,aes(x=x.beg,xend=x.end,y=y.beg,yend=y.end),lineend="round")
    pp<-pp+scale_y_continuous("",labels=NULL)
    if(is_x_date)
        pp<-pp+scale_x_date("Time",labels=date_format("%Y"),minor_breaks="1 year")
        else
            pp<-pp+scale_x_continuous("Time")
            
            if(tip_labels){
                pp<-pp+geom_text(data=df_tip,aes(x=x,y=y,label=label),hjust=0)
            }
    if(!is.null(var_tip_colour)){
       pp<-pp+geom_point(data=df_tip,aes_string(x="x",y="y",colour=var_tip_colour))
    }
      
   #print(size)
   # change sizes so that it allows for 40% of space for the tree
   total_width_df=sum(size[-1])
   phylo_width=sum(size[-1])*0.4
   size<-(size[-1])/(phylo_width+total_width_df)
   
   #print(size)
   #do.call("grid.arrange",c(p, ncol=3))
   grid.newpage()
   vp_nb<-length(p)+1
   pushViewport(viewport(layout=grid.layout(1,vp_nb,widths=plotsize)))
   print(pp, vp=viewport(layout.pos.row=1,layout.pos.col=1))
   for (i in 1:length(p)){
     print(p[[i]], vp=viewport(layout.pos.row=1,layout.pos.col=i+1))
   }
}


