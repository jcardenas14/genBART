library(ggplot2)
#library(fastcluster)
#library(qusage)
#library(NMF)
#library(pca3d)
#library(shinyjs)

enableBookmarking(store = "server")
load("moduleinfo.Rdata")
load("moduleinfo_rna.RData")
moduleinfo2<-qusage::read.gmt("BaylorModules.gmt")
x1<-unlist(moduleinfo2)
x2<-rep(names(moduleinfo2),times=lapply(moduleinfo2,length))
x2<-gsub("_",".",x2)
moduleinfo2<-data.frame(SYMBOL=x1,Module=x2)
names(moduleinfo)[6] <- "affy"
names(moduleinfo)[8] <- "Module_V3"
module_annotations <- read.csv("v2_annotated_module_list.csv", header = T)
moduleinfo$Modulev2_Annotation <- module_annotations[match(moduleinfo$Module, module_annotations$Module), 2]
moduleinfo_rna$Modulev2_Annotation <- module_annotations[match(moduleinfo_rna$Module, module_annotations$Module), 2]

get_all_tables <- function(dat, my_i){
  nam <- names(dat)
  index1 <- grep("p.Value",nam,fixed=T)
  index2 <- grep("log10",nam,fixed=T)
  p.nams <- nam[setdiff(index1,index2)]
  nams <- gsub("p.Value.for.Estimate.of.", "", p.nams)
  sc_n <- paste0("sc_", my_i)
  assign(sc_n, data.frame(dat[, c(1, 2)], dat[, grep(nams[my_i], names(dat))]))
  write.csv(get(paste0("sc_", my_i)), file = paste0(nams[my_i], ".csv"))
}

helpPopup <- function(title, content,
                      placement=c('right', 'top', 'left', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual')) {
  tagList(
    singleton(
      tags$head(
        tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })")
      )
    ),
    tags$a(
      href = "#", class = "btn btn-mini", `data-toggle` = "popover",
      title = title, `data-content` = content, `data-animation` = TRUE,
      `data-placement` = match.arg(placement, several.ok=TRUE)[1],
      `data-trigger` = match.arg(trigger, several.ok=TRUE)[1],
      
      tags$i(class="fa fa-question-circle")
    )
  )
}

infoPopup <- function(title, content,
                      placement=c('right', 'top', 'left', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual')) {
  tagList(
    singleton(
      tags$head(
        tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })")
      )
    ),
    tags$a(
      href = "#", class = "btn btn-mini", `data-toggle` = "popover",
      title = title, `data-content` = content, `data-animation` = TRUE,
      `data-placement` = match.arg(placement, several.ok=TRUE)[1],
      `data-trigger` = match.arg(trigger, several.ok=TRUE)[1],
      
      tags$i(class="fa fa-info-circle")
    )
  )
}

source("HeatMap2.txt")

mycorrection<-function(dat,alpha,correction){
  
  k<-1:length(dat)
  index<-alpha*k/length(dat)
  sort.dat<-dat[order(dat)]
  
  if(correction=="FDR") maxk<-ifelse( length(which(sort.dat<=index))==0,0,max(which(sort.dat<=index)))
  if(correction=="BONF") maxk<-ifelse( length(which(sort.dat<=(alpha/length(sort.dat))))==0,0,max(which(sort.dat<=(alpha/length(sort.dat)))))
  if(correction=="RAW") maxk<-max(which(sort.dat<=alpha))
  
  percentage<-maxk/length(k)
  return(percentage)
}

flowlistmaker <- function(dat, var_name, alpha = 0.05, method = "FDR"){
  
  var_name_1 <- paste0("P.Value for ", var_name)
  pvalues <- dat[, var_name_1]
  sorted_p <- pvalues[order(pvalues)]
  index3 <- grep(var_name, names(dat), fixed = TRUE)
  
  #Fixing bug if contrast are written very similarly where
  #grep can't tell the difference        
  allpvalnames_ind <- grep("P.Value", names(dat)[index3])
  pvalnames <- names(dat)[index3][allpvalnames_ind]
  #log10pval_index <- grep("Sign_neg_log10", names(dat)[index3])
  #ind<-setdiff(allpvalnames_ind,log10pval_index)
  #pvalnames <- names(dat)[index3][ind]
  n_pvalnames <- length(pvalnames)
  
  ############################
  newdat <- data.frame(dat[,1], dat[, index3], p.adjust(pvalues,"fdr"), p.adjust(pvalues, "bonferroni"))
  names(newdat) <- c("Flow_Variable", names(dat)[index3], paste0(var_name,"FDR"), paste0(var_name,"BONF"))
  
  sort.dat <- newdat[order(dat[, var_name_1]),]
  if(method=="FDR"){    #k<-1:(dim(dat)[1])
                        #index<-alpha*k/dim(dat)[1]
                        maxk<-ifelse( length(which(p.adjust(sorted_p, "fdr")<=alpha))==0,0,max(which(p.adjust(sorted_p, "fdr")<=alpha)))
                        if(maxk!=0){flow_list<-sort.dat[1:maxk,]}
                        if(maxk==0){flow_list<- "No Variables Present"}       
  }
  if(method=="Raw"){       maxp <- ifelse( length(which(sorted_p<=alpha))==0,0,max(which(sorted_p<=alpha)))
                           if(maxp!=0){flow_list<-sort.dat[1:maxp,]}
                           if(maxp==0){flow_list<- "No Variables Present"}
  }
  if(method=="Bonferroni"){
                              maxp<-ifelse( length(which(p.adjust(sorted_p, "bonferroni")<=alpha))==0,0,max(which(p.adjust(sorted_p, "bonferroni")<=alpha)))
                              if(maxp!=0){flow_list<-sort.dat[1:maxp,]}
                              if(maxp==0){flow_list<- "No Variables Present"}
  }
  return(data.frame(Flow_Variable = flow_list))
}

movetolast <- function(dat, move) {
  dat[c(setdiff(names(dat), move), move)]
}

genelistmaker<-function(dat,var_name,alpha=.05,method="FDR",module_merge=FALSE){
  pvalues<-dat[,var_name]
  sorted_p<-pvalues[order(pvalues)]
  index3<-grep(substring(var_name,24),names(dat),fixed=T)
  
  #Fixing bug if contrast are written very similarly where
  #grep can't tell the difference        
  allpvalnames_ind<-grep("p.Val",names(dat)[index3])
  allpvalnames<-names(dat)[index3][allpvalnames_ind]
  log10pval_index<-grep("log10",names(dat)[index3])
  ind<-setdiff(allpvalnames_ind,log10pval_index)
  pvalnames<-names(dat)[index3][ind]
  n_pvalnames<-length(pvalnames)
  
  if(n_pvalnames>1){div<-length(index3)/n_pvalnames
  indicator1<-nchar(pvalnames)  
  trueval<-nchar(var_name)
  trueindex<-which(indicator1==trueval)
  newindex<-c()
  for(i in 1:div){
    newindex[i]<-index3[(trueindex+n_pvalnames*(i-1))]     }
  index3<-newindex  
  }
  ############################
  newdat<-cbind(dat$PROBE_ID,dat$SYMBOL,dat[,index3],p.adjust(pvalues,"fdr"),p.adjust(pvalues,"bonferroni"))
  names(newdat)<-c("PROBE_ID","SYMBOL",names(dat)[index3],paste(var_name,"FDR",sep=""),paste(var_name,"BONF",sep=""))
  
  sort.dat<-newdat[order(dat[,var_name]),]
  if(method=="FDR"){       #k<-1:(dim(dat)[1])
                           #index<-alpha*k/dim(dat)[1]
                           maxk<-ifelse( length(which(p.adjust(sorted_p,"fdr")<=alpha))==0,0,max(which(p.adjust(sorted_p,"fdr")<=alpha)))
                           if(maxk!=0){PROBE_ID<-sort.dat[1:maxk,]}
                           if(maxk==0){PROBE_ID<- "No Genes Present"}       
  }
  if(method=="Raw"){       maxp<- ifelse(length(which(sorted_p<=alpha)) == 0,0,max(which(sorted_p<=alpha)))
                           if(maxp != 0){PROBE_ID<-sort.dat[1:maxp,]}
                           if(maxp == 0){PROBE_ID<- "No Genes Present"}
  }
  if(method=="Bonferroni"){     maxp<-ifelse( length(which(p.adjust(sorted_p,"bonferroni")<=alpha))==0,0,max(which(p.adjust(sorted_p,"bonferroni")<=alpha)))
                                if(maxp!=0){PROBE_ID<-sort.dat[1:maxp,]}
                                if(maxp==0){PROBE_ID<- "No Genes Present"}
  }
  if(module_merge==TRUE){ PROBE_ID<-merge(PROBE_ID,moduleinfo,by="PROBE_ID",all.x = TRUE)
  index4<- which( !(names(PROBE_ID)%in%c("PROBE_ID","SYMBOL","Module","Module_V3","Modulev2_Annotation","Modulev3_Annotation")))   
  name2<-names(PROBE_ID)[index4]
  PROBE_ID<-cbind(PROBE_ID[,c("PROBE_ID","SYMBOL","Module","Modulev2_Annotation")],PROBE_ID[,index4])    
  names(PROBE_ID)<-c("PROBE_ID","SYMBOL","Module","Modulev2_Annotation",name2)
  PROBE_ID<-PROBE_ID[order(PROBE_ID[,var_name]),]
  }
  PROBE_ID <- data.frame(PROBE_ID)
  return(PROBE_ID)
}



getColor = function(x) {

 colorRampPalette(RColorBrewer::brewer.pal(7, x)[7:1], interpolate = "spline")
 
}

setPalettes = function(n) {
  pal = list("Spectral", "PuOr", "PiYG","RdGy","PRGn","RdBu","RdYlBu")
  sel = pal[1:n]
  lapply(sel, getColor)
}

#BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral


getLevels = function(x) {
  as.character(levels(x))
}


data.manipulate<-function(exp,des,basevariable,baselevel,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=T){
  index.samples<-which(names(exp) %in% des$columnname)

  if(format=="Probes"){ 
            index.ps<-which(names(exp) %in% c("PROBE_ID","SYMBOL"))
                               }
  if(format=="Modules"){
            index.ps<-which(names(exp) %in% c("Module"))
                                }
exp<-exp[,c(index.ps,index.samples)]

if(allsamples==FALSE){
if (longitudinal==TRUE){
            index.basevar<-which(names(des)==basevariable)
            base<-des[which(des[,index.basevar]==baselevel),]
		all_time_base_rm<-des[which(des[,index.basevar]!=baselevel),]

            index.subjectid<-which(names(des)==subjects)
		base_subject_ids<-unique(base[,index.subjectid])
		
            all_time_base_rm_subject_ids<-unique(all_time_base_rm[,index.subjectid])

		index.include<-base_subject_ids%in%all_time_base_rm_subject_ids

		#"complete_ids" is a list of subject ids that have a "basevariable" of "baselevel" as well as other "baselevel" observations

		complete_ids<-base_subject_ids[index.include]
         
		columnname<-names(exp)[-index.ps]

get.norm<-function(id,mynames,index.sid,index.basevar,baselevel,keep=TRUE){ data1<-des[which(des[,index.sid]==id),mynames]
				base_name<-data1[which(data1[,3]==baselevel),"columnname"]
                        if(keep==TRUE){
                        base_norm1<-exp[,as.character(data1$columnname)]-exp[,as.character(base_name)]
                        }
                        if(keep==FALSE){
                        
                        base_norm1<-exp[,setdiff(as.character(data1$columnname),as.character(base_name))]-exp[,as.character(base_name)]
                        if(length(setdiff(as.character(data1$columnname),as.character(base_name)))==1){
                                                  base_norm1=data.frame(base_norm1)
                                                  names(base_norm1)<-as.character(data1$columnname[which(as.character(data1$columnname)!=as.character(base_name))])}
				}
				return(base_norm1)
			     }


            base_norm<-get.norm(id=complete_ids[1],mynames=c(subjects,"columnname",basevariable),index.sid=index.subjectid,index.basevar=index.basevar,baselevel=baselevel,keep=keepbase)
		for(i in 2:(length(complete_ids))){base_norm<-cbind(base_norm,get.norm(complete_ids[i],c(subjects,"columnname",basevariable),index.subjectid,index.basevar,baselevel,keep=keepbase))}

		design_base_norm<-des[(des[,index.subjectid]%in%complete_ids),]
            if(keepbase==FALSE){design_base_norm<-design_base_norm[which(design_base_norm[,index.basevar]!=baselevel),]}
                        }



if (longitudinal==FALSE){
		index.basevar<-which(names(des)==basevariable)
            base<-des[which(des[,index.basevar]==baselevel),]
            base_rm<-des[which(des[,index.basevar]!=baselevel),]
            base_names<-base[,"columnname"]
            base_rm_names<-base_rm[,"columnname"]
             
get.norm2<-function(baseline,rest,keep=TRUE){index.base<-which(names(exp)%in%baseline)
                        index.rest<-which(names(exp) %in%rest )
				if(keep==FALSE){
				normdat<-exp[,index.rest]-apply(exp[,index.base],1,mean)
				#names(normdat)<-names(exp)[index.rest]
				}
                        if(keep==TRUE){
				#normdat<-exp[,-c(1,2)]-apply(exp[,index.base],1,mean)
                        normdat<-exp[,c(index.base,index.rest)]-apply(exp[,index.base],1,mean)
                        #names(normadat)<-names(exp[,-c(1,2)])  
                        }
                        return(normdat)
		           }
		base_norm<-get.norm2(base_names,base_rm_names,keep=keepbase)
            
            if(keepbase==FALSE){design_base_norm<-base_rm}
            if(keepbase==TRUE){design_base_norm<-des}
			}	
    }

if(allsamples==TRUE){
    myfactor<-as.vector(apply(exp[,-(1:2)],1,mean))
    base_norm<-exp[,-(1:2)]-myfactor
    design_base_norm<-des
    }
if(lg2==TRUE){base_norm=log(base_norm,2)}

#Data is normalized with an updated design file.  Normalized exp file still doesn't have
#Probe and Symbol attached will do in a moment.

#Now we rearrange the columns based on "ordernames".
# if(length(ordernames)>0){
#            data1<-design_base_norm[,c("Columnname",ordernames)]
#            inside<-paste("data1$",ordernames,sep="")
#            column_index<-eval(parse(text=paste("order(",paste(inside,collapse=","),")",sep="")))
#	      sorted.data<-base_norm[,column_index]}
#        if(length(ordernames)==0){sorted.data=base_norm}
		
#final_norm_shuff<-cbind(exp[,index.ps],sorted.data)
#names(final_norm_shuff)<-c(names(exp[,index.ps]),names(sorted.data))

#final_norm_shuff<-as.matrix(sorted.data)
 final_norm_shuff<-as.matrix(base_norm)
 colnames(final_norm_shuff)<-names(base_norm)
 if(format=="Probes"){rownames(final_norm_shuff)<-exp$SYMBOL}
 if(format=="Modules"){rownames(final_norm_shuff)<-exp$Module}


return(list(heatexp=final_norm_shuff,heatdes=design_base_norm))
}


moduleinfo1<-read.table("module_v2only_complete.txt",header=T,sep=",")
modnames<-unique(moduleinfo1$Module)
y<-unlist(strsplit(as.character(modnames),".",fixed=TRUE))
y<-matrix(y,byrow=T,nrow=260,ncol=2)
sortmodnames<-as.character(modnames)[order(y[,1],y[,2])]
Modulelist<-list(as.character(moduleinfo1$PROBE_ID[which(moduleinfo1$Module==sortmodnames[1])]),as.character(moduleinfo1$PROBE_ID[which(moduleinfo1$Module==sortmodnames[2])]))
for(i in 3:length(sortmodnames)){Modulelist<-c(Modulelist,list(as.character(moduleinfo1$PROBE_ID[which(moduleinfo1$Module==sortmodnames[i])])))}
names(Modulelist)<-sortmodnames

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  #library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



