library(stringr)
library(Hmisc)
library(ggplot2)
library(clarify)
library(tidyverse)
library(RColorBrewer)
library(gtools)
library(doParallel)
library(foreach)

# note identical in both OC(C) - renamed incase of future edits to avoid namespace collision 
OneRegressionComb<-function(data, treats="", DV="", model=c("linear", "logit"), 
                  controls="",subgroups=""){
  
    if(controls[1]!=""){
        RHS<-paste0(controls, collapse="+")
        RHS<-paste0(treats,"+",RHS)
    }else{
        RHS<-treats
    }

    form<-as.formula(paste0(DV, "~", RHS))
  
    if(subgroups[1]!=""){
        groups<-paste0(subgroups, collapse="+")
        textform<-gsub("treats", paste0(DV,"~","treats*",groups),RHS)
        form<-as.formula(as.character(textform))
    }

    if(model=="linear"){
        reg<-lm(form, data=data)
    } else if(model=="logit"){
        reg<-glm(form, data=data, family = "binomial")
    }else stop("Please specify model. Options are 'linear' or 'logit'.")           
    
    return(reg)
}

# generates for one subgroup
gen_subgroup_bootcombs<-function(subtreatboots,sg,comb="comb"){
  
  if(comb=="comb"){
    combs<-combinations(n=ncol(subtreatboots), r=2) 
  } else if(comb=="perm"){
    combs<-permutations(n=ncol(subtreatboots), r=2) 
  }
  
  bootcombs<-list()
  cols_<-list() # quick fix, names do not "stick" when converted to data.frame
  for(i in 1:nrow(combs)){
    # gets the difference
    bootcombs[[i]]<-subtreatboots[combs[i,]][,1]-subtreatboots[combs[i,]][,2] 
    # creates the name
    t <- subtreatboots[combs[i,]] %>% 
      colnames %>% 
      sapply(function(x) str_match_all(x,"\\.\\.\\.\\.[a-zA-Z0-9_]*")) %>% array %>% unlist %>% gsub("[.]","",.) %>% 
      gsub("treat_", "",.)
    names(bootcombs)[i] <- paste(
      sg,
      paste0(
        t[c(1,3)],
        collapse=" - "
        ),
      sep=": "
      )
    cols_[i] <- names(bootcombs)[i] # quick fix, names do not "stick" when converted to data.frame
  }

  bootcombs <- bootcombs %>% data.frame
  colnames(bootcombs) <-cols_ # quick fix, names do not "stick" when converted to data.frame

  return(bootcombs)
}

# across all subgroups
gen_allbootcombs<-function(treatboots,data,subgroups="",comb="comb"){

    if(subgroups[1]!=""){
      allbootcombs <- list()
      
      for(sg in unique(data[subgroups])[[1]]){
        
        allbootcombs <- dplyr::select(treatboots,matches(paste0("\\.",sg,".$"))) %>% gen_subgroup_bootcombs(.,sg,comb) %>% append(allbootcombs,.)
      }
    }  
    else{
            if(comb=="comb"){
                    combs<-combinations(n=ncol(treatboots), r=2) 
            } else if(comb=="perm"){
                    combs<-permutations(n=ncol(treatboots), r=2) 
            }
            
    allbootcombs<-list()
      
    for(i in 1:nrow(combs)){
      allbootcombs[[i]]<-treatboots[combs[i,]][,1]-treatboots[combs[i,]][,2]
      names(allbootcombs)[i] <-treatboots[combs[i,]] %>% 
        colnames %>% 
        strsplit("[.]") %>% 
        sapply(function(x) x[length(x)]) %>% 
        paste0(collapse=" - ") %>% 
        gsub("treat_", "",.)
      }
    }
    
    return(allbootcombs)

}

OneSimulateComb<-function(data, reg, treats="", subgroups="",sims=10000,controlgroup=T, comb="comb", n_workers=""){

  treatlist<-list()
  treatlist[[treats]]<-unlist(unique(data[treats]))
  
  if(subgroups[1]!=""){
    for(i in 1:length(subgroups)){
      treatlist[[1+i]]<-unlist(unique(data[subgroups[i]]))
      names(treatlist)[[1+i]]<-subgroups[i]
    }
  }
  
  # parallel option containing core simulation code 
  if(n_workers[1]!=""){
    registerDoParallel(n_workers)
    est<-foreach(i=1:4)%dopar%{
      s<-clarify::sim(reg, n=ceiling(sims/n_workers))
      est<-sim_setx(s, x = treatlist,
                    verbose = TRUE)
      return(est)
    }
    est<-est %>% lapply(data.frame) %>% do.call(rbind,.)
  } else {
    s<-clarify::sim(reg, n=sims)
    est <- sim_setx(s, x = treatlist,
                    verbose = TRUE)
  }
    
  #get estimated bootstraps
  boots<-est %>% data.frame
  
  #Only Keep Treatments
  if(controlgroup==F){
    treatboots<-boots %>% dplyr::select(!contains("control"))
  }else{
    treatboots<-boots
  }
  
  allbootcombs<-gen_allbootcombs(treatboots,data,subgroups=subgroups,comb=comb)

  return(allbootcombs)
}

OneSummarizeComb<-function(allbootcombs, treats="", DV="", subgroups="",order=F){

    #get confidence intervals for treatment effects
    CI95<-sapply(allbootcombs, quantile, c(.025, .5, .975))
    CI85<-sapply(allbootcombs, quantile, c(.075, .5, .925))
    
    #clean 95% CIs
    CI95<-CI95 %>% t %>% data.frame %>% 
      mutate(treat=rownames(.)) %>% mutate(treat=gsub("[.]", " - ",treat))
    rownames(CI95)<-NULL
    colnames(CI95)<-c("min", "med", "max", "treat")
    CI95$treat<-CI95$treat %>% gsub("treat....tg_", "",.) %>% gsub("[.]", "",.) 
    CI95$conf<-"95% Interval"
    
    #clean 85% CIs
    CI85<-CI85 %>% t %>% data.frame %>% 
      mutate(treat=rownames(.)) %>% mutate(treat=gsub("[.]", " - ",treat))
    rownames(CI85)<-NULL
    colnames(CI85)<-c("min", "med", "max", "treat")
    CI85$treat<-CI85$treat %>% gsub("treat....tg_", "",.) %>% gsub("[.]", "",.)
    CI85$conf<-"85% Interval"
    
    #bind both ranges
    allest<-rbind(CI95, CI85)
    
    if(subgroups[1]!=""){

      allest[ncol(allest)+1]<-allest$treat %>% gsub(treats, "",.) %>% 
        sapply(function(x) strsplit(x, ": ")[[1]][1])
      names(allest)[ncol(allest)]<-subgroups
      
      allest$treat<-allest$treat %>% gsub(treats, "",.) %>% 
        sapply(function(x) strsplit(x, ": ")[[1]][2])
    }
    
    # modified to return different symbols for 3 significance categories
    highconf<-allest %>% filter(conf=="95% Interval")
    highsignif<-(highconf$max<0 | highconf$min>0)
    
    lowconf<-allest %>% filter(conf=="85% Interval")
    lowsignif<-(lowconf$max<0 | lowconf$min>0)
    
    signif<-ifelse(
      highsignif,
      "Significant at 95%",
      ifelse(
        lowsignif,
        "Significant at 85%",
        "Not significant"
      )
    )
    
    allest$signif_dif<-signif
    allest$DV<-DV
  
    return(allest)
}

OnePlotComb<-function(allest, treats="", DV="", model=c("linear", "logit"), 
                  controls="",subgroups="",sims=10000){
    
    
    # ANDY COMMENT: this ylim section was not working so I commented it out here and below
    # adj_scalar = 1.1
    # maxy<-allest$max %>% max*adj_scalar
    # miny<-allest$min %>% min*adj_scalar
    # medy<-c(miny, maxy) %>% mean
    # medy<-medy
    
    if(subgroups[1]==""){
      #plot
      gg<-ggplot(allest, aes(x=treat, ymin=min, y=med, ymax=max))+
        geom_pointrange(lwd=2,size=1,  aes(col=conf, pch=signif_dif))+
        geom_hline(yintercept=0, col="red", lty=2)+
        coord_flip()+
        theme_minimal()+
        scale_colour_brewer(palette ="Paired")+
        labs(y="Difference", 
             x="Treatments", 
             col="Confidence Interval",
             pch="Significance",
             title= DV,
             caption=paste0("function:", "OneClarifyComb", "; model: ", model, "; treats: ", treats, ";\n controls: ", paste0(controls, collapse=", "), "; subgroups: ", subgroups, "; sims: ", sims))
      # +
      #   ylim(miny, maxy)   
    } else{
      
      facetform<-as.formula(paste0(".~", subgroups))
      
      gg<-ggplot(allest, aes(x=treat, ymin=min, y=med, ymax=max))+
        geom_pointrange(lwd=2,size=1,  aes(col=conf, pch=signif_dif))+
        facet_wrap(facetform)+
        geom_hline(yintercept=0, col="red", lty=2)+
        coord_flip()+
        theme_minimal()+
        scale_colour_brewer(palette ="Paired")+
        labs(y="Difference", 
             x="Treatments", 
             col="Confidence Interval",
             pch="Significance",
             title= DV,
             caption=paste0("function:", "OneClarifyComb", "; model: ", model, "; treats: ", treats, ";\n controls: ", paste0(controls, collapse=", "), "; subgroups: ", subgroups, "; sims: ", sims)
             )
      # +
      #   ylim(miny, maxy)   
    }
    return(gg)
}

# cleans inputs by factorizing and replacing apostrophe, space, dash, slash, and pipe
MiniTidy<-function(obj,coerce=FALSE){
  
  if(class(obj)=="character" | coerce){
    obj <- obj %>% as.factor %>% 
      gsub("'","",.) %>% gsub(" |-|/|\\|","_",.) 
  }
  
  return(obj)
}

TidyWrapper<-function(data,treats,subgroups,controls){
  
  for(input in list(treats,subgroups)){
    if(input!=""){
      data[[input]]<-MiniTidy(data[[input]],coerce=TRUE)
    }
  }
  if (length(controls) > 1){
    for(control in controls){
      data[[control]]<-MiniTidy(data[[control]])
    }
  }
  else{
    if (controls != ""){
      data[[controls]]<-MiniTidy(data[[controls]])
    }
  }
  
  return(data)
}

OneClarifyComb<-function(data, treats="", DV="", model=c("linear", "logit"), 
                     controls="",subgroups="",sims=10000, order=FALSE,
                     controlgroup=TRUE, comb="comb", n_workers=""){
  
        start<-Sys.time()
  
        # cleans input variables to avoid picky errors
        data<-TidyWrapper(data,treats,subgroups,controls)
        # perform the regression
        reg<-OneRegressionComb(data,treats,DV,model,controls,subgroups)
        # run simulations and compile the combinations
        allbootcombs<-OneSimulateComb(data,reg,treats,subgroups,sims,controlgroup,comb,n_workers)
        # get summary statistics from the simulations 
        allest<-OneSummarizeComb(allbootcombs,treats,DV,subgroups)
        # plot the results of the simulations 
        gg<-OnePlotComb(allest,treats,DV,model,controls,subgroups,sims)
        #detail specification so people can easily report what they did
        specification<-c("treats", "DV", "model", "controls", "subgroups","sims", "controlgroup","function",
                         treats, DV, model, paste0(controls, collapse=", "), subgroups, sims, controlgroup,"OneClarifyComb") %>% matrix(ncol=2, nrow=8)
        colnames(specification)<-c("argument", "user-specified")
        
        out<-list("reg"=reg, "allest"=allest, "specification"=specification, "gg"=gg)
        
        end<-Sys.time()
        
        return(out)
}


