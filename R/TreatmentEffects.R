
RegModel<-function(data, treats="", DV="", model=c("linear", "logit"),
                  controls="",subgroups=""){

    if(controls[1]!=""){
        RHS<-paste0(controls, collapse="+")
        RHS<-paste0(treats,"+",RHS)
    }else{
        RHS<-treats
    }

    if(subgroups[1]!=""){
        groups<-paste0(subgroups, collapse="+")
        RHS<-gsub(treats, paste0(treats, "*",groups),RHS)
    }

    form<-as.formula(paste0(DV, "~", RHS))

    if(model=="linear"){
        reg<-lm(form, data=data)
    } else if(model=="logit"){
        reg<-glm(form, data=data, family = "binomial")
    }else stop("Please specify model. Options are 'linear' or 'logit'.")

    return(reg)
}

SimModel<-function(data, reg, treats="", subgroups="",sims=10000,controlgroup=T, comb="comb"){
        treatlist<-list()
        treatlist[[treats]]<-data[treats] %>% unlist %>% unique

        if(subgroups[1]!=""){
                for(i in 1:length(subgroups)){
                        treatlist[[1+i]]<-data[subgroups[i]] %>% unlist %>% unique
                        names(treatlist)[[1+i]]<-subgroups[i]
                }
        }

        #Run simulations
        s<-clarify::sim(reg, n=sims)
        est <- sim_setx(s, x = treatlist,
                        verbose = TRUE)

        #get estimated bootstraps
        boots<-est %>% data.frame

        return(boots)
}

FirstDifference<-function(boots, comb, treats, subgroups){
        FD<-list()

        if(subgroups[1]==""){
                if(comb=="comb"){
                        combs<-combinations(n=ncol(boots), r=2)
                } else if(comb=="perm"){
                        combs<-permutations(n=ncol(boots), r=2)
                }

                for(i in 1:nrow(combs)){
                        FD[[i]]<-boots[combs[i,1]]- boots[combs[i,2]]
                        names(FD[[i]])<- paste0(names(boots[combs[i,1]]),
                                                      " - ",
                                                      names(boots[combs[i,2]])) %>%
                                gsub(treats, "",.) %>%
                                gsub("[.]","",.)
                }
                FD<-do.call(cbind, FD)
        }

        if(subgroups[1]!=""){
                sublevel<-data[subgroups] %>% unlist %>% unique
                for(n in 1:length(sublevel)){
                        levelboots<-select(boots, contains(sublevel[n]))

                        if(comb=="comb"){
                                combs<-combinations(n=ncol(levelboots), r=2)
                        } else if(comb=="perm"){
                                combs<-permutations(n=ncol(levelboots), r=2)
                        }

                        levelFD<-list()
                        for(i in 1:nrow(combs)){
                                levelFD[[i]]<-levelboots[combs[i,1]]- levelboots[combs[i,2]]
                                names(levelFD[[i]])<- paste0(names(levelboots[combs[i,1]]),
                                                                   " - ",
                                                                   names(levelboots[combs[i,2]])) %>%
                                        gsub(treats, "",.) %>%
                                        gsub("[.]","",.) %>%
                                        gsub(subgroups, "", .) %>%
                                        gsub(sublevel[n], "",.)
                        }
                        levelFD<-do.call(cbind, levelFD)
                        names(levelFD)<-paste0(sublevel, ": ",names(levelFD))
                        FD[[n]]<-levelFD
                }
                FD<-do.call(cbind, FD)
        }
        return(FD)
}

SummarizeData<-function(FD, treats, DV, subgroups){

        #get confidence intervals for treatment effects
        CI95<-sapply(FD, quantile, c(.025, .5, .975)) %>% t
        CI85<-sapply(FD, quantile, c(.075, .5, .925)) %>% t

        if(subgroups!=""){
                ph<-strsplit(rownames(CI95), "[:]") %>% data.frame %>% t %>% trimws()
                CI95<-cbind(data.frame(CI95), data.frame(ph))
                CI95$confint<-"95% Interval"
                names(CI95)<-c("min", "med", "max", "group", "treat", "conf")
                highsignif<-(CI95$max<0 | CI95$min>0)

                ph<-strsplit(rownames(CI85), "[:]") %>% data.frame %>% t %>% trimws()
                CI85<-cbind(data.frame(CI85), data.frame(ph))
                CI85$confint<-"85% Interval"
                names(CI85)<-c("min", "med", "max", "group", "treat", "conf")
                lowsignif<-(CI85$max<0 | CI85$min>0)
        } else{
                ph<-rownames(CI95) %>% trimws()
                CI95<-cbind(data.frame(CI95), data.frame(ph))
                CI95$confint<-"95% Interval"
                names(CI95)<-c("min", "med", "max",  "treat", "conf")
                highsignif<-(CI95$max<0 | CI95$min>0)

                ph<-rownames(CI85) %>% trimws()
                CI85<-cbind(data.frame(CI85), data.frame(ph))
                CI85$confint<-"85% Interval"
                names(CI85)<-c("min", "med", "max",  "treat", "conf")
                lowsignif<-(CI85$max<0 | CI85$min>0)

        }
        allest<-rbind(CI95, CI85)

        signif<-ifelse(highsignif,
                       "Significant at 95%",
                        ifelse(lowsignif,"Significant at 85%","Not significant"))

        allest$signif_dif<-signif
        allest$DV<-DV

        return(allest)
}

PlotEffects<-function(allest, treats="", DV="", model=c("linear", "logit"),
                      controls="",subgroups="",sims=10000){
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
        }
        return(gg)
}

TreatmentEffects<-function(data, treats="", DV="", model=c("linear", "logit"),
                           controls="",subgroups="",sims=10000, order=FALSE,
                           comb="perm"){
        library(stringr)
        library(Hmisc)
        library(ggplot2)
        library(clarify)
        library(tidyverse)
        library(RColorBrewer)
        library(gtools)
        library(doParallel)
        library(foreach)

        start<-Sys.time()
        # perform the regression
        reg<-RegModel(data,treats,DV,model,controls,subgroups)
        boots<-SimModel(data, reg, treats, subgroups,sims,comb)
        FD<-FirstDifference(boots, comb, treats, subgroups)
        clean_df<-SummarizeData(FD, treats, DV, subgroups)
        # plot the results of the simulations
        gg<-PlotEffects(clean_df,treats,DV,model,controls,subgroups,sims)
        #detail specification so people can easily report what they did
        specification<-c("treats", "DV", "model", "controls", "subgroups","sims", "controlgroup","function",
                         treats, DV, model, paste0(controls, collapse=", "), subgroups, sims, controlgroup,"OneClarifyComb") %>% matrix(ncol=2, nrow=8)
        colnames(specification)<-c("argument", "user-specified")

        out<-list("reg"=reg, "clean_df"=clean_df, "specification"=specification, "gg"=gg)

        end<-Sys.time()

        return(out)
}
