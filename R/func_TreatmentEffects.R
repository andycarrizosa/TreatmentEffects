
RegModel<-function(data, treats, DV, model,
                  controls,subgroups){

        if(controls[1]!=""){
                #set up model with controls
                RHS<-paste0(controls, collapse="+")
                RHS<-paste0(treats,"+",RHS)
        }else{
                #set up model without controls
                RHS<-treats
        }

        if(subgroups[1]!=""){
                #include interactions b/w treatment and subgroups if applicable
                groups<-paste0(subgroups, collapse="+")
                RHS<-gsub(treats, paste0(treats, "*",groups),RHS)
        }

        # Turn text into a formula
        form<-as.formula(paste0(DV, "~", RHS))

        if(model=="linear"){
                reg<-lm(form, data=data) # Run linear model if specified
        } else if(model=="logit"){
                reg<-glm(form, data=data, family = "binomial") # Run logit model if specified
        }else stop("Please specify model. Options are 'linear' or 'logit'.") #otherwise return error

    return(reg)
}

SimModel<-function(data, reg, treats, subgroups,sims, comb, clust){
        # Set up all conditions to vary for simulations
        # We will vary treatment conditions and if applicable we will also vary subgroup conditions
        treatlist<-list()
        treatlist[[treats]]<-data[treats] %>% unlist %>% unique

        if(subgroups[1]!=""){
                for(i in 1:length(subgroups)){
                        treatlist[[1+i]]<-data[subgroups[i]] %>% unlist %>% unique
                        names(treatlist)[[1+i]]<-subgroups[i]
                }
        }

        #Run simulations for treatment and if applicable also subgroup conditions
        if(clust<2){
                #run simulations in a single cluster
                s<-clarify::sim(reg, n=sims)
                est <- sim_setx(s, x = treatlist,
                                verbose = TRUE)
                boots<-est %>% data.frame

        }else{
                #run simulations as a parallel process through multiple clusters to speed estimation
                cl<-makeCluster(clust)
                registerDoParallel(cl)
                options(cores=clust)
                minisims<-ceiling(sims/clust)

                est<-foreach(i=1:length(clust))%dopar%{
                        library(clarify)
                        s<-clarify::sim(reg, n=minisims)
                        est <- sim_setx(s, x = treatlist,
                                        verbose = TRUE)
                        return(est)
                }

                boots<-lapply(est, data.frame) %>%
                        do.call(rbind, .)
                stopCluster(cl)
        }
        return(boots)
}

FirstDifference<-function(data,boots, comb, treats, subgroups){
        FD<-list()
        if(subgroups[1]==""){
                # take first difference between all estimated bootstrap columns to get differences between all treatments
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
                # first subset by each subgroup category, and then perform first differences of all treatments within each category
                sublevel<-data[subgroups] %>% unlist %>% unique %>% gsub(" ", ".",.) %>%
                        gsub("[-]", ".",.) %>% gsub("[+]",".",.)
                for(n in 1:length(sublevel)){
                        levelboots<-select(boots, contains(sublevel[n],ignore.case=F))

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
                                        gsub(treats, "",.)  %>%
                                        gsub(subgroups, "", .) %>%
                                        gsub(sublevel[n], "",.)%>%
                                        gsub("[.]","",.)
                        }
                        levelFD<-do.call(cbind, levelFD)
                        names(levelFD)<-paste0(sublevel[n], ": ",names(levelFD))
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
                # clean data and create a variable for each subgroup of the subgroup analysis
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
                #clean data without subgroup variable (not applicable)
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
        clean_df<-rbind(CI95, CI85)

        signif<-ifelse(highsignif,
                       "Significant at 95%",
                        ifelse(lowsignif,"Significant at 85%","Not significant"))

        clean_df$signif_dif<-signif
        clean_df$DV<-DV

        return(clean_df)
}

PlotEffects<-function(clean_df, treats, DV, model,
                      controls,subgroups,sims){
        if(subgroups[1]==""){
                #plot figure if subgroups do not exist
                gg<-ggplot(clean_df, aes(x=treat, ymin=min, y=med, ymax=max))+
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
                # plot figure if subgroups do exist
                gg<-ggplot(clean_df, aes(x=treat, ymin=min, y=med, ymax=max))+
                        geom_pointrange(lwd=2,size=1,  aes(col=conf, pch=signif_dif))+
                        facet_wrap(.~group)+
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
                           controls="",subgroups="",sims=10000,
                           comb="perm", clust=1){
        library(stringr)
        library(Hmisc)
        library(ggplot2)
        library(clarify)
        library(tidyverse)
        library(RColorBrewer)
        library(gtools)
        library(doParallel)
        library(foreach)

        reg<-RegModel(data,treats,DV,model,controls,subgroups) # run the regression model
        boots<-SimModel(data, reg, treats, subgroups,sims,comb, clust) #perform the simulations
        FD<-FirstDifference(data,boots, comb, treats, subgroups) # get differences between bootstraps of simulations
        clean_df<-SummarizeData(FD, treats, DV, subgroups) #summarize by gettintg quantiles and median of each bootstrapped first-difference
        gg<-PlotEffects(clean_df,treats,DV,model,controls,subgroups,sims) #plot simplified data

        # Return summary information about how the model was specified
        specification<-c("treats", "DV", "model", "controls", "subgroups","sims",
                         treats, DV, model, paste0(controls, collapse=", "), subgroups, sims) %>% matrix(ncol=2, nrow=6)
        colnames(specification)<-c("argument", "user-specified")

        out<-list("reg"=reg, "clean_df"=clean_df, "specification"=specification, "gg"=gg)

        return(out)
        end<-Sys.time()
}
