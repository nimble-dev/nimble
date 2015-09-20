make_pages<-function(mods,allModels,ncol,dir){
    curDir <- getwd()
    outputDir <- file.path(dir, 'html')
    if(!file.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
    setwd(outputDir)
    on.exit(setwd(curDir))
    make_main(allModels)
    for(i in 1:length(mods)){
### create timing table 
        vars=mods[[i]][[1]]
        time=mods[[i]][[2]]
        
        df=data.frame()
        for(j in 1:length(dimnames(vars)[[3]])){
            temp=as.data.frame(vars[,,j])
            temp$size_time=temp$effectiveSize/time$time
            temp$method=rownames(temp)
            temp$var=rep(dimnames(vars)[[3]][j],length(rownames(temp))) 
            df=rbind(df,temp)
        }
### create factors to order facet plots
        df$var=factor(df$var,levels=unique(as.vector(factor(df$var))))
        
        img_h=floor(dim(vars)[3]*4.5/3)
### effective size faceted plots
        p=ggplot(df,aes(x=method,y=size_time,fill=method))+
            geom_bar(position=position_dodge(),stat='identity')+
                ggtitle("Effective Size per Time")+
                    facet_wrap(~ var,ncol=ncol,scales='free')+
                        ggsave(paste(allModels[i],'_effectiveSize.jpg',sep=''), height = img_h, width = 12,limitsize=F)
        
### raw data scatter plot
        p<-ggplot(df,aes(x=method,y=mean))+
            geom_point(aes(colour=method,size=1))+
                ggtitle("Posterior Statistics")+
                    guides(size=F,colour=F)+
                        geom_point(aes(x=method,y=median,size=1),shape=4)+
                            facet_wrap(~ var,ncol=ncol,scales='free')+
                                geom_errorbar(aes(ymax=CI95_upp,ymin=CI95_low),width=.25)+
                                    ggsave(paste(allModels[i],'_rawData.jpg',sep=""),height=img_h,width=12,limitsize=F)
                                        #alphabetically?
        
                                        #creates example html
        make_example_html(allModels[i],time,img_h)
    }
    invisible(NULL)
}
