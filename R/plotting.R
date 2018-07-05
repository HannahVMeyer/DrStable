plotStability <- function(stability, colorPalette="YlGn", saveTo=NULL) {
    p_pass <- ggplot(stability$pass_long, 
                     aes(x=as.factor(threshold), y=componentsPass))
    p_pass <- p_pass +
        geom_boxplot(outlier.colour=NA) + 
        geom_jitter(width = 0.2, size=0.5) +
        labs(x = "Correlation threshold for filtering", 
             y = "Number of components") +
        geom_text(data = stability$medians_pass, aes(x = as.factor(threshold), 
                                           y = componentsPass, 
                                           label= componentsPass), 
                  nudge_x =-0.7, color="red") +
        theme_bw()
    
    p_maxcor <- ggplot(stability$maxcor_long, 
                       aes(x=as.factor(component), y=abs(correlation), 
                                        fill=as.factor(threshold)))
    p_maxcor <- p_maxcor + 
        geom_boxplot(outlier.colour=NA) +
        scale_fill_brewer(type="seq", palette=colorPalette) +
        guides(fill=guide_legend(title="Threshold")) +
        labs(x = "Components", y = "Spearman correlation") +
        theme_bw() +
        theme(axis.text.x = element_text(angle=90))
    
    if (!is.null(saveTo)) {
        ggsave(filename=paste(saveTo, "/CorrelationComponentsPassThreshold_", 
                              method, ".pdf", sep=""), 
           plot=p_pass, height=12, width=16)
        ggsave(filename=paste(saveTo, "/CorrelationThresholdPerComponent_", 
                              method, ".pdf", sep=""), 
           plot=p_maxcor, height=12, width=16)
    }
    return(list(pass=p_pass, corr=p_maxcor))
}
    
