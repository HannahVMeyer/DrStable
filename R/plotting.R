plotStability <- function(corpass, maxcor, medians_pass, color) {
    p_pass <- ggplot(corpass, aes(x=as.factor(threshold), y=value))
    p_pass <- p_pass +
        geom_boxplot(outlier.colour=NA) + 
        geom_jitter(width = 0.2, size=0.5) +
        labs(x = "Correlation threshold for filtering", 
             y = "Number of components") +
        geom_text(data = medians_pass, aes(x = as.factor(threshold), 
                                           y = value, label= value), 
                  nudge_x =-0.7, color="red") +
        theme_bw()
    
    p_maxcor <- ggplot(maxcor, aes(x=as.factor(component), y=abs(value), 
                                        fill=as.factor(threshold)))
    p_maxcor <- p_maxcor + 
        geom_boxplot(outlier.colour=NA) +
        scale_fill_manual(
            values=color[(16-length(levels(as.factor(maxcor$threshold)))):16]) +
        guides(fill=guide_legend(title="Threshold")) +
        labs(x = "Components", y = "Spearman correlation") +
        theme_bw() +
        theme(axis.text.x = element_text(angle=90))
    
    ggsave(filename=paste(outdir, "CorrelationComponentsPass_", method, 
                          ".pdf", sep=""), 
           plot=p_pass, height=12, width=16)
    ggsave(filename=paste(outdir, "CorrelationComponentsMaxcor_", method, 
                          ".pdf", sep=""), 
           plot=p_maxcor, height=12, width=16)
}
    
