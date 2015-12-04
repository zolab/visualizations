####################################################
### Multiple histogram function
####################################################

mhist <- function(dt, names, results, outcome, file) {

    # plotting multiple alive:death comparisons per page   
    dt[, outcome := as.factor(eval(parse(text=outcome)))]  # outcome var
    color_palette <- colorRampPalette(c('#003399', '#FF9900'))(n = 2)
    pct.round <- function(x) .bincode(x, sort(quantile(x, (0:100)/100)), include.lowest = TRUE)

     nm.input <- unlist(as.list(unique(dt[, names, with = F])))
    nm.unique <- unique(unlist(lapply(nm.input, function(x) {
                split <- strsplit(x, "_")
                last <- lapply(split, function(y) y[length(y)])
                unique(unlist(last))
    })))
    
    nm <- list()
    nm_count <- list()
    nm.new <-  paste0("_",nm.unique)

    for (i in 1:length(nm.unique)) {
        # ordering by name (rather than type)
        nm[[i]] <- grep(paste0(nm.new[[i]], "$"), nm.input, value = T)
        nm_count[[i]] = length(nm[[i]])
    }
    
    nm <- unlist(nm)

   plots <- lapply(nm, function(x) {
            toplot <- copy(dt)[get(names) == x]
            toplot[, result := eval(parse(text=results))]  # assign result to new col
           toplot[, pctile := pct.round(result)]  # result pctile
           toplot[pctile > 99, result := max(toplot[pctile <= 99, result])]  # set top 1 pctile to max top 2 pctile
           # toplot[pctile < 2, result := min(toplot[pctile >= 2, result])]  # set top 1 pctile to max top 2 pctile
            hlabs.d <- 
            ggplot(data = toplot, aes(x = result, fill = outcome)) + 
            geom_histogram(alpha=.7, position='identity', colour='black', 
                                aes(y=(..density..))) +
            scale_fill_manual(values=color_palette, guide = F) +
            theme(plot.background = element_blank(), panel.background = element_blank(), axis.text=element_text(size=8), axis.title=element_text(size=10, face="bold"), axis.text.y = element_blank(), axis.ticks = element_blank()) +
            xlab(x) +
            ylab("")
            
            return(hlabs.d)
    })

    temp1 <- 1
    temp2 <- nm_count[[1]]
    iter <- length(nm.unique) - 1

    pdf(file, width = 12)
    for (k in 1:iter) {
        onepage <- plots[temp1:temp2]
        do.call(grid.arrange,onepage)
        temp1 <- temp2 + 1
        temp2 <- temp2 + nm_count[[k+1]]
     }
    
    dev.off()

}
####################################################
