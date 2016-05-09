# wrappers for commonly used visualization tools

sim_hist <- function(counts, observed, xlab = "Simulated Outcomes", main = "Comparing Observation to Simulation", col = "cyan"){

  # makes a histogram of simulation outcomes and plots observed counts as dotted black line

  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  left_buffer = min(min(counts), observed)
  breaks = seq(-0.5, max(max(counts), observed) + 0.5 + left_buffer/4, 1)
  h = hist(counts, xlab=xlab, main=main, breaks = breaks, col=col, xaxt="n")
  axis(side=1, at = breaks+0.5, labels=breaks + 0.5)
  abline(v=observed, col="black", lty=3, lw=5)
}

density_hist_plot <- function(df, xfactor, fillfactor, title, xlabel, legend_title = "Data Set", binwidth=1){
  # xfactor and fillfactor should be factors of df
  ggplot(df, aes_string(x = xfactor, fill = fillfactor)) +
    geom_bar(binwidth = 1, position = "dodge", aes(y = ..density..)) +
    geom_density(alpha = 0.50, kernel = "gaussian") +
    ggtitle(title) +
    xlab(xlabel) +
    ylab("Density") +
    guides(fill = guide_legend(override.aes = list(colour = NULL))) +
    scale_fill_discrete(name = legend_title) +
    theme(plot.title = element_text(size = 16),
          axis.text = element_text(size = 14))
}

poisson_bar_ggplot <- function(p, obs, main){

  p = p[1:(obs+51)]
  counts = unlist(mapply(function(count, times) rep(count, times), seq(0, obs+50), round(1000000*p)))
  sim_hist(counts, obs, main = main)

}
