if (!require("ggplot2")){
  install.packages('ggplot2')
}
if (!require("svglite")){
  install.packages('svglite')
}

library(ggplot2)
library(svglite)

#' Given a Formatted Distance to TSS TSV file Create Stacked Bargraphs
#' Author: Benjamin Tudor Price
#' 
#' === Version History ===
#' 
#'  1.0 - 09/14/2023
#'    - Written using code provided in Clayton's tutorial.
#' 
#' === Usage ===
#' 
#' $ Rscript make_distance_to_tss_barplot.R <path_to_tss_info_file> <out.image.svg> identifiers...
#' 
#' Where each identifier maps to exactly 5 entries of the group column in the
#' input. These will form one bar on the barplot, in the order they are given in identifiers.
#' 
#' === Usage Format Example ===
#' $ cat tss_info.tsv
#' 
#' Group	DistanceToTSS	Peak.Fraction
#' lost.D81WT_vs_D81A	0-100	0.0256248687250578
#' lost.D81WT_vs_D81A	100-1000	0.0593362738920395
#' lost.D81WT_vs_D81A	1000-10000	0.210670027305188
#' lost.D81WT_vs_D81A	10000-100000	0.454316320100819
#' lost.D81WT_vs_D81A	gt100000	0.250052509976896
#' gained.D81WT_vs_D81A	0-100	0.0216487870880448
#' gained.D81WT_vs_D81A	100-1000	0.0581811152991205
#' gained.D81WT_vs_D81A	1000-10000	0.233207693051126
#' gained.D81WT_vs_D81A	10000-100000	0.482941915531072
#' gained.D81WT_vs_D81A	gt100000	0.204020489030637
#' 
#' $ Rscript make_distance_to_tss_barplot.R tss_info.tsv out.svg gained lost
#' 
#' === Pipeline Flow ===
#' 
#' peaks bedfiles (often subset in an interesting way) -> $soft/addNearestGeneToBED.pl 
#' -> $soft/getDistanceToTSSmatrix.pl -> make_distance_to_tss_barplot.R -> barplots
#' 
#' 

# Parse args and check inputs.
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3){
  exit('Expected at least 3 arguments.')
}

tss_tsv_file <-  args[1]
out_svg <- args[2]
identifiers <- args[3:length(args)]

if (! file.exists(tss_tsv_file)){
  exit(paste0(c(tss_tsv_file, 'does not exist.'), collapse=' '))
}

x <- read.table(tss_tsv_file, sep="\t", header=T)
x_sorted <- x[-(1:nrow(x)),]

# Map identifier to rows.
for (identifier in identifiers){
  matching <- grepl(identifier, x$Group)
  
  # Make sure that the identifier maps to exactly 5 rows.
  if (sum(matching) > 5){
    stop(paste0(c('Identifier:', identifier, 'maps to too many columns.'), 
                collapse = ' '))
  } else if (sum(matching) < 5){
    stop(paste0(c('Identifier:', identifier, 'maps to too few columns.'), 
                collapse = ' '))
  }
  
  # Substitute New Name for old name, add to sorted dataframe.
  x$Group[matching] <- identifier
  x_sorted <- rbind(x_sorted, x[matching,])
}
x <- x_sorted

# Clayton's Code
x$Group <- factor(x$Group, levels = identifiers)
x$DistanceToTSS <- factor(x$DistanceToTSS, levels = c("gt100000","10000-100000","1000-10000","100-1000","0-100"))
plt <- ggplot(data = x, aes(x = Group, y = Peak.Fraction, fill = DistanceToTSS)) + 
  geom_bar(stat = "identity") +  
  scale_fill_manual(values=c("red4", "red3", "gray70","royalblue3","royalblue4")) +
  theme(axis.line = element_line(colour = "black"))
ggsave(out_svg, plt)



