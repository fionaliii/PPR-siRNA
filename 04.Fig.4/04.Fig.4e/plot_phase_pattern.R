# Load necessary libraries

library(dplyr)
library(ggplot2)
library(gridExtra)

# Define the Plot_Score function for data processing and phase pattern visualization
generate_phase_pattern <- function(count_file, score_file, region, pattern = 21, osname) {
  
  # determine the plot region
  plot_start <- region %>% strsplit("[:-]") %>% .[[1]] %>% .[2] %>% as.numeric
  plot_end <- region %>% strsplit("[:-]") %>% .[[1]] %>% .[3] %>% as.numeric

  plot_site <- seq(plot_start, plot_end, by = 1)
  pattern_temp <- data.frame(Coordinate = plot_site, Abundance = 0)
 
  # organize abundance data
  fsize <- file.info(count_file)$size
  if (fsize > 0) {
    plot_data <- read.delim2(count_file, header = FALSE, sep = " ")
    plot_data <- apply(plot_data[, c(2, 3)], c(1, 2), as.numeric)
    plot_data <- data.frame(Coordinate = plot_data[, 1], Abundance = plot_data[, 2])

    # identify the site with highest abundant reads
    biggest_score <- max(plot_data$Abundance)
    phase_at <- plot_data$Coordinate[which(plot_data$Abundance %in% biggest_score)][1]
    
    match_ids <- match(plot_data$Coordinate, plot_site)
    pattern_temp[match_ids[!is.na(match_ids)], "Abundance"] <- plot_data$Abundance[!is.na(match_ids)]
  } else {
    phase_at <- as.numeric(pattern_temp[1, 1])
  }

  pattern_temp <- data.frame(pattern_temp)
  pattern_y_value <- max(pattern_temp$Abundance)
 
  # define the y-axis limits 
  if (pattern_y_value > 10) {
    dis <- round((pattern_y_value + 4) / 5)
    pattern_y_break <- round(seq(0, dis * 5, by = dis))
  } else {
    pattern_y_break <- round(seq(0, pattern_y_value, by = 1))
  }
  
  # identify phase registered main site
  diff_start <- phase_at - plot_start
  if (diff_start < pattern) {
    pre_start <- plot_start
  } else {
    pre_start <- seq(phase_at, plot_start, by = -(pattern))
  }
  
  diff_end <- plot_end - phase_at
  if (diff_end < pattern) {
    pro_site <- plot_end
  } else {
    pro_site <- seq(phase_at, plot_end, by = pattern)
  }
  
  x_break <- sort(unique(c(pre_start, pro_site)))

  # organize phasing score data
  score_temp <- data.frame(Coordinate = plot_site, Abundance = 0)
  fsize <- file.info(score_file)$size

  if (fsize > 0) {
    score_data <- read.delim2(score_file, header = TRUE, sep = "\t")
    score_data <- apply(score_data[, c(2, 12)], c(1, 2), as.numeric)
    score_data <- data.frame(Coordinate = score_data[, 1], Abundance = score_data[, 2])
    match_score_ids <- match(score_data$Coordinate, plot_site)
    score_temp$Abundance[match_score_ids[!is.na(match_score_ids)]] <- score_data$Abundance[!is.na(match_score_ids)]
  }
	
  score_temp = data.frame(score_temp)
  score_y_value = max(score_temp$Abundance)
  if(score_y_value > 10){
    dis = round((score_y_value + 4) / 5)
    score_y_break = round(seq(0, dis * 5, by = dis)) 
  }else{
    score_y_break = round(seq(0, score_y_value, by = 1))
  }
  score_temp$Class = "grey"
   
  # identify multiple phase pattern
  maxat = score_temp$Coordinate[max(score_temp$Abundance) == score_temp$Abundance][1]
  pstart = seq(maxat, min(score_temp$Coordinate), by = -(pattern))
  pend = seq(maxat, max(score_temp$Coordinate), by = (pattern))
  pstartend = sort(unique(c(pstart, pend)))
  pstartend = unlist(sapply(pstartend, function(x){
    aa = c((x - 2 ) : (x + 2))
    return(aa)
  }))
  pstartend = unique(pstartend)
  score_temp$Class[match(pstartend, score_temp$Coordinate)] = "red"
    
  pstartend = sort(unique(c(pstart, pend)))
  pstartend = unlist(sapply(pstartend, function(x){
    aa = c((x - 10 - 2) : (x - 10 + 2))
    return(aa)
  }))
  pstartend = unique(pstartend)
  score_temp$Class[match(pstartend, score_temp$Coordinate)] = "blue"
    
  pstartend = sort(unique(c(pstart, pend)))
  pstartend = unlist(sapply(pstartend, function(x){
    aa = c((x - 5 - 2 ) : (x - 5 + 2 ))
    return(aa)
  }))
  pstartend = unique(pstartend)
  score_temp$Class[match(pstartend, score_temp$Coordinate)] = "orange"
  
  # plot abundance data  
  tmp = data.frame(pattern_temp)
  p1 <- ggplot() +
         geom_line(data = tmp, aes(x = Coordinate, y = Abundance), color = "#05017B", size = 0.5)+ 
	 theme_minimal()+
	 xlab("")+
	 ylab("Abundance")+
	 scale_x_continuous(breaks= x_break, limits = c(plot_start, plot_end), expand=c(0, 0))+
	 scale_y_continuous(breaks= pattern_y_break, limits = c(0, max(pattern_y_break)))+
	 theme(plot.margin = grid::unit(c(2, 3, 2, 3),'cm'))+
	 theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5, colour = "black"))+
         theme(axis.text.y = element_text(size = 10, colour = "black"))+
	 theme(axis.line = element_line(size = 0.5, colour = "#000000"))+
	 theme(panel.grid.major.x = element_line(colour = "grey90", size = 0.5))+
	 theme(panel.grid.major.y = element_line(colour = "#E2BB76", size = 0.5))+
	 theme(panel.grid.minor = element_line(colour = "#FFFFFF", size = 0.5)) 
   
  # define phase score patterns with different colors
  tmp = data.frame(score_temp)
  tmp_red = tmp_blue = tmp_grey = tmp_orange = tmp
  tmp_red$Abundance [ !tmp_red$Class %in% "red" ] = 0
  tmp_red$Class [ !tmp_red$Class %in% "red" ] = "red"
  tmp_blue$Abundance [ !tmp_blue$Class %in% "blue" ] = 0
  tmp_blue$Class [ !tmp_blue$Class %in% "blue" ] = "blue"
  tmp_orange$Abundance [ !tmp_orange$Class %in% "orange" ] = 0
  tmp_orange$Class [ !tmp_orange$Class %in% "orange" ] = "orange"
  tmp_grey$Abundance [ !tmp_grey$Class %in% c("grey") ] = 0 
  tmp_grey$Class [ !tmp_grey$Class %in% c("grey") ] = "grey"
  tmp = rbind(tmp_red, tmp_blue)
  tmp = rbind(tmp, tmp_orange)
  tmp = rbind(tmp, tmp_grey)
     
  tmp$Abundance = as.numeric(tmp$Abundance)
  tmp$Coordinate = as.numeric(tmp$Coordinate)
  tmp$Class = factor(tmp$Class, levels = c("red", "blue", "orange", "grey"))

  # plot phasing pattern
  p2 <- ggplot() +
         geom_line(data = tmp, aes(x = Coordinate, y = Abundance, group = Class, colour = Class), size = 0.5)+
         xlab("")+
         ylab("Phasing score")+
         scale_color_manual(values=c("red" = 'red',"blue" = 'blue',"orange" = '#AC28F6',"grey" = '#7EC5B8')) + #F0BFC6
         scale_x_continuous(breaks = x_break, limits = c(plot_start, plot_end), expand = c(0, 0))+
         scale_y_continuous(breaks = score_y_break, limits = c(0, max(score_y_break)), expand = c(0, 0))+
         theme_minimal()+
         theme(plot.margin = grid::unit(c(2, 3, 2, 3), 'cm'))+
         theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5, colour = "black"))+
         theme(axis.text.y = element_text(size = 10, colour = "black"))+
         theme(axis.line = element_line(size = 0.5, colour = "#000000"))+
         theme(panel.grid.major.x = element_line(colour = "grey90", size = 0.5)) + ##E2BB76
         theme(panel.grid.major.y = element_line(colour = "#E2BB76", size = 0.5)) +
         theme(panel.grid.minor = element_line(colour = "#FFFFFF", size = 0.5)) +
         theme(legend.position = "none")


   write.table(tmp,"./test.txt") 
   save_pdf = paste(osname,"_pattern_score.v2.pdf", sep="")
   ggsave(file = save_pdf, arrangeGrob(p1, p2), width = 12, height = 6.5)
	
}	

#####################################################################################

# Usage function
print_usage <- function() {
  cat("Usage:\n")
  cat("  Rscript plot_phase_pattern.R <inbam> <locus> <phase> <fname>\n\n")
  cat("Arguments:\n")
  cat("  <inbam>: Input BAM file.\n")
  cat("  <locus>: Genomic locus or region of interest (e.g., chr1:1000-2000).\n")
  cat("  <phase>: Phase value as a numeric parameter.\n")
  cat("  <fname>: Output file name prefix.\n\n")
  cat("Description:\n")
  cat("  This script analyzes the phasing score within a specific genome locus. It takes the following parameters:\n")
  cat("  - <inbam>: Input BAM file containing sRNA data.\n")
  cat("  - <locus>: Genomic locus or region of interest in the format 'chr:start-end'.\n")
  cat("  - <phase>: Specifying the registered phase length.\n")
  cat("  - <fname>: Prefix for the output file names.\n\n")
  cat("Example:\n")
  cat("  Rscript plot_phase_pattern.R input.bam chr1:1000-2000 21 output_prefix\n")
}

# Check for command line arguments
options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  cat("Error: Incorrect number of arguments.\n")
  print_usage()
  quit(status = 1)
}

# set parameters
options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
inbam <- args[1]
locus <- args[2]
phase <- as.numeric(args[3])
fname <- args[4]

# extract read count
cat("STEP1. Extract read count\n")
command <- paste("samtools view ", inbam, " ", locus, "| awk -F'_' 'length($1) == 21' | awk -F'[\t_]' '{if($5 == 16){print $6\"_\"$7+2\"\t\"1/$3};if($5 == 0){print $6\"_\"$7\"\t\"1/$3}}' | awk '{sum[$1]+=$2}END{for(key in sum){print key,sum[key]}}' | awk -F'_' '{print $1\" \"$2}' | sort +0 -1 +1n -2 > locus_sRNA_rc.txt", sep = "")
system(command)

# organize count data
count <- read.delim2("locus_sRNA_rc.txt", sep = " ", header = FALSE)
plot_start <- locus %>% strsplit("[:-]") %>% .[[1]] %>% .[2] %>% as.numeric
plot_end <- locus %>% strsplit("[:-]") %>% .[[1]] %>% .[3] %>% as.numeric
data <- data.frame(Chr = count[1, 1], Position = seq(plot_start, plot_end, by = 1), Count = 0)
mm <- match(count[, 2], data$Position)
data$Count[mm[!is.na(mm)]] <- count[, 3][which(!is.na(mm))]
write.table(data, "locus_sRNA_rc.revised.txt", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

# calculate phasing score
cat("STEP2. Calculate phasing score\n")
command <- paste("perl phasing_analysis.pl locus_sRNA_rc.revised.txt ", phase, sep = "")
system(command)

# plot phase pattern
cat("STEP3. Plot phase pattern\n")
count_file <- "locus_sRNA_rc.txt"
score_file <- "phase_score.txt"
generate_phase_pattern(count_file, score_file, region = locus, phase, fname)



