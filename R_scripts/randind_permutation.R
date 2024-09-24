library(rstudioapi)

# Get the path of the script's directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set working directory one level above the script's directory
setwd(dirname(script_dir))


# Define the directory path
res_path <- "results/R_scripts_plots/"
data_path <- "data/processed_data/"

# Check if the directory exists; if not, create it
if (!file.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
  cat("Directory created:", res_path, "\n")
} else {
  cat("Directory already exists:", res_path, "\n")
}



library(ggplot2)


randindexdata = read.csv(paste0(data_path,'randindexpermu.csv'), sep = '\t', row.names = 1)


# from https://rdrr.io/github/fruciano/GeometricMorphometricsMix/src/R/adjRand_test.R
adjRand_test=function(A, B, perm=1000) {
  if (length(A)!=length(B)) { stop("A and B should have the same length") }
  # Make sure that the two groups of partitions have the same length
  
  ARIo=mclust::adjustedRandIndex(A, B)
  # Observed adjusted Rand index
  print(ARIo)
  
  Aperm=lapply(seq_len(perm), function(X) sample(A, length(A), replace=FALSE))
  Bperm=lapply(seq_len(perm), function(X) sample(B, length(B), replace=FALSE))
  # Generate permuted samples
  
  ARIperm=unlist(lapply(seq_len(perm),
                        function(i) mclust::adjustedRandIndex(Aperm[[i]], Bperm[[i]])))
  # Compute adjusted Rand index for the permuted samples
  
  m=mean(ARIperm); v=var(ARIperm)
  # compute mean and variance of the permuted samples
  
  NARI=(ARIperm-m)/sqrt(v)
  # compute NARI (normalized ARI)
  
  NARIo=(ARIo-m)/sqrt(v)
  # compute observed NARI
  
  p_value=(length(which(NARI>NARIo))+1)/(perm+1)
  # Compute p value as proportion of permuted NARI larger than the observed
  print(p_value)
  
  return(ARIperm)
}

set.seed(42)
values <- adjRand_test(randindexdata$to, randindexdata$from)
# Convert the array to a data frame and calculate the frequency of each value
df <- as.data.frame(table(values))

# Rename the columns for clarity
colnames(df) <- c("Value", "Frequency")
# Convert Value to numeric
df$Value <- as.numeric(as.character(df$Value))
df$AboveThreshold <- ifelse(df$Value > 0.0163, "Above", "Below")

# Create the bar plot
bar_plot <- ggplot(df, aes(x = Value, y = Frequency, fill = AboveThreshold)) +
  geom_bar(stat = "identity") +   # Use 'identity' to plot actual values
  scale_fill_manual(values = c("Below" = "gray", "Above" = "red")) +
  labs(title = "", x = "Adjusted Rand Index", y = "Frequency") +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_vline(xintercept = 0.0163, linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = 0.0163, linetype = "dashed", color = "red", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) +
  annotate("text", x = 0.05, y = 25, label = "Observed Adjusted Rand index = 0.0163\nP-value = 0.1309", angle = 0, vjust = 1.5, color = "red", size = 5) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    axis.text.x = element_text(size = 14),   # Increase x-axis text size
    axis.text.y = element_text(size = 12),    # Increase y-axis text size
    panel.grid.major = element_blank(),      # Remove major gridlines
    panel.grid.minor = element_blank(),      # Remove minor gridlines
  )

ggsave(paste0(res_path,'randindexpermutations.png'), bar_plot, scale = 1, width = 10, height = 8)
