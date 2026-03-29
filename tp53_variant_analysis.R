# TP53 variant analysis - alignment scores and visualisation
# MSc project | IBAB Bangalore
# used R for the first time properly here, still getting used to ggplot

library(ggplot2)

# manually entered alignment scores from blast output
# identity % values taken from blastp results for each variant vs wildtype

variants <- c("R175H", "R248W", "R248Q", "R273H", "R273C",
              "R282W", "G245S", "V143A", "R249S", "Y220C")

alignment_scores <- c(312, 289, 301, 275, 285, 263, 318, 245, 271, 258)

identity_pct <- c(88.2, 85.4, 87.1, 83.9, 84.7,
                  81.3, 89.5, 79.6, 82.5, 80.2)

# putting it all in a dataframe
df <- data.frame(variants, alignment_scores, identity_pct)

# quick check
print(df)

# basic summary
cat("mean alignment score:", mean(alignment_scores), "\n")
cat("lowest score:", min(alignment_scores), "->", variants[which.min(alignment_scores)], "\n")
cat("highest score:", max(alignment_scores), "->", variants[which.max(alignment_scores)], "\n")

# bar plot of alignment scores
# tried a few colour options, went with this one
ggplot(df, aes(x = reorder(variants, alignment_scores), y = alignment_scores)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "BLASTp Alignment Scores of TP53 Missense Variants",
       x = "Variant",
       y = "Alignment Score") +
  theme_minimal()

ggsave("tp53_alignment_scores.png", width = 7, height = 5)
cat("plot saved\n")

# scatter plot - identity % for each variant
ggplot(df, aes(x = variants, y = identity_pct)) +
  geom_point(size = 3, colour = "darkred") +
  geom_hline(yintercept = mean(identity_pct), linetype = "dashed", colour = "grey50") +
  labs(title = "Sequence Identity % - TP53 Variants vs Wild Type",
       x = "Variant",
       y = "Identity (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("tp53_identity_scatter.png", width = 7, height = 5)
cat("scatter plot saved\n")

# R248W and R273H had notably lower identity scores
# consistent with what the literature says about these being
# highly destabilising mutations in the DNA binding domain
