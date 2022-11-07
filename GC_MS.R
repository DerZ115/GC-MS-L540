library(openxlsx)
library(tidyverse)
library(data.table)
library(rstatix)
library(ggplot2)
library(ggpubr)

source("lib_parse.R") # Load functions

integ_file <- "data/Integration_report_L_540.xlsx"

xlsx_sheets <- getSheetNames(integ_file)

results = list()
all_compounds <- c()

dir.create("merged", showWarnings=F)

for (i in 1:length(xlsx_sheets)) {
  integ_result <- read.xlsx(integ_file, sheet=i, startRow=4, colNames=T, cols=c(2,5))
  search_file <- paste0("data/", xlsx_sheets[i], ".txt")
  search_result <- parse_libsearch(search_file)
  full_result <- merge(integ_result, search_result, by="Ret.Time")
  full_result <- full_result[full_result$Quality >= 70, ]
  results[[i]] <- full_result
  all_compounds <- c(all_compounds, results[[i]]$CAS)
  
  write.csv(results[[i]], file=paste0("merged/", xlsx_sheets[i], ".csv"), row.names=F)
}

unique_compounds <- as.data.frame(table(all_compounds), stringsAsFactors = FALSE)
colnames(unique_compounds) <- c("cas", "count")
unique_compounds[c("name", "category", "weight_orig", "weight_deriv")] <- NA

for (i in 1:length(unique_compounds$cas)) {
  cas <- unique_compounds$cas[i]
  for (df in results) {
    if (cas %in% df$CAS) {
       unique_compounds$name[i] <- df$Compound[df$CAS == cas][1]
       break
    }
  }
}

names(results) <- xlsx_sheets

### Uncomment if running for the first time, to generate a new compound file
# write.csv(unique_compounds, "compounds.csv", row.names=F)

unique_compounds <- read_csv("compounds.csv", col_types = "ciccdd")

res_sum <- list()

for (df in results) {
  df <- inner_join(df, unique_compounds, by=c("CAS" = "cas"))
  std_idx <- which(df$name == "Cinnamic Acid")
  std <- df[std_idx,]
  df <- df[-std_idx,]
  
  std$Concentration_deriv <- 2.5e-3 * 5e-6 * std$weight_deriv * 1e9 # ng per 100 uL
  std$Concentration_orig <- 2.5e-3 * 5e-6 * std$weight_orig * 1e9
  area_per_ng <- std$Area / std$Concentration_deriv
  
  df$Concentration_deriv <- df$Area / area_per_ng
  df$Concentration_orig <- df$Concentration_deriv * df$weight_orig / df$weight_deriv
  
  res_sum[[length(res_sum)+1]] <- df %>% 
    group_by(name) %>%
    summarise(Concentration = sum(Concentration_orig))
  
}

names(res_sum) <- xlsx_sheets

result <- res_sum %>% 
  reduce(full_join, by="name") %>% column_to_rownames(., var = "name")

colnames(result) <- xlsx_sheets

counts <- deframe(read.csv("data/cell_counts.csv", header = F))

for (col in colnames(result)) {
  result[,col] <- result[,col] * 1e6 / counts[substring(col, 1, nchar(col)-1)]
}

result <- result[rowSums(is.na(result)) <= 12,]
result <- as.data.frame(t(result))
result[is.na(result)] <- 0

rows <- rownames(result)
result$Treatment <- substring(rows, 1, nchar(rows)-2)


result.C_E <- result %>%
  filter(Treatment %in% c("C", "E")) %>%
  pivot_longer(-Treatment, names_to = "Compound", values_to = "Concentration")

t.test.C_E <- result.C_E %>%
  group_by(Compound) %>%
  t_test(Concentration ~ Treatment) %>%
  add_significance()

ttest.top6.C_E <- t.test.C_E %>% 
  arrange(p) %>%
  head(6) %>%
  add_xy_position(x = "Treatment", scales = "free_y")

top6.C_E <- result.C_E %>%
  filter(Compound %in% ttest.top6.C_E$Compound)

p <- ggboxplot(
  top6.C_E, x="Treatment", y="Concentration",
  legend = "none", ggtheme = theme_pubr(border = T)
) + 
  facet_wrap(~factor(Compound, levels = ttest.top6.C_E$Compound), scales = "free_y") + 
  stat_pvalue_manual(ttest.top6.C_E, label = "p = {p}", vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(y=bquote("Concentration (ng per "~10^6~" cells)")) +
  theme(strip.text = element_text(face = "bold"))

p


result.C_R <- result %>%
  filter(Treatment %in% c("C", "R")) %>%
  pivot_longer(-Treatment, names_to = "Compound", values_to = "Concentration")

t.test.C_R <- result.C_R %>%
  group_by(Compound) %>%
  t_test(Concentration ~ Treatment) %>%
  add_significance()

ttest.top6.C_R <- t.test.C_R %>% 
  arrange(p) %>%
  head(6) %>%
  add_xy_position(x = "Treatment", scales = "free_y")

top6.C_R <- result.C_R %>%
  filter(Compound %in% ttest.top6.C_R$Compound)

p <- ggboxplot(
  top6.C_R, x="Treatment", y="Concentration",
  legend = "none", ggtheme = theme_pubr(border = T)
) + 
  facet_wrap(~factor(Compound, levels = ttest.top6.C_R$Compound), scales = "free_y") + 
  stat_pvalue_manual(ttest.top6.C_R, label = "p = {p}", vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(y=bquote("Concentration (ng per "~10^6~" cells)")) +
  theme(strip.text = element_text(face = "bold"))

p


result.C_RE <- result %>%
  filter(Treatment %in% c("C", "RE")) %>%
  pivot_longer(-Treatment, names_to = "Compound", values_to = "Concentration")

t.test.C_RE <- result.C_RE %>%
  group_by(Compound) %>%
  t_test(Concentration ~ Treatment) %>%
  add_significance()

ttest.top6.C_RE <- t.test.C_RE %>% 
  arrange(p) %>%
  head(6) %>%
  add_xy_position(x = "Treatment", scales = "free_y")

top6.C_RE <- result.C_RE %>%
  filter(Compound %in% ttest.top6.C_RE$Compound)

p <- ggboxplot(
  top6.C_RE, x="Treatment", y="Concentration",
  legend = "none", ggtheme = theme_pubr(border = T)
) + 
  facet_wrap(~factor(Compound, levels = ttest.top6.C_RE$Compound), scales = "free_y") + 
  stat_pvalue_manual(ttest.top6.C_RE, label = "p = {p}", vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(y=bquote("Concentration (ng per "~10^6~" cells)")) +
  theme(strip.text = element_text(face = "bold"))

p


result.E_RE <- result %>%
  filter(Treatment %in% c("E", "RE")) %>%
  pivot_longer(-Treatment, names_to = "Compound", values_to = "Concentration")

t.test.E_RE <- result.E_RE %>%
  group_by(Compound) %>%
  t_test(Concentration ~ Treatment) %>%
  add_significance()

ttest.top6.E_RE <- t.test.E_RE %>% 
  arrange(p) %>%
  head(6) %>%
  add_xy_position(x = "Treatment", scales = "free_y")

top6.E_RE <- result.E_RE %>%
  filter(Compound %in% ttest.top6.E_RE$Compound)

p <- ggboxplot(
  top6.E_RE, x="Treatment", y="Concentration",
  legend = "none", ggtheme = theme_pubr(border = T)
) + 
  facet_wrap(~factor(Compound, levels = ttest.top6.E_RE$Compound), scales = "free_y") + 
  stat_pvalue_manual(ttest.top6.E_RE, label = "p = {p}", vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(y=bquote("Concentration (ng per "~10^6~" cells)")) +
  theme(strip.text = element_text(face = "bold"))

p
