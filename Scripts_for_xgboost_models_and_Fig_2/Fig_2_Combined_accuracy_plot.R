
library(tidyverse)
library(data.table)
library(readxl)
library(ggbeeswarm)
library(cowplot)


# xgboost fin test set metrics
test_metrics <- fread('results/xgboost_test_metrics.tsv')
test_metrics <- filter(test_metrics, 
                       Metric %in% c("Accuracy","AccuracyLower","Sensitivity","Specificity",
                                     "Pos.Pred.Value","Neg.Pred.Value","F1","Prevalence","Balanced.Accuracy"))
test_metrics$Metric <- str_replace_all(test_metrics$Metric, 'Balanced.Accuracy', 'Balanced accuracy')
test_metrics$Metric <- str_replace_all(test_metrics$Metric, 'AccuracyLower', 'Accuracy lower')
test_metrics$Metric <- str_replace_all(test_metrics$Metric, 'Pos.Pred.Value', 'PPV')
test_metrics$Metric <- str_replace_all(test_metrics$Metric, 'Neg.Pred.Value', 'NPV')

colnames(test_metrics)[2] <- 'Antigen'
test_metrics <- test_metrics %>% dplyr::select(-Antigen) %>% group_by(., Metric) %>% 
  summarise(Mmean=mean(value)) %>% left_join(., test_metrics) %>% arrange(Mmean)
test_metrics$name <- factor(test_metrics$Metric, levels = test_metrics$Metric %>% unique)

# rf fin test set metrics
tmp <- readRDS('data/Conf_mat_test.rds')
tmp.names <- tmp %>% names %>% .[-c(1)]
tmp <- c(tmp[[1]][1:4], map(2:length(tmp), function(i) tmp[[i]][[2]]) )
names(tmp)[-c(1:4)] <- tmp.names 

test_metrics.2 <- map2(1:length(tmp), names(tmp), function(i, y) {
  res <- tmp[[i]]
  data.frame(res$overall[c(1,3)] %>% t, res$byClass %>% t, Antigen=y)
}) %>% do.call(rbind, .)
test_metrics.2 %>% dim
test_metrics.2 <- pivot_longer(test_metrics.2, 1:13)
test_metrics.2 %>% head
colnames(test_metrics.2)[2] <- 'Metric'
test_metrics.2 <- filter(test_metrics.2, 
                         Metric %in% c("Accuracy","AccuracyLower","Sensitivity","Specificity",
                                       "Pos.Pred.Value","Neg.Pred.Value","F1","Prevalence","Balanced.Accuracy"))
test_metrics.2$Metric <- str_replace_all(test_metrics.2$Metric, 'Balanced.Accuracy', 'Balanced accuracy')
test_metrics.2$Metric <- str_replace_all(test_metrics.2$Metric, 'AccuracyLower', 'Accuracy lower')
test_metrics.2$Metric <- str_replace_all(test_metrics.2$Metric, 'Pos.Pred.Value', 'PPV')
test_metrics.2$Metric <- str_replace_all(test_metrics.2$Metric, 'Neg.Pred.Value', 'NPV')
# test_metrics.2 <- test_metrics.2 %>% dplyr::select(-Antigen) %>% group_by(., Metric) %>% 
#   summarise(Mmean=mean(value)) %>% left_join(., test_metrics.2) %>% arrange(Mmean)
# test_metrics.2$Metric <- factor(test_metrics.2$Metric, levels = test_metrics.2$Metric %>% unique)
test_metrics.2$name <- factor(test_metrics.2$Metric, levels = test_metrics.2$Metric %>% unique)

# fin model in Danish full data
test.metrics.FD <- read_xlsx('data/suppl4.xlsx', sheet = 1)[, -1]
# dan model in Danish test data
test.metrics.DDtest <- read_xlsx('data/suppl4.xlsx', sheet = 2)[, -1]

colnames(test.metrics.FD)[c(1,4,5)] <- colnames(test.metrics.DDtest)[c(1,4,5)] <- c('Antigen', 'PPV', 'NPV')
test.metrics.FD <- test.metrics.FD %>% pivot_longer(2:6)
test.metrics.DDtest <- test.metrics.DDtest %>% pivot_longer(2:6)



# combine data for plotting
dat <- rbind(data.frame(test_metrics.2 %>% dplyr::select(c(Antigen, name, value)), Model='Finnish random forest\nFinnish test set'),
             data.frame(test_metrics %>% dplyr::select(c(Antigen, name, value)), Model='Finnish gradient boosting\nFinnish test set'),
             data.frame(test.metrics.FD %>% dplyr::select(c(Antigen, name, value)), Model='Finnish random forest\nDanish full set'),
             data.frame(test.metrics.DDtest %>% dplyr::select(c(Antigen, name, value)), Model='Danish random forest\nDanish test set'))
fwrite(dat, 'data/BG_test_datasets.tsv', sep='\t')
common.antigens <- map(dat$Model %>% unique, function(x) filter(dat, Model==x)$Antigen) %>% Reduce(intersect, .)
dat <- filter(dat, name %in% c("Sensitivity","Specificity","PPV","NPV","Balanced accuracy"),
              Antigen %in% common.antigens)
dat$name <- factor(dat$name)
dat$name <- dat$name %>% droplevels()
dat$value <- dat$value %>% as.numeric()
dat$Model <- factor(dat$Model, levels = dat$Model %>% unique())
  
# swarm plot
p.rf.xg <- dat %>% ggplot(aes(value, name, color=value)) +
  # geom_jitter(alpha=0.6, size = 1.5, height = 0.3) +
  geom_quasirandom(groupOnX = F, width = 0.4, alpha=0.6, size = 1.5, position = 'pseudorandom',
                   bandwidth = 0.8, nbins = NULL, method = 'pseudorandom', dodge.width = 10) +
  # geom_beeswarm(groupOnX = F, width = 0.1, alpha=0.6, size = 2) +
  # xlab('Metric value') + ylab('Metric name') +
  xlab('') + ylab('') +
  # geom_vline(xintercept = 0.95, linewidth = 0.3, color = 'grey') +
  scale_x_continuous(breaks = seq(0, 1, by=0.25),
                     labels = c('0', '0.25', '0.50', '0.75', '1')) +
  scale_color_gradient(high="green4", low="firebrick2") +
  facet_wrap(~ Model, ncol = 2) +
  theme_minimal() +
  theme(strip.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_blank(),
        # axis.ticks.length.x = unit(100, 'points'),
        axis.ticks.x = element_line(linewidth = 0.3),
        legend.position = "none",
        title = element_text(size = 12),
        panel.border = element_rect(fill = NA),
        panel.spacing = unit(0.6, "lines")) 

p.rf.xg <- ggdraw(p.rf.xg) + 
  draw_label("A", x = 0.27,  y = 0.95, size=13) + 
  draw_label("B", x = 0.64,  y = 0.95, size=13) +
  draw_label("C", x = 0.27,  y = 0.5,  size=13) +
  draw_label("D", x = 0.64,  y = 0.5,  size=13) 

jpeg('./results/rf_xg_accuracy_metrics_test_summary.jpeg', 
     width=6, height=4, res=600, units='in')
p.rf.xg
dev.off()



