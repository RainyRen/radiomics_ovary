library(tidyverse)
library(caret)
library(pROC)
library(glmnet)
library(DMwR)
library(rmda)
library(ggpubr)
library(ModelGood)
library(rms)
library(mRMRe)
library(DescTools)
library(Publish)
library(pheatmap)

# # ========== Environment settings ==========
# # 移除现有的环境变量
rm(list=ls())
setwd("C:/Workstation/RProject/duan_code")
getwd()
# # =========================================================

source('assist_file.R')

set.seed(1234)

clinics_column <- 16

fpath <- choose.files()

dt_all <- read.csv(fpath)

dt_all$Label <- factor(dt_all$Label, ordered = T)

dt <- dt_all[, -c(2:clinics_column)]

if(!is_empty(nearZeroVar(dt)))
{
  dt <- dt[, -nearZeroVar(dt)]
}


# clinical data

dt_cli <- dt_all[, c(1:clinics_column)]

acc_val_train <- numeric(0)
acc_val_test <- numeric(0)

all_idx <- list()
train_list <- list()
test_list <- list()
cvfit_list <- list()
fit_list <- list()
out_list <- list()

for(numIt  in c(1:10))
{
  s_pre <- preProcess(dt, method = 'medianImpute')
  dt <- predict(s_pre, dt)
  idx_train <- createDataPartition(dt$Label, p = 0.7, list = F)
  dt_train <- dt[idx_train, ]
  dt_test <- dt[-idx_train, ]
  
  # browser()
  
  all_idx[[numIt]] <- idx_train
  step_pre <- preProcess(dt_train, method = c('center', 'scale'))
  dt_train_pre <- predict(step_pre, dt_train) %>% as.data.frame
  
  dt_train_pre0 <- dt_train_pre
  # dt_train_pre <- SMOTE(Label~., data = dt_train_pre)
  dt_test_pre <- predict(step_pre, dt_test)
  
  dt_mrmr <- mRMR.data(dt_train_pre)
  f_sel <- mRMR.classic(data = dt_mrmr, target_indices = c(1), feature_count = 20)
  
  # browser()
  
  dt_train_pre <- select(dt_train_pre, c('Label', featureNames(f_sel)[unlist(solutions(f_sel))]))
  dt_train_pre0 <- select(dt_train_pre0, c('Label', featureNames(f_sel)[unlist(solutions(f_sel))]))
  dt_test_pre <- select(dt_test_pre, c('Label', featureNames(f_sel)[unlist(solutions(f_sel))]))
  
  x <- as.matrix(dt_train_pre[, -1])
  y <- dt_train_pre$Label
  
  cv.fit <- cv.glmnet(x, y, family = 'binomial')
  fit <- glmnet(x, y, family = 'binomial')
  
  
  train_list[[numIt]] <- dt_train_pre0
  test_list[[numIt]] <- dt_test_pre
  cvfit_list[[numIt]] <- cv.fit
  fit_list[[numIt]] <- fit
  
  # browser()
  pre_res_test <- as.vector(predict(fit, newx = as.matrix(dt_test_pre[, -1]), s = cv.fit$lambda.min))
  roc_res_test <- pROC::roc(dt_test_pre$Label, pre_res_test)
  
  
  
  pre_res_train <- as.vector(predict(fit, newx = x, s = cv.fit$lambda.min))
  roc_res_train <- pROC::roc(dt_train_pre$Label, pre_res_train)
  
  
  dir_sign_test <- roc_res_test$direction
  dir_sign_train <- roc_res_train$direction
  
  if(dir_sign_test == dir_sign_train)
  {
    acc_val_test <- c(acc_val_test, pROC::auc(roc_res_test))
    acc_val_train <- c(acc_val_train, pROC::auc(roc_res_train))
  }
  else
  {
    acc_val_test <- c(acc_val_test, 0)
    acc_val_train <- c(acc_val_train, 0)
  }
}

idx_vec <- c(1:length(acc_val_test))
idx_vec <- idx_vec[acc_val_train > acc_val_test]
acc_val <- acc_val_test[acc_val_train > acc_val_test]
init_idx <- which.max(acc_val)
sel_idx <- idx_vec[init_idx]

idx_train <- all_idx[[sel_idx]]
grp_info <- tibble(Label = dt$Label, Group = 'Test')
grp_info$Radscore <- 0
grp_info$Group[idx_train] <- 'Training'


dt_train_final <- train_list[[sel_idx]]
dt_test_final <- test_list[[sel_idx]]
cvfit <- cvfit_list[[sel_idx]]
fit <- fit_list[[sel_idx]]

s = cvfit$lambda.min


pre_res_test <- as.vector(predict(fit, newx = as.matrix(dt_test_final[, -1]), s = s))
pre_res_test_prob <- as.vector(predict(fit, newx = as.matrix(dt_test_final[, -1]), s = s, 
                                       type = 'link'))
roc_res_test <- pROC::roc(dt_test_final$Label, pre_res_test, ci = T)

out_res_test <- ifelse(pre_res_test > coords(roc_res_test, x = 'best')[1], 1, 0)
conf_mat_test <- confusionMatrix(as.factor(out_res_test),as.factor(dt_test_final$Label))
rec_test <- c(conf_mat_test$overall[c(1, 3, 4)], conf_mat_test$byClass[c(1:4)])

pre_res_train <- as.vector(predict(fit, newx = as.matrix(dt_train_final[, -1]), s = s))
pre_res_train_prob <- as.vector(predict(fit, newx = as.matrix(dt_train_final[, -1]), s = s, 
                                        type = 'link'))
roc_res_train <- pROC::roc(dt_train_final$Label, pre_res_train, ci = T)

out_res_train <- ifelse(pre_res_train > coords(roc_res_train, x = 'best')[1], 1, 0)
conf_mat_train <- confusionMatrix(as.factor(out_res_train), as.factor(dt_train_final$Label))
rec_train <- c(conf_mat_train$overall[c(1, 3, 4)], conf_mat_train$byClass[c(1:4)])

rec_rad <- data.frame(rbind(rec_train, rec_test), row.names = c('Train', 'Test'))

write.csv(rec_rad, file = 'res_radiomics.csv')

grp_info$Radscore[idx_train] <- pre_res_train
grp_info$Radscore[-idx_train] <- pre_res_test

write_csv(grp_info, 'group_info.csv')

cutoff_radiomics <- coords(roc_res_train, x = 'best')

## rad score
dt_final_test <- tibble(Label = dt_test_final$Label, rad_score = pre_res_test)
dt_final_arr <- arrange(dt_final_test, rad_score)
dt_final_arr$x <- 1:nrow(dt_final_arr)

dt_final_train <- tibble(Label = dt_train_final$Label, rad_score = pre_res_train)

p_train <- ggboxplot(x = 'Label', y = 'rad_score', data = dt_final_train,
                     add = 'jitter', color = 'Label', palette = 'jco') + 
  ylim(-3, 3) + 
  stat_compare_means(method = 'wilcox.test') + 
  geom_hline(yintercept = coords(roc_res_train, x = 'best')[1]) + theme_bw()
p_test <- ggboxplot(x = 'Label', y = 'rad_score', data = dt_final_test,
                    add = 'jitter', color = 'Label', palette = 'jco') +
  ylim(-3, 3) +
  stat_compare_means(method = 'wilcox.test') + 
  geom_hline(yintercept = coords(roc_res_train, x = 'best')[1]) + theme_bw()



# p_rad <- ggplot(aes(x = x, y = rad_score), data = dt_final_arr)
# p_rad <- p_rad + geom_col(aes(fill = Label)) + labs(x = '', y = 'Rad Score') + 
#   theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())



coefs <- coefficients(fit, s = s)
useful_feature <- unlist(coefs@Dimnames)[coefs@i + 1]
useful_feature <- useful_feature[-1]

dt_coef <- data.frame(Feature = useful_feature, Coef = coefs@x[-1])
dt_coef <- arrange(dt_coef, desc(Coef))
dt_coef$Feature <- factor(dt_coef$Feature, 
                          levels = as.character(dt_coef$Feature))

p_coef <- ggplot(aes(x = Feature, y = Coef), data = dt_coef)
p_coef <- p_coef + geom_col(fill = 'blue', width = 0.7) + coord_flip() + 
  theme_bw() + ylab('Coefficients')

final_data_test <- add_column(dt_test_final, radscore = pre_res_test)
final_data_test <- select(final_data_test, c('Label', useful_feature, 'radscore'))
write_csv(final_data_test, path = 'dataset_test.csv')
final_data_train <- add_column(dt_train_final, radscore = pre_res_train)
final_data_train <- select(final_data_train, c('Label', useful_feature, 'radscore'))
write_csv(final_data_train, path = 'dataset_train.csv')



fit_train <- glm(Label~rad_score, data = dt_final_train, family = 'binomial')

dt_dca_train <- dt_final_train
dt_dca_train$Label <- as.numeric(dt_dca_train$Label) - 1
dca_curve <- decision_curve(Label~rad_score, data = dt_dca_train)

radscore <- paste('Radscore = ', paste(round(coefs@x[-1], 3), 
                                       useful_feature, sep = '*', collapse = '+'), '+', round(coefs@x[1], 3))
print(radscore)
write_file(radscore, path = 'radscore.txt')

# figure1
oldpar <- par(mfrow = c(2, 1))
plot(cvfit)
# figure2
plot(fit, s = s, xvar = 'lambda')
abline(v = log(cvfit$lambda.min), lty = 2)
par(oldpar)
# figure3

oldpar <- par(mfrow = c(1, 2))
plot(roc_res_train, print.auc = T, 
     print.auc.pattern = 'AUC: %.2f(%.2f-%.2f)', legacy.axes = T)

# figure4
plot(roc_res_test, print.auc = T, legacy.axes = T,
     print.auc.pattern = 'AUC: %.2f(%.2f-%.2f)')

par(oldpar)
# figure5
ggarrange(p_train, p_test, ncol = 2)



#figure 6
p_coef + theme(axis.text.y = element_text(size = 12))

# clinics analysis

idx_train <- all_idx[[sel_idx]]
dt_cli_train <- dt_cli[idx_train, ]
dt_cli_test <- dt_cli[-idx_train, ]

p_list <- lapply(dt_cli_train[, -1], inter_test, y = dt_cli_train$Label)

#p_list_test <- lapply(dt_cli_test[, -1], inter_test, y = dt_cli_test$Label)


sel_name <- c(names(which(p_list < 0.1)))
dt_cli_train_1 <- select(dt_cli_train, c('Label', sel_name))

res_ulogit_list <- lapply(colnames(dt_cli_train_1)[-1], ulogit_test, dt = dt_cli_train_1)
res_ulogit <- bind_rows(res_ulogit_list)
res_ulogit <- res_ulogit[-seq(from = 1, to = nrow(res_ulogit), by = 2), ]
res_ulogit$Var <- colnames(dt_cli_train_1)[-1]

res_ulogit_sel <- filter(res_ulogit, p_val < 0.1)
cli_name <- res_ulogit_sel$Var

dt_cli_train_2 <- select(dt_cli_train_1, c('Label', cli_name))
dt_cli_test_2 <- select(dt_cli_test, c('Label', cli_name))
  
write.csv(res_ulogit, file = 'ulogit_cli.csv')

log_fit <- glm(Label~., data = dt_cli_train_2, family = binomial)
vif_val <- vif(log_fit)
vif_max <- max(vif_val)
vif_max_idx <- which.max(vif_val)
while(vif_max > 10)
{
  dt_cli_train_i <- select(dt_cli_train_2, -c(names(vif_max_idx)))
  log_fit <- glm(Label~., data = dt_cli_train_i, family = binomial)
  vif_val <- vif(log_fit)
  vif_max <- max(vif_val)
  vif_max_idx <- which.max(vif_val)
}
log_fit_final <- step(log_fit)

output_mlogit(log_fit_final, 'mlogit_cli.csv')

dt_cli_train_final <- model.frame(log_fit_final, data = dt_cli_train)
dt_cli_test_final <- model.frame(log_fit_final, data = dt_cli_test)

dt_combined_train <- dt_cli_train_final
dt_combined_train$rad_score <- pre_res_train_prob

dt_combined_test <- dt_cli_test_final
dt_combined_test$rad_score <- pre_res_test_prob


com_form <- paste('Label', paste(colnames(dt_combined_train)[-1], collapse = '+'), sep = '~') %>% as.formula


mod_com_final <- glm(com_form, data = dt_combined_train, family = 'binomial')

output_mlogit(mod_com_final, 'mlogit_com.csv')

dt_combined_train_final <- model.frame(mod_com_final, data = dt_combined_train)
dt_combined_test_final <- model.frame(mod_com_final, data = dt_combined_test)

com_form <- paste('Label', paste(colnames(dt_combined_train_final)[-1], collapse = '+'), sep = '~') %>% as.formula

dt_nom <- filter(dt_combined_train_final, rad_score > -5 & rad_score < 10)

ddist_train_com <- datadist(dt_nom)
options(datadist = 'ddist_train_com')
mod_train <- lrm(com_form, 
                 data = dt_nom, x = TRUE, y = TRUE)
nom_com <- nomogram(mod_train, lp = F, fun = plogis, fun.at = c(0.1, 0.4, 0.9), 
                    funlabel = 'Risk')
plot(nom_com)


mod_test <- lrm(com_form, 
                data = dt_combined_test_final, x = TRUE, y = TRUE)

oldpar <- par(mfrow = c(1, 2))
cal_train <- calPlot2(mod_train, data = dt_combined_train_final, 
                      legend = F, col = 	'#FF6EB4', lty = 2)
cal_test <- calPlot2(mod_train, data = dt_combined_test_final,
                     legend = F, col = 	'#FF6EB4', lty = 2)
par(oldpar)


cli_form <- paste('Label', paste(colnames(dt_combined_train_final)[-c(1, ncol(dt_combined_train_final))], collapse = '+'), 
                  sep = '~') %>% as.formula

HosmerLemeshowTest(cal_train$Frame$lrm, cal_train$Frame$jack)
HosmerLemeshowTest(cal_test$Frame$lrm, cal_test$Frame$jack)


res_train_final <- predict(mod_com_final, newdata = dt_combined_train)
res_roc_com_train <- pROC::roc(dt_combined_train$Label, res_train_final, 
                               ci = T)
cutoff <- coords(res_roc_com_train, x = 'best')[1]
res_train_final_bin <- as.factor(ifelse(res_train_final > cutoff, 1, 0))
res_conf_com_train <- confusionMatrix(dt_combined_train$Label, 
                                      res_train_final_bin, positive = '1')
res_test_final <- predict(mod_train, newdata = dt_combined_test)
res_roc_com_test <- pROC::roc(dt_combined_test$Label, res_test_final,
                              ci = T)
cutoff <- coords(res_roc_com_test, x = 'best')[1]
res_test_final_bin <- as.factor(ifelse(res_test_final > cutoff, 1, 0))
res_conf_com_test <- confusionMatrix(dt_combined_test$Label, 
                                     res_test_final_bin, positive = '1')


rec_train_com <- c(res_conf_com_train$overall[c(1, 3, 4)], res_conf_com_train$byClass[c(1:4)])
rec_test_com <- c(res_conf_com_test$overall[c(1, 3, 4)], res_conf_com_test$byClass[c(1:4)])

rec_all <- data.frame(rbind(rec_train_com, rec_test_com), row.names = c('Train', 'Test'))
write.csv(rec_all, file = 'res_combined.csv')


dt_dca <- dt_combined_train_final
dt_dca$Label <- ifelse(dt_dca$Label == '0', 0, 1)
dca1 <- decision_curve(com_form, 
                       data = dt_dca)
dca2 <- decision_curve(cli_form, 
                       data = dt_dca)

plot_decision_curve(list(dca1, dca2), confidence.intervals = F, 
                    col = c('red', 'green', 'blue', 'black'), 
                    curve.names = c('With Radscore', 'Without Radscore'),
                    legend.position = 'topright', cost.benefits = FALSE)

cli_mod <- glm(cli_form, data = dt_combined_train_final, family = 'binomial')

res_cli_train <- predict(cli_mod, newdata = dt_combined_train_final, type = 'link')
roc_cli_train <- pROC::roc(dt_combined_train_final$Label, res_cli_train,ci=T)
cutoff_cli <- coords(roc_cli_train, x = 'best')[1]

res_cli_test <- predict(cli_mod, newdata = dt_combined_test_final, type = 'link')
roc_cli_test <- pROC::roc(dt_combined_test_final$Label, res_cli_test,ci=T)


res_cli_train_bin <- as.factor(ifelse(res_cli_train > cutoff_cli, 1, 0))
res_cli_test_bin <- as.factor(ifelse(res_cli_test > cutoff_cli, 1, 0))



conf_mat_cli_test <- confusionMatrix(dt_combined_test_final$Label, 
                                     res_cli_test_bin, positive = '1')

conf_mat_cli_train <- confusionMatrix(dt_combined_train_final$Label,
                                      res_cli_train_bin, positive = '1')


rec_train_cli <- c(conf_mat_cli_train$overall[c(1, 3, 4)], conf_mat_cli_train$byClass[c(1:4)])
rec_test_cli <- c(conf_mat_cli_test$overall[c(1, 3, 4)], conf_mat_cli_test$byClass[c(1:4)])

rec_cli <- data.frame(rbind(rec_train_cli, rec_test_cli), row.names = c('Train', 'Test'))
write.csv(rec_cli, file = 'res_clinics.csv')


oldpar <- par(mfrow = c(1, 2))
plot(res_roc_com_train, print.auc = T, print.auc.pattern = 'AUC: %.2f (%.2f - %.2f)',
     legacy.axes = T, col = 'red')
plot(roc_res_train, print.auc = T, print.auc.pattern = 'AUC: %.2f (%.2f - %.2f)',
     add = T, col = 'blue', print.auc.y = 0.45)
plot(roc_cli_train, print.auc =T, print.auc.pattern = 'AUC: %.2f (%.2f - %.2f)', 
     add = T, col = 'green', print.auc.y = 0.4)
legend(x = 0.3, y = 0.2, legend = c('Combined', 'Radiomics', 'Clinics'), 
       col = c('red', 'blue', 'green'), lty = 1)


plot(res_roc_com_test, print.auc = T, print.auc.pattern = 'AUC: %.2f (%.2f - %.2f)',
     legacy.axes = T, col = 'red', print.auc.y = 0.5)
plot(roc_res_test, print.auc = T, print.auc.pattern = 'AUC: %.2f (%.2f - %.2f)',
     add = T, col = 'blue', print.auc.y = 0.45)
plot(roc_cli_test, print.auc =T, print.auc.pattern = 'AUC: %.2f (%.2f - %.2f)', 
     add = T, col = 'green', print.auc.y = 0.4)

legend(x = 0.3, y = 0.2, legend = c('Combined', 'Radiomics', 'Clinics'), 
       col = c('red', 'blue', 'green'), lty = 1)
par(oldpar)

rec_final <- bind_rows(rec_cli, rec_rad, rec_all)
rec_final$Group <- rep(c('Trainig', 'Test'), times = 3)
rec_final$Model <- rep(c('Clinics', 'Radiomics', 'Nomogram'), times = c(2, 2, 2))

#roc.test(res_roc_com_train,res_roc_com_test)

#roc.test(roc_cli_train,roc_cli_test)

#roc.test(roc_res_train,roc_res_test)



rmarkdown::render('results_final.Rmd')
