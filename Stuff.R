library(ggplot2)
library(gdata)

setwd("./Project")
CpGs <- read.csv("34_CpGs_10x_with_regression_and_CpGs_info_present_in3.csv", sep = "\t")
CpGs$p_fdr <- p.adjust(CpGs$TCDD_P, method = "fdr", n = length(CpGs$TCDD_P))
write.csv(CpGs, "34_CpGs_10x_with_regression_and_CpGs_info_present_in3_fdr.csv", row.names = FALSE, col.names = FALSE)
CpGs_sign <- CpGs[CpGs$p_fdr < 0.05,]
CpGs_1 <- CpGs[,2:35]
write.csv(CpGs_1, "CpGs.csv", row.names = FALSE, col.names = FALSE)
CpGs_11 <- read.csv("CpGs.csv", sep = ",")
CpGs_11$n <- count.fields("CpGs.csv", sep = ",", quote = "-", skip = 1)
?count.fields
colnames(CpGs) <- names
write.csv(CpGs, "34_CpGs_10x_with_regression_and_CpGs_info_sign_2.csv", sep = "\t", row.names = FALSE)


smoke <- read.csv("smoking_34_subjects.csv",header=TRUE, sep = ";")
smoke$smoking_last6months <- ifelse(smoke$LASTONE == 5, 0, smoke$SMOKEDAY)
smoke$smoking_last6months[is.na(smoke$smoking_last6months)] <- 0
smoke$TCDD <- IDs$nstq2378d
smoke$smoking_last6months_reduced <- smoke$smoking_last6months
smoke$smoking_last6months_reduced <- ifelse (smoke$smoking_last6months_reduced <3 &  smoke$smoking_last6months_reduced !=0, 1, smoke$smoking_last6months_reduced)
smoke$smoking_last6months_reduced <- ifelse (smoke$smoking_last6months_reduced >= 3, 2, smoke$smoking_last6months_reduced)

write.csv(smoke, "Smoke_and_TCDD.csv", sep = "\t", row.names = FALSE)

fit_smoke_tcdd <- lm (smoke$smoking_last6months ~ smoke$TCDD, data = smoke)
summary(fit_smoke_tcdd)
layout(matrix(1,1,1)) # optional 4 graphs/page 
plot(fit_smoke_tcdd)


####

setwd("./Project")
CpGs_smoking <- read.csv("34_10x_CpG~Smoking.csv", sep = "\t")

CpGs_smoking$p_fdr_3 <- p.adjust(CpGs_smoking$TCDD_P_3, method = "fdr", n = length(CpGs_smoking$TCDD_P_3))
CpGs_smoking$p_fdr_6 <- p.adjust(CpGs_smoking$TCDD_P_6, method = "fdr", n = length(CpGs_smoking$TCDD_P_6))
write.csv(CpGs_smoking, "34_10x_CpG~Smoking.csv", row.names = FALSE, col.names = FALSE)

CpGs_sign <- CpGs[CpGs$p_fdr < 0.05,]

####
#smoke6 and TCDD conc

setwd("./Project")
CpGs_smoking6_TCDD <- read.csv("34_10x_CpG~TCDD_and_Smoke6_cat.csv", sep = "\t")
kek <- CpGs_smoking6_TCDD
CpGs_smoking6_TCDD <- kek


CpGs_smoking6_TCDD$B_Smoke_T0 <- kek$B_TCDD
CpGs_smoking6_TCDD$B_Smoke_T1 <- kek$B_Smoke_T0
CpGs_smoking6_TCDD$B_Smoke_T2 <- kek$B_Smoke_T1
CpGs_smoking6_TCDD$B_Smoke_T3 <- kek$B_smoke_T2
CpGs_smoking6_TCDD$B_Smoke_T4 <- kek$B_smoke_T3
CpGs_smoking6_TCDD$B_TCDD <- kek$B_smoke_T4
CpGs_smoking6_TCDD$P_smoke_T2 <- NULL
CpGs_smoking6_TCDD$P_TCDD <- NULL
CpGs_smoking6_TCDD$P_smoke_T0 <- NULL
CpGs_smoking6_TCDD$p_smoke_T2 <- NULL
CpGs_smoking6_TCDD$P_smoke_T4 <- NULL
CpGs_smoking6_TCDD$P_smoke_T1 <- NULL
CpGs_smoking6_TCDD$P_Smoke_T0 <- kek$P_TCDD
CpGs_smoking6_TCDD$P_Smoke_T1 <- kek$P_smoke_T0
CpGs_smoking6_TCDD$P_Smoke_T2 <- kek$P_smoke_T1
CpGs_smoking6_TCDD$P_Smoke_T3 <- kek$p_smoke_T2
CpGs_smoking6_TCDD$P_Smoke_T4 <- kek$P_smoke_T3
CpGs_smoking6_TCDD$P_TCDD <- kek$P_smoke_T4

CpGs_smoking6_TCDD$P_summ_fdr <- p.adjust(CpGs_smoking6_TCDD$P_summ, method = "fdr", n = length(CpGs_smoking6_TCDD$P_summ))
CpGs_smoking6_TCDD$P_summ_fdr_1 <- p.adjust(CpGs_smoking6_TCDD$P_summ, method = "fdr", n = 2611773)


write.csv(CpGs_smoking6_TCDD, "34_10x_CpG~TCDD_and_Smoke6_cat.csv", row.names = FALSE, col.names = FALSE)
CpGs_sign_sm6_TCDD <- CpGs_smoking6_TCDD[CpGs_smoking6_TCDD$P_summ_fdr <= 0.05,]
CpGs_sign_sm6_TCDD <- CpGs_sign_sm6_TCDD[abs(CpGs_sign_sm6_TCDD$R2) >= 0.7,]
CpGs_sign_sm6_TCDD <- CpGs_sign_sm6_TCDD[complete.cases(CpGs_sign_sm6_TCDD), ]

CpGs_sign_sm6_TCDD_2 <- CpGs_smoking6_TCDD[CpGs_smoking6_TCDD$P_summ_fdr_1 <= 0.05,]
CpGs_sign_sm6_TCDD_2 <- CpGs_sign_sm6_TCDD_2[abs(CpGs_sign_sm6_TCDD_2$R2) >= 0.7,]
CpGs_sign_sm6_TCDD_2 <- CpGs_sign_sm6_TCDD_2[complete.cases(CpGs_sign_sm6_TCDD_2), ]
write.csv(CpGs_sign_sm6_TCDD, "34_10x_CpG~TCDD_and_Smoke6_sign.csv", row.names = FALSE, col.names = FALSE)

###
#smoke3_cat
smoke3 <- read.csv("34_10x_CpG~Smoke3_cat.csv", sep = "\t")
smoke3$P_summ_fdr <- p.adjust(smoke3$P_summ, method = "fdr", n = length(smoke3$P_summ))
write.csv(smoke6, "34_10x_CpG~Smoke3_cat.csv", row.names = FALSE, col.names = FALSE)

smoke3$P_summ_fdr_1 <- p.adjust(smoke3$P_summ, method = "fdr", n = 2611773)
smoke3_sign <- smoke3[smoke3$P_summ_fdr <= 0.05,]
smoke3_sign <- smoke3_sign[abs(smoke3_sign$R2) >= 0.7,]
smoke3_sign <- smoke3_sign[complete.cases(smoke3_sign), ]

smoke3_sign_2 <- smoke3[smoke3$P_summ_fdr_1 <= 0.05,]
smoke3_sign_2 <- smoke3_sign_2[abs(smoke3_sign_2$R2) >= 0.7,]
smoke3_sign_2 <- smoke3_sign_2[complete.cases(smoke3_sign_2), ]

write.csv(smoke3_sign, "34_10x_CpG~Smoke3_cat.csv_sign.csv", row.names = FALSE, col.names = FALSE)

#smoke6_cat
smoke6 <- read.csv("34_10x_CpG~Smoke6_cat.csv", sep = "\t")
smoke6$P_summ_fdr <- p.adjust(smoke6$P_summ, method = "fdr", n = length(smoke6$P_summ))
write.csv(smoke6, "34_10x_CpG~Smoke6_cat.csv", row.names = FALSE, col.names = FALSE)

smoke6$P_summ_fdr_1 <- p.adjust(smoke6$P_summ, method = "fdr", n = 2611773)
smoke6_sign <- smoke6[smoke6$P_summ_fdr <= 0.05,]
smoke6_sign <- smoke6_sign[abs(smoke6_sign$R2) >= 0.7,]
smoke6_sign <- smoke6_sign[complete.cases(smoke6_sign), ]

smoke6_sign_2 <- smoke6[smoke6$P_summ_fdr_1 <= 0.05,]
smoke6_sign_2 <- smoke6_sign_2[complete.cases(smoke6_sign_2), ]
smoke6_sign_2 <- smoke6_sign_2[abs(smoke6_sign_2$R2) >= 0.7,]
smoke6_sign_2 <- smoke6_sign_2[complete.cases(smoke6_sign_2), ]

write.csv(smoke6_sign, "34_10x_CpG~Smoke6_cat.csv_sign.csv", row.names = FALSE, col.names = FALSE)

#terciles_TCDD
terc <- read.csv("34_10x_CpG~TCDD_terc_cat.csv", sep = "\t")
terc$P_summ_fdr <- p.adjust(terc$P_summ, method = "fdr", n = length(terc$P_summ))
write.csv(terc, "34_10x_CpG~TCDD_terc_cat.csv", row.names = FALSE, col.names = FALSE)

terc$P_summ_fdr_1 <- p.adjust(terc$P_summ, method = "fdr", n = 2611773)
terc_sign <- terc[terc$P_summ_fdr <= 0.05,]
terc_sign <- terc_sign[abs(terc_sign$R2) >= 0.7,]
terc_sign <- terc_sign[complete.cases(terc_sign), ]

terc_sign_2 <- terc[terc$P_summ_fdr_1 <= 0.05,]
terc_sign_2 <- terc_sign_2[abs(terc_sign_2$R2) >= 0.7,]
terc_sign_2 <- terc_sign_2[complete.cases(terc_sign_2), ]

write.csv(terc_sign, "34_10x_CpG~TCDD_terc_cat_sign.csv", row.names = FALSE, col.names = FALSE)

####terciles TCDD + smoke3

terc_smoke <- read.csv("34_10x_CpG~TCDD_terc_and_smoke3_cat.csv", sep = "\t")
terc_smoke$P_summ_fdr <- p.adjust(terc_smoke$P_summ, method = "fdr", n = length(terc_smoke$P_summ))
write.csv(terc_smoke, "34_10x_CpG~TCDD_terc_and_smoke3_cat.csv", row.names = FALSE, col.names = FALSE)

terc_smoke$P_summ_fdr_1 <- p.adjust(terc_smoke$P_summ, method = "fdr", n = 2611773)
terc_smoke_sign <- terc_smoke[terc_smoke$P_summ_fdr <= 0.05,]
terc_smoke_sign <- terc_smoke_sign[abs(terc_smoke_sign$R2) >= 0.7,]
terc_smoke_sign <- terc_smoke_sign[complete.cases(terc_smoke_sign), ]

terc_smoke_sign_2 <- terc_smoke[terc_smoke$P_summ_fdr_1 <= 0.05,]
terc_smoke_sign_2 <- terc_smoke_sign_2[abs(terc_smoke_sign_2$R2) >= 0.7,]
terc_smoke_sign_2 <- terc_smoke_sign_2[complete.cases(terc_smoke_sign_2), ]

write.csv(terc_smoke_sign, "34_10x_CpG~TCDD_terc_and_smoke3_cat_sign.csv", row.names = FALSE, col.names = FALSE)


###
setwd("./Project")
CpG_terciles <- read.csv("34_10x_Terc~CpG.csv", sep = "\t")
CpG_terciles$p_fdr <- p.adjust(CpG_terciles$TCDD_P, method = "fdr", n = length(CpG_terciles$TCDD_P))
write.csv(CpG_terciles, "34_10x_Terc~CpG.csv", row.names = FALSE, col.names = FALSE)

CpG_terciles_sign <- read.csv("34_10x_Terc~CpG_sign.csv", sep = "\t")
CpG_conc_sign <- read.csv("TCDD~CpG_sign.csv", sep = "\t")
CpGs_diox_sign_merged <- df <- merge(CpG_terciles_sign, CpG_conc_sign, by.x = 'P', by.y = 'P')

#### histograms
CpG4 <- CpGs_1$X44sp[CpGs_1$X44sp != "-"]
CpG4 <- as.character(CpG4)
CpG4 <- as.numeric(CpG4)

hist(CpG4, breaks = 100, ylim = c(0, 20000), main="ID = 44sp", xlab="B-value", col = "blue")

CpG478 <- CpGs_1$X478_up1[CpGs_1$X478_up1 != "-"]
CpG478 <- as.character(CpG478)
CpG478 <- as.numeric(CpG478)

hist(CpG478, breaks = 100, ylim = c(0, 20000), main="ID = 478_up1", xlab="B-value", col = "blue")

#ALL
CpG__ <- CpGs_1[CpGs_1 != "-"]
CpG__ <- as.character(CpG__)
CpG__ <- as.numeric(CpG__)

hist(CpG__, breaks = 100,  ylim = c(0, 700000), main="ALL", xlab="B-value", col = "blue")

#terciles
CpG_1st <- CpGs_1[which(IDs$TCDD_quartiles == 1)]
CpG_1st <- CpG_1st[CpG_1st != "-"]
CpG_1st <- as.character(CpG_1st)
CpG_1st <- as.numeric(CpG_1st)

CpG_2nd <- CpGs_1[which(IDs$TCDD_quartiles == 2)]
CpG_2nd  <- CpG_2nd[CpG_2nd  != "-"]
CpG_2nd <- as.character(CpG_2nd)
CpG_2nd  <- as.numeric(CpG_2nd )

CpG_3rd <- CpGs_1[which(IDs$TCDD_quartiles == 3)]
CpG_3rd  <- CpG_3rd[CpG_3rd != "-"]
CpG_3rd <- as.character(CpG_3rd)
CpG_3rd <- as.numeric(CpG_3rd)

hist(CpG_1st, breaks = 100,  ylim = c(0, 210000),main="First tercile", xlab="B-value", col = "blue")
hist(CpG_2nd, breaks = 100,  ylim = c(0, 210000), main="Second tercile", xlab="B-value", col = "blue")
hist(CpG_3rd, breaks = 100, ylim = c(0, 210000), main="Third tercile", xlab="B-value", col = "blue")

hist(CpG_1st, breaks = 100,  ylim = c(0, 7000000),main="First tercile", xlab="B-value", col = "blue")
hist(CpG_2nd, breaks = 100,  ylim = c(0, 7000000), main="Second tercile", xlab="B-value", col = "blue")
hist(CpG_3rd, breaks = 100, ylim = c(0, 7000000), main="Third tercile", xlab="B-value", col = "blue")

####

setwd("./Project")

TCDD_CpG_genes <- read.csv("TCDD~CpG_sign_genes_uniq.csv", header = FALSE, sep = "\t")

colnames(TCDD_CpG_genes) <- c("chr", "CpG_pos", "TCDD_R2", "type", "start_pos", "end_pos", "gene_id", "gene_type", "gene_name")
write.csv(TCDD_CpG_genes, "TCDD~CpG_sign_genes_uniq.csv", row.names = FALSE)


#### smoke and TCDD regression
smoke_TCDD <- read.csv("Smoke_and_TCDD.csv", header = TRUE, sep = ",")
TCDDs <- read.csv("RRBS_TCDD_34_selected_subjects.csv", header = TRUE, sep = ";")
smoke_TCDD$TCDD_terc <- TCDDs$TCDD_quartiles
str(smoke_TCDD)
smoke_TCDD$TCDD_terc <- smoke_TCDD$TCDD_terc
smoke_TCDD$smoking_last6months_reduced <- as.factor(smoke_TCDD$smoking_last6months_reduced)
smoke_TCDD$smoking_last6months <- as.factor(smoke_TCDD$smoking_last6months)

lin_model <- lm(smoking_last6months_reduced ~ TCDD_terc, data = smoke_TCDD)
summary(lin_model)

plot(smoke_TCDD$TCDD_terc, smoke_TCDD$smoking_last6months_reduced, 
     pch = 16, cex = 1.3, col = "blue", main = "Regression smoke3~TCDD_terc",
     xlab = "TCDD_terc", ylab = "Smoking_3")
abline(lin_model)

library(ggplot2)
plot <- ggplot(smoke_TCDD, aes(x=TCDD_terc, y=smoking_last6months_reduced, colour=factor(TCDD_terc)))
plot + stat_smooth(method=lm, fullrange=FALSE) + geom_point()


install.packages("nnet")
library("nnet")

model_2 <- multinom(smoking_last6months_reduced ~ TCDD_terc,data=smoke_TCDD, family = "binomial")
summary(model_2)
z <- summary(model_2)$coefficients/summary(model_2)$standard.errors
z
# 2-tailed z test
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

model_3 <- multinom(smoking_last6months ~ TCDD_terc,family=binomial(link='logit'),data=smoke_TCDD)
summary(model_3)
z <- summary(model_3)$coefficients/summary(model_3)$standard.errors
z
# 2-tailed z test
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

model_4 <- multinom(smoking_last6months ~ TCDD,family=binomial(link='logit'),data=smoke_TCDD)
summary(model_4)
z <- summary(model_4)$coefficients/summary(model_4)$standard.errors
z
# 2-tailed z test
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p
summary(model)

model$rsquared
model_lin <- lm(smoking_last6months_reduced ~ TCDD_terc,data=smoke_TCDD)
summary(model_lin)


####compare two models

TCDD <- read.csv("TCDD~CpG_sign.csv", header = TRUE, sep = "\t")
smoke3_sign$P <- gsub('[\',)\\(]', '', smoke3_sign$P)
merged_smoke <- merge(smoke6_sign, smoke3_sign, by.x = "P", by.y = "P")
merged_TCDD <- merge(TCDD, terc_sign, by.x = "P", by.y = "P")
merged_TCDD_Smoke6 <- merge(TCDD, smoke6_sign, by.x = "P", by.y = "P")
merged_smoke6_smTCDD <- merge(smoke6_sign, CpGs_sign_sm6_TCDD, by.x = "P", by.y = "P")
merged_TCDD_smTCDD <- merge(TCDD, CpGs_sign_sm6_TCDD, by.x = "P", by.y = "P")


terc_smoke_sign_1$P <- gsub( "\\(", "", as.character(terc_smoke_sign_1$P))
test_df <- merge(terc_smoke_sign_1, CpGs_sign_sm6_TCDD, by.x = "P", by.y = "P")

#some plots
smoke_TCDD$TCDD_terc <- as.factor(smoke_TCDD$TCDD_terc)
hist(smoke_TCDD$TCDD, col = smoke_TCDD$TCDD_terc, breaks = 5, main = "TCDD_terciles")

ggplot(smoke_TCDD, aes(x = smoke_TCDD$TCDD, fill = smoke_TCDD$TCDD_terc)) + geom_histogram(position =
                                                         "dodge")

ggplot(data=smoke_TCDD, aes(smoke_TCDD$TCDD, fill = TCDD_terc)) + 
  geom_histogram(binwidth = 1.7, position = "identity", alpha = 0.7, color = "black")+
  ggtitle("Histogram of TCDD concentration") +
  labs(x="TCDD concentration", y="Count")+
  scale_fill_manual(values=c("lightgreen", "lightblue", "pink"))+
  theme_classic()

smoke_TCDD$smoking_last6months_reduced <- as.factor(smoke_TCDD$smoking_last6months_reduced)
ggplot(data=smoke_TCDD, aes(smoke_TCDD$TCDD, fill = smoking_last6months_reduced)) + 
  geom_histogram(binwidth = 1,  alpha = 0.7, color = "black")+
  ggtitle("Histogram of TCDD concentration") +
  labs(x="TCDD concentration", y="Count")+
  #scale_fill_manual(values=c("lightgreen", "lightblue", "pink"))+
  theme_classic()

hist(as.numeric(as.character(smoke_TCDD$smoking_last6months)), breaks = 5, main = "Smoking last 6 months", col = "pink", 
     xlab = "Smoke6")

ggplot(data=smoke_TCDD, aes(smoking_last6months)) + 
  geom_bar(alpha = 0.7, color = "black", fill = "pink")+
  ggtitle("Smoking in the last 6 months") +
  labs(x="Smoke6", y="Count")+
  #scale_fill_manual(values=c("lightgreen", "lightblue", "pink"))+
  theme_classic()
ggplot(data=smoke_TCDD, aes(smoking_last6months_reduced)) + 
  geom_bar(alpha = 0.7, color = "black", fill = "pink")+
  ggtitle("Smoking last 6 months") +
  labs(x="Smoke6", y="Count")+
  #scale_fill_manual(values=c("lightgreen", "lightblue", "pink"))+
  theme_classic()
