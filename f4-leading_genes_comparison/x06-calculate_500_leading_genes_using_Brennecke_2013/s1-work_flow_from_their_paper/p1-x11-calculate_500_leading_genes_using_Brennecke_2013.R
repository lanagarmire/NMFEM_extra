stop()
q('no')

library(ggplot2)
library(dplyr)
library(statmod)

epitSE2 <- readRDS('epitSE2.RDS')

head(assay(epitSE2))
colData(epitSE2)

sf <- colMeans(assays(epitSE2)$counts) / colMeans(assays(epitSE2)$normalized)

rw <- assays(epitSE2)$normalized
rwl <- zx::log_trans(rw)


##### cv plot & analysis #####
#==== for their sample data ====#
# rw <- read.table('nmeth.2645-S10.csv', sep=',', header=T, row.names=1)
# sample(rownames(rw), 20)
# rw <- rw[grep('^ENSG', row.names(rw)),]
# all(apply(rw, 1, function(x)any(x>0)))
# rw <- rw[apply(rw, 1, function(x)any(x>0)),]
#===============================#

m_vs_cv2 <- data.frame(gene=rownames(rw), m=apply(rw, 1, mean), v =apply(rw, 1, var)) %>%
  mutate(cv2=v/m^2)

m_vs_cv2 %>% 
  arrange(cv2) %>%
  DataFrame()

min_mean_for_fit <- quantile(m_vs_cv2[m_vs_cv2$cv2 > 5,]$m, 0.95)
m_vs_cv2_f <- m_vs_cv2 %>% filter(m > min_mean_for_fit)

fit <- glmgam.fit(data.frame(a0=1, a1tilde=1/m_vs_cv2_f$m), m_vs_cv2_f$cv2)

xi <- mean(1/sf)

a0 <- fit$coefficients['a0']
a1 <- fit$coefficients['a1tilde'] - xi

df <- ncol(rw) - 1

fitted_curve_ggdat <- data.frame(x=10^seq(-1, 2, length.out=100)) %>%
  mutate(y=(xi+a1)/x + a0) %>%
  mutate(ylow=y*qchisq(0.999, df)/df) %>%
  mutate(yhigh=y*qchisq(0.001, df)/df)

psia1theta <- mean(1/sf) + a1
minBiolDisp <- .5^2
m <- ncol(rw)
cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
testDenom <- ( m_vs_cv2$m * psia1theta + m_vs_cv2$m^2 * cv2th ) / ( 1 + cv2th/m )
p <- 1 - pchisq( m_vs_cv2$v * (m-1) / testDenom, m-1 )
padj <- p.adjust( p, "BH" )
sig <- padj < .1
table( sig )

#==== plot with curves ====#
ggplot() +
  geom_point(data=m_vs_cv2, aes(x=m, y=cv2), alpha=0.1) +
  geom_point(data=m_vs_cv2[sig,], aes(x=m, y=cv2), color='blue', alpha=0.1) +
  scale_x_log10(name='average normalized expression counts') +
  scale_y_log10(name='squared coefficient of variance', limits=c(1e-01,1e+02)) +
  geom_line(data=fitted_curve_ggdat, aes(x=x, y=y)) +
  geom_line(data=fitted_curve_ggdat, aes(x=x, y=ylow), linetype='dashed') +
  geom_line(data=fitted_curve_ggdat, aes(x=x, y=yhigh), linetype='dashed')
#==========================#

epitSE2_B2013 <- epitSE2[sig, ]

saveRDS(epitSE2_B2013, '../s3-RDS_after_B2013/epitSE2_B2013.RDS')

leadingGenesB2013 <- rownames(rw)[sig]

saveRDS(leadingGenesB2013, '../s2-leading_genes/leadingGenesB2013.RDS')
###############################


