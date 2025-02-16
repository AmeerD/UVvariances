library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(tibble)
library(sampling)
library(survey)
library(purrr)
library(kableExtra)
library(latex2exp)

load("ACS Analysis/WAdat.rda")
load("ACS Analysis/ACSresults.rda")

WAdat <- WAdat %>% arrange(puma) %>% filter(income > 0) %>% 
  mutate(log.income = log(income), idx=row_number()) %>% 
  group_by(puma) %>% mutate(nc = n()) %>% ungroup

robs <- mean(WAdat$income)
lobs <- mean(WAdat$log.income)

vardf <- as_tibble(vardf)
estdf <- as_tibble(estdf)

## Create variance ratios plots
empdf <- estdf %>% 
  group_by(sample.size, design, outcome) %>% 
  summarise(empirical = var(value))

vrdf <- vardf %>% 
  left_join(empdf) %>%
  mutate(vratio = value/empirical,
         outcome = ifelse(outcome == "log", "log-income", "income")) %>%
  filter(!grepl("2b", method))

vrplot <- function(dat, d, trunc=NULL) {
  if (d == "stratified") {
    dat <- dat %>% mutate(sample.size = as.factor(pmax(ceiling(sample.size/length(unique(WAdat$puma))))))
  } else if (d == "clustered") {
    dat <- dat %>% mutate(sample.size = factor(sample.size, levels=c(30,50,100,250,500), 
                                               labels=c("15(2)", "25(2)", "50(2)", "50(5)", "50(10)")))
  } else {
    dat <- dat %>% mutate(sample.size = as.factor(sample.size))
  }
  
  alpha <- 0.5
  
  if (!is.null(trunc)) {
    dat <- dat %>%
      group_by(sample.size, method, design, outcome) %>%
      filter(vratio <= quantile(vratio, 1-trunc))
    alpha <- 0.5/(1-trunc)
  }
  
  dat %>%
    filter(design == d) %>% 
    mutate(method = case_when(
             method == "asymptotic" ~ "Asymptotic",
             method == "Hdecomp" ~ "H-decomposition",
             method == "IJ" ~ "Infinitesimal jackknife",
             TRUE ~ method
           )) %>% 
    ggplot(aes(x=sample.size, y=vratio, colour=method)) +
    geom_hline(yintercept=1) +
    geom_violin(draw_quantiles = c(alpha)) +
    facet_wrap(~outcome, scales="free") +
    xlab(ifelse(d == "stratified", "Samples per strata", 
                ifelse(d == "clustered", "Number of clusters (samples per cluster)", 
                       ifelse(d == "srswor", "Sample size", "Expected sample size")))) + 
    ylab("Variance ratio") +
    theme(legend.position="right")
}

# for (d in unique(vrdf$design)) {
#   ggsave(paste0("ACS Analysis/Figures/vr", d, ".pdf"), vrplot(vrdf, d), height=4, width=8)
#   ggsave(paste0("ACS Analysis/Figures/vr", d, "_0.05.pdf"), vrplot(vrdf, d, 0.05), height=4, width=8)
# }

lp <- '
AABB
CCDD
EEFF
'

tc <- 0.05

vrtrunc <- wrap_plots(
  vrplot(vrdf, "bernoulli", tc) + ggtitle("Bernoulli"),
  vrplot(vrdf, "poisson", tc) + ggtitle("Poisson"),
  vrplot(vrdf, "srswor", tc) + ggtitle("SRSWOR"),
  vrplot(vrdf, "stratified", tc) + ggtitle("Stratified"),
  vrplot(vrdf, "clustered", tc) + ggtitle("Two-stage cluster sampling"),
  axis_titles="collect_y"
  ) + guide_area() + plot_layout(design=lp, guides='collect') & theme(legend.position="right")
vrtrunc

ggsave("ACS Analysis/Figures/vrtrunc.pdf", vrtrunc, height=8, width=9)

tc <- NULL

vrtrunc <- wrap_plots(
  vrplot(vrdf, "bernoulli", tc) + ggtitle("Bernoulli"),
  vrplot(vrdf, "poisson", tc) + ggtitle("Poisson"),
  vrplot(vrdf, "srswor", tc) + ggtitle("SRSWOR"),
  vrplot(vrdf, "stratified", tc) + ggtitle("Stratified"),
  vrplot(vrdf, "clustered", tc) + ggtitle("Two-stage cluster sampling"),
  axis_titles="collect_y"
) + guide_area() + plot_layout(design=lp, guides='collect') & theme(legend.position="right")
vrtrunc

ggsave("ACS Analysis/Figures/vrsupp.pdf", vrtrunc, height=8, width=9)

## Create a coverage table
q80 <- qnorm(0.8 + (1-0.8)/2)
q90 <- qnorm(0.9 + (1-0.9)/2)
q95 <- qnorm(0.95 + (1-0.95)/2)

covdf <- estdf %>%
  left_join(vardf %>% select(-runtime) %>% pivot_wider(names_from=method, values_from=value)) %>%
  left_join(empdf) %>%
  select(-ends_with("2b")) %>% 
  pivot_longer(asymptotic:empirical, names_to="method", values_to="variance") %>% 
  mutate(
    cov80 = ifelse(outcome == "real", (robs >= (value - sqrt(variance)*q80)) & (robs <= (value + sqrt(variance)*q80)),
                   (lobs >= (value - sqrt(variance)*q80)) & (lobs <= (value + sqrt(variance)*q80))),
    cov90 = ifelse(outcome == "real", (robs >= (value - sqrt(variance)*q90)) & (robs <= (value + sqrt(variance)*q90)),
                   (lobs >= (value - sqrt(variance)*q90)) & (lobs <= (value + sqrt(variance)*q90))),
    cov95 = ifelse(outcome == "real", (robs >= (value - sqrt(variance)*q95)) & (robs <= (value + sqrt(variance)*q95)),
                  (lobs >= (value - sqrt(variance)*q95)) & (lobs <= (value + sqrt(variance)*q95))),
  ) %>%
  group_by(sample.size, design, outcome, method) %>%
  summarise(
    cov80 = sum(cov80)/n(),
    cov90 = sum(cov90)/n(),
    cov95 = sum(cov95)/n()
  ) %>% 
  pivot_longer(starts_with("cov"), names_to="level") %>%
  ungroup %>% 
  mutate(level = gsub("cov", "", level))
  

covplot <- function(dat, d) {
  if (d == "stratified") {
    dat <- dat %>% mutate(sample.size = as.factor(pmax(ceiling(sample.size/length(unique(WAdat$puma))))))
  } else if (d == "clustered") {
    dat <- dat %>% mutate(sample.size = factor(sample.size, levels=c(30,50,100,250,500), 
                                               labels=c("15(2)", "25(2)", "50(2)", "50(5)", "50(10)")))
  } else {
    dat <- dat %>% mutate(sample.size = as.factor(sample.size))
  }
  
  dat %>%
    filter(design == d) %>% 
    mutate(method = case_when(
      method == "asymptotic" ~ "Asymptotic",
      method == "Hdecomp" ~ "H-decomposition",
      method == "IJ" ~ "Infinitesimal jackknife",
      method == "empirical" ~ "Empirical",
      TRUE ~ method
    )) %>% 
    mutate(level = paste0("alpha == ", level)) %>%
    ggplot(aes(x=sample.size, y=value, group=method, colour=method)) +
    geom_hline(yintercept=c(0,1)) +
    geom_line() + 
    geom_hline(aes(yintercept=as.numeric(gsub("alpha == ", "", level))/100), linetype=2) +
    facet_wrap(~level, nrow=1, labeller=label_parsed) +
    xlab(ifelse(d == "stratified", "Samples per strata", 
                ifelse(d == "clustered", "Number of clusters (samples per cluster)", 
                       ifelse(d == "srswor", "Sample size", "Expected sample size")))) + 
    ylab("Coverage probability") +
    theme(legend.position="bottom")
}

lp = '
AAABBB
CCCDDD
EEEEFF
'

covlog <- wrap_plots(
  covplot(covdf %>% filter(outcome == "log"), "bernoulli") + ggtitle("Bernoulli"),
  covplot(covdf %>% filter(outcome == "log"), "poisson") + ggtitle("Poisson"),
  covplot(covdf %>% filter(outcome == "log"), "srswor") + ggtitle("SRSWOR"),
  covplot(covdf %>% filter(outcome == "log"), "stratified") + ggtitle("Stratified"),
  covplot(covdf %>% filter(outcome == "log"), "clustered") + ggtitle("Two-stage cluster sampling")
) + guide_area() + plot_layout(design=lp, guides='collect', axis_titles="collect_y") & theme(legend.position="right")
covlog

ggsave("ACS Analysis/Figures/covlog.pdf", covlog, height=8, width=9)

covreal <- wrap_plots(
  covplot(covdf %>% filter(outcome == "real"), "bernoulli") + ggtitle("Bernoulli"),
  covplot(covdf %>% filter(outcome == "real"), "poisson") + ggtitle("Poisson"),
  covplot(covdf %>% filter(outcome == "real"), "srswor") + ggtitle("SRSWOR"),
  covplot(covdf %>% filter(outcome == "real"), "stratified") + ggtitle("Stratified"),
  covplot(covdf %>% filter(outcome == "real"), "clustered") + ggtitle("Two-stage cluster sampling")
) + guide_area() + plot_layout(design=lp, guides='collect', axis_titles="collect_y") & theme(legend.position="right")
covreal

ggsave("ACS Analysis/Figures/covreal.pdf", covreal, height=8, width=9)

## Create a plot of the ratio between the estimators with and without tau2b
tau2b <- vardf %>% 
  filter(grepl("Hdecomp", method)) %>%
  select(-runtime) %>%
  pivot_wider(names_from=method, values_from=value) %>% 
  mutate(ratio = Hdecomp/Hdecomp.2b,
         outcome = ifelse(outcome == "log", "log-income", "income"),
         sample.size=as.factor(sample.size),
         design = case_when(
           design == "bernoulli" ~ "Bernoulli",
           design == "poisson" ~ "Poisson",
           design == "srswor" ~ "SRSWOR",
           design == "stratified" ~ "Stratified",
           TRUE ~ "Two-stage cluster sampling"
         )) %>%
  group_by(sample.size, design, outcome) %>%
  summarise(lower = min(ratio), upper = max(ratio)) %>%
  ggplot(aes(x=sample.size)) +
  geom_hline(yintercept=1) +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="seagreen") +
  facet_grid(outcome~design, scales="free_x") +
  ylim(c(0.9998,1.0002)) +
  theme(text = element_text(size = 8)) +
  xlab("Sample Size") + ylab("Ratio between (14) and (15)") 
tau2b
ggsave("ACS Analysis/Figures/tau2b.pdf", tau2b, height=4, width=8)


## Run Shapiro-Wilk normality tests on each set of GREG estimates
swplot <- estdf %>% 
  group_by(sample.size, design, outcome) %>% 
  summarise(log.sw=log(shapiro.test(value)$p.value)) %>%
  mutate(sample.size=as.factor(sample.size),
         design = case_when(
           design == "bernoulli" ~ "Bernoulli",
           design == "poisson" ~ "Poisson",
           design == "srswor" ~ "SRSWOR",
           design == "stratified" ~ "Stratified",
           TRUE ~ "Two-stage cluster sampling"
         ),
         outcome = ifelse(outcome == "log", "log-income", "income")) %>%
  ggplot(aes(x=sample.size, y=log.sw, colour=design)) +
  geom_point() +
  geom_hline(yintercept=log(0.05)) +
  facet_wrap(~outcome) +
  xlab("Sample Size") + ylab("P-value (log scale)") +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(nrow=2))

ggsave("ACS Analysis/Figures/swtest.pdf", swplot, height=4, width=8)


## Runtime results
vardf %>%
  filter(method != "Hdecomp.2b") %>%
  group_by(sample.size, method) %>%
  summarise(mi=min(runtime), ma=max(runtime)) %>%
  mutate(rt = paste0(round(mi, digits=3), "-", round(ma, digits=3)), .keep="unused") %>%
  pivot_wider(names_from=method, values_from=rt) %>% 
  select(sample.size, asymptotic, Hdecomp, IJ) %>%
  kable(col.names=c("Sample Size", "Asymptotic", "H-decomposition", "Infinitesimal Jackknife"),
        format="latex")










