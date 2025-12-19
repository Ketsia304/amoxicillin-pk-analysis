# ============================================================
# Amoxicillin Non-Compartmental Pharmacokinetic Analysis
# ============================================================

# ---- Libraries ----
library(tidyverse)
library(DescTools)
library(PKNCA)
library(kableExtra)

# ---- Data Import ----
amoxicillin <- read_csv('amoxicillin.csv')

# ---- Summary Characteristics by Individual ----
pk_individual <- amoxicillin %>%
  group_by(ID)%>%
  summarise(
    cmin=min(conc),
    cmax=max(conc),
    tmin=time[which.min(conc)],
    tmax=time[which.max(conc)],
    median_conc=median(conc),
    iqr_conc=IQR(conc),
    sd_conc= sd(conc),
    .groups = 'drop'
    )

pk_population <- pk_individual %>%
  summarise(
    mean_cmin = mean(cmin),
    mean_cmax = mean(cmax),
    mean_sd   = mean(sd_conc)
  )

# ---- Concentration–Time Plots ----
p_linear <- ggplot(amoxicillin, aes(time,conc, group=ID))+
  geom_point()+
  labs(
    title = "Amoxicillin Concentration-Time Profiles", 
    x= "Time (hours)", 
    y= "Concentration (mg/L)")+
  geom_line(alpha=0.7)+
  theme_bw()+
  facet_wrap(~ID)

p_linear_grouped <- ggplot(amoxicillin, aes(x=time, y=conc, group=ID))+
  geom_point(aes(colour=ID))+
  labs(title = "Concentration vs Time graph", x= "Time(hrs)", y= "concentration(mg)")+
  geom_line(aes(colour=ID))+
  theme_bw()

p_log <- ggplot(amoxicillin, aes(time, conc, group=ID, colour = ID))+
  geom_point()+
  scale_y_log10()+
  labs(
    title = "Log-Scale Concentration-Time Profiles", 
    x= "Time (hours)", 
    y= "Log Concentration (mg/L)")+
  geom_line()+
  theme_bw()

ggsave("amox_ct_linear.pdf", p_linear, width = 7, height = 5)
ggsave('amox_ct_linear_grouped.pdf', p_linear_grouped, width = 7, height = 5)
ggsave("amox_ct_log.pdf", p_log, width = 7, height = 5)

# ---- AUC (0–8h) ----
AUC_by_ID <- amoxicillin %>%
  group_by(ID) %>%
  summarise(
    AUC0_8=AUC(time,conc,method="trapezoid"),
    .groups = 'drop')

AUC_summary <- AUC_by_ID %>%
  summarise(
    Mean_AUC=mean(AUC0_8), 
    SD_AUC=SD(AUC0_8))

AUC_table <- AUC_by_ID %>%
  kbl(booktabs=TRUE,caption="AUC (0-8h) by Individual") %>%
  kable_styling(bootstrap_options = "striped")

save_kable(
  AUC_table,
  file = "auc_by_individual.html")
# ---- PKNCA Analysis ----
DoseAmox <- amoxicillin %>%
  subset(time==0) %>%
  mutate(dose=1000)

amox_conc <- PKNCAconc(amoxicillin, conc~time|ID)
amox_dose <- PKNCAdose(DoseAmox, dose~time|ID)

amox_pkdata <- PKNCAdata(
  amox_conc,
  amox_dose,
  options = list(auc.method = "linear")
)

amox_nca <- pk.nca(amox_pkdata)

summary(amox_nca)

nca_results <- amox_nca$result

rm(dose)
