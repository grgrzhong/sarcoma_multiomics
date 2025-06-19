library("survminer")
source(here::here("scripts/lib/study_lib.R"))

require("survival")
fit <- survfit(Surv(time, status) ~ sex, data = lung)
print(fit)
ggsurvplot(fit, data = lung, pval = TRUE)
ggsurvplot(fit, data = lung, censor.shape="|", censor.size = 4)

ggsurvplot(
    fit,
    data = lung,
    size = 1, # change line size
    palette =
        c("#E7B800", "#2E9FDF"), # custom color palettes
    conf.int = TRUE, # Add confidence interval
    pval = TRUE, # Add p-value
    risk.table = TRUE, # Add risk table
    risk.table.col = "strata", # Risk table color by groups
    legend.labs =
        c("Male", "Female"), # Change legend labels
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw() # Change ggplot2 theme
)

ggsurvplot(
    fit, # survfit object with calculated statistics.
    data = lung, # data used to fit survival curves.
    risk.table = TRUE, # show risk table.
    pval = TRUE, # show p-value of log-rank test.
    conf.int = TRUE, # show confidence intervals for
    # point estimates of survival curves.
    xlim = c(0, 500), # present narrower X axis, but not affect
    # survival estimates.
    xlab = "Time in days", # customize X axis label.
    break.time.by = 100, # break X axis in time intervals by 500.
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
)
