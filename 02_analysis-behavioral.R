# Environment -------------------------------------------------------------

rm(list = ls())

# Packages ----------------------------------------------------------------

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggeffects)
library(effects)
library(codingMatrices)
library(broom.mixed)
library(flextable)
library(ftExtra)
library(ggeffects)
library(gghalves)

devtools::load_all()

select <- dplyr::select

# Functions ---------------------------------------------------------------

gg_theme <- function(base_size = 20){
    theme_minimal(base_size = base_size)
}

custom_tab <- function(data){
    data |> 
        flextable() |> 
        autofit() |> 
        theme_vanilla() |> 
        colformat_double(digits = 3)
}

exclude_ran <- function(data){
    cols <- c("effect", "term", "estimate")
    data[data$effect == "ran_pars", !names(data) %in% cols] <- NA
    data
}

# Importing Data ----------------------------------------------------------

dat <- read_rds("data/clean/behavioral/dat_clean.rds")

# Subjects of the CDA analysis --------------------------------------------

cda_subj <- c("sj_0003", "sj_0005", "sj_0007", "sj_0008", "sj_0010", "sj_0012", "sj_0014", "sj_0015", "sj_0017", "sj_0018", "sj_0019", "sj_0020", "sj_0022", "sj_0024")

cda_subj <- parse_number(cda_subj)

# MODELS ------------------------------------------------------------------

## acc ~ pas --------------------------------------------------------------

datx <- dat |> 
  filter(trial_type == "valid") |> 
  group_by(subject) |> 
  mutate(trial = 1:n()) |> 
  ungroup()

contrasts(dat$pasf_sdiff) <- codingMatrices::contr.diff(4)

fit_acc_pas <- glmer(test_acc ~ pasf_sdiff + (1|subject),
                     data = datx,
                     family = binomial(link = "logit"),
                     contrasts = list(pasf_sdiff = codingMatrices::contr.diff(4)),
                     control = glmerControl(optimizer="bobyqa"))

fit_acc_pas_cda <- glmer(test_acc ~ pasf_sdiff + (1|subject),
                     data = filter(datx, subject %in% cda_subj),
                     family = binomial(link = "logit"),
                     contrasts = list(pasf_sdiff = codingMatrices::contr.diff(4)),
                     control = glmerControl(optimizer="bobyqa"))

## sdt model pas ----------------------------------------------------------

# using the glm to compute the dprime the intepretation of the criterion
# is flipped compared to the psycho::dprime() function
# negative --> tendency to say NO
# positive --> tendency to say YES

dat$is_signal <- factor(ifelse(dat$is_signal, 1, 0))
dat$say_signal <- ifelse(dat$say_signal, 1, 0)

fit_sdt_pas <- glmer(say_signal ~ is_signal + (is_signal|subject),
                     data = dat,
                     family = binomial(link = "probit"),
                     contrasts = list(is_signal = contr.sum(2)/2*-1))

summary(fit_sdt_pas)

# TABLES ------------------------------------------------------------------

table_fit_acc_pas <- fit_acc_pas |> 
    tidy(conf.int = TRUE, exponentiate = TRUE) |> 
    select(effect, term, starts_with("estimate"), everything(), -group) |> 
    mutate(
        term = case_when(
            term == "(Intercept)" ~ "PAS 1 vs 50%",
            term == "pasf_sdiff2-1" ~ "PAS 2 vs PAS 1",
            term == "pasf_sdiff3-2" ~ "PAS 3 vs PAS 2",
            term == "pasf_sdiff4-3" ~ "PAS 4 vs PAS 3",
            term == "sd__(Intercept)" ~ "$\\sigma_{id}$"
        ),
        p.value = ifelse(p.value < 0.001, "$p < 0.001$", sprintf("$p = %s$", round(p.value, 3))),
        ci = sprintf("[%s, %s]", round(conf.low, 3), round(conf.high, 3)),
        estimate = ifelse(term == "PAS 1 vs 50%", sprintf("%.3f ($\\hat p_{correct} = %.3f$)", estimate, plogis(log(estimate))), sprintf("%.3f", estimate))) |> 
    exclude_ran() |> 
    select(term, estimate, std.error, ci, statistic, p.value) |> 
    custom_tab() |> 
    set_header_labels(values = list(
        term = "Parameter",
        estimate = "$\\beta$",
        std.error = "SE",
        statistic = "z",
        p.value = "*p*",
        ci = "95% CI"
    )) |> 
    align(part = "all", align = "center") |> 
    colformat_md(part = "all")

table_fit_acc_pas_cda <- fit_acc_pas_cda |> 
    tidy(conf.int = TRUE, exponentiate = TRUE) |> 
    select(effect, term, starts_with("estimate"), everything(), -group) |> 
    mutate(
        term = case_when(
            term == "(Intercept)" ~ "PAS 1 vs 50%",
            term == "pasf_sdiff2-1" ~ "PAS 2 vs PAS 1",
            term == "pasf_sdiff3-2" ~ "PAS 3 vs PAS 2",
            term == "pasf_sdiff4-3" ~ "PAS 4 vs PAS 3",
            term == "sd__(Intercept)" ~ "$\\sigma_{id}$"
        ),
        p.value = ifelse(p.value < 0.001, "$p < 0.001$", sprintf("$p = %s$", round(p.value, 3))),
        ci = sprintf("[%s, %s]", round(conf.low, 3), round(conf.high, 3)),
        estimate = ifelse(term == "PAS 1 vs 50%", sprintf("%.3f ($\\hat p_{correct} = %.3f$)", estimate, plogis(log(estimate))), sprintf("%.3f", estimate))) |> 
    exclude_ran() |> 
    select(term, estimate, std.error, ci, statistic, p.value) |> 
    custom_tab() |> 
    set_header_labels(values = list(
        term = "Parameter",
        estimate = "$\\text{exp}(\\beta)$",
        std.error = "SE",
        statistic = "z",
        p.value = "*p*",
        ci = "95% CI"
    )) |> 
    align(part = "all", align = "center") |> 
    colformat_md(part = "all")


# PLOTS -------------------------------------------------------------------

## pas ~ trial_type -------------------------------------------------------

plot_pas_trial_type_bar <- dat |> 
    mutate(isaware = ifelse(pas == 1, "Unaware", "Aware"),
           trial_type = str_to_title(trial_type)) |> 
    ggplot(aes(x = pas)) +
    geom_bar(aes(fill = isaware), color = "black") +
    facet_wrap(~trial_type) +
    gg_theme(base_size = 20) +
    theme(#axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
    ylab("Trials") +
    scale_fill_manual(values = c("#E74C3C", "#3498DB"))

## pas ~ trial_type (average) ---------------------------------------------

plot_pas_trial_type_box <- dat |> 
    mutate(trial_type = str_to_title(trial_type)) |> 
    group_by(subject, trial_type) |> 
    mutate(tot = n()) |> 
    group_by(subject, trial_type, pas, tot) |> 
    summarise(n = n()) |> 
    mutate(p = n/tot,
           isaware = ifelse(pas == 1, "Unaware", "Aware")) |> 
    ungroup() |> 
    ggplot(aes(x = factor(pas), y = n, fill = isaware)) +
    geom_boxplot() +
    facet_grid(~trial_type) +
    gg_theme() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_fill_manual(values = c("#E74C3C", "#3498DB"))

plot_pas_trial_type <- cowplot::plot_grid(plot_pas_trial_type_bar, 
                                          plot_pas_trial_type_box, 
                                          nrow = 1, align = "h")

## acc ~ pas --------------------------------------------------------------

dat_fit_acc_pas_agg <- ggeffect(fit_acc_pas)[[1]] |> 
    data.frame() |> 
    mutate(x = paste("PAS", x))

dat_fit_acc_plot <- fit_acc_pas@frame |> 
    group_by(subject, pasf_sdiff) |> 
    summarise(acc = mean(test_acc),
              n = n()) |> 
    filter(acc > 0.2) |> 
    mutate(pasf = paste("PAS", pasf_sdiff))

plot_acc_pas <- ggplot() +
    geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.7) +
    geom_half_point(data = dat_fit_acc_plot,
                    aes(x = pasf, y = acc, fill = pasf),
                    size = 4,
                    position = position_jitter(width = 0.01, height = 0, seed = 2022),
                    colour="black",
                    pch=21) +
    ylim(c(0.3, 1)) +
    geom_pointrange(data = dat_fit_acc_pas_agg,
                    aes(x = x, 
                        y = predicted, 
                        ymin = conf.low, 
                        ymax = conf.high),
                    size = 0.8) +
    cowplot::theme_minimal_grid(font_size = 18) +
    theme(axis.title.x = element_blank(),
          legend.position = "none") +
    ylab("Accuracy")

## Saving ---------------------------------------------------------------

# FIGURE 2

ggsave("figures/figure-2.pdf", plot_pas_trial_type, width = 10, height = 5)
ggsave("figures/figure-2.png", plot_pas_trial_type, width = 8, height = 5, dpi = 300)

# FIGURE 3

ggsave("figures/figure-3.pdf", plot_acc_pas, width = 8, height = 5)
ggsave("figures/figure-3.png", plot_acc_pas, width = 8, height = 5, dpi = 300)

# TABLE 1
save_as_docx(table_fit_acc_pas, path = "tables/table1.docx")

# TABLE 3
save_as_docx(table_fit_acc_pas_cda, path = "tables/table3.docx")
