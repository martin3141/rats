# install/load required R libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(spant, reshape2, ggpubr)

# change the working directory to the source file location
# when "Sourcing" from the RStudio GUI
if (Sys.getenv("RSTUDIO") == "1") {
  this.dir <- dirname(parent.frame(2)$ofile)
  setwd(this.dir)
}

# if figure dir doesn't exist, export plots to the current working dir
fig_dir <- file.path(dirname(getwd()), "figures")
if (!file.exists(fig_dir)) fig_dir <- getwd()

# ggplot2 theme
theme_set(theme_bw())

# set random num generator seed
set.seed(1)

# simulate an ideal spectrum that looks like normal brain
ideal_spec <- sim_brain_1h(type = "normal_v1") %>% lb(5) %>% td2fd

# normalise the spectral intensity to NAA peak height
ideal_spec <- ideal_spec / as.numeric(peak_info(ideal_spec)$height)

# number of spectra to simulate
num_spec <- 512

# SNR range to investigate
SNRS <- c(2.5, 5.0, 7.5, 10.0, 15.0, 20.0, 25.0)

################################################################################
# Frequency drift test
################################################################################
force_sim1 = TRUE
if (!exists("sim1") | force_sim1) {
  sim1 <- TRUE
  runs <- length(SNRS)
  rats_freq_sd <- vector(mode = "numeric", length = runs)
  td_freq_sd <- vector(mode = "numeric", length = runs)
  rats_phase_sd <- vector(mode = "numeric", length = runs)
  td_phase_sd <- vector(mode = "numeric", length = runs)
  
  for (n in 1:runs) {
    # ideal unshifted data
    sim_data_0 <- ideal_spec %>% rep_dyn(num_spec) %>% 
      + sim_noise(1 / (SNRS[n]), dyns = num_spec)
    
    # define frequency adjustments for shift only example
    shifts <- seq(0, 10, length.out = num_spec)
    
    sim_data_1 <- sim_data_0 %>% shift(shifts, units = "hz")
    
    # correct using rats
    rats_res <- rats(sim_data_1, get_dyns(sim_data_1, 1))
    
    # correct using td method
    td_res <- tdsr(sim_data_1, get_dyns(sim_data_1, 1))
    #td_res <- tdsr(sim_data_1, xlim = c(15, -5))
    
    # create a dataframe of true and estimated shift values for SNR = 5
    if (n == 2) {
      shifts_df_snr <- data.frame(shifts, RATS = as.numeric(rats_res$shifts),
                              TDSR = as.numeric(td_res$shifts))
    
      shifts_df_snr <- melt(shifts_df_snr, "shifts", variable.name = "Method",
                        value.name = "estimate")
    
      shifts_df_snr$err <- shifts_df_snr$shifts - shifts_df_snr$estimate
    }
   
    rats_freq_sd[n]   <- sd(shifts - rats_res$shifts)
    td_freq_sd[n] <- sd(shifts - td_res$shifts)
    
    rats_phase_sd[n]   <- sd(rats_res$phases)
    td_phase_sd[n] <- sd(td_res$phases)
  }
}

freq_df <- data.frame(SNR = SNRS, RATS = rats_freq_sd, TDSR = td_freq_sd)
freq_df <- melt(freq_df, "SNR", variable.name = "Method", value.name = "error")

phase_df <- data.frame(SNR = SNRS, RATS = rats_phase_sd, TDSR = td_phase_sd)
phase_df <- melt(phase_df, "SNR", variable.name = "Method", value.name = "error")

freq <- ggline(freq_df, "SNR", "error", color = "Method", palette = "jco",
       ylab = "Frequency error (Hz)", numeric.x.axis = T) 

phase <- ggline(phase_df, "SNR", "error", color = "Method", palette = "jco",
       ylab = "Phase error (degrees)", numeric.x.axis = T)

lab_font <- list(size = 20, color = "black",
                 face = "bold", family = NULL)


shift_scatter <- ggscatter(shifts_df_snr, x = "shifts", y = "estimate",
                 col = "Method", ylab = "Estimated frequency shift (Hz)",
                 xlab = "True frequency shift (Hz)", palette = "jco", shape = 1) +
                 geom_abline(col = "black")

# FIGURE 1
setEPS()
postscript(file.path(fig_dir, "figure1.eps"), width = 10, height = 3.5)
fig1 <- ggarrange(freq, phase, shift_scatter, labels = "auto", 
                  font.label = lab_font, ncol = 3)
print(fig1)
dev.off()

################################################################################
# Frequency drift and random phase test
################################################################################
set.seed(2)
force_sim2 = TRUE
if (!exists("sim2") | force_sim2) {
  sim2 <- TRUE
  runs <- length(SNRS)
  rats_freq_sd <- vector(mode = "numeric", length = runs)
  td_freq_sd <- vector(mode = "numeric", length = runs)
  rats_phase_sd <- vector(mode = "numeric", length = runs)
  td_phase_sd <- vector(mode = "numeric", length = runs)
  
  for (n in 1:runs) {
    # ideal unshifted data
    sim_data_0 <- ideal_spec %>% rep_dyn(num_spec) %>% 
      + sim_noise(1 / (SNRS[n]), dyns = num_spec)
    
    # define frequency adjustments for shift only example
    shifts <- seq(0, 10, length.out = num_spec)
    
    sim_data_1 <- sim_data_0 %>% shift(shifts, units = "hz")
    
    phases <- c(0, runif(num_spec - 1, -180, 180))
    sim_data_2 <- sim_data_1 %>% phase(phases)
    
    # correct using rats
    rats_res <- rats(sim_data_2, get_dyns(sim_data_2, 1))
    
    # correct using td method
    td_res <- tdsr(sim_data_2, get_dyns(sim_data_2, 1))
    #td_res <- tdsr(sim_data_2, xlim = c(15, -5))
    
    rats_freq_sd[n]   <- sd(shifts - rats_res$shifts)
    td_freq_sd[n] <- sd(shifts - td_res$shifts)
    
    rats_phase_err <- phases - rats_res$phases
    # correct for multiples of 2pi
    rats_phase_err <- rats_phase_err - round(rats_phase_err / 180) * 180
    rats_phase_sd[n]   <- sd(rats_phase_err)
    
    td_phase_err <- phases - td_res$phases
    # correct for multiples of 2pi
    td_phase_err <- td_phase_err - round(td_phase_err / 180) * 180
    td_phase_sd[n] <- sd(td_phase_err)
  }
}

freq_df <- data.frame(SNR = SNRS, RATS = rats_freq_sd, TDSR = td_freq_sd)
freq_df <- melt(freq_df, "SNR", variable.name = "Method", value.name = "error")

phase_df <- data.frame(SNR = SNRS, RATS = rats_phase_sd, TDSR = td_phase_sd)
phase_df <- melt(phase_df, "SNR", variable.name = "Method", value.name = "error")

freq <- ggline(freq_df, "SNR", "error", color = "Method", palette = "jco",
       ylab = "Frequency error (Hz)", numeric.x.axis = T) 

phase <- ggline(phase_df, "SNR", "error", color = "Method", palette = "jco",
       ylab = "Phase error (degrees)", numeric.x.axis = T)

lab_font <- list(size = 20, color = "black",
                 face = "bold", family = NULL)

# FIGURE 2
setEPS()
postscript(file.path(fig_dir, "figure2.eps"), width = 6.67, height = 3.5)
fig2 <- ggarrange(freq, phase, labels = "auto", 
                  font.label = lab_font, ncol = 2)
print(fig2)
dev.off()

#lw <- 0.35
setEPS()
postscript(file.path(fig_dir, "figure_s1.eps"), width = 6.0, height = 9.0)
stackplot(get_dyns(sim_data_2, seq(from = 1, to = num_spec, by = 64)),
          xlim = c(4, 0.5), y_offset = 210)
dev.off()

################################################################################
# Unstable residual water test
################################################################################

# simulate a water resonance
water_ideal <- sim_resonances(freq = 4.65, amp = 0.5, lw = 10)

set.seed(3)
water_phases <- c(0, runif(num_spec - 1, -180, 180))
water_sig <- water_ideal %>% rep_dyn(num_spec) %>% phase(water_phases)

# ideal unshifted data
SNR <- 15
sim_data_0 <- ideal_spec %>% rep_dyn(num_spec) %>% 
  + sim_noise(1 / SNR, dyns = num_spec)

# define frequency adjustments for shift only example
shifts <- seq(0, 10, length.out = num_spec)

sim_data_1 <- sim_data_0 %>% shift(shifts, units = "hz")

set.seed(4)
phases <- c(0, runif(num_spec - 1, -180, 180))
sim_data_2 <- sim_data_1 %>% phase(phases)

sim_data_3 <- sim_data_2 + water_sig

# correct using rats
rats_res <- rats(sim_data_3, get_dyns(sim_data_3, 1))

# correct using td method
td_res <- tdsr(sim_data_3, get_dyns(sim_data_3, 1))

# create a dataframe of true and estimated shift values
shifts_df <- data.frame(shifts, RATS = as.numeric(rats_res$shifts),
                        TDSR = as.numeric(td_res$shifts))

shifts_df <- melt(shifts_df, "shifts", variable.name = "Method",
                  value.name = "estimate")

shifts_df$err <- shifts_df$shifts - shifts_df$estimate

shift_box <- ggboxplot(shifts_df, x = "Method", y = "err",
                 ylab = "True frequency shift - estimate (hz)") + 
                 theme(plot.margin = unit(c(2,1,1,2), "lines"))

# create a dataframe of true and estimated phase values
phase_df <- data.frame(phases, RATS = as.numeric(rats_res$phases),
                 TDSR = as.numeric(td_res$phases))

phase_df <- melt(phase_df, "phases", variable.name = "Method",
                  value.name = "estimate")

phase_df$err <- phase_df$phases - phase_df$estimate
phase_df$err <- phase_df$err - round(phase_df$err / 180) * 180

phase_box <- ggboxplot(phase_df, x = "Method", y = "err",
                 ylab = "True phase - estimate (degrees)") + 
                 theme(plot.margin = unit(c(2,1,1,2), "lines"))

# FIGURE 3
setEPS()
postscript(file.path(fig_dir, "figure3.eps"), width = 8, height = 4)
fig3 <- ggarrange(shift_box, phase_box, labels = "auto", font.label = lab_font)
print(fig3)
dev.off()

# FIGURE 4
lw <- 0.35
subset <- seq(from = 1, to = 512, length.out = 32)
setEPS()
postscript(file.path(fig_dir, "figure4.eps"), width = 6, height = 2.5)

par(mfrow = c(1, 3))
stackplot(get_dyns(sim_data_3, subset), xlim = c(4, 0.5), y_offset = 0,
          restore_def_par = FALSE, lwd = lw)

mtext("a", at = 4, padj = 0.5, font = 2)

stackplot(get_dyns(td_res$corrected, subset), xlim = c(4, 0.5), y_offset = 0,
          restore_def_par = FALSE, lwd = lw)

mtext("b", at = 4, padj = 0.5, font = 2)

stackplot(get_dyns(rats_res$corrected, subset), xlim = c(4, 0.5), y_offset = 0,
          restore_def_par = FALSE, lwd = lw)

mtext("c", at = 4, padj = 0.5, font = 2)
dev.off()

# moderate baseline issue

num_spec <- 32

# simulate a water resonance
water_ideal <- sim_resonances(freq = 4.65, amp = 0.2, lw = 10)

set.seed(3)
water_phases <- c(0, runif(num_spec - 1, -180, 180))
water_sig <- water_ideal %>% rep_dyn(num_spec) %>% phase(water_phases)

# ideal unshifted data
SNR <- 15
sim_data_0 <- ideal_spec %>% rep_dyn(num_spec) %>% 
  + sim_noise(1 / SNR, dyns = num_spec)

# define frequency adjustments for shift only example
shifts <- seq(0, 10, length.out = num_spec)

sim_data_1 <- sim_data_0 %>% shift(shifts, units = "hz")

set.seed(4)
phases <- c(0, runif(num_spec - 1, -180, 180))
sim_data_2 <- sim_data_1 %>% phase(phases)

sim_data_3 <- sim_data_2 + water_sig

# correct using rats
rats_res <- rats(sim_data_3, get_dyns(sim_data_3, 1))

# correct using td method
td_res <- tdsr(sim_data_3, get_dyns(sim_data_3, 1))

# create a dataframe of true and estimated shift values
shifts_df <- data.frame(shifts, RATS = as.numeric(rats_res$shifts),
                        TDSR = as.numeric(td_res$shifts))

shifts_df <- melt(shifts_df, "shifts", variable.name = "Method",
                  value.name = "estimate")

shifts_df$err <- shifts_df$shifts - shifts_df$estimate

shift_box <- ggboxplot(shifts_df, x = "Method", y = "err",
                 ylab = "True frequency shift - estimate (hz)") + 
                 theme(plot.margin = unit(c(2,1,1,2), "lines"))

# create a dataframe of true and estimated phase values
phase_df <- data.frame(phases, RATS = as.numeric(rats_res$phases),
                 TDSR = as.numeric(td_res$phases))

phase_df <- melt(phase_df, "phases", variable.name = "Method",
                  value.name = "estimate")

phase_df$err <- phase_df$phases - phase_df$estimate
phase_df$err <- phase_df$err - round(phase_df$err / 180) * 180

phase_box <- ggboxplot(phase_df, x = "Method", y = "err",
                 ylab = "True phase - estimate (degrees)") + 
                 theme(plot.margin = unit(c(2,1,1,2), "lines"))

# FIGURE S1
lw <- 0.35
setEPS()
postscript(file.path(fig_dir, "figure_s2.eps"), width = 6, height = 2.5)

par(mfrow = c(1, 3))
stackplot(sim_data_3, xlim = c(4, 0.5), y_offset = 0,
          restore_def_par = FALSE, lwd = lw)

mtext("a", at = 4, padj = 0.5, font = 2)

stackplot(td_res$corrected, xlim = c(4, 0.5), y_offset = 0,
          restore_def_par = FALSE, lwd = lw)

mtext("b", at = 4, padj = 0.5, font = 2)

stackplot(rats_res$corrected, xlim = c(4, 0.5), y_offset = 0,
          restore_def_par = FALSE, lwd = lw)

mtext("c", at = 4, padj = 0.5, font = 2)
dev.off()

################################################################################
# GSH MEGA-PRESS example pipeline
################################################################################

# load raw_data object
load("gsh_mpress_data.Rda")

# assume all edit off scans make up the first half of the dynamics followed
# by the edit on scans
proc <- interleave_dyns(raw_data)

# split ed-on ed-off data
odd_data <- get_odd_dyns(proc)
even_data <- -get_even_dyns(proc) # undo the phase flip
median_spec <- median_dyns(append_dyns(odd_data, even_data))

# find the scan closest to the median
median_odd <-  median_dyns(odd_data)
median_even <- median_dyns(even_data)
odd_diff <- odd_data - rep_dyn(median_odd, dyns(odd_data))
even_diff <- even_data - rep_dyn(median_even, dyns(even_data))
odd_min <- which.min(int_spec(odd_diff, xlim = c(4, 0.5), mode = "mod"))
even_min <- which.min(int_spec(even_diff, xlim = c(4, 0.5), mode = "mod"))
odd_ref <- get_dyns(odd_data, odd_min)
even_ref <- get_dyns(even_data, even_min)

# align data to the dynamic closest to the median scans (rats)
odd_data_rats_corr_res <-  rats(odd_data,  odd_ref)
even_data_rats_corr_res <- rats(even_data, even_ref)
odd_data_rats_corr <- odd_data_rats_corr_res$corrected
even_data_rats_corr <- even_data_rats_corr_res$corrected

# calcuate the mean of the corrected data (rats)
mean_even_rats <- mean_dyns(even_data_rats_corr)
mean_odd_rats  <- mean_dyns(odd_data_rats_corr)

# align data to the dynamics closest to the median scans (tdsr)
odd_data_td_corr_res <-  tdsr(odd_data,  odd_ref)
even_data_td_corr_res <- tdsr(even_data, even_ref)
odd_data_td_corr <- odd_data_td_corr_res$corrected
even_data_td_corr <- even_data_td_corr_res$corrected

# calcuate the mean of the corrected data (tdsr)
mean_even_td <- mean_dyns(even_data_td_corr)
mean_odd_td  <- mean_dyns(odd_data_td_corr)

# even_mean_corr_rats <- rats(mean_even_rats, mean_odd_rats, xlim = c(2.2, 1.8), p_deg = -1)$corrected
even_mean_corr_rats <- rats(mean_even_rats, mean_odd_rats, xlim = c(2.2, 1.8))$corrected
even_mean_corr_td <- tdsr(mean_even_td, mean_odd_td,
                                xlim = c(2.2, 1.8))$corrected

uncorr_ed <- append_dyns(odd_data, -even_data) %>% mean_dyns()
uncorr_ed <- uncorr_ed * 2 # to maintain the same scaling as corrected data

# subtract edited scan and apply line-broading and zero-filling

rats_ed <- mean_odd_rats - even_mean_corr_rats
rats_ed_orig <- rats_ed
rats_ed <- rats_ed %>% lb(3.0) %>% zf(4)
td_ed <- mean_odd_td - even_mean_corr_td
td_ed_orig <- td_ed
td_ed <- td_ed %>% lb(3.0) %>% zf(4)
uncorr_ed_orig <- uncorr_ed
uncorr_ed <- uncorr_ed %>% lb(3.0) %>% zf(4)

################################################################################
# GSH fitting (TARQUIN needs to be on the system path for this to work...)
# version 4.3.11 used to generate paper results
################################################################################

gsh <- get_uncoupled_mol("GSH", 2.95, "1H", 1, 8, 0)
p275 <- get_uncoupled_mol("P275", 2.74, "1H", 1, 2, 0)
p270 <- get_uncoupled_mol("P270", 2.68, "1H", -1, 2, 0)
p261 <- get_uncoupled_mol("P261", 2.62, "1H", 1, 2, 0)
p259 <- get_uncoupled_mol("P259", 2.58, "1H", 1, 2, 0)
p251 <- get_uncoupled_mol("P251", 2.51, "1H", -1, 2, 0)
p245 <- get_uncoupled_mol("P245", 2.45, "1H", 1, 2, 0)
p235 <- get_uncoupled_mol("P236", 2.37, "1H", 1, 2, 0)

mol_list <- list(gsh, p275, p270, p261, p259, p251, p245, p235)

basis <- sim_basis(mol_list)

opts <- "--rescale_lcm_basis false --auto_ref false --auto_phase false --start_pnt 10 --max_iters 100"

rats_ed_res <- fit_mrs(rats_ed_orig, method = "TARQUIN", basis = basis,
                       opts = opts)

td_ed_res <- fit_mrs(td_ed_orig, method = "TARQUIN", basis = basis, opts = opts)

uncorr_ed_res <- fit_mrs(uncorr_ed_orig, method = "TARQUIN", basis = basis,
                         opts = opts)

print(uncorr_ed_res$res_tab$GSH / rats_ed_res$res_tab$GSH * 100)
print(td_ed_res$res_tab$GSH / rats_ed_res$res_tab$GSH * 100)

# FIGURE 5
setEPS()
postscript(file.path(fig_dir, "figure5.eps"), width = 8, height = 8)

# six part figure
par(mfrow = c(2, 3))
mar <- c(3, 1.5, 1, 1)

plot(uncorr_ed, xlim = c(4, 1.0), restore_def_par = FALSE)
points(2, 0, cex = 10, col = "red")
mtext("a", at = 4, padj = 0.5, cex = 1.5, font = 2)

plot(td_ed, xlim = c(4, 1.0), restore_def_par = FALSE)
points(2, 0, cex = 10, col = "red")
mtext("b", at = 4, padj = 0.5, cex = 1.5, font = 2)

plot(rats_ed, xlim = c(4, 1.0), restore_def_par = FALSE)
points(2, 0, cex = 10, col = "red")
mtext("c", at = 4, padj = 0.5, cex = 1.5, font = 2)

plot(uncorr_ed_res, xlim = c(3.3,2), restore_def_par = FALSE, mar = mar)
mtext("d", at = 3.37, padj = 0.9, cex = 1.5, font = 2)

plot(td_ed_res, xlim = c(3.3,2), restore_def_par = FALSE, mar = mar)
mtext("e", at = 3.37, padj = 0.9, cex = 1.5, font = 2)

plot(rats_ed_res, xlim = c(3.3,2), restore_def_par = FALSE, mar = mar)
mtext("f", at = 3.37, padj = 0.9, cex = 1.5, font = 2)
dev.off()

################################################################################
# quick example to estimate the typical time needed to correct 128 averages
################################################################################

sim_data_0 <- ideal_spec %>% rep_dyn(128) %>% 
  + sim_noise(1 / (10), dyns = 128)

# define frequency adjustments for shift only example
shifts <- seq(0, 10, length.out = 128)

sim_data_1 <- sim_data_0 %>% shift(shifts, units = "hz")

# correct using rats
system.time(rats(sim_data_1))

# correct using td method
system.time(tdsr(sim_data_1))
