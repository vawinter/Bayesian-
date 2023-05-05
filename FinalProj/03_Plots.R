# making plots for final proj

library(MCMCvis)
library(knitr)
library(kableExtra)

# Load in data
mr.ss.re.ab <- readRDS("../IPM/Data/ipm.rds")

# Abundance ----
png(file="../IPM/Datavis/abundance.png",
     width=12, height=12, units="in", res=100)
par(mfrow = c(2, 2))
# adult
abundance <- as.data.frame(mr.ss.re.ab$mean$abun.ad.wmu)
row.names(abundance) <- c("2020", "2021", "2022")
colnames(abundance) <- unique(full$WMU)

# boxplot
boxplot(abundance, xlab = "WMU", ylab = "Mean Nad",
        main = "Adult abundance per WMU",
        ylim = c(1, 4), col = "darkolivegreen2")

# mcmc plot
# MCMCplot(mr.ss.re.ab, params = 'abun.ad.wmu')

# line plot
matplot(abundance, type = "l", ylab = "Mean Nad",
        main = "Adult abundance across WMUs", xlim = c(1.0, 3.5),
        ylim = c(0, 4), xlab = "Year", 
        xaxt = 'n')
# Display legend
legend("right", legend = colnames(abundance), col = 1:3, lty = 1:3)
# Adding custom second axis
axis(side = 1,
     at = c(1.0, 2.0, 3.0),
     #line = 2,
     lwd = 0,
     c(2020, 2021, 2022),
     col = "steelblue")

# juvenile
abundanceJ <- as.data.frame(mr.ss.re.ab$mean$abun.juv.wmu)
row.names(abundanceJ) <- c("2020", "2021", "2022")
colnames(abundanceJ) <- unique(full$WMU)

# boxplot
boxplot(abundanceJ, xlab = "WMU", ylab = "Mean Njuv",
        main = "Juvenile abundance per WMU",
        ylim = c(1, 1.5), col = "darkorange")

# matplot
matplot(abundanceJ, type = "l", ylab = "Njuv abundance",
        main = "Juvenile abundance across WMUs", 
        xlim = c(1.0, 3.5), xaxt = 'n', xlab = "Year",
        ylim = c(0, 2))
# Display legend
legend("right", legend = colnames(abundance), col = 1:3, lty = 1:3)
# Adding custom second axis
axis(side = 1,
     at = c(1.0, 2.0, 3.0),
    # line = 2,
     lwd = 0,
     c(2020, 2021, 2022),
     col = "steelblue")
dev.off()

# Harvest ----
png(file="../IPM/Datavis/harvest.png",
     width=12, height=12, units="in", res=100)

par(mfrow = c(2, 2))
# adult
harvest <- as.data.frame(mr.ss.re.ab$mean$h.ad.wmu)
row.names(harvest) <- c("2020", "2021", "2022")
colnames(harvest) <- unique(full$WMU)

boxplot(harvest, xlab = "WMU", ylab = "Mean Had",
        main = "Adult harvest per WMU",
        col = "darkolivegreen2", ylim = c(0.1, 0.5))

#matplot(t.abun, type = "l")
matplot(harvest, type = "l", ylab = "Mean H ad",
        main = "Adult harvest across WMUs", xlim = c(1.0, 3.5),
        ylim = c(0,1), xlab = "Year",
        xaxt = 'n')
# Display legend
legend("right", legend = colnames(abundance), col = 1:3, lty = 1:3)
# Adding custom second axis
axis(side = 1,
     at = c(1.0, 2.0, 3.0),
    # line = 2,
     lwd = 0,
     c(2020, 2021, 2022),
     col = "steelblue")

# Juvenile harvest
harvestJ <- as.data.frame(mr.ss.re.ab$mean$h.juv.wmu)
row.names(harvestJ) <- c("2020", "2021", "2022")
colnames(harvestJ) <- unique(full$WMU)

# boxplot
boxplot(harvest, xlab = "WMU", ylab = "Mean Had",
        main = "Juvenile harvest per WMU",
        col = "darkorange", ylim = c(0.1, 0.5))

# matplot
matplot(harvestJ, type = "l", ylab = "Mean H juv",
        main = "Juvenile harvest across WMUs", xlim = c(1.0, 3.5),
        ylim = c(0, 1), xaxt = 'n', xlab = "Year")

# Display legend
legend("right", legend = colnames(harvestJ), col = 1:3, lty = 1:3)
# Adding custom second axis
axis(side = 1,
     at = c(1.0, 2.0, 3.0),
     #line = 1,
     lwd = 0,
     c(2020, 2021, 2022),
     col = "steelblue")
dev.off()

# Survival ----
png(file="../IPM/Datavis/survival.png",
     width=12, height=12, units="in", res=100)
par(mfrow = c(2, 2))
## adult
surv <- as.data.frame(mr.ss.re.ab$mean$s.ad.wmu)
row.names(surv) <- c("2020", "2021", "2022")
colnames(surv) <- unique(full$WMU)


boxplot(surv, xlab = "WMU", ylab = expression(paste("Mean ", phi, "ad")),
        main = "Adult survival per WMU",
        col = "darkolivegreen2")

matplot(surv, type = "l", ylab = expression(paste("Mean ", phi, "ad")),
        main = "Adult survival across WMUs", xlim = c(1.0, 3.5),
        ylim = c(0, 1), xaxt = 'n', xlab = "Year")
# Display legend
legend("right", legend = colnames(abundance), col = 1:3, lty = 1:3)
# Adding custom second axis
axis(side = 1,
     at = c(1.0, 2.0, 3.0),
     #line = 1,
     lwd = 0,
     c(2020, 2021, 2022),
     col = "steelblue")

## Juvenile
survj <- as.data.frame(mr.ss.re.ab$mean$s.juv.wmu)
row.names(survj) <- c("2020", "2021", "2022")
colnames(survj) <- unique(full$WMU)


boxplot(survj, xlab = "WMU", ylab = expression(paste("Mean ", phi, "juv")),
        main = "Juvenile survival per WMU",
        col = "darkorange")

matplot(survj, type = "l", ylab = expression(paste("Mean ", phi, "juv")),
        main = "Juvenile survival across WMUs", xlim = c(1.0, 3.5),
        ylim = c(0, 1), xaxt = 'n', xlab = "Year")
# Display legend
legend("right", legend = colnames(survj), col = 1:3, lty = 1:3)
# Adding custom second axis
axis(side = 1,
     at = c(1.0, 2.0, 3.0),
     #line = 1,
     lwd = 0,
     c(2020, 2021, 2022),
     col = "steelblue")
dev.off()



# Print summary of results
table(cat("Mean of Juvenile harvest (State):", mean(mr.ss.re.ab$sims.list$h.juv), "\n"),
cat("Mean of Juvenile harvest (WMU):", mean(mr.ss.re.ab$sims.list$h.juv.wmu), "\n"),
cat("Mean of Adult harvest (State):", mean(mr.ss.re.ab$sims.list$h.ad), "\n"),
cat("Mean of Adult harvest (WMU):", mean(mr.ss.re.ab$sims.list$h.ad.wmu), "\n"),
cat("Mean of Juvenile survival (State):", mean(mr.ss.re.ab$sims.list$s.jv), "\n"),
cat("Mean of Juvenile survival (WMU):", mean(mr.ss.re.ab$sims.list$s.juv.wmu), "\n"),
cat("Mean of Adult survival (State):", mean(mr.ss.re.ab$sims.list$s.ad), "\n"),
cat("Mean of Adult survival (WMU):", mean(mr.ss.re.ab$sims.list$s.ad.wmu), "\n"),
cat("Mean of Juvenile abundance (State):", mean(mr.ss.re.ab$sims.list$abun.juv), "\n"),
cat("Mean of Juvenile abundance (WMU):", mean(mr.ss.re.ab$sims.list$abun.juv.wmu), "\n"),
cat("Mean of Adult abundance (State):", mean(mr.ss.re.ab$sims.list$abun.ad), "\n"),
cat("Mean of Adult abundance (WMU):", mean(mr.ss.re.ab$sims.list$abun.ad.wmu), "\n"),
cat("Mean of Reruitment (State):", mean(mr.ss.re.ab$sims.list$P), "\n"),
cat("Mean of Reruitment (WMU):", mean(mr.ss.re.ab$sims.list$P.wmu), "\n"),
cat("95% credible interval Juvenile harvest (State):", quantile(mr.ss.re.ab$sims.list$h.juv, c(0.025, 0.975)), "\n"),
cat("95% credible interval Juvenile harvest (WMU):", quantile(mr.ss.re.ab$sims.list$h.juv.wmu, c(0.025, 0.975)), "\n"),
cat("95% credible interval for Adult harvest (State):", quantile(mr.ss.re.ab$sims.list$h.ad, c(0.025, 0.975)), "\n"),
cat("95% credible interval for Adult harvest (WMU):", quantile(mr.ss.re.ab$sims.list$h.ad.wmu, c(0.025, 0.975)), "\n"),
cat("95% credible interval for Juvenile survial (State):", quantile(mr.ss.re.ab$sims.list$s.jv, c(0.025, 0.975)), "\n"),
cat("95% credible interval for Juvenile survival (WMU):", quantile(mr.ss.re.ab$sims.list$s.juv.wmu, c(0.025, 0.975)), "\n"),
cat("95% credible interval for Adult survival (State):", quantile(mr.ss.re.ab$sims.list$s.ad, c(0.025, 0.975)), "\n"),
cat("95% credible interval for Adult survival (WMU):", quantile(mr.ss.re.ab$sims.list$s.ad.wmu, c(0.025, 0.975)), "\n"),
cat("95% credible interval for Juvenile abundance (State):", quantile(mr.ss.re.ab$sims.list$abun.juv, c(0.025, 0.975)), "\n"),
cat("95% credible interval for Juvenile abundance (WMU):", quantile(mr.ss.re.ab$sims.list$abun.juv.wmu, c(0.025, 0.975)), "\n"),
cat("95% credible interval for Adult abundance (State):", quantile(mr.ss.re.ab$sims.list$abun.ad, c(0.025, 0.975)), "\n"),
cat("95% credible interval for Adult abundance (WMU):", quantile(mr.ss.re.ab$sims.list$abun.ad.wmu, c(0.025, 0.975)), "\n"))


# get summary from MCMC
x <- MCMCsummary(mr.ss.re.ab, params = c('abun.ad','abun.juv','s.jv', 's.ad',
                                    'h.ad', 'h.juv'), round = 2)

MCMCtrace(mr.ss.re.ab, params = c('abun.ad','abun.juv','s.jv', 's.ad',
                                    'h.ad', 'h.juv', 'beta', 'gamma'), pdf = T)

# create table for overleaf
sum.table <- data.frame(x) %>% 
  knitr::kable(booktabs = T,
               escape = F,
               caption = "MCMC summary of abundance (abun), survival (s), and harvest (h),
               per age class (ad, juv), and year. Mean and standard deviation are reported
               along with reduction statistic statistic (Rhat) and effective size (n.eff). Credible intervals
               of 2.5 to 97.5 are reported as well.",
               format = "latex",
               align = c("lrccccccccc")) %>% 
  column_spec(2,bold=T,latex_column_spec = ">{\\\\color{black}}c") %>% 
  collapse_rows(columns = 2, latex_hline = "major",valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header", "striped"), full_width = FALSE)

writeClipboard(sum.table)

# Done!
