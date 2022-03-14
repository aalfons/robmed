# load packages
library("dplyr")
library("ggplot2")
library("robmed")
library("scales")

# control parameters
tuning.psi <- reg_control()$tuning.psi
col <- hue_pal()(4)[c(1, 3)]
lty <- c(2, 1)
line_size <- 2/3

# grid of x-values
x <- seq(-4, 4, length.out = 100)

# labels to be used for methods and panels
methods <- c("OLS", "MM")
functions <- c("Loss", "Weight")

# loss functions
OLS_loss <- data.frame(x = x,
                       y = x^2,
                       Method = factor(methods[1], levels = methods),
                       Function = factor(functions[1], levels = functions))
robust_loss <- data.frame(x = x,
                          y = Mpsi(x, tuning.psi, "bisquare", deriv = -1),
                          Method = factor(methods[2], levels = methods),
                          Function = factor(functions[1], levels = functions))

# weight functions
OLS_weight <- data.frame(x = x,
                         y = 1,
                         Method = factor(methods[1], levels = methods),
                         Function = factor(functions[2], levels = functions))
robust_weight <- data.frame(x = x,
                            y = Mwgt(x, tuning.psi, "bisquare"),
                            Method = factor(methods[2], levels = methods),
                            Function = factor(functions[2], levels = functions))

# data frame for plot
df <- rbind(OLS_loss, robust_loss, OLS_weight, robust_weight) %>% filter(y <= 4)

# write plot to file
pdf(file = "vignettes/static-figures/loss.pdf", width = 8.5, height = 3.6)
ggplot() +
  geom_line(aes(x = x, y = y, color = Method, linetype = Method),
            data = df, size = line_size) +
  labs(x = NULL, y = NULL) +
  scale_linetype_manual(values = lty) +
  facet_wrap(~ Function, scales = "free_y") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 13),
                     legend.key = element_rect(color = NA),
                     legend.key.width = unit(1.5, "line"),
                     legend.text = element_text(size = 12),
                     legend.title = element_text(size = 13),
                     strip.text = element_text(size = 12))
dev.off()

