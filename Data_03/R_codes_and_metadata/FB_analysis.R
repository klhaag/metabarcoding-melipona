library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)

crop.data <- read.csv("supplementation_bees.csv", header = TRUE, colClasses = c("factor", "factor", "factor", "factor", "numeric"))
normality <- shapiro.test(crop.data$FB)
one.way1 <- aov(FB ~ Pertreatment, data = crop.data)
summary(one.way1)
one.way2 <- aov(FB ~ Period, data = crop.data)
summary(one.way2)
one.way3 <- aov(FB ~ Colony, data = crop.data)
summary(one.way3)
two.way <- aov(FB ~ Period + Pertreatment, data = crop.data)
summary(two.way)
interaction <- aov(FB ~ Period*Pertreatment, data = crop.data)
summary(interaction)

crop.data %>%
  select(Colony, Treatment, Pertreatment, Period, FB) %>%
  pivot_longer(cols=c(Colony, Treatment, Pertreatment, Period), names_to="characteristic", values_to="value") %>%
  drop_na() %>%
  nest(data = -characteristic) %>%
  mutate(tests = map(data, ~tidy(kruskal.test(FB ~ value, data=.x)))) %>%
  unnest(cols=tests) %>%
  select(-data) %>%
  mutate(p.value.adj = p.adjust(p.value, method="BH"))

pairwise.wilcox.test(g=crop.data$Pertreatment, x=crop.data$FB, p.adjust.method="BH")

ggplot(crop.data, aes(x=Treatment, y=FB, color=Period, fill=Period)) +
  geom_boxplot(alpha=0.9, outlier.colour=NA) +
  geom_jitter(shape=19)+
  scale_color_manual(name=NULL,
                     values=c("darkblue","darkred"),
                     breaks=c("A","D"),
                     labels=c("Before","After")) +
  scale_fill_manual(name=NULL,
                    values=c("lightblue","pink"),
                    breaks=c("A","D"),
                    labels=c("Before","After")) +
  scale_x_discrete(limits=c("AS","S","C"),
                   labels=c("AP+SA","SA","CT")) +
  labs(title=NULL,
       x=NULL,
       y="Fat Body Index") +
  theme_classic()