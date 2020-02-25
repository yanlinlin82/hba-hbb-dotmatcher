library(tidyverse)

window_size <- 3
threshold <- 11

blosum <- read.table("BLOSUM62", check.names = FALSE)
hba <- readLines("HBA.fa")[-1] %>% paste(collapse = "")
hbb <- readLines("HBB.fa")[-1] %>% paste(collapse = "")

sliding_window <- function(x, size) {
  sapply(1:(nchar(x)-(size-1)), function(i) substr(x, i, i+(size-1)))
}
hba_win <- sliding_window(hba, window_size)
hbb_win <- sliding_window(hbb, window_size)

calc_score <- function(a, b) {
  a <- strsplit(a, "")
  b <- strsplit(b, "")
  sapply(seq_along(a), function(i) {
           sum(sapply(seq_along(a[[i]]), function(j) {
                    blosum[a[[i]][j],b[[i]][j]]
                  }))
         })
}
a <- expand.grid(a = seq_along(hba_win),
                 b = seq_along(hbb_win),
                 stringsAsFactors = FALSE) %>%
  as_tibble %>%
  mutate(hba = hba_win[a],
         hbb = hbb_win[b],
         score = calc_score(hba, hbb))

g <- a %>%
  ggplot() +
  geom_abline(color = "blue", size = .5, alpha = .5) +
  geom_rect(aes(xmin = a - 1, xmax = a, ymin = b - 1, ymax = b, fill = (score >= threshold))) +
  scale_x_continuous(limits = c(0,length(hba_win)), breaks = seq(from = 0, to = length(hba_win), by = 10)) +
  scale_y_continuous(limits = c(0,length(hbb_win)), breaks = seq(from = 0, to = length(hbb_win), by = 10)) +
  scale_fill_manual(values = c(NA, "black")) +
  labs(x = "HBA", y = "HBB") +
  guides(fill = FALSE) +
  theme_bw() +
  theme(title = element_text(size = 40),
        text = element_text(size = 30))
g %>% ggsave(filename = "HBA-vs-HBB.png", width = 15, height = 15)

g <- a %>%
  ggplot() +
  geom_abline(color = "blue", size = .5, alpha = .5) +
  geom_rect(aes(xmin = b - 1, xmax = b, ymin = a - 1, ymax = a, fill = (score >= threshold))) +
  scale_x_continuous(limits = c(0,length(hbb_win)), breaks = (0:15)*10) +
  scale_y_continuous(limits = c(0,length(hba_win)), breaks = (0:15)*10) +
  scale_fill_manual(values = c(NA, "black")) +
  labs(x = "HBB", y = "HBA") +
  guides(fill = FALSE) +
  theme_bw() +
  theme(title = element_text(size = 40),
        text = element_text(size = 30))
g %>% ggsave(filename = "HBB-vs-HBA.png", width = 15, height = 15)
