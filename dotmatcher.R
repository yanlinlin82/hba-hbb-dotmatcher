library(tidyverse)

dotmatcher <- function(a_seq, b_seq, blosum, window_size, threshold) {
  sliding_window <- function(x, size) {
    sapply(1:(nchar(x)-(size-1)), function(i) substr(x, i, i+(size-1)))
  }
  a_win <- sliding_window(a_seq, window_size)
  b_win <- sliding_window(b_seq, window_size)

  calc_score <- function(a, b) {
    a <- strsplit(a, "")
    b <- strsplit(b, "")
    sapply(seq_along(a), function(i) {
             sum(sapply(seq_along(a[[i]]), function(j) {
                      blosum[a[[i]][j],b[[i]][j]]
                    }))
           })
  }
  a <- expand.grid(a_index = seq_along(a_win),
                   b_index = seq_along(b_win),
                   stringsAsFactors = FALSE) %>%
    as_tibble %>%
    mutate(a_seg = a_win[a_index],
           b_seg = b_win[b_index],
           score = calc_score(a_seg, b_seg))

  a2 <- lapply(window_size:1, function(d) {
                 a %>%
                   filter(score >= threshold) %>%
                   mutate(a_index = a_index + d - 1,
                          b_index = b_index + d - 1,
                          delta = d - 1)
          }) %>% do.call("rbind", .) %>% as_tibble

  g <- a2 %>%
    ggplot() +
    geom_abline(color = "blue", size = .5, alpha = .5) +
    geom_rect(aes(xmin = a_index - 1, xmax = a_index, ymin = b_index - 1, ymax = b_index, fill = (delta > 0))) +
    scale_x_continuous(limits = c(0,nchar(a_seq)), breaks = seq(from = 0, to = nchar(a_seq), by = 10)) +
    scale_y_continuous(limits = c(0,nchar(b_seq)), breaks = seq(from = 0, to = nchar(b_seq), by = 10)) +
    scale_fill_manual(values = c("red", "black")) +
    labs(x = as.character(substitute(a_seq)),
         y = as.character(substitute(b_seq))) +
    guides(fill = FALSE) +
    theme_bw() +
    theme(title = element_text(size = 40),
          text = element_text(size = 30))
  return(g)
}

window_size <- 3
threshold <- 11

blosum <- read.table("BLOSUM62", check.names = FALSE)
hba <- readLines("HBA.fa")[-1] %>% paste(collapse = "")
hbb <- readLines("HBB.fa")[-1] %>% paste(collapse = "")

g <- dotmatcher(hba, hbb, blosum, window_size, threshold)
g %>% ggsave(filename = "HBA-vs-HBB.png", width = 15, height = 15)

g <- dotmatcher(hbb, hba, blosum, window_size, threshold)
g %>% ggsave(filename = "HBB-vs-HBA.png", width = 15, height = 15)
