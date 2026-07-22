load("./results/RData/sampling_complete.RData")
load("./results/RData/sampling_complete_fungi.RData")
load("./results/RData/sampling_complete_plants.RData")

p ##sampling completeness bacteria
p_fung ##sampling completeness fungi
p_pl ##sampling completeness plants

#combine into FIgS1

library(patchwork)


# add panel tags A) and B)
combined_sampl_compl <- p + p_fung + p_pl + plot_annotation(tag_levels = 'A')


ggsave(
  "results/combined_sampl_compl.png",
  combined_sampl_compl,
  width = 15,
  height = 6,
  dpi = 300
)


save(combined_sampl_compl, file = "results/RData/combined_sampl_compl.RData")
