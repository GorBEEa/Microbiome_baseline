## =========================================================
##  PACKAGES
## =========================================================
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(scales)
library(ggrepel)
library(patchwork)

## =========================================================
##  COMMON SETUP
## =========================================================

# ensure vertex groups are set on g_all
if (is.null(V(g_all)$group)) {
  V(g_all)$group <- ifelse(
    grepl("^Plants__", V(g_all)$name),
    "Plants",
    ifelse(
      grepl("^Bacteria__", V(g_all)$name),
      "Bacteria",
      ifelse(grepl("^Fungi__", V(g_all)$name), "Fungi", "Other")
    )
  )
}

# palettes
domain_cols <- c(
  Plants = "#02401B",
  Bacteria = "#D8B70A",
  Fungi = "#A2A475",
  Other = "grey80"
)

edge_cols <- c(positive = "#81A88D", negative = "#972D15")

# ---- helpers for taxonomy labels ----
.pick_lowest <- function(df, ranks = c("Genus", "Family", "Order")) {
  tmp <- df[, intersect(ranks, names(df)), drop = FALSE]
  tmp[] <- lapply(tmp, function(x) {
    x[is.na(x) | x == ""] <- NA
    x
  })
  lowest_name <- dplyr::coalesce(!!!tmp)
  rr <- as.matrix(tmp)
  lowest_rank <- apply(rr, 1, function(v) {
    idx <- which(!is.na(v))
    if (length(idx) == 0) NA_character_ else colnames(rr)[idx[1]]
  })
  tibble(tax_label = lowest_name, tax_rank = lowest_rank)
}

annotate_domain <- function(
  df_dom,
  tax,
  key_tax,
  ranks = c("Genus", "Family", "Order")
) {
  if (nrow(df_dom) == 0) {
    return(
      df_dom %>% mutate(tax_label = NA_character_, tax_rank = NA_character_)
    )
  }
  tax_small <- tax %>%
    distinct(.data[[key_tax]], .keep_all = TRUE) %>%
    select(all_of(c(key_tax, ranks[ranks %in% names(tax)])))
  out <- df_dom %>%
    left_join(tax_small, by = c("taxon" = key_tax))
  picked <- .pick_lowest(out, ranks = ranks)
  out$tax_label <- picked$tax_label
  out$tax_rank <- picked$tax_rank
  out
}

## =========================================================
##  1. FULL TRIPARTITE GRAPH VERSION (NO LABELS)
## =========================================================

min_r_plot <- 0.15 # show only edges with |r| >= this

g_full_tbl <- as_tbl_graph(g_all) %>%
  activate(edges) %>%
  filter(abs(r) >= min_r_plot) %>%
  activate(nodes) %>%
  mutate(
    group = ifelse(
      grepl("^Plants__", name),
      "Plants",
      ifelse(
        grepl("^Bacteria__", name),
        "Bacteria",
        ifelse(grepl("^Fungi__", name), "Fungi", "Other")
      )
    ),
    degree = centrality_degree()
  )

set.seed(123)
p_full <- ggraph(g_full_tbl, layout = "fr") +
  geom_edge_link(
    aes(width = abs(r), alpha = abs(r), colour = sign),
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(size = degree, fill = group),
    shape = 21,
    colour = "black"
  ) +
  scale_fill_manual(values = domain_cols, name = "Domain") +
  scale_edge_colour_manual(values = edge_cols) +
  scale_edge_width(range = c(0.2, 1.8)) +
  scale_edge_alpha(range = c(0.1, 0.7)) +
  scale_size_continuous(range = c(2, 10), guide = "none") +
  theme_graph() +
  ggtitle("Tripartite cross-domain co-occurrence network")

## =========================================================
##  2. ZOOM: STRONGEST POSITIVE CO-OCCURRENCES (P–B–F)
## =========================================================

top_edges <- 40 # number of strongest positive edges to show

edges_all <- as_data_frame(g_all, "edges")
verts_all <- as_data_frame(g_all, "vertices")

# take strongest positive edges overall
edges_pos_top <- edges_all %>%
  filter(sign == "positive") %>%
  arrange(desc(r)) %>%
  slice_head(n = top_edges)

verts_zoom <- verts_all %>%
  filter(name %in% c(edges_pos_top$from, edges_pos_top$to))

g_zoom <- graph_from_data_frame(
  edges_pos_top,
  directed = FALSE,
  vertices = verts_zoom
)

# ensure group on zoom graph
if (is.null(V(g_zoom)$group)) {
  V(g_zoom)$group <- ifelse(
    grepl("^Plants__", V(g_zoom)$name),
    "Plants",
    ifelse(
      grepl("^Bacteria__", V(g_zoom)$name),
      "Bacteria",
      ifelse(grepl("^Fungi__", V(g_zoom)$name), "Fungi", "Other")
    )
  )
}

# vertex table for annotation
verts_zoom_df <- as_data_frame(g_zoom, "vertices") %>%
  mutate(
    domain = case_when(
      grepl("^Plants__", name) ~ "Plants",
      grepl("^Bacteria__", name) ~ "Bacteria",
      grepl("^Fungi__", name) ~ "Fungi",
      TRUE ~ "Other"
    ),
    taxon = sub("^[A-Za-z]+__", "", name)
  )

# annotate each domain with taxonomy
df_pl <- verts_zoom_df %>% filter(domain == "Plants")
df_bac <- verts_zoom_df %>% filter(domain == "Bacteria")
df_fun <- verts_zoom_df %>% filter(domain == "Fungi")

df_pl_a <- annotate_domain(df_pl, tax_pln, "ASV_code")
df_bac_a <- annotate_domain(df_bac, tax_bac, "ASV_code")
df_fun_a <- annotate_domain(df_fun, tax_fun, "ASV_code")

labels_zoom <- bind_rows(df_pl_a, df_bac_a, df_fun_a) %>%
  select(name, domain, tax_label, tax_rank)

# tidygraph + merge labels
g_zoom_tbl <- as_tbl_graph(g_zoom) %>%
  activate(nodes) %>%
  left_join(labels_zoom, by = "name") %>%
  mutate(
    group = ifelse(
      grepl("^Plants__", name),
      "Plants",
      ifelse(
        grepl("^Bacteria__", name),
        "Bacteria",
        ifelse(grepl("^Fungi__", name), "Fungi", "Other")
      )
    )
  )

# compute node strength (sum of |r|)
strength_vec <- igraph::strength(g_zoom, weights = abs(E(g_zoom)$r))

g_zoom_tbl <- g_zoom_tbl %>%
  mutate(
    strength = strength_vec,
    label_plot = ifelse(is.na(tax_label), name, tax_label)
  )

# build layout explicitly so x,y are available for ggrepel
set.seed(123)
lay_zoom <- create_layout(g_zoom_tbl, layout = "fr")

p_zoom <- ggplot(lay_zoom) +
  geom_edge_link(aes(
    x = x,
    y = y,
    xend = xend,
    yend = yend,
    width = abs(r),
    alpha = abs(r),
    colour = sign
  )) +
  geom_point(
    aes(x = x, y = y, size = strength, fill = group),
    shape = 21,
    colour = "black"
  ) +
  geom_text_repel(
    aes(x = x, y = y, label = label_plot),
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.size = 0.2,
    max.overlaps = Inf
  ) +
  scale_fill_manual(values = domain_cols, name = "Domain") +
  scale_edge_colour_manual(values = edge_cols) +
  scale_edge_width(range = c(0.4, 2.2), guide = "none") +
  scale_edge_alpha(range = c(0.4, 0.9), guide = "none") +
  scale_size_continuous(range = c(4, 11), guide = "none") +
  theme_graph() +
  ggtitle("Zoom: strongest positive co-occurrences\n(Plants–Bacteria–Fungi)")

## =========================================================
##  3. COMBINED FIGURE
## =========================================================

combined <- p_full + p_zoom + plot_layout(widths = c(1.1, 1))

combined

ggsave(
  "results/tripartite_full_plus_zoom_positive.png",
  combined,
  width = 14,
  height = 6,
  dpi = 400
)
