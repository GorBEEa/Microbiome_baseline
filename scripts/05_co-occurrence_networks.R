
load("./results/RData/data_co-occurrence.RData")


####
bac<-df_asv_t
fun<-df_fung_t
pln<-df_asv_pl_t

tax_bac <- tax_table.cl    

####REMOVE CHLOROPLASTS FROM BACTERIA DATA

# ---- 1) Identify chloroplast ASVs in the bacteria taxonomy ----
# Choose the ASV id column name in your taxonomy table:
key_bac <- if ("ASV_code" %in% names(tax_bac)) "ASV_code" else
  if ("ASV"      %in% names(tax_bac)) "ASV"      else
    stop("Could not find ASV id column in tax_bac (tried 'ASV_code' and 'ASV').")

# Taxonomy columns to scan for the word "Chloroplast"
tax_cols <- intersect(c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Lineage","Taxon"),
                      names(tax_bac))


tax_bac[[key_bac]] <- trimws(as.character(tax_bac[[key_bac]]))
colnames(bac)      <- trimws(colnames(bac))


# TRUE if any rank mentions "chloroplast" (case-insensitive)
is_chloro <- apply(tax_bac[, tax_cols, drop = FALSE], 1, function(x)
  any(grepl("chloroplast", x, ignore.case = TRUE))
)

chloro_asvs <- trimws(unique(tax_bac[[key_bac]][is_chloro]))



# ---- 2) Remove those ASVs from the bacteria count matrix ----
drop_cols <- intersect(colnames(bac), chloro_asvs)
bac <- bac[, !(colnames(bac) %in% drop_cols), drop = FALSE]


# ---- 3) (Recommended) Keep taxonomy in sync (remove chloroplast rows) ----
tax_bac <- tax_bac[!is_chloro, , drop = FALSE]

# ---- 4) Quick summary ----
message(sprintf("Removed %d chloroplast ASVs from bacteria (%.1f%% of features).",
                length(drop_cols),
                ifelse(ncol(df_asv_t) > 0, 100*length(drop_cols)/ncol(df_asv_t), 0)))

stopifnot(!any(colnames(bac) %in% chloro_asvs))
message(sprintf("Removed %d chloroplast ASVs from bacteria.", length(drop_cols)))



# 0) Normalize sample IDs (keep only the part before the first "_")
normalize_sample_ids <- function(x) sub("(_.*)$", "", x)

rownames(bac) <- normalize_sample_ids(rownames(bac))
rownames(fun) <- normalize_sample_ids(rownames(fun))
rownames(pln) <- normalize_sample_ids(rownames(pln))

# 1) Collapse replicates that now share the same ID (sum counts)
collapse_by_id <- function(M) {
  M <- as.matrix(M); storage.mode(M) <- "numeric"
  rowsum(M, group = rownames(M), reorder = FALSE)
}
bac <- collapse_by_id(bac)
fun <- collapse_by_id(fun)
pln <- collapse_by_id(pln)



# Packages
library(igraph)

# ---------- helpers ----------
subset_samples_min_reads <- function(X, min_reads = 300) {
  X[rowSums(X) >= min_reads, , drop = FALSE]
}

select_top_taxa <- function(X, top_n = 100) {
  if (ncol(X) <= top_n) return(X)
  keep <- order(colSums(X), decreasing = TRUE)[seq_len(top_n)]
  X[, keep, drop = FALSE]
}

clr_transform <- function(X, pseudo = 1) {
  # X: samples x taxa
  Y <- X + pseudo
  gm <- exp(rowMeans(log(Y)))     # geometric mean per sample
  sweep(log(Y), 1, log(gm), FUN = "-")
}

# Build cross-domain network (pairwise)
build_cross_network5 <- function(X, Y, group1 = "G1", group2 = "G2",
                                 top_n = 100, min_reads_X = 300, min_reads_Y = 300,
                                 min_r = 0.1, pseudo = 1) {
  
  stopifnot(!is.null(rownames(X)), !is.null(rownames(Y)))
  
  # 1) Align samples
  common <- intersect(rownames(X), rownames(Y))
  if (length(common) == 0) stop("No shared samples between datasets.")
  X <- as.matrix(X[common, , drop = FALSE]); storage.mode(X) <- "numeric"
  Y <- as.matrix(Y[common, , drop = FALSE]); storage.mode(Y) <- "numeric"
  
  # 2) Per-dataset read thresholds
  keep <- (rowSums(X) >= min_reads_X) & (rowSums(Y) >= min_reads_Y)
  if (!any(keep)) {
    stop(sprintf("No common samples pass thresholds: X(>= %d)=%d, Y(>= %d)=%d.",
                 min_reads_X, sum(rowSums(X) >= min_reads_X),
                 min_reads_Y, sum(rowSums(Y) >= min_reads_Y)))
  }
  X <- X[keep, , drop = FALSE]
  Y <- Y[keep, , drop = FALSE]
  
  # 3) Select top features (by overall abundance)
  select_top <- function(M, k) {
    if (ncol(M) <= k) return(M)
    M[, order(colSums(M), decreasing = TRUE)[seq_len(k)], drop = FALSE]
  }
  X <- select_top(X, top_n)
  Y <- select_top(Y, top_n)
  if (ncol(X) + ncol(Y) < 2) stop("Not enough features after filtering.")
  
  # 4) Ensure unique, non-NA colnames across domains
  fix_names <- function(nm, prefix) {
    if (is.null(nm)) nm <- paste0("feat", seq_len(ncol(X)))
    nm <- paste0(prefix, "__", nm)
    make.unique(nm, sep = "_")
  }
  colnames(X) <- fix_names(colnames(X), group1)
  colnames(Y) <- fix_names(colnames(Y), group2)
  
  # 5) CLR on combined matrix
  Z <- cbind(X, Y) + pseudo
  gm <- exp(rowMeans(log(Z)))
  Z_clr <- sweep(log(Z), 1, log(gm), "-")
  
  # Drop zero-variance features
  sds <- apply(Z_clr, 2, sd)
  Z_clr <- Z_clr[, sds > 0, drop = FALSE]
  if (ncol(Z_clr) < 2) stop("All features have zero variance after CLR/filtering.")
  
  # 6) Pearson correlations
  R <- cor(Z_clr, method = "pearson", use = "pairwise.complete.obs")
  cn <- colnames(R)
  
  # 7) Edge list from upper triangle
  idx <- which(abs(R) > min_r & upper.tri(R), arr.ind = TRUE)
  if (nrow(idx) == 0) {
    # Return empty edge list but keep vertices
    groups <- ifelse(startsWith(cn, paste0(group1, "__")), group1, group2)
    verts  <- data.frame(name = cn, group = groups, stringsAsFactors = FALSE)
    g <- igraph::graph_from_data_frame(data.frame(from=character(), to=character()),
                                       directed = FALSE, vertices = verts)
    return(list(graph = g,
                edges = data.frame(from=character(), to=character(), r=numeric(),
                                   sign=character(), from_group=character(), to_group=character()),
                R = R, samples_used = rownames(Z_clr)))
  }
  
  edges <- data.frame(
    from = cn[idx[,1]],
    to   = cn[idx[,2]],
    r    = R[idx],
    stringsAsFactors = FALSE
  )
  
  # keep only cross-domain edges
  groups <- ifelse(startsWith(cn, paste0(group1, "__")), group1, group2)
  names(groups) <- cn
  edges$from_group <- groups[edges$from]
  edges$to_group   <- groups[edges$to]
  edges <- edges[edges$from_group != edges$to_group, , drop = FALSE]
  
  if (nrow(edges) == 0) {
    verts  <- data.frame(name = cn, group = groups, stringsAsFactors = FALSE)
    g <- igraph::graph_from_data_frame(data.frame(from=character(), to=character()),
                                       directed = FALSE, vertices = verts)
    return(list(graph = g, edges = edges, R = R, samples_used = rownames(Z_clr)))
  }
  
  edges$sign <- ifelse(edges$r > 0, "positive", "negative")
  
  # 8) Graph from edge list (robust to dimnames issues)
  verts <- data.frame(name = cn, group = groups, stringsAsFactors = FALSE)
  g <- igraph::graph_from_data_frame(edges[, c("from","to")], directed = FALSE, vertices = verts)
  igraph::E(g)$r    <- edges$r
  igraph::E(g)$sign <- edges$sign
  
  list(graph = g, edges = edges, R = R, samples_used = rownames(Z_clr))
}


# Bacteria–Fungi (unchanged thresholds)
bf <- build_cross_network5(bac, fun,
                           group1 = "Bacteria", group2 = "Fungi",
                           top_n = 100, min_reads_X = 300, min_reads_Y = 300, min_r = 0.1)

# Plants–Bacteria (lower plant threshold if needed)
pb <- build_cross_network5(pln, bac,
                           group1 = "Plants", group2 = "Bacteria",
                           top_n = 100, min_reads_X = 100, min_reads_Y = 300, min_r = 0.1)

# Plants–Fungi
pf <- build_cross_network5(pln, fun,
                           group1 = "Plants", group2 = "Fungi",
                           top_n = 100, min_reads_X = 100, min_reads_Y = 300, min_r = 0.1)



library(dplyr)
library(igraph)

edges_all <- dplyr::bind_rows(
  bf$edges %>% mutate(pair = "Bacteria–Fungi"),
  pb$edges %>% mutate(pair = "Plants–Bacteria"),
  pf$edges %>% mutate(pair = "Plants–Fungi")
) %>% dplyr::distinct(from, to, .keep_all = TRUE)

get_group <- function(nm) {
  ifelse(grepl("^Plants__", nm), "Plants",
         ifelse(grepl("^Bacteria__", nm), "Bacteria",
                ifelse(grepl("^Fungi__", nm), "Fungi", "Other")))
}
verts <- data.frame(name = unique(c(edges_all$from, edges_all$to)))
verts$group <- get_group(verts$name)

g_all <- igraph::graph_from_data_frame(edges_all[, c("from","to")],
                                       directed = FALSE, vertices = verts)
E(g_all)$r    <- edges_all$r
E(g_all)$sign <- edges_all$sign






# 4) Quick sanity checks
length(intersect(rownames(pln), rownames(bac)))  # shared samples plants-bacteria
length(intersect(rownames(pln), rownames(fun)))  # shared samples plants-fungi


# Quick plots (basic)
plot(bf$graph, vertex.size = 5, vertex.label = NA)
plot(pb$graph, vertex.size = 5, vertex.label = NA)
plot(pf$graph, vertex.size = 5, vertex.label = NA)


# Save edge tables
# write.csv(bf$edges, "edges_bacteria_fungi.csv", row.names = FALSE)
# write.csv(pb$edges, "edges_plants_bacteria.csv", row.names = FALSE)
# write.csv(pf$edges, "edges_plants_fungi.csv", row.names = FALSE)




# deps
# install.packages(c("igraph","ggraph","tidygraph","dplyr","scales"))
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)

# ========= Paletas y helper =========
# ============ Paletas y helper ============
# Formas "rellenables": 21=círculo, 22=cuadrado, 24=triángulo
.node_shapes <- c(Bacteria = 21, Fungi = 22, Plants = 24)
.node_fills  <- c(Bacteria = "#D8B70A", Fungi = "#A2A475", Plants = "#02401B")
.edge_cols   <- c(positive = "#81A88D", negative = "#972D15")

prep_graph <- function(g){
  if (is.null(igraph::V(g)$group)) {
    igraph::V(g)$group <- ifelse(grepl("^Plants__",   igraph::V(g)$name), "Plants",
                                 ifelse(grepl("^Bacteria__", igraph::V(g)$name), "Bacteria",
                                        ifelse(grepl("^Fungi__",    igraph::V(g)$name), "Fungi", "Other")))
  }
  igraph::V(g)$deg <- igraph::degree(g)
  g
}
# ============ TRIPARTITO ============
plot_tripartite_v4 <- function(g, title = "Cross-domain co-occurrence") {
  g <- prep_graph(g)
  
  # x por dominio; y por rango de grado dentro del dominio
  xmap <- c("Plants" = 0, "Bacteria" = 1, "Fungi" = 2)
  vd <- data.frame(name = igraph::V(g)$name,
                   group = igraph::V(g)$group,
                   deg   = igraph::V(g)$deg)
  vd <- dplyr::group_by(vd, group) |>
    dplyr::arrange(dplyr::desc(deg), name, .by_group = TRUE) |>
    dplyr::mutate(rank = dplyr::row_number(),
                  y = scales::rescale(-(rank), to = c(-1, 1))) |>
    dplyr::ungroup()
  
  ord <- match(igraph::V(g)$name, vd$name)
  x <- xmap[vd$group[ord]]
  y <- vd$y[ord]
  
  p <- ggraph::ggraph(g, layout = "manual", x = x, y = y)
  
  # Aristas (bezier si disponible, si no rectas)
  if ("geom_edge_bezier" %in% getNamespaceExports("ggraph")) {
    p <- p + ggraph::geom_edge_bezier(
      ggplot2::aes(color = sign,
                   width = pmax(0.2, abs(r)),
                   alpha = pmin(1, scales::rescale(abs(r)))),
      curvature = 0.2, show.legend = TRUE
    )
  } else {
    p <- p + ggraph::geom_edge_link(
      ggplot2::aes(color = sign,
                   width = pmax(0.2, abs(r)),
                   alpha = pmin(1, scales::rescale(abs(r)))),
      show.legend = TRUE
    )
  }
  
  # Nodos: forma+relleno por dominio (contorno negro fino)
  p <- p + ggraph::geom_node_point(
    ggplot2::aes(shape = group, fill = group, size = pmax(1, deg)),
    color = "black", stroke = 0.3, show.legend = TRUE
  )
  
  # Escalas de arista
  if ("scale_edge_width_continuous" %in% getNamespaceExports("ggraph")) {
    p <- p + ggraph::scale_edge_width_continuous(
      name = "Pearson coefficient of correlation (|r|)", range = c(0.2, 1.6))
  } else {
    p <- p + ggraph::scale_edge_width(name = "Pearson coefficient of correlation (|r|)",
                                      range = c(0.2, 1.6))
  }
  if ("scale_edge_alpha_continuous" %in% getNamespaceExports("ggraph")) {
    p <- p + ggraph::scale_edge_alpha_continuous(guide = "none",
                                                 range = c(0.25, 0.9))
  } else {
    p <- p + ggraph::scale_edge_alpha(guide = "none", range = c(0.25, 0.9))
  }
  p <- p + ggraph::scale_edge_colour_manual(values = .edge_cols, guide = "none")
  
  # Leyenda de dominio con formas visibles (usamos la leyenda de 'fill' y la sobreescribimos con shapes)
  domain_order <- c("Plants","Bacteria","Fungi")
  p +
    ggplot2::scale_shape_manual(values = .node_shapes, guide = "none") +
    ggplot2::scale_fill_manual(
      values  = .node_fills,
      breaks  = domain_order,
      name    = "domain",
      guide   = ggplot2::guide_legend(
        override.aes = list(
          shape  = unname(.node_shapes[domain_order]),
          fill   = unname(.node_fills[domain_order]),
          colour = "black",
          size   = 4
        )
      )
    ) +
    ggplot2::scale_size_continuous(guide = "none") +
    ggraph::theme_graph() +
    ggplot2::labs(title = title)
}
# ============ BIPARTITO ============
plot_pair_network_v2 <- function(g, left_group, right_group, title = NULL) {
  g <- prep_graph(g)
  g <- igraph::induced_subgraph(g, vids = igraph::V(g)[igraph::V(g)$group %in% c(left_group, right_group)])
  igraph::V(g)$type <- igraph::V(g)$group == right_group
  
  L <- igraph::layout_as_bipartite(g, types = igraph::V(g)$type)
  x <- L[,1]; y <- L[,2]
  
  p <- ggraph::ggraph(g, layout = "manual", x = x, y = y) +
    ggraph::geom_edge_link(
      ggplot2::aes(width = pmax(0.2, abs(r)),
                   alpha = pmin(1, scales::rescale(abs(r))),
                   color = sign),
      show.legend = TRUE
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(shape = group, fill = group, size = pmax(1, deg)),
      color = "black", stroke = 0.3, show.legend = TRUE
    ) +
    ggraph::scale_edge_width_continuous(name = "Pearson coefficient of correlation (|r|)",
                                        range = c(0.2, 1.5)) +
    ggraph::scale_edge_alpha_continuous(guide = "none", range = c(0.3, 0.9)) +
    ggraph::scale_edge_colour_manual(values = .edge_cols, guide = "none")
  
  # Leyenda de dominio con formas: solo los dos presentes, en orden fijo
  domain_order <- c("Plants","Bacteria","Fungi")
  present_groups <- domain_order[domain_order %in% unique(igraph::V(g)$group)]
  
  p +
    ggplot2::scale_shape_manual(values = .node_shapes, guide = "none") +
    ggplot2::scale_fill_manual(
      values = .node_fills,
      breaks = present_groups,
      name   = "domain",
      guide  = ggplot2::guide_legend(
        override.aes = list(
          shape  = unname(.node_shapes[present_groups]),
          fill   = unname(.node_fills[present_groups]),
          colour = "black",
          size   = 4
        )
      )
    ) +
    ggplot2::scale_size_continuous(guide = "none") +
    ggraph::theme_graph() +
    ggplot2::labs(title = if (is.null(title)) paste(left_group, "–", right_group) else title)
}


p_tri<-plot_tripartite_v4(g_all, title = "Plants | Bacteria | Fungi cross-domain co-occurrence")
p_bf <- plot_pair_network_v2(g_all, "Bacteria", "Fungi",    "Bacteria–Fungi")
p_pb <- plot_pair_network_v2(g_all, "Plants",   "Bacteria", "Plants–Bacteria")
p_pf <- plot_pair_network_v2(g_all, "Plants",   "Fungi",    "Plants–Fungi")

p_tri
p_bf 
p_pb 
p_pf



ggsave(
  "results/tripartite_coex_plot.png",
  plot = p_tri,
  width = 14,
  height = 5,
  dpi = 400
)



ggsave(
  "results/bipartite_coex_plot_bac_fun.png",
  plot = p_bf,
  width = 9,
  height = 5,
  dpi = 400
)



ggsave(
  "results/bipartite_coex_plot_plant_bac.png",
  plot = p_pb,
  width = 9,
  height = 5,
  dpi = 400
)



ggsave(
  "results/bipartite_coex_plot_plant_fun.png",
  plot = p_pf,
  width = 9,
  height = 5,
  dpi = 400
)


save(
  p_tri, p_bf, 
  p_pb, 
  p_pf,
  file = "./results/RData/tripartite_bipartite_coex_plots.RData"
)


# print to device (e.g., grid.arrange) or save
# ggsave("tripartite_crossdomain.png", width = 8, height = 5, dpi = 300)
# ggsave("panel_bacteria_fungi.png",   p_bf, width = 6, height = 4, dpi = 300)
# ggsave("panel_plants_bacteria.png",  p_pb, width = 6, height = 4, dpi = 300)
# ggsave("panel_plants_fungi.png",     p_pf, width = 6, height = 4, dpi = 300)



library(dplyr)
# Edges/vertices from the graph you plotted
edges <- igraph::as_data_frame(g_all, what = "edges")
verts <- igraph::as_data_frame(g_all, what = "vertices")

# Make sure groups exist on vertices (Plants/Bacteria/Fungi)
if (is.null(verts$group)) {
  verts$group <- ifelse(grepl("^Plants__",   verts$name), "Plants",
                        ifelse(grepl("^Bacteria__", verts$name), "Bacteria",
                               ifelse(grepl("^Fungi__",    verts$name), "Fungi", "Other")))
}

# Attach groups to edges and build pair label
edges <- edges %>%
  mutate(from_group = verts$group[match(from, verts$name)],
         to_group   = verts$group[match(to,   verts$name)]) %>%
  filter(from_group != to_group) %>%
  mutate(pair = ifelse(from_group < to_group,
                       paste(from_group, to_group, sep = "–"),
                       paste(to_group, from_group, sep = "–")))

# 1) Nodes per domain
nodes_tbl <- verts %>% count(group, name = "n_nodes")

# 2) How many possible cross-domain pairs (for densities)
nP <- sum(verts$group == "Plants")
nB <- sum(verts$group == "Bacteria")
nF <- sum(verts$group == "Fungi")
possible_tbl <- tibble(
  pair = c("Bacteria–Fungi","Plants–Bacteria","Plants–Fungi"),
  n_possible = c(nB*nF, nP*nB, nP*nF)
)

# 3) Edge summary by pair
pair_summary <- edges %>%
  group_by(pair) %>%
  summarise(
    n_edges   = n(),
    pos       = sum(sign == "positive", na.rm = TRUE),
    neg       = sum(sign == "negative", na.rm = TRUE),
    pct_pos   = 100 * pos / n_edges,
    median_r  = median(r, na.rm = TRUE),
    mean_absr = mean(abs(r), na.rm = TRUE),
    max_absr  = max(abs(r), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(possible_tbl, by = "pair") %>%
  mutate(density = n_edges / n_possible)

# 4) Overall stats
overall <- summarise(edges,
                     total_edges = n(),
                     pct_positive_overall = 100 * mean(sign == "positive"),
                     median_r_overall     = median(r, na.rm = TRUE),
                     mean_absr_overall    = mean(abs(r), na.rm = TRUE)
)

# 5) Hubs (top-5 by degree and by strength |r|)
verts$degree  <- igraph::degree(g_all)
verts$strength <- igraph::strength(g_all, weights = abs(igraph::E(g_all)$r))
top_degree <- verts %>% group_by(group) %>% arrange(desc(degree))  %>% slice_head(n = 5)
top_strength <- verts %>% group_by(group) %>% arrange(desc(strength)) %>% slice_head(n = 5)

# 6) Optional: community structure (weighted)
cl <- igraph::cluster_louvain(g_all, weights = abs(igraph::E(g_all)$r))
mod_val <- igraph::modularity(cl)
comm_sizes <- sizes(cl)


data.frame(
  pair = c("Bacteria–Fungi","Plants–Bacteria","Plants–Fungi"),
  n_samples = c(length(unique(bf$samples_used)),
                length(unique(pb$samples_used)),
                length(unique(pf$samples_used)))
)




    library(igraph)
  
  # Make sure vertex groups exist (from your prefixed names)
  .ensure_groups <- function(g){
    if (is.null(V(g)$group)) {
      V(g)$group <- ifelse(grepl("^Plants__",   V(g)$name), "Plants",
                           ifelse(grepl("^Bacteria__",V(g)$name), "Bacteria",
                                  ifelse(grepl("^Fungi__",   V(g)$name), "Fungi", "Other")))
    }
    g
  }
  
  # Build sentence + Table Sy
  hub_summary <- function(g, per_domain = 3) {
    g <- .ensure_groups(g)
    
    # Degree and weighted strength (use |r| if present)
    w <- if ("r" %in% edge_attr_names(g)) abs(E(g)$r) else rep(1, igraph::gsize(g))
    V(g)$degree   <- igraph::degree(g, mode = "all")
    V(g)$strength <- igraph::strength(g, weights = w)
    
    strip_prefix <- function(nm) sub("^[A-Za-z]+__", "", nm)
    doms <- intersect(c("Plants","Bacteria","Fungi"), unique(V(g)$group))
    
    # Top by degree & by strength per domain
    picks <- lapply(doms, function(d){
      vids <- V(g)[group == d]
      k <- min(per_domain, length(vids))
      ord_deg <- order(V(g)$degree[vids],   decreasing = TRUE)
      ord_str <- order(V(g)$strength[vids], decreasing = TRUE)
      list(
        top_degree_names   = strip_prefix(V(g)$name[vids][ord_deg][seq_len(k)]),
        top_strength_names = strip_prefix(V(g)$name[vids][ord_str][seq_len(k)])
      )
    })
    names(picks) <- doms
    
    # Sentence for the manuscript
    parts <- mapply(function(d, x) paste0(d, ": ", paste(x$top_degree_names, collapse = ", ")),
                    doms, picks, USE.NAMES = FALSE)
    sentence <- paste0(
      "Hub taxa (top-degree) included ",
      paste(parts, collapse = "; "),
      "; results were consistent when ranking by weighted strength (Table Sy)."
    )
    
    # Table Sy: top 5 by weighted strength per domain
    table_sy <- do.call(rbind, lapply(doms, function(d){
      vids <- V(g)[group == d]
      k <- min(5, length(vids))
      ord <- order(V(g)$strength[vids], decreasing = TRUE)[seq_len(k)]
      data.frame(
        domain   = d,
        taxon    = strip_prefix(V(g)$name[vids][ord]),
        degree   = V(g)$degree[vids][ord],
        strength = V(g)$strength[vids][ord],
        stringsAsFactors = FALSE
      )
    }))
    
    list(sentence = sentence, table_sy = table_sy)
  }
  
  # ---- Use it on your graph g_all ----
  res <- hub_summary(g_all, per_domain = 3)
  cat(res$sentence, "\n\n")
  print(res$table_sy)
  
  
  df.hub<-res$table_sy
  
  ##code to join this table with actual names of those ASVs per group
  
  library(dplyr)
  
  # --- helper: pick first non-missing among Genus, Family, Order, and record the rank ---
  .pick_lowest <- function(df, ranks = c("Genus","Family","Order")) {
    # treat "" as NA
    tmp <- df[, intersect(ranks, names(df)), drop = FALSE]
    tmp[] <- lapply(tmp, function(x) { x[is.na(x) | x == ""] <- NA; x })
    
    # chosen name
    lowest_name <- do.call(coalesce, tmp)
    
    # chosen rank (row-wise first non-NA)
    rr <- matrix(unlist(tmp, use.names = FALSE), nrow = nrow(tmp))
    colnames(rr) <- colnames(tmp)
    lowest_rank <- apply(rr, 1, function(v) {
      idx <- which(!is.na(v))
      if (length(idx) == 0) NA_character_ else colnames(rr)[idx[1]]
    })
    
    tibble(tax_label = lowest_name, tax_rank = lowest_rank)
  }
  
  # --- INPUTS you should provide (one taxonomy table per domain) ---
  # BACTERIA taxonomy: columns at least ASV_code, Genus, Family, Order
  # FUNGI    taxonomy: columns at least ASV_code, Genus, Family, Order
  # PLANTS   taxonomy: column with plant feature code (e.g., 'Feature' or 'ASV_code' if applicable),
  #                    plus Genus/Family/Order if you have them (or whatever lowest ranks you track)
  
  # Example object names (rename to yours):
     # tax_bac <- tax_table.cl        # your bacteria taxonomy
     tax_fun <- tax_table.cl_fung        # your fungi taxonomy
     
     ##plant taxonomy
     
     ps.gbp23_pl <- readRDS("data/gbp23.plant.decontam.0.5.RDS")
     count_tab.cl_pl <- otu_table(ps.gbp23_pl)
     count_tab.cl_pl <- as.data.frame(count_tab.cl_pl)
     sample_data_tab.cl_pl <- as(sample_data(ps.gbp23_pl), "data.frame")
     tax_pln <- as.data.frame(tax_table(ps.gbp23_pl))
     tax_pln$ASV_code<-row.names(tax_pln)
     
      # And their join keys:
  key_bac <- "ASV_code"   # in tax_bac
  key_fun <- "ASV_code"   # in tax_fun
  key_pln <- "ASV_code"    # e.g., the column that holds 'feat1', 'feat47', ... in tax_pln
  
  # --- annotate each domain separately, then bind back together ---
  annotate_domain <- function(df_dom, tax, key_tax, ranks = c("Genus","Family","Order")) {
    if (nrow(df_dom) == 0) return(df_dom %>% mutate(tax_label = NA_character_, tax_rank = NA_character_))
    
    tax_small <- tax %>%
      distinct(.data[[key_tax]], .keep_all = TRUE) %>%
      select(all_of(c(key_tax, ranks[ranks %in% names(tax)])))
    
    out <- df_dom %>%
      left_join(tax_small, by = c("taxon" = key_tax))
    
    picked <- .pick_lowest(out, ranks = ranks)
    out$tax_label <- picked$tax_label
    out$tax_rank  <- picked$tax_rank
    out
  }
  
  # Split df.hub
  df_pl  <- df.hub %>% filter(domain == "Plants")
  df_bac <- df.hub %>% filter(domain == "Bacteria")
  df_fun <- df.hub %>% filter(domain == "Fungi")
  
  # Annotate (REPLACE tax_* objects and key_* with your actual names)
  df_pl_a  <- annotate_domain(df_pl,  tax_pln, key_pln, ranks = c("Genus","Family","Order"))
  df_bac_a <- annotate_domain(df_bac, tax_bac, key_bac, ranks = c("Genus","Family","Order"))
  df_fun_a <- annotate_domain(df_fun, tax_fun, key_fun, ranks = c("Genus","Family","Order"))
  
  # Combine
  df.hub_annot <- bind_rows(df_pl_a, df_bac_a, df_fun_a)
  
  # Optional: fill unresolved with a placeholder
  df.hub_annot <- df.hub_annot %>%
    mutate(tax_label = ifelse(is.na(tax_label), "unclassified", tax_label),
           tax_rank  = ifelse(is.na(tax_rank),  "unclassified", tax_rank))
  
  # Inspect
  df.hub_annot %>% arrange(domain, desc(strength)) %>% head(15)
  
  

  # Optionally save:
   write.csv(  df.hub_annot, "results/Table_Sy_top_strength_hubs.csv", row.names = FALSE)
  
  
  
  
  
#reporting results
  

nP #retained plant taxa, 
nB #retained bacterial ASVs 
nF #retained fungal ASVs after filtering. 

#The cross-domain network contained 
overall$total_edges #associations 


#(median r = 

round(overall$median_r_overall, 2)

#mean ∣∣r∣ = 
round(overall$mean_absr_overall, 2))

#of which 
round(overall$pct_positive_overall,1)

#% were positive. By pair, we detected 

pair_summary$n_edges[pair_summary$pair=="Plants–Bacteria"] #Plants–Bacteria links 
round(pair_summary$pct_pos[ pair_summary$pair=="Plants–Bacteria"],1)#% positive; … Plants–Fungi, and r … Bacteria–Fungi (Table Sx). 

    
    