# Functions used in this RMarkdown file are used in the project. 
# To load the functions below read in this file ("source("microbiome_functions.Rmd")")

# Define a function to check and install packages
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) install.packages(new_packages)
}

# Summarise the number of samples, taxa and reads in a given ps object
ps_stats <- function(ps, stats, step_name) {
  stats <- add_row(stats, step = step_name, 
                   samples = nsamples(ps), 
                   taxa = ntaxa(ps), 
                   reads = sum(otu_table(ps)))
  return(stats)
}

# Filter unwanted taxa 
create_taxa_filter <- function(taxa) {
  # remove Chloroplast/Mitochondria/Archae/Eukaryota/no kingdom-level annotation
  # see Park et al. (Raes). Microbiome. 2019. / Winkler et al. Cell. 2020.
  taxa_filter <- (str_detect(replace_na(taxa[, "Family"], ""), "Mitochondria") |
                    str_detect(replace_na(taxa[, "Order"], ""), "Chloroplast") |
                    str_detect(replace_na(taxa[, "Kingdom"], ""), "Archae|Eukaryota") |
                    is.na(taxa[, "Kingdom"]))
  return(taxa_filter)
}

# Generate taxa names
create_taxa_names <- function(tax_table) {
  tax_names <- apply(tax_table, 1, function(x) {
    if(any(names(na.omit(x)) %in% "ASV")) {
      name <- str_c(tail(na.omit(x), 2), collapse = "_")
    } else { 
      name <- tail(na.omit(x), 1)
    }
    return(name)
  })
  
  tax_names %<>% 
    str_replace_all("-| ", "_") %>%
    str_remove_all("\\[|\\]|\\(|\\)")
  
  return(str_c(tax_names, 1:nrow(tax_table), sep = "_"))
}

# Create a phyloseq object manually
create_phylo_manual <- function(seqtab, taxa, meta) {
  ps <- phyloseq(otu_table(t(seqtab), taxa_are_rows=T), 
                 sample_data(meta), # TODO add meta
                 tax_table(taxa))
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- create_taxa_names(taxa)
  tax_table(ps) <- cbind(tax_table(ps), ASV = taxa_names(ps))
  return(ps)
}

# Create a filtered phyloseq object.
create_ps <- function(seqtab, taxa, meta) {
  #order mapping/seqtab
  seqtab_ordered <- seqtab[match(rownames(meta), rownames(seqtab)), ]
  #filter taxa
  taxa_filter <- create_taxa_filter(taxa)
  taxa_filtered <- taxa[!taxa_filter, ]
  seqtab_filtered <- seqtab_ordered[, !taxa_filter]
  n_taxa_excl <- nrow(taxa) - nrow(taxa_filtered)
  print(glue::glue("ASVs beloning to the phylogenetic groups 'Mitochondria', 'Chloroplast', 'Archae' or 'Eukaryota' were
                   excluded (n = {n_taxa_excl} ASVs), resulting in a total of {scales::comma(nrow(taxa))} ASVs after exclusion."))
  ps <- create_phylo_manual(seqtab_filtered, taxa_filtered, meta)
  return(ps)
}

# Creates an unfiltered phyloseq object.
create_crude_ps <- function(seqtab, species) {
  ps <- phyloseq(otu_table(t(seqtab), taxa_are_rows = T), tax_table(species %>% as.matrix))
  return(ps)
}
 
# Removes the 7th taxonomic rank.
remove_c7_taxa <- function(all_ps) {
  x <- vector("list")
  for(i in all_ps) {
    tax_table(i) <- tax_table(i)[,c(1,2,3,4,5,6,8)]
    x[[length(x) + 1]] <- i}
  return(x)}

# Merges crude phyloseq objects.
merge_ps_crude <- function(dada2_files) {
  all_ps <- dada2_files$ps_crude
  names(all_ps) <- dada2_files$name
  list2env(all_ps, envir = environment())
  ps_merged <-  do.call(merge_phyloseq, map(names(all_ps), as.symbol))
  return(ps_merged)
}

# Reorders taxa by abundance.
taxa_sum_reorder <- function(ps) {
  index <- order(taxa_sums(ps), decreasing = T)
  tax_table(ps) <- tax_table(ps)[index,]
  otu_table(ps) <- otu_table(ps)[,index]
  return(ps)
}

# To visualize qPCR and picogreen across Miseq runs
plot_sample_data <- function(df, x, y) {
  df %>%
    ggplot(aes(x = x, y = y, color=niche)) +
    geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.6, shape = 16) + 
    labs(color = "Sample type") +
    geom_boxplot(alpha = 0.4, outlier.shape = NA, fill = "white", show.legend = FALSE) + 
    theme(legend.position="top") + 
    scale_color_manual(values = c("burlywood3", "cadetblue3", "coral3", "coral4", 
                                  "darkolivegreen3", "darkolivegreen4", "#FFC857", "azure3")) + 
    guides(color = guide_legend(override.aes = list(alpha=1, size = 4, shape = 15), nrow = 1))
}

# Separates species into multiple columns.
separate_species <- function(ps) {
  seqtab <- as(otu_table(ps), "matrix")
  taxa <- as.data.frame(tax_table(ps)) %>% separate(Species, c('species1', 'species2', "species3"), sep = "/", remove = F) %>% 
    mutate(Species_one = ifelse(is.na(species2), species1, NA),
           Species_one = ifelse(is.na(Genus), NA, Species_one),
           Species_one = ifelse(is.na(Species_one), Genus, paste0(Genus, "_", Species_one)))
  taxa2 <- subset(taxa, select = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species_one")) %>% as.matrix()
  
  taxa2[taxa2 == "Prevotella_7"] <- "Prevotella"
  
  meta2 <- sample_data(ps) %>% as.data.frame()
  ps <- create_ps(seqtab, taxa2, meta2) 
  
  return(ps)
}

# Generate rarefaction curves for the samples.
rarecurve <- function (x, step = 1, sample, xlab = "Sample Size", 
                       ylab = "Species", 
                       label = TRUE, col, lty, tidy = FALSE, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  minobs <- min(x[x > 0])
  if (minobs > 1) 
    warning(gettextf("most observed count data have counts 1, 
                     but smallest count is %d", minobs))
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i], use.names = FALSE)
    }
    drop(suppressWarnings(rarefy(x[i, ], n)))
  })
  if (tidy) {
    len <- sapply(out, length)
    nm <- rownames(x)
    df <- data.frame(Site = factor(rep(nm, len), levels = nm), 
                     Sample = unlist(lapply(out, attr, which = "Subsample")), 
                     Species = unlist(out))
    return(df)
  }
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample)
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    abline(h = rare, lwd = 0.5)
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) {
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}

# The sample data (bacterial density) was checked using this function
plot_sample_data <- function(df, x, y) {
df %>% ggplot(aes(x = x, y = y, color=niche)) +
  geom_point(position=position_jitterdodge(), size = 2.5, 
             alpha = 0.6, shape = 16) + 
  labs(color = "Sample type") +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, fill = "white", 
               show.legend = FALSE) + 
  theme(legend.position="none") + 
  scale_color_manual(values = c("burlywood3", "cadetblue3", 
                                "coral4", "#FFC857", "darkolivegreen4")) + 
  guides(color = guide_legend(override.aes = list(alpha=1, size = 4, 
                                                  shape = 15), nrow = 1))}

# Check the composition of the different samples types
create_ordered_bar <- function(ps, type, n) {
  
  ps_RA <- ps %>%
    prune_samples(sample_data(ps)$niche == type, .) %>%
    to_RA()
  
  otu_RA_m <- as(otu_table(ps_RA), "matrix")
  
  bc <- vegdist(t(otu_RA_m), "bray")
  hc <- hclust(bc, method = "average")
  
  ps_RA %>%
    prep_bar(n = n) %>%
    mutate(sample_id = fct_relevel(sample_id, hc$labels[hc$order])) %>%
    create_bar(df_topn = ., n = n, ncol_legend = 1, name_legend = "ASV") +
    coord_flip() +
    theme(legend.position = "right") + 
    scale_fill_manual(values = c("#F0F0F0", paletteer_d("pals::stepped")[2:(n+1)])) +
    labs(fill = "ASV")
} 

# Plot contaminants from decontam
plot_contaminants <- function(data) {
  ggplot(data, aes(x = as.character(Miseq.Run), y = value)) +
    facet_wrap(~ ASV, scales = "free_y", ncol = 5) +
    theme(strip.text.x = element_markdown(), legend.position = "top") +
    geom_jitter(aes(color = factor(niche)), size = 1, alpha = 0.4, shape = 16) +
    labs(colour = "Niche") +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4, 
                                                     shape = 15))) +
    geom_boxplot(aes(color = factor(niche)), alpha = 0.4, 
                 outlier.shape = NA, fill = "white") +
    guides(col = "none") +
    scale_colour_manual(values = c("burlywood4", "cadetblue4")) +
    labs(y = "Relative abundance", x = "Isolation-run") +
    scale_y_continuous(trans = 'log10',
                       breaks = trans_breaks('log10', function(x) 10^x),
                       labels = trans_format('log10', math_format(10^.x)))
}


# Create excel file from gt table
export_to_excel <- function(table) {
  table %>%
    gtsummary::as_tibble() %>% 
    rename_with(~str_remove_all(., "\\**"))
} 

# Remove rows with NAs from the dataset
NA_remove <- function(ord, meta, vars_all) { 
  # Calculate the initial number of rows
  initial_rows <- nrow(meta)
  # Remove missing values for all variables in vars_all
  for (var in vars_all) {
    meta <- subset(meta, !is.na(meta[var]))
  }
  # Calculate the final number of rows
  final_rows <- nrow(meta)
  # Calculate the number of rows removed
  rows_removed <- initial_rows - final_rows
  # Subset the distance matrix and metadata based on remaining samples
  ord_sub <- ord[rownames(ord) %in% meta$sample_id, colnames(ord) %in% meta$sample_id]
  ord_sub <- as.dist(ord_sub)
  return(list(ord_sub = ord_sub, m_sub = meta, rows_removed = rows_removed))
}


# PERMANOVA functions 
beta_df <- function(ord, meta, var) {
  m_sub <- subset(meta, !is.na(meta[[var]]))
  ord_sub <- ord[rownames(ord) %in% m_sub$sample_id, colnames(ord) %in% m_sub$sample_id]
  list(ord_sub = as.dist(ord_sub), m_sub = m_sub)
}

do_permanova <- function(ord, meta, var) {
  adonis2(as.formula(glue("ord ~ {var}")), permutations = 999, data = meta) %>%
    broom::tidy() %>%
    subset(term == var) %>%
    mutate(
      label_p.value = case_when(
        p.value < 0.05 ~ "***",
        p.value < 0.1 ~ "**",
        p.value < 0.25 ~ "*",
        TRUE ~ NA_character_
      )
    )
}

do_permanova_com <- function(ord, meta, var) {
  data <- beta_df(ord, meta, var)
  do_permanova(data$ord_sub, data$m_sub, var)
}
