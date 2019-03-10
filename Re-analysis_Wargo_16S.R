# Preparation ------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)

library(pacman)
p_load(tidyverse, readxl, dada2, phyloseq, vegan, ggplot2, ggpubr, VennDiagram)

# Set palette ------------------------------------------------------------
# 74 distinctive colors in R
if (F) {
  library(RColorBrewer)
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  distinctive_colors <- unlist(mapply(brewer.pal, 
                                      qual_col_pals$maxcolors, 
                                      rownames(qual_col_pals))) %>% 
    .[-c(2)]
}

# Random Colors
if (F) {
  library(randomcoloR)
  palette <- distinctColorPalette(n)
}

# DADA2 pipeline ----------------------------------------------------------
if (F) {
  
  # Remove filtered folder before re-run dada2 !!!
  
  # Filtering -------------------------------------------------------------
  # the directory containing demultiplexed fastq.gz files
  path <- "/home/yeguanhua/Wargo/PRJEB22894/ERR2162225/demux"
  # Filtered files go into the filtered file
  filtpath <- file.path(path, "filtered")
  fns <- list.files(path, pattern="fastq.gz")
  
  out <- filterAndTrim(file.path(path,fns), file.path(filtpath,fns), 
                       truncLen=240,  maxEE=5, truncQ=2, 
                       rm.phix=TRUE, compress=TRUE, verbose=TRUE, 
                       multithread=TRUE)
  #out_print <- out %>% as.data.frame() %>% select(reads.in, reads.out)
  #out_print
  
  # Sample inference ------------------------------------------------------
  # File parsing
  filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE)
  sample.names <- dir("/home/yeguanhua/Wargo/PRJEB22894/ERR2162225/demux/filtered/") %>% 
    str_replace("\\.fastq\\.gz", "")
  names(filts) <- sample.names
  # Learn error rates
  set.seed(100)
  err <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE)
  # Infer sequence variants(won't dereplicate)
  dds <- vector("list", length(sample.names))
  names(dds) <- sample.names
  for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    derep <- derepFastq(filts[[sam]])
    dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
  }
  # Construct sequence table
  seqtab <- makeSequenceTable(dds)
  # Check the demension of sequence table
  #dim(seqtab)
  
  # Remove Chimeras -----------------------------------------------------
  seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                      multithread=TRUE)
  # Check the demension
  #dim(seqtab_nochim)
  # Inspect sequence remove rate:
  #sum(seqtab_nochim)/sum(seqtab)
  
  # Assign taxonomy -----------------------------------------------------
  tax_silva <- assignTaxonomy(seqtab, 
            "/home/yeguanhua/database/silva_nr_v132_train_set.fa.gz", 
            multithread=TRUE)
  tax_gg <- assignTaxonomy(seqtab, 
            "/home/yeguanhua/database/gg_13_8_train_set_97.fa.gz", 
            multithread=TRUE)
  
  # Track Reads -----------------------------------------------------------
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dds, getN), rowSums(seqtab_nochim))
  colnames(track) <- c("input", "filtered", "dereplicated", "nonchim")
  rownames(track) <- sample.names
  track_plot <- rownames_to_column(as.data.frame(track), var = "subject_id")
  track_plot$subject_id <- factor(track_plot$subject_id, 
                                  levels = track_plot$subject_id)
  track_plot <- gather(track_plot, 2:ncol(track_plot), 
                       key = "stages", value = "reads")
  track_plot$stages <- factor(track_plot$stages, 
                              levels = c("input", "filtered", "dereplicated", "nonchim"))
  ggplot(track_plot, 
         aes(x = stages, y = reads, color = subject_id)) + 
    scale_color_manual(values = distinctive_colors) + 
    geom_point() + 
    geom_line(aes(group = subject_id))
  
  # Write to disk
  #write.csv(seqtab_nochim, 
  #  file = "/home/yeguanhua/Wargo/PRJEB22894/ERR2162225/demux/filtered/wargo_seqtab_nochim.csv", 
  #  quote = F)
  #write.csv(tax_silva, 
  #  file = "/home/yeguanhua/Wargo/PRJEB22894/ERR2162225/demux/filtered/tax_silva.csv", 
  #  quote = F)
  #write.csv(tax_gg, 
  #  file = "/home/yeguanhua/Wargo/PRJEB22894/ERR2162225/demux/filtered/tax_gg.csv", 
  #  quote = F)
  
  dada2_res <- list()
  dada2_res$seq_tab <- seqtab_nochim
  dada2_res$tax_tab <- tax_silva
  dada2_res$reads_track <- track
}

# DADA2 results to OTU table ----------------------------------------------
if (F) {
  # merge reads to genus (by Jiangwei)
  read_to_genus=function(taxa,otu){
    map_data=taxa[colnames(otu),] %>% as.data.frame() %>% 
      mutate(id=paste(#.$Kingdom,
                      #.$Phylum,
                      #.$Class,
                      #.$Order,
                      .$Family,
                      .$Genus,
                      sep = " "),names=colnames(otu))
    otu_genus=NULL
    for (i in 1:length(unique(map_data$id))){
      if(is.na(str_match(unique(map_data$id)[i],"NA"))){
        #otu_merge %>% mutate(unique(map_data$id))[i]=sum(.[,filter(map_data,id==unique(map_data$id))[i])]))
        x=unique(filter(map_data,id==unique(map_data$id)[i])$names)
        if(length(x)==1){otu_genus[[unique(map_data$id)[i]]]=otu[,x]}
        else{
          otu_genus[[unique(map_data$id)[i]]]=unique(filter(map_data,id==unique(map_data$id)[i])$names) %>% 
            otu[,.] %>% apply(.,1,sum) %>% as.numeric()
        }
      }
    }
    otu_genus=as.data.frame(otu_genus)
  }
  otu_genus2 <- read_to_genus(taxa = tax_silva, otu = seqtab_nochim)
  
  # export DADA2 sequences to fasta file
  seqs2fasta = function(seqtab, out_path) {
    seqtab.t = as.data.frame(t(seqtab))
    seqs = row.names(seqtab.t)
    row.names(seqtab.t) = paste0("OTU", 1:nrow(seqtab.t))
    seqs = as.list(seqs)
    seqinr::write.fasta(sequences = seqs, names = row.names(seqtab.t), 
                        file.out = out_path)
  }
  # must have a exist output path
  seqs2fasta(seqtab = seqtab_nochim, 
             out_path = "/home/yeguanhua/Wargo/PRJEB22894/ERR2162225/demux/filtered/seqs.fasta")
  
  # Construct OTU table (by Yip)
  construct_otu_table <- function(seq, tax, level = "all", lefse = FALSE) {
    
    # Set options
    options(stringsAsFactors = FALSE)
    
    seq_tab <- t(seq) %>% as.data.frame()
    tax_tab <- as.data.frame(tax)
    for (i in rownames(tax_tab)) {
      if (tax_tab[rownames(tax_tab) == i, ] %>% is.na() %>% all()) {
        tax_tab <- subset(tax_tab, rownames(tax_tab) != i)
        seq_tab <- subset(seq_tab, rownames(seq_tab) != i)
      }
    }
    if (level == "all") {
      for (i in rownames(tax_tab)) {
        tax_tab[rownames(tax_tab) == i, 'Taxonomy'] <- tax[rownames(tax) == i,] %>% 
          na.omit() %>% str_c(collapse = "|")
      }
      tax_tab <- subset(tax_tab, select = Taxonomy)
    } else {
      tax_tab <- tax_tab %>% as.data.frame() %>% select(level)
      for (i in rownames(tax_tab)) {
        if (tax_tab[rownames(tax_tab) == i, ] %>% is.na()) {
          tax_tab <- subset(tax_tab, rownames(tax_tab) != i)
          seq_tab <- subset(seq_tab, rownames(seq_tab) != i)
        }
      }
      if (lefse == TRUE) {
        tax <- tax %>% as.data.frame() %>% select(1:which(colnames(tax) == level))
        for (i in rownames(tax_tab)) {
          tax_tab[rownames(tax_tab) == i, level] <- tax[rownames(tax) == i,] %>% 
            na.omit() %>% str_c(collapse = "|")
        }
      }
    }
    otu_tab <- left_join(rownames_to_column(seq_tab), rownames_to_column(tax_tab), 
                         by = "rowname") %>% .[,-1]
    otu_tab <- t(otu_tab)
    colnames(otu_tab) <- otu_tab[nrow(otu_tab),]
    otu_tab <- otu_tab[-(nrow(otu_tab)),]
    otu_tab0 <- apply(otu_tab, 2, as.numeric)
    rownames(otu_tab0) <- rownames(otu_tab)
    otu_tab <- rowsum(t(otu_tab0), group = colnames(otu_tab0)) %>% 
      as.data.frame()
    return(otu_tab)
  }
  otu_silva_all <- construct_otu_table(seq = seqtab_nochim, tax = tax_silva)
  otu_silva_Kingdom <- construct_otu_table(seq = seqtab_nochim, tax = tax_silva, 
                                         level = "Kingdom", lefse = TRUE)
  otu_silva_Genus <- construct_otu_table(seqtab_nochim, tax_silva, "Genus")
}

# Taxonomy comparison with original results -------------------------------
if (F) {
  # Load original taxonomy classification
  tax_origin <- read_excel(path = "/home/yeguanhua/Wargo/original_taxonomic_classification_top50.xlsx", sheet = 2)
  # Extract top 50 genus of GreenGene database
  gg_Genus <- rowSums(otu_gg_Genus) %>% as.data.frame()
  colnames(gg_Genus) <- "abundance"
  gg_Genus <- gg_Genus %>% rownames_to_column() %>% arrange(desc(abundance))
  gg_Genus <- gg_Genus[-2,]
  gg_Genus <- gg_Genus$rowname %>% .[1:length(unique(tax_origin$gg_Genus_origin))] %>% str_replace("g__", "")
  # Extract top 50 genus of Silva database
  silva_Genus <- rowSums(otu_silva_Genus) %>% as.data.frame()
  colnames(silva_Genus) <- "abundance"
  silva_Genus <- silva_Genus %>% rownames_to_column() %>% arrange(desc(abundance))
  silva_Genus <- silva_Genus$rowname %>% .[1:length(unique(tax_origin$silva_Genus_origin))]
  # Venn Diagram
  Genus4venn <- cbind(gg_Genus, silva_Genus)
  colnames(Genus4venn) <- c("gg_Genus", "silva_Genus")
  Genus4venn <- as.data.frame(Genus4venn)
  # Comparison of Genus of GreenGene database
  gg_Genus_venn <- venn.diagram(list("Original" = tax_origin$gg_Genus_origin, 
                                     "Re-analysis" = Genus4venn$gg_Genus), 
                                filename = NULL, fill = rainbow(2), 
                                cat.cex = 1, cex = 2.5, 
                                main = "Comparison of Genus overlap", 
                                sub = "GreenGene database")
  # Comparison of Genus of Silva database
  silva_Genus_venn <- venn.diagram(list("Original" = tax_origin$silva_Genus_origin, 
                                        "Re-analysis" = Genus4venn$silva_Genus), 
                                   filename = NULL, fill = rainbow(2), 
                                   cat.cex = 1, cex = 2.5, 
                                   main = "Comparison of Genus overlap", 
                                   sub = "Silva database")
  
  grid::grid.draw(silva_Genus_venn)
}

# OTU distribution --------------------------------------------------------
# plot OTU disdistribution via RAM::group.OTU
if (F) {
  #BiocManager::install(pkgs=c("RAM"))
  #install.packages("RAM")
  library(RAM)
  pdf("/home/yeguanhua/Wargo/otu_distribution_all.pdf")
  # taxonomy must be from Kingdom to Genus
  group.OTU(otu = otu_grouped, otuIDs = rownames(otu_grouped), 
            meta = metadata_wargo, meta.factor = "phenotype")
  dev.off()
}

# Genus level
if (F) {
  # Construct a long table for otu distribution
  otu_silva_genus <- construct_otu_table(seq = seqtab_nochim, tax = tax_silva, level = "Genus")
  
  # Construct otu for ditribution plot
  otu_distribution_genus <- otu_silva_genus %>% t() %>% as.data.frame() %>% rownames_to_column()
  
  colnames(otu_distribution_genus)[1] <- "subject_id"
  
  # Turn OTU to ggplot_type
  otu_distribution_genus <- gather(otu_distribution_genus, 
                                   2:(ncol(otu_distribution_genus)-1), 
                                   key = "Genus", value = "genus_abundance")
  
  # Convert abundance to log10 and remove infinite value
  otu_distribution_genus <- mutate(otu_distribution_genus, 
                                   `log10(genus_abundance)` = log10(genus_abundance)) %>% 
    filter(!is.infinite(.$`log10(genus_abundance)`))
  
  median_value <- otu_distribution_genus %>% 
    group_by(subject_id) %>% 
    dplyr::summarise(median_value = median(`log10(genus_abundance)`)) %>% 
    dplyr::arrange(median_value) %>% 
    left_join(as.data.frame(rownames_to_column(phenotype_wargo)), 
              by = c("subject_id" = "rowname"))
  
  otu_distribution_genus$subject_id <- factor(otu_distribution_genus$subject_id, 
                                              levels = median_value$subject_id)
  
  # 可以用arrange进行排序
  #otu_distribution_genus <- left_join(otu_distribution_genus, 
  #                                    median_value, 
  #                                    by = "subject_id") %>% 
  #  arrange(median_value)
  
  ggboxplot(otu_distribution_genus, 
            x = "subject_id", 
            y = "log10(genus_abundance)", 
            add = "jitter", 
            add.params = list(size = 0.5), 
            outlier.shape = NA, 
            fill = "subject_id", 
            palette = distinctive_colors) + 
    # phenotype的顺序需要和subject_id对应！！！
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     color = otu_distribution_genus$phenotype))
}

# Stacked bar plot -------------------------------------------------------
# phyloseq::plot_bar
if (F) {
  seqtab_percent <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
  ps2 <- phyloseq(tax_table(tax_silva), sample_data(metadata_wargo), 
                  otu_table(otu_percent, taxa_are_rows = F), phy_tree(fitGTR$tree))
  #colourCount = length(unique(metadata_wargo$subject_id))
  #getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  plot_bar(ps2, x = "Sample", y = "Abundance", fill = "Order") + 
    # set filling color
    #scale_fill_manual(values = getPalette(colourCount)) + 
    scale_fill_manual(values = distinctive_colors) + 
    #scale_fill_manual(values = distinctColorPalette(38)) + 
    # set lines-around-bars' color
    geom_bar(color = "transparent", stat = "identity")
}

# manually plot with ggplot
if (F) {
  seqtab_percent <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
  otu_silva_order_percent <- construct_otu_table(seq = seqtab_percent, 
                                                 tax = tax_silva, 
                                                 level = "Order")
  # Construct data frame for stacked bar plot
  stacked_bar_plot <- otu_silva_order_percent %>% t() %>% as.data.frame() %>% 
    # Order columns in descending order（从大到小）
    .[,order(colSums(.), decreasing = TRUE)] %>% 
    .[order(.[,1], decreasing = TRUE),] %>% 
    rownames_to_column()
  colnames(stacked_bar_plot)[1] <- "subject_id"
  # Join phenotype
  stacked_bar_plot <- left_join(stacked_bar_plot, 
                                as.data.frame(rownames_to_column(phenotype_wargo)), 
                                by = c("subject_id" = "rowname"))
  # Add levels to subject_id
  stacked_bar_plot$subject_id <- factor(stacked_bar_plot$subject_id, 
                                        levels = stacked_bar_plot$subject_id)
  # Turn data frame to ggplot_type
  stacked_bar_plot <- gather(stacked_bar_plot, 
                             colnames(stacked_bar_plot)[2:(ncol(stacked_bar_plot)-1)], 
                             key = "Order", value = "abundance")
  # Bar plot (Black names are R, red names are NR)
  ggplot(stacked_bar_plot, aes(x = subject_id, y = abundance)) + 
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     color = stacked_bar_plot$phenotype)) + 
    scale_fill_manual(values = distinctive_colors) + 
    geom_bar(mapping = aes(fill = Order), position = "fill", stat = "identity")
}

# Matson plot
if (F) {
  mt_phenotype <- matson_metadata %>% select(Run, subject_id, Response)
  mt_phenotype$Response <- factor(mt_phenotype$Response, 
                                  levels = c("Responder", "NonResponder"))
  
  mt_seqtab_percent <- t(apply(otu_mt_uniform, 1, function(x) x / sum(x))) %>% 
    as.data.frame()
  
  mt_otu_order <- construct_otu_table(seq = mt_seqtab_percent, tax = mt_taxon, 
                                      level = "Order")
  
  mt_stacked_bar <- mt_otu_order %>% t() %>% as.data.frame() %>% 
    # Order columns in descending order（从大到小）
    .[,order(colSums(.), decreasing = TRUE)] %>% 
    .[order(.[,1], decreasing = TRUE),] %>% 
    rownames_to_column()
  
  mt_stacked_bar$rowname <- str_replace(mt_stacked_bar$rowname, 
                                        "(SRR.+)_1.fastq.gz", "\\1")
  mt_stacked_bar <- left_join(mt_stacked_bar, 
                              mt_phenotype, 
                              by = c("rowname" = "Run")) %>% .[,-1]
  
  mt_stacked_bar$subject_id <- factor(mt_stacked_bar$subject_id, 
                                      levels = mt_stacked_bar$subject_id)
  
  mt_stacked_bar <- gather(mt_stacked_bar, 
                           colnames(mt_stacked_bar)[1:(ncol(mt_stacked_bar)-2)], 
                           key = "Order", value = "abundance")
  
  ggplot(mt_stacked_bar, aes(x = subject_id, y = abundance)) + 
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     color = mt_stacked_bar$Response)) + 
    scale_fill_manual(values = distinctive_colors) + 
    geom_bar(mapping = aes(fill = Order), position = "fill", stat = "identity")
}

# Rarefaction plot -------------------------------------------------------
if (F) {
  rarecurve(seqtab_nochim, 
            step = 40, 
            sample = min(rowSums(seqtab_nochim)), 
            col = c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink"), 
            lty = c("solid", "dashed", "longdash", "dotdash"), 
            label = F)
}

# Construct the phylogenetic tree ----------------------------------------
# Method 1: Biocunductor workflow
if (F) {
  seqs <- getSequences(seqtab_nochim)
  names(seqs) <- seqs #This propagates to the tip labels of the tree.(标签一直延伸到进化树顶端)
  
  # Version 1: using msa package
  #BiocManager::install(pkgs=c("msa"))
  #library(msa)
  #mult <- msa(seqs, method="ClustalW", type="dna", order="input")
  
  # Version 2: performing a multiple-alignment using the DECIPHER package
  #BiocManager::install(pkgs=c("DECIPHER"))
  library(DECIPHER)
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
  # The phangorn R package is then used to construct a phylogenetic tree
  #BiocManager::install(pkgs=c("phangorn"))
  library(phangorn)
  #phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab)) #v1
  phang.align <- phyDat(as(alignment, "matrix"), type="DNA")  #v2
  dm <- dist.ml(phang.align)
  treeNJ <- NJ(dm)  #Note: tip order != sequence order 
  fit = pml(treeNJ, data=phang.align)
  fitGTR <- update(fit, k=4, inv=0.2)
  start_time <- as.POSIXct(Sys.time())
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
                      rearrangement = "stochastic", 
                      control = pml.control(trace = 0))
  end_time <- as.POSIXct(Sys.time())
  # optim.pml run time
  difftime(end_time, start_time, units = "secs")
  detach("package:phangorn", unload=TRUE)
}

# Method 2: Qiime2
if (F) {
  #cd /home/yeguanhua/Wargo
  #sh tree.sh
  #rooted_tree <- read_tree(treefile = "/home/yeguanhua/Wargo/PRJEB22894/ERR2162225/demux/filtered/rooted-tree/tree/data/tree.nwk")
}

# Method 3: Use RAxML to contruct a maximum likelihood tree in R
if (F) {
  # Deploy RAxML in server
  # Download RAxML from https://github.com/stamatak/standard-RAxML
  #unzip standard-RAxML-master.zip
  #cd standard-RAxML-master
  #make -f Makefile.AVX.PTHREADS.gcc
  
  setwd("~/Wargo/PRJEB22894/ERR2162225/demux/filtered")
  seqs <- getSequences(seqtab_nochim)
  names(seqs) <- seqs
  # Performing a multiple-alignment using the DECIPHER R package
  #BiocManager::install(pkgs=c("DECIPHER"))
  library(DECIPHER)
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
  #BiocManager::install(pkgs=c("phangorn"))
  library(phangorn)
  phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
  write.phyDat(phang.align, 
               file="/home/yeguanhua/Wargo/PRJEB22894/ERR2162225/demux/filtered/alignment.fasta", 
               format="fasta")
  # Use RAxML to contruct a maximum likelihood tree in R
  #BiocManager::install(pkgs=c("ips"))
  library(ips)
  alignment <- read.dna(file = "/home/yeguanhua/Wargo/PRJEB22894/ERR2162225/demux/filtered/alignment.fasta", 
                        format="fasta", 
                        as.matrix=TRUE)
  fitGTR <- raxml(alignment, m = "GTRGAMMAIX", f = "a", p = 1234, x = 2345, N = 100, threads = 30, 
                  exec = "/home/yeguanhua/standard-RAxML-master/raxmlHPC-PTHREADS-AVX")
  
  # RAxML can't finish construct maximum likelihood tree after 13 hours using 30 cores.
  # Not uesable.
}

# Import data into phyloseq ----------------------------------------------
if (F) {
  # Read metadata
  wargo_metadata <- read_excel(path = "/home/yeguanhua/Wargo/Wargo_metadata.xlsx", 
                               sheet = 1)
  # remove Wargo.116774(has the least reads count)
  # remove Wargo.121585(mixed with Wargo.133270)
  wargo_metadata <- wargo_metadata[-c(3, 10),] %>% arrange(subject_id)
  rownames(wargo_metadata) <- wargo_metadata$subject_id
  
  otu <- seqtab_nochim[wargo_metadata$subject_id,]
  
  # use phyloseq to plot alpha diversity
  #wargo_sample_data <- sample_data(wargo_metadata)
  #wargo_otu_table <- otu_table(otu, taxa_are_rows = F)
  #phylo1 <- phyloseq(wargo_sample_data, wargo_otu_table)
  ps <- phyloseq(tax_table(tax_silva), sample_data(wargo_metadata), 
                 otu_table(otu, taxa_are_rows = F), 
                 phy_tree(fitGTR$tree))
}

# Taxonomic Filtering-----------------------------------------------------
# Phylum（门）
if (F) {
  # Extract table, number of features for each phyla from tax_silva
  table(tax_table(ps)[, "Phylum"], exclude = NULL)
  
  # Phylum of NA are probably artifacts and should be remove
  ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
  
  # Compute prevalence of each feature, store as data.frame
  prevdf <-  apply(X = otu_table(ps), 
                 # if taxa_are_rows=yes, MARGIN=1(apply FUN on every row)
                 # else MARGIN=2(apply FUN on every column)
                 MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2), 
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf <-  data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps), 
                      tax_table(ps))
  # Check phylum of low prevalence
  plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence), 
                                                    sum(df1$Prevalence))})
  # Manually filter Phylum of low prevalence 
  ps = subset_taxa(ps, !Phylum %in% c("Deferribacteres", "Incertae_Sedis"))
  
  ggplot(prevdf, aes(x = TotalAbundance, 
                     y = Prevalence / nsamples(ps), 
                     color = Phylum)) + 
    # Include a guess for parameter
    geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + 
    geom_point(size = 2, alpha = 0.7) + scale_x_log10() + 
    xlab("Total Abundance") + 
    ylab("Prevalence [Frac. Samples]") + facet_wrap( ~ Phylum) + 
    theme(legend.position="none")
}

if (F) {
  ## Taxonomic Filtering:
  # Remove Genus of NA
  ps <- subset_taxa(ps, !is.na(Genus) & !Genus %in% c("", "uncharacterized"))
  # Explore feature prevalence in the dataset:
  # Compute prevalence of each feature, store as data.frame
  prevdf <-  apply(X = otu_table(ps), 
                   MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2), 
                   FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf <- data.frame(Prevalence = prevdf, 
                       TotalAbundance = taxa_sums(ps), 
                       tax_table(ps))
  # Compute the total and average prevalences of the features in each Genus and filter Genus less than 1%
  filter_Genus <- plyr::ddply(prevdf, 
                              "Genus", 
                              function(df1){cbind(mean(df1$Prevalence), 
                                                  sum(df1$Prevalence))}) %>% 
    arrange(desc(`2`)) %>% filter(., `2`>(sum(.$`2`)*0.01))
  # Filter Genus
  ps <-  subset_taxa(ps, !Genus %in% filter_Genus$Genus)
  ## Taxa prevalence versus total counts
  # Subset to the remaining phyla
  prevdf <- subset(prevdf, Genus %in% get_taxa_unique(ps, "Genus"))
  ggplot(prevdf, aes(x = TotalAbundance, y = Prevalence / nsamples(ps),color = Genus)) + 
    # Include a guess for parameter
    geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + 
    geom_point(size = 2, alpha = 0.7) + scale_x_log10() + xlab("Total Abundance") + 
    ylab("Prevalence [Frac. Samples]") + facet_wrap( ~ Phylum) + 
    theme(legend.position="none")
  # In the Genus prevalence plot figure, every point is a different taxa.
  
  # After filtering taxonomy, alpha and beta diversity will be changed, so don't use it if not sure.
}

# Alpha diversity --------------------------------------------------------
if (F) {
  # Inverse Simpson diversity
  InvSimpson_alpha <- plot_richness(ps, x = "phenotype", measures = "InvSimpson")
  # Wilcox-test p-value (Mann-Whitney U test)
  InvSimpson_w_p <- wilcox.test(InvSimpson_alpha$data$value ~ phenotype_4_wilcox)$p.value
  # Box plot
  ggboxplot(InvSimpson_alpha$data, 
            x = "phenotype", 
            y = "value", 
            add = "jitter", 
            add.params = list(size = 3), 
            color = "phenotype", 
            outlier.shape = NA, 
            palette = c("#00AFBB", "#FC4E07")) + 
    ylab("InvSimpson Diversity") + 
    annotate("text", x = 2, y = 1, 
             label = paste("wilcox p-value = ", round(InvSimpson_w_p, 4), sep = ""), size = 3)
}

# Beta diversity ---------------------------------------------------------
# weighted-unifrac
if (F) {
  # Calculate PCoA by Bray-Curtis or weighted-Unifrac(phylogenetic tree required)
  #pcoa <- cmdscale(vegdist(otu, method = "bray"), k = 2, eig = T)
  pcoa <- cmdscale(phyloseq::distance(ps, method = "wunifrac"), k = 2, eig = T)
  
  pcoa_points <- as.data.frame(pcoa$points)
  
  colnames(pcoa_points) <- c("pc1", "pc2")
  
  data_ord <- merge(pcoa_points, wargo_metadata, by = "row.names")
  
  ggplot(data = data_ord, aes(x = pc1, y = pc2)) + 
    geom_point(aes(color = phenotype), size = 5) + 
    xlab(paste("PC1 ", round(100*as.numeric(pcoa$eig[1]/sum(pcoa$eig)), 2), 
               "%", sep = " ")) + 
    ylab(paste("PC2 ", round(100*as.numeric(pcoa$eig[2]/sum(pcoa$eig)), 2), 
               "%", sep = " "))
  
  # permanova test for Bray-Curtis(by Jiang Wei)
  #permanova_test=function(d,meta,feature,metric){
  #  id=which(!is.na(meta[,feature]) & meta[,feature]!="Others")
  #  res=adonis(d[id,] ~ meta[id,feature],data=meta,permutations = 999,method=metric)
  #  return(res)
  #}
  #perm <- permanova_test(d = otu, meta = medatada.wargo, feature = "phenotype", 
  #                       metric = "bray")
}

# Bray Curtis
if (F) {
  # Bray Curtis
  bc_pcoa <- cmdscale(vegdist(otu, method = "bray"), k = 2, eig = T)
  bc_mds1 <- as.data.frame(bc_pcoa$points)
  colnames(bc_mds1) <- c("pc1", "pc2")
  bc_data_ord <- merge(bc_mds1, metadata_wargo, by = "row.names")
  # Permanova test p-value
  bc_permanova <- adonis(otu ~ metadata_wargo$phenotype, 
                         data = metadata_wargo, 
                         permutations = 999, 
                         method = "bray")
  bc_permanova_p <- bc_permanova$aov.tab$`Pr(>F)`[1]
  # Beta diversity plots
  ggplot(data = bc_data_ord, aes(x = pc1, y = pc2)) + 
    geom_point(aes(color = phenotype), size = 5) + 
    xlab(paste("PC1 ", round(100*as.numeric(bc_pcoa$eig[1]/sum(bc_pcoa$eig)), 2), 
               "%", sep = " ")) + 
    ylab(paste("PC2 ", round(100*as.numeric(bc_pcoa$eig[2]/sum(bc_pcoa$eig)), 2), 
               "%", sep = " ")) + 
    annotate("text",
             x=max(bc_data_ord$pc1),
             y=min(bc_data_ord$pc2),
             label=paste("Permannova Test: p = ",bc_permanova_p,sep=""),
             hjust=1,vjust=-1,size=5)
}

# NMDS (Non-metric Multidimensional scaling) : 群落分异分析 --------------
if (F) {
  library(MASS)
  library(plotly)
  
  otu_dis <- vegdist(otu_table(ps), method = "mountford")
  
  otu_nmds <- metaMDS(otu_table(ps), distance = "mountford")
  
  stressplot(otu_nmds, otu_dis)
  
  otu_nmds_points <- merge(otu_nmds$points, phenotype_wargo, by = "row.names")
  
  ggplot(otu_nmds_points, aes(MDS1, MDS2, col = phenotype)) + 
    geom_point() + 
    theme(axis.title = element_text(size = 14), 
          axis.text = element_text(size = 12), 
          legend.title = element_blank()) + 
    scale_color_manual(values=c("#00AFBB", "#FC4E07"))
  
  #plot(otu_nmds)  #Red crosses are sequences(3211 tatal), white points are samples(40 tatal).
  #plot(otu_nmds,type="t", display="sites", choices = c(1,2))
  #orditorp(otu_nmds,display="species",col="red",air=0.01)
  #orditorp(otu_nmds,display="sites",cex=1.25,air=0.01)
  
}

# ANOSIM (Analysis of Group Similarities): 差异显著性分析（R和NR）--------
if (F) {
  ANOSIM_otu <- left_join(rownames_to_column(as.data.frame(otu)), 
                          rownames_to_column(phenotype_wargo))
  ANOSIM <- anosim(vegdist(otu), ANOSIM_otu$phenotype)
  plot(ANOSIM)
}

# PICRUSt: 通路注释 ------------------------------------------------------
# Export sequences to fasta file
if (F) {
  seqs2fasta = function(seqtab, out_path) {
    seqtab.t = as.data.frame(t(seqtab))
    seqs = row.names(seqtab.t)
    row.names(seqtab.t) = paste0("OTU", 1:nrow(seqtab.t))
    seqs = as.list(seqs)
    seqinr::write.fasta(sequences = seqs, names = row.names(seqtab.t), file.out = out_path)
  }
  seqs2fasta(seqtab = seqtab_nochim, out_path = "/home/yeguanhua/Wargo/PRJEB22894/ERR2162225/demux/filtered/seqs.fasta")
}

# BLAST otu_id by Shell commands
if (F) {
  system("cd /home/yeguanhua/Wargo/PICRUSt")
  system("/home/qinbingcai/work/bin/pip/metagenomics/v1/download/ncbi-blast/makeblastdb -in /home/DataShare/Database/Microbiome/marker_ref/gg_13_8_otus/rep_set/99_otus.fasta -input_type fasta -dbtype nucl -max_file_sz 1GB -out subject")
  system("/home/qinbingcai/work/bin/pip/metagenomics/v1/download/ncbi-blast/blastn -db subject -query seqs.fasta -out otu_blast.tsv -outfmt 6 -num_threads 10 -evalue 1e-5 -max_target_seqs 1")
}

# Construct PICRUSt_otu
if (F) {
  otu_id <- read_tsv(file = "/home/yeguanhua/Wargo/otu_blast.tsv", 
                     col_names = F) %>% select(X1, X2)
  
  PICRUSt_otu = as.data.frame(t(seqtab_nochim))
  row.names(PICRUSt_otu) = paste0("OTU", 1:nrow(PICRUSt_otu))
  PICRUSt_otu <- rownames_to_column(PICRUSt_otu)
  PICRUSt_otu <- left_join(otu_id, PICRUSt_otu, by = c("X1" = "rowname"))
  PICRUSt_otu <- PICRUSt_otu[,-1] %>% t()
  colnames(PICRUSt_otu) <- PICRUSt_otu[1,]
  PICRUSt_otu <- PICRUSt_otu[-1,]
  PICRUSt_otu <- rowsum(t(PICRUSt_otu), group = colnames(PICRUSt_otu)) %>% 
    as.data.frame()
}

# PICURSt in R
if (F) {
  #install.packages("themetagenomics")
  #BiocManager::install(pkgs=c("themetagenomics"))
  library(themetagenomics)
  
  ref <- "/home/yeguanhua/Wargo/themetagenomics_data"
  
  PICRUSt_result <- picrust(PICRUSt_otu, rows_are_taxa=TRUE, reference='gg_ko', 
                            reference_path=ref, cn_normalize=TRUE, 
                            sample_normalize=FALSE, drop=TRUE)
  
  PICRUSt_result$fxn_table[1:5,1:5]
  names(PICRUSt_result$fxn_meta)
  head(PICRUSt_result$fxn_meta$KEGG_Description)
  
  write.csv(PICRUSt_result$fxn_table, file = "PICRUSt_result_fxn_table.csv", row.names = TRUE)
  write.csv(PICRUSt_result$fxn_meta, file = "PICRUSt_result_fxn_meta.csv", row.names = TRUE)
  write.csv(PICRUSt_result$method_meta, file = "PICRUSt_result_method_meta.csv", row.names = TRUE)
  # manually combine 3 results to 1 file
}

# BJ & Matson PICRUSt
if (F) {
  library(themetagenomics)
  ref <- "/home/yeguanhua/Wargo/themetagenomics_data"
  
  source("/home/jiangwei/xbiome/tumor/tumor/response_classifier/R/pre_processing.R")
  # BJCaner results
  load("/home/jiangwei/xbiome/tumor/tumor/response_classifier/bj_otu_meta.RData")
  meta_info=as.data.frame(bjcancer_metadata_clean)
  row.names(meta_info)=meta_info$样本编号
  otu=otu_merge[meta_info$样本编号,]
  otu_bj=otu[,names(which(apply(otu,2,sum)>20))]
  colnames(otu_bj)=str_sub(colnames(otu_bj),1,204)
  otu_bj_uniform=t(rowsum(t(otu_bj),group=colnames(otu_bj)))
  #set.seed(1)
  bj_taxon=assignTaxonomy(otu_bj_uniform,
                          "/home/jiangwei/xbiome/tumor/pipeline-tuning/data/silva_nr_v132_train_set.fa.gz",
                          multithread = TRUE)
  # 保存
  save(bj_taxon, otu_bj_uniform, meta_info, file = "/home/yeguanhua/Wargo/PICRUSt/bj_res.rdata")
  # 导出序列文件
  seqs2fasta(seqtab = as.data.frame(otu_bj_uniform), 
             out_path = "/home/yeguanhua/Wargo/bj_seqs.fasta")
  # sequences blast
  system("/home/qinbingcai/work/bin/pip/metagenomics/v1/download/ncbi-blast/blastn -db subject -query bj_seqs.fasta -out bj_otu_blast.tsv -outfmt 6 -num_threads 10 -evalue 1e-5 -max_target_seqs 1")
  # Construct PICRUSt OTU
  bj_otu_id <- read_tsv(file = "/home/yeguanhua/Wargo/PICRUSt/bj_otu_blast.tsv", 
                        col_names = F) %>% select(X1, X2)
  bj_PICRUSt_otu = as.data.frame(t(otu_bj_uniform))
  row.names(bj_PICRUSt_otu) = paste0("OTU", 1:nrow(bj_PICRUSt_otu))
  bj_PICRUSt_otu <- rownames_to_column(bj_PICRUSt_otu)
  bj_PICRUSt_otu <- left_join(bj_otu_id, bj_PICRUSt_otu, by = c("X1" = "rowname"))
  bj_PICRUSt_otu <- bj_PICRUSt_otu[,-1] %>% t()
  colnames(bj_PICRUSt_otu) <- bj_PICRUSt_otu[1,]
  bj_PICRUSt_otu <- bj_PICRUSt_otu[-1,]
  bj_PICRUSt_otu <- rowsum(t(bj_PICRUSt_otu), group = colnames(bj_PICRUSt_otu)) %>% 
    as.data.frame()
  # PICRUSt
  bj_PICRUSt_result <- picrust(bj_PICRUSt_otu, rows_are_taxa=TRUE, reference='gg_ko', 
                               reference_path=ref, cn_normalize=TRUE, 
                               sample_normalize=FALSE, drop=TRUE)
  write.csv(bj_PICRUSt_result$fxn_table, 
            file = "/home/yeguanhua/Wargo/PICRUSt/bj_PICRUSt_result_fxn_table.csv", 
            row.names = TRUE)
  write.csv(bj_PICRUSt_result$fxn_meta, 
            file = "/home/yeguanhua/Wargo/PICRUSt/bj_PICRUSt_result_fxn_meta.csv", 
            row.names = TRUE)
  write.csv(bj_PICRUSt_result$method_meta, 
            file = "/home/yeguanhua/Wargo/PICRUSt/bj_PICRUSt_result_method_meta.csv", 
            row.names = TRUE)
  
  # Matson results
  load("/home/jiangwei/xbiome/tumor/tumor/response_classifier/Rdata/mt.Rdata")
  matson_metadata <- read_delim("/home/jiangwei/xbiome/tumor/tumor/response_classifier/data/HS16S_metadata.tsv", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
  otu_mt_prep=pre_process(mt_seqRes_double_exp)
  otu_mt=otu_mt_prep$Abd[,names(which(apply(otu_mt_prep$Abd,2,sum)>20))]
  colnames(otu_mt)=str_sub(colnames(otu_mt),1,204)
  otu_mt_uniform=t(rowsum(t(otu_mt),group=colnames(otu_mt)))
  otu_mt_rel=sweep(otu_mt_uniform,1,rowSums(otu_mt_uniform),'/')
  #set.seed(1)
  mt_taxon=assignTaxonomy(otu_mt_uniform,
                          "/home/jiangwei/xbiome/tumor/pipeline-tuning/data/silva_nr_v132_train_set.fa.gz",
                          multithread = TRUE)
  save(mt_taxon, otu_mt_uniform, matson_metadata, file="/home/yeguanhua/Wargo/PICRUSt/mt_res.rdata")
  # 导出序列文件
  seqs2fasta(seqtab = as.data.frame(otu_mt_uniform), 
             out_path = "/home/yeguanhua/Wargo/mt_seqs.fasta")
  # sequences blast
  system("/home/qinbingcai/work/bin/pip/metagenomics/v1/download/ncbi-blast/blastn -db subject -query mt_seqs.fasta -out mt_otu_blast.tsv -outfmt 6 -num_threads 10 -evalue 1e-5 -max_target_seqs 1")
  # Construct PICRUSt OTU
  mt_otu_id <- read_tsv(file = "/home/yeguanhua/Wargo/PICRUSt/mt_otu_blast.tsv", 
                        col_names = F) %>% select(X1, X2)
  mt_PICRUSt_otu = as.data.frame(t(otu_mt_uniform))
  row.names(mt_PICRUSt_otu) = paste0("OTU", 1:nrow(mt_PICRUSt_otu))
  mt_PICRUSt_otu <- rownames_to_column(mt_PICRUSt_otu)
  mt_PICRUSt_otu <- left_join(mt_otu_id, mt_PICRUSt_otu, by = c("X1" = "rowname"))
  mt_PICRUSt_otu <- mt_PICRUSt_otu[,-1] %>% t()
  colnames(mt_PICRUSt_otu) <- mt_PICRUSt_otu[1,]
  mt_PICRUSt_otu <- mt_PICRUSt_otu[-1,]
  mt_PICRUSt_otu <- rowsum(t(mt_PICRUSt_otu), group = colnames(mt_PICRUSt_otu)) %>% 
    as.data.frame()
  # PICRUSt
  mt_PICRUSt_result <- picrust(mt_PICRUSt_otu, rows_are_taxa=TRUE, reference='gg_ko', 
                               reference_path=ref, cn_normalize=TRUE, 
                               sample_normalize=FALSE, drop=TRUE)
  write.csv(mt_PICRUSt_result$fxn_table, 
            file = "/home/yeguanhua/Wargo/PICRUSt/mt_PICRUSt_result_fxn_table.csv", 
            row.names = TRUE)
  write.csv(mt_PICRUSt_result$fxn_meta, 
            file = "/home/yeguanhua/Wargo/PICRUSt/mt_PICRUSt_result_fxn_meta.csv", 
            row.names = TRUE)
  write.csv(mt_PICRUSt_result$method_meta, 
            file = "/home/yeguanhua/Wargo/PICRUSt/mt_PICRUSt_result_method_meta.csv", 
            row.names = TRUE)
}

# LEfSe ------------------------------------------------------------------
if (F) {
  lefse_input <- rbind(construct_otu_table(seqtab_nochim, tax_silva, "Kingdom", lefse = TRUE), 
                       construct_otu_table(seqtab_nochim, tax_silva, "Phylum", lefse = TRUE), 
                       construct_otu_table(seqtab_nochim, tax_silva, "Class", lefse = TRUE), 
                       construct_otu_table(seqtab_nochim, tax_silva, "Order", lefse = TRUE), 
                       construct_otu_table(seqtab_nochim, tax_silva, "Family", lefse = TRUE), 
                       construct_otu_table(seqtab_nochim, tax_silva, "Genus", lefse = TRUE)) %>% 
    as.data.frame()
  
  # convert to abundance
  lefse_input <- sweep(lefse_input, 2, colSums(lefse_input), '/') %>% rownames_to_column()
  
  lefse_phenotype <- phenotype_wargo %>% rownames_to_column() %>% t() %>% as.data.frame() %>% 
    rownames_to_column()
  
  lefse_phenotype[1,1] <- "subject_id"
  
  colnames(lefse_phenotype) <- colnames(lefse_input)
  
  lefse_phenotype <- lefse_phenotype[nrow(lefse_phenotype):1,]
  
  lefse_input <- rbind(lefse_phenotype, lefse_input)
  
  write_tsv(lefse_input, path = "/home/yeguanhua/Wargo/LEfSe/lefse_inpiut.tsv", col_names = FALSE)
}

# Matson
if (F) {
  # replace Run names with subject_id
  mt_subject_id <- select(matson_metadata, Run, subject_id)
  rownames(otu_mt_uniform) <- str_replace(rownames(otu_mt_uniform), "(.+)_1.fastq.gz", "\\1")
  for(i in 1:nrow(otu_mt_uniform)) {
    rownames(otu_mt_uniform)[i] <- mt_phenotype[mt_subject_id$Run == rownames(otu_mt_uniform)[i],'subject_id']
  }
  
  mt_lefse_input <- rbind(construct_otu_table(otu_mt_uniform, mt_taxon, "Kingdom", lefse = TRUE), 
                          construct_otu_table(otu_mt_uniform, mt_taxon, "Phylum", lefse = TRUE), 
                          construct_otu_table(otu_mt_uniform, mt_taxon, "Class", lefse = TRUE), 
                          construct_otu_table(otu_mt_uniform, mt_taxon, "Order", lefse = TRUE), 
                          construct_otu_table(otu_mt_uniform, mt_taxon, "Family", lefse = TRUE), 
                          construct_otu_table(otu_mt_uniform, mt_taxon, "Genus", lefse = TRUE)) %>% 
    as.data.frame()
  mt_lefse_input <- sweep(mt_lefse_input, 2, colSums(mt_lefse_input), '/') %>% 
    rownames_to_column()
  
  # metadata
  mt_response <- select(matson_metadata, subject_id, Response)
  mt_lefse_response <- mt_lefse_response %>% t() %>% as.data.frame() %>% rownames_to_column()
  mt_lefse_response <- mt_lefse_response[nrow(mt_lefse_response):1,]
  colnames(mt_lefse_response) <- colnames(mt_lefse_input)
  
  # combine
  mt_lefse_input <- rbind(mt_lefse_response, mt_lefse_input)
  
  write_tsv(mt_lefse_input, path = "/home/yeguanhua/Wargo/Matson/mt_lefse_inpiut.tsv", col_names = FALSE)
  
  # in Linux
  system("/home/yeguanhua/Wargo/Matson")  # the path that has lefse_input.tsv and run.sh
  system("docker run -v $PWD:/test -u root:root -t -i lefse:v2")
  # go to /test
  system("cd ..")
  system("cd ..")
  system("cd test/")
  system("sh run.sh")
  # done
}

# Log2 Fold change -------------------------------------------------------
# Manually
if (F) {
  #a <- 1:21
  #b <- 21:1
  #library(gtools)
  #f <- foldchange(a,b)
  #cbind(a,b,f)
  
  #library(gtools)
  fold_change <- function(seq, tax, tax_level) {
    fc_tax_level <- construct_otu_table(seq = seq, tax = tax, level = tax_level)
    for (i in colnames(fc_tax_level)) {
      colnames(fc_tax_level)[colnames(fc_tax_level)==i] <- phenotype_wargo[i, "phenotype"] %>% 
        as.character()
    }
    fc_tax_level <- rowsum(t(fc_tax_level), group = colnames(fc_tax_level)) %>% 
      t() %>% as.data.frame()
    fc_tax_level <- fc_tax_level[,c("R", "NR")]
    # 将0替换成0.001
    fc_tax_level[fc_tax_level == 0] <- 0.001
    fc <- foldchange(fc_tax_level$R, fc_tax_level$NR) %>% as.data.frame()
    rownames(fc) <- rownames(fc_tax_level)
    colnames(fc) <- "fold_change"
    fc <- cbind(fc_tax_level, fc)
    return(fc)
  }
  
  fc <- fold_change(seq=seqtab_nochim, tax = tax_silva, tax_level = "Genus")
  
  # Mann-Whitney-Wilcoxon Test计算p-value
  #wilcox.test(genus ~ phenotype)
}

# DESeq2 (Genus level)
if (F) {
  
  library(DESeq2)
  
  deseq2_dds <- DESeqDataSetFromMatrix(countData = otu_silva_genus, 
                                       colData = phenotype_wargo, 
                                       design = ~ phenotype)

  deseq2_dds <- DESeq(deseq2_dds)
  
  deseq2_res <- results(deseq2_dds, addMLE = FALSE)
  
  # order our results table by the smallest p-value
  deseq2_res <- deseq2_res[order(deseq2_res$pvalue),]
  
  #res_significant <- subset(res, pvalue <= 0.05)
  
  #plotMA(res, alpha=0.05, ylim=c(-5, 5))
  
  #plot(res$log2FoldChange, -log10(res$pvalue), log="x", 
  #     xlab="log2(FC)", 
  #     ylab=expression(-log[10](pvalue)))
  
  deseq2_plot <- deseq2_res %>% as.data.frame()
  colnames(deseq2_plot) <- colnames(deseq2_res)
  rownames(deseq2_plot) <- rownames(deseq2_res)
  
  ggplot(data = deseq2_plot, 
         aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(alpha = 0.5, size = 1.75) +
    theme(legend.position = "none") +
    xlim(c(-5, 5)) + ylim(c(0, 2)) +
    labs(x = expression(log[2](FC)), y = expression(-log[10](P))) +
    theme(axis.title.x=element_text(size=20), 
          axis.text.x=element_text(size=15)) +
    theme(axis.title.y=element_text(size=20),
          axis.text.y=element_text(size=15))
  
  # 与LEfSe的结果进行对比：
  # 提取pvalue < 0.5, |log2FC| > 1 的结果:
  # 保留rownames
  deseq2_res2 <- deseq2_plot[which(deseq2_plot$pvalue < 0.5 & abs(deseq2_plot$log2FoldChange) > 1), ] %>% 
    rownames_to_column()
  # 无法保留rownames
  #filter(deseq2_plot, abs(deseq2_plot$log2FoldChange)>1 & deseq2_plot$pvalue<0.5)
  
  lefse_res <- read_xlsx(path = "/home/yeguanhua/Wargo/LEfSe/lefse_res.xlsx", col_names = FALSE)
  lefse_res <- lefse_res[which(!is.na(lefse_res$X__3)),]
  for (i in 1:nrow(lefse_res)) {
    lefse_res[i,1] <- str_match_all(lefse_res$X__1[i], "\\.(\\w+)") %>% unlist() %>% .[length(.)]
  }
  
  foldchange_venn <- venn.diagram(list("Log2FC" = deseq2_res2$rowname, "LEfSe" = lefse_res$X__1), 
                                  filename = NULL, fill = c("#00AFBB", "#FC4E07"), cat.cex = 1, 
                                  cex = 2.5, main = "Taxonomy overlap between Log2FC and LEfSe")
  grid::grid.draw(foldchange_venn)
}

# Matson (Genus level)
if (F) {
  
  mt_genus <- construct_otu_table(seq = otu_mt_uniform, tax = mt_taxon, level = "Genus")
  
  library(DESeq2)
  
  mt_deseq2_dds <- DESeqDataSetFromMatrix(countData = mt_genus, 
                                          colData = mt_response, 
                                          design = ~ Response)
  
  mt_deseq2_dds <- DESeq(mt_deseq2_dds)
  
  mt_deseq2_dds <- results(mt_deseq2_dds, addMLE = FALSE)
  
  # order our results table by the smallest p-value
  mt_deseq2_dds <- mt_deseq2_dds[order(mt_deseq2_dds$pvalue),]
  
  write.csv(mt_deseq2_dds, file = "/home/yeguanhua/Wargo/Matson/mt_log2fc.csv", 
            quote = FALSE)
  
  mt_deseq2_plot <- mt_deseq2_dds %>% as.data.frame()
  colnames(mt_deseq2_plot) <- colnames(mt_deseq2_dds)
  rownames(mt_deseq2_plot) <- rownames(mt_deseq2_dds)
  
  ggplot(data = mt_deseq2_plot, 
         aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(alpha = 0.5, size = 1.75) +
    theme(legend.position = "none") +
    xlim(c(-5, 5)) + ylim(c(0, 2)) +
    labs(x = expression(log[2](FC)), y = expression(-log[10](P))) +
    theme(axis.title.x=element_text(size=20), 
          axis.text.x=element_text(size=15)) +
    theme(axis.title.y=element_text(size=20),
          axis.text.y=element_text(size=15))
}

# Hierarchical multiple testing ------------------------------------------
# Using sequences table
if (F) {
  #BiocManager::install(pkgs = c("structSSI"))
  library(structSSI)
  pslog <- transform_sample_counts(ps, function(x) log(1 + x))
  ps_dds <- phyloseq_to_deseq2(ps, ~ phenotype)
  # Add function to provide error:
  # every gene contains at least one zero, cannot compute log geometric means
  gm_mean = function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
  geoMeans = apply(counts(ps_dds), 1, gm_mean)
  
  ps_dds <- estimateSizeFactors(ps_dds, geoMeans = geoMeans)
  ps_dds <- estimateDispersions(ps_dds)
  abund <- getVarianceStabilizedData(ps_dds)
  short_names <- substr(rownames(abund), 1, 5) %>% make.names(unique = TRUE)
  rownames(abund) <- short_names
  el <- phy_tree(pslog)$edge
  el0 <- el[nrow(el):1, ]
  el_names <- c(short_names, seq_len(phy_tree(pslog)$Nnode))
  el[, 1] <- el_names[el0[, 1]]
  el[, 2] <- el_names[as.numeric(el0[, 2])]
  unadj_p <- treePValues(el, abund, sample_data(pslog)$phenotype)
  hfdr_res <- hFDR.adjust(unadj_p, el, .75)
  abund_sums <- rbind(data.frame(sum = colSums(abund), 
                                 sample = colnames(abund), 
                                 type = "DESeq2"), 
                      data.frame(sum = rowSums(otu_table(pslog)), 
                                 sample = rownames(otu_table(pslog)), 
                                 type = "log(1 + x)"))
  ggplot(abund_sums) + 
    geom_histogram(aes(x = sum), binwidth = 20) + 
    facet_grid(type ~ .) + 
    xlab("Total abundance within sample")
  plot(hfdr_res, height = 5000)  #opens in a browser
}

# Using other otu table
if (F) {
  dds_genus <- DESeqDataSetFromMatrix(countData = otu_silva_Genus, 
                                      colData = phenotype_wargo, 
                                      design = ~ phenotype)
  dds_genus <- DESeq(dds_genus)
  ps_genus <- phyloseq(tax_table(tax_silva), sample_data(metadata_wargo), 
                       otu_table(otu_silva_Genus, taxa_are_rows = T), 
                       phy_tree(fitGTR$tree))
  pslog <- transform_sample_counts(ps_genus, function(x) log(1 + x))
  dds_genus <- estimateSizeFactors(dds_genus)
  dds_genus <- estimateDispersions(dds_genus)
  abund <- getVarianceStabilizedData(dds_genus)
  #short_names <- substr(rownames(abund), 1, 5) %>% make.names(unique = TRUE)
  #rownames(abund) <- short_names
  el <- phy_tree(pslog)$edge
  #el0 <- el[nrow(el):1, ]
  #el_names <- c(short_names, seq_len(phy_tree(pslog)$Nnode))
  #el[, 1] <- el_names[el0[, 1]]
  #el[, 2] <- el_names[as.numeric(el0[, 2])]
  unadj_p <- treePValues(el, abund, sample_data(pslog)$phenotype)
  hfdr_res <- hFDR.adjust(unadj_p, el, .75)
  abund_sums <- rbind(data.frame(sum = colSums(abund), 
                                 sample = colnames(abund), 
                                 type = "DESeq2"), 
                      data.frame(sum = rowSums(otu_table(pslog)), 
                                 sample = rownames(otu_table(pslog)), 
                                 type = "log(1 + x)"))
  ggplot(abund_sums) + 
    geom_histogram(aes(x = sum), binwidth = 20) + 
    facet_grid(type ~ .) + 
    xlab("Total abundance within sample")
}

# DPCoA ------------------------------------------------------------------
if (F) {
  pslog <- transform_sample_counts(ps, function(x) log(1 + x))
  
  out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
  
  evals <- out.dpcoa.log$eig
  
  # DPCoA for samples
  plot_ordination(pslog, 
                  out.dpcoa.log, 
                  type = "samples", 
                  color = "phenotype") + 
    coord_fixed(sqrt(evals[2] / evals[1])) + 
    scale_color_manual(values=c("#00AFBB", "#FC4E07"))
  
  # DPCoA for taxonomy
  plot_ordination(pslog, 
                  out.dpcoa.log, 
                  type = "species", 
                  color = "Order") + 
    coord_fixed(sqrt(evals[2] / evals[1])) + 
    scale_color_manual(values = distinctive_colors)
}

# K-means ----------------------------------------------------------------
if (F) {
  # Construct percentage OTU table
  otu_silva_order_percent <- construct_otu_table(seq = seqtab_percent, 
                                                 tax = tax_silva, 
                                                 level = "Order")
  # Construct K-means input table
  # Rows are subject_id, columns are Order
  km_otu <- otu_silva_order_percent %>% t() %>% as.data.frame() %>% 
    .[,order(colSums(.), decreasing = TRUE)] %>% rownames_to_column()
  colnames(km_otu)[1] <- "subject_id"
  # Join phenotype
  km_otu <- left_join(km_otu, rownames_to_column(as.data.frame(phenotype_wargo)), by = c("subject_id" = "rowname"))
  # Calculate k-means
  km <- kmeans(km_otu[,2:(ncol(km_otu)-1)], centers = 2)
  # Plot k-means on the two most abundant taxonomy, labeled by R vs NR
  km_plot <- km_otu[,-c(4:(ncol(km_otu)-1))] %>% cbind(as.factor(km$cluster))
  colnames(km_plot)[ncol(km_plot)] <- "cluster"
  ggplot(km_plot, aes(Bacteroidales, Clostridiales)) + 
    geom_point(mapping = aes(color = phenotype, shape = cluster), size = 4) + 
    scale_color_manual(values=c("#00AFBB", "#FC4E07"))
}
