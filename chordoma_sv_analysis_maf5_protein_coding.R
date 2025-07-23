#!/usr/bin/env Rscript
# Author: Dongjing Wu
# Chordoma SV Analysis Following Gong et al. 2025 Nature Communications
# Complete workflow from VCF to potentially pathogenic SVs
# Using MAF-filtered VCF and protein-coding genes only
# VERSION: CORRECTED for SVTYPE misclassification bug

# Set up environment ----------------------------------------------------------
cat("Setting up environment...\n")
setwd("/data/wud7/sv_analysis/MAF_filtering")

# Check and install required packages
if (!require("StructuralVariantUtil")) {
  if (!require("devtools")) install.packages("devtools")
  devtools::install_github("tgong1/StructuralVariantUtil")
}
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")

# Load required libraries
library(StructuralVariantUtil)
library(dplyr)
library(ggplot2)

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

# Define file paths - UPDATED FOR MAF FILTERING AND PROTEIN-CODING GENES
vcf_file <- "/data/wud7/sv_analysis/MAF_filtering/MAF_filter_VCF/chordoma.maf5.vcf"
gene_bed_file <- "/data/wud7/sv_analysis/MAF_filtering/gene_annotations/hg19_protein_coding_genes.bed"
bedtools_path <- "/usr/local/apps/bedtools/2.31.1/bin/bedtools"

# Step 1: Load gene annotations -----------------------------------------------
cat("\n=== STEP 1: Loading Protein-Coding Gene Annotations ===\n")

# Load pre-processed protein-coding genes
gene_bed_simple <- read.table(gene_bed_file, 
                             header = TRUE, 
                             sep = "\t", 
                             stringsAsFactors = FALSE)

cat("Loaded", nrow(gene_bed_simple), "protein-coding genes\n")
cat("Chromosome distribution:\n")
print(table(gene_bed_simple$chrom))

# Step 2: Process VCF and classify SVs ----------------------------------------
cat("\n=== STEP 2: Processing MAF-Filtered VCF File ===\n")

# Check VCF file
if (!file.exists(vcf_file)) {
  stop("VCF file not found!")
}

file_info <- file.info(vcf_file)
cat("VCF file:", basename(vcf_file), "\n")
cat("File size:", round(file_info$size / 1024^2, 2), "MB\n")

# Convert VCF to dataframe
cat("\nConverting VCF to dataframe...\n")
df_vcf <- vcf_to_dataframe(vcf_file)
cat("Total variants in MAF-filtered VCF:", nrow(df_vcf), "\n")

# Show original SVTYPE distribution from VCF
cat("\nOriginal SVTYPE distribution from VCF INFO field:\n")
if("INFO_SVTYPE" %in% names(df_vcf)) {
  print(table(df_vcf$INFO_SVTYPE))
}

# Classify SVs
cat("\nClassifying structural variants...\n")
sv_bedpe <- simple_SVTYPE_classification(SV_data = vcf_file, 
                                        caller_name = "chordoma_sv")

# Show potentially misclassified distribution
cat("\nSV Type Distribution after simple_SVTYPE_classification (POTENTIALLY INCORRECT):\n")
print(table(sv_bedpe$SVTYPE))

# FIX THE SVTYPE MISCLASSIFICATION BUG
cat("\n=== FIXING SVTYPE MISCLASSIFICATION BUG ===\n")
cat("The simple_SVTYPE_classification function incorrectly classifies many variants as INS.\n")
cat("Applying correction based on original VCF INFO field...\n")

# Check if we have INFO_SVTYPE from vcf_to_dataframe
if("INFO_SVTYPE" %in% names(df_vcf)) {
  # Check that row counts match
  if(nrow(df_vcf) == nrow(sv_bedpe)) {
    # Save original misclassified SVTYPE for reference
    sv_bedpe$SVTYPE_misclassified <- sv_bedpe$SVTYPE
    
    # Apply correct SVTYPE from VCF INFO field
    sv_bedpe$SVTYPE <- df_vcf$INFO_SVTYPE
    
    # Convert BND to TRA (standard bedpe convention)
    sv_bedpe$SVTYPE[sv_bedpe$SVTYPE == "BND"] <- "TRA"
    
    cat("\nCorrected SVTYPE distribution:\n")
    print(table(sv_bedpe$SVTYPE))
    
    # Show correction summary
    mismatches <- sum(sv_bedpe$SVTYPE != sv_bedpe$SVTYPE_misclassified)
    cat("\nCorrected", mismatches, "misclassified variants out of", nrow(sv_bedpe), "\n")
    
    # Show what types were misclassified
    if(mismatches > 0) {
      cat("\nMisclassification summary:\n")
      mismatch_table <- table(
        Original = sv_bedpe$SVTYPE[sv_bedpe$SVTYPE != sv_bedpe$SVTYPE_misclassified],
        Misclassified = sv_bedpe$SVTYPE_misclassified[sv_bedpe$SVTYPE != sv_bedpe$SVTYPE_misclassified]
      )
      print(mismatch_table)
    }
  } else {
    stop("Row count mismatch between VCF and bedpe - cannot apply correction!")
  }
} else {
  warning("INFO_SVTYPE not found in VCF data - cannot correct misclassification!")
}

# Show final corrected SV distribution
sv_summary <- table(sv_bedpe$SVTYPE)
cat("\n=== FINAL CORRECTED SV Type Distribution (MAF ≤ 5%) ===\n")
print(sv_summary)
cat("Total SVs after MAF filtering:", sum(sv_summary), "\n")

# Calculate percentages
cat("\nPercentages:\n")
for(type in names(sv_summary)) {
  cat(sprintf("%s: %d (%.1f%%)\n", type, sv_summary[type], 
              100 * sv_summary[type] / sum(sv_summary)))
}

# Save corrected classified SVs
write.table(sv_bedpe, "results/chordoma_sv_classified_maf5_CORRECTED.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Step 3: Annotate SVs with protein-coding genes -----------------------------
cat("\n=== STEP 3: Protein-Coding Gene Annotation ===\n")

# Custom annotation function
annotate_svs_custom <- function(sv_data, gene_bed) {
  sv_data$pos1_overlap_gene <- NA
  sv_data$pos2_overlap_gene <- NA
  
  gene_file <- tempfile(fileext = ".bed")
  write.table(gene_bed[, c("chrom", "start", "end", "gene_name")], 
              gene_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  
  annotate_breakpoints <- function(chroms, starts, ends, ids) {
    valid <- !is.na(chroms) & !is.na(starts) & !is.na(ends)
    if (sum(valid) == 0) return(rep(NA, length(chroms)))
    
    sv_file <- tempfile(fileext = ".bed")
    sv_df <- data.frame(
      chrom = chroms[valid],
      start = starts[valid],
      end = ends[valid],
      id = ids[valid]
    )
    write.table(sv_df, sv_file, sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = FALSE)
    
    out_file <- tempfile(fileext = ".bed")
    cmd <- paste(bedtools_path, "intersect -a", sv_file, "-b", gene_file, "-wa -wb >", out_file)
    system(cmd, ignore.stderr = TRUE)
    
    result <- rep(NA, length(chroms))
    if (file.info(out_file)$size > 0) {
      overlaps <- read.table(out_file, sep = "\t", stringsAsFactors = FALSE)
      colnames(overlaps) <- c("sv_chr", "sv_start", "sv_end", "sv_id",
                             "gene_chr", "gene_start", "gene_end", "gene_name")
      
      gene_list <- aggregate(gene_name ~ sv_id, data = overlaps, 
                           FUN = function(x) paste(unique(x), collapse = ";"))
      
      for (i in 1:nrow(gene_list)) {
        idx <- which(ids == gene_list$sv_id[i])
        if (length(idx) > 0) {
          result[idx] <- gene_list$gene_name[i]
        }
      }
    }
    
    file.remove(sv_file, out_file)
    return(result)
  }
  
  cat("Annotating breakpoint 1...\n")
  sv_data$pos1_overlap_gene <- annotate_breakpoints(
    sv_data$chrom1, sv_data$start1, sv_data$end1, sv_data$ID
  )
  
  cat("Annotating breakpoint 2...\n")
  sv_data$pos2_overlap_gene <- annotate_breakpoints(
    sv_data$chrom2, sv_data$start2, sv_data$end2, sv_data$ID_mate
  )
  
  file.remove(gene_file)
  
  # Create gene fusion summary
  fusions <- sv_data[!is.na(sv_data$pos1_overlap_gene) & 
                     !is.na(sv_data$pos2_overlap_gene) &
                     sv_data$pos1_overlap_gene != sv_data$pos2_overlap_gene, ]
  
  if (nrow(fusions) > 0) {
    gene_fusions <- data.frame(
      gene1 = fusions$pos1_overlap_gene,
      gene2 = fusions$pos2_overlap_gene,
      sv_type = fusions$SVTYPE,
      chrom1 = fusions$chrom1,
      pos1 = fusions$start1,
      chrom2 = fusions$chrom2,
      pos2 = fusions$start2,
      stringsAsFactors = FALSE
    )
  } else {
    gene_fusions <- data.frame(
      gene1 = character(0),
      gene2 = character(0),
      sv_type = character(0),
      chrom1 = character(0),
      pos1 = integer(0),
      chrom2 = character(0),
      pos2 = integer(0),
      stringsAsFactors = FALSE
    )
  }
  
  return(list(sv_annotated = sv_data, gene_fusions = gene_fusions))
}

# Run annotation
cat("\nRunning gene annotation with protein-coding genes only...\n")
annotation_results <- annotate_svs_custom(sv_bedpe, gene_bed_simple)

sv_annotated <- annotation_results$sv_annotated
gene_fusions <- annotation_results$gene_fusions

# Summary
svs_with_genes <- sum(!is.na(sv_annotated$pos1_overlap_gene) | 
                      !is.na(sv_annotated$pos2_overlap_gene))

cat("\n=== Gene Annotation Summary ===\n")
cat("Total SVs (MAF ≤ 5%):", nrow(sv_annotated), "\n")
cat("SVs overlapping protein-coding genes:", svs_with_genes, 
    "(", round(100 * svs_with_genes / nrow(sv_annotated), 1), "%)\n")
cat("Potential gene fusions:", nrow(gene_fusions), "\n")

# Calculate SV type distribution for gene-overlapping SVs
svs_with_genes_df <- sv_annotated[!is.na(sv_annotated$pos1_overlap_gene) | 
                                   !is.na(sv_annotated$pos2_overlap_gene), ]

cat("\n=== SV Type Distribution for Gene-Overlapping SVs ===\n")
gene_sv_type_dist <- table(svs_with_genes_df$SVTYPE)
cat("Total SVs overlapping genes:", nrow(svs_with_genes_df), "\n")
print(gene_sv_type_dist)
cat("\nPercentages of gene-overlapping SVs:\n")
for(type in names(gene_sv_type_dist)) {
  cat(sprintf("%s: %d (%.1f%%)\n", type, gene_sv_type_dist[type], 
              100 * gene_sv_type_dist[type] / sum(gene_sv_type_dist)))
}

# Save results
write.table(sv_annotated, "results/chordoma_sv_gene_annotated_maf5_protein_coding_CORRECTED.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

if (nrow(gene_fusions) > 0) {
  write.table(gene_fusions, "results/chordoma_gene_fusions_maf5_protein_coding.tsv", 
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# Step 4: Categorize functional impacts ---------------------------------------
cat("\n=== STEP 4: Functional Impact Categorization ===\n")

sv_annotated$gene_impact_type <- NA
sv_annotated$affected_genes <- NA

# Count impact types as we go
impact_counts <- list(DEL_pLoF = 0, DUP_IED = 0, DUP_CG = 0, INS_pLoF = 0, INV_pLoF = 0, TRA_pLoF = 0)

for (i in 1:nrow(sv_annotated)) {
  sv <- sv_annotated[i, ]
  
  genes <- unique(c(
    if (!is.na(sv$pos1_overlap_gene)) strsplit(sv$pos1_overlap_gene, ";")[[1]] else NULL,
    if (!is.na(sv$pos2_overlap_gene)) strsplit(sv$pos2_overlap_gene, ";")[[1]] else NULL
  ))
  
  if (length(genes) > 0) {
    sv_annotated$affected_genes[i] <- paste(unique(genes), collapse = ";")
    
    if (sv$SVTYPE == "DEL") {
      sv_annotated$gene_impact_type[i] <- "pLoF"
      impact_counts$DEL_pLoF <- impact_counts$DEL_pLoF + 1
      
    } else if (sv$SVTYPE == "DUP") {
      genes1 <- if (!is.na(sv$pos1_overlap_gene)) strsplit(sv$pos1_overlap_gene, ";")[[1]] else character(0)
      genes2 <- if (!is.na(sv$pos2_overlap_gene)) strsplit(sv$pos2_overlap_gene, ";")[[1]] else character(0)
      
      if (length(genes1) > 0 && length(genes2) > 0 && any(genes1 %in% genes2)) {
        sv_annotated$gene_impact_type[i] <- "IED"
        impact_counts$DUP_IED <- impact_counts$DUP_IED + 1
      } else {
        sv_annotated$gene_impact_type[i] <- "CG"
        impact_counts$DUP_CG <- impact_counts$DUP_CG + 1
      }
      
    } else if (sv$SVTYPE %in% c("INV", "TRA")) {
      sv_annotated$gene_impact_type[i] <- "pLoF"
      if(sv$SVTYPE == "INV") {
        impact_counts$INV_pLoF <- impact_counts$INV_pLoF + 1
      } else {
        impact_counts$TRA_pLoF <- impact_counts$TRA_pLoF + 1
      }
      
    } else if (sv$SVTYPE == "INS") {
      sv_annotated$gene_impact_type[i] <- "pLoF"
      impact_counts$INS_pLoF <- impact_counts$INS_pLoF + 1
    }
  }
}

# Show detailed impact classification
cat("\nDetailed impact classification:\n")
for(name in names(impact_counts)) {
  if(impact_counts[[name]] > 0) {
    cat(sprintf("  %s: %d\n", name, impact_counts[[name]]))
  }
}

# Filter for gene-disruptive SVs
gene_disruptive_svs <- sv_annotated[
  !is.na(sv_annotated$gene_impact_type) & 
  sv_annotated$gene_impact_type %in% c("pLoF", "CG", "IED"), 
]

cat("\nGene-disruptive SVs (protein-coding genes only):", nrow(gene_disruptive_svs), "\n")
cat("\nImpact type summary:\n")
print(table(gene_disruptive_svs$gene_impact_type))

cat("\nSV types in gene-disruptive set:\n")
print(table(gene_disruptive_svs$SVTYPE))

# Save results
write.table(sv_annotated, "results/chordoma_sv_with_functional_impact_maf5_protein_coding_CORRECTED.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_disruptive_svs, "results/chordoma_gene_disruptive_svs_final_maf5_protein_coding_CORRECTED.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Step 5: Generate visualizations ---------------------------------------------
cat("\n=== STEP 5: Creating Visualizations ===\n")

# 1. SV type distribution - overall
sv_type_df <- as.data.frame(table(sv_bedpe$SVTYPE))
colnames(sv_type_df) <- c("SV_Type", "Count")

p1 <- ggplot(sv_type_df, aes(x = reorder(SV_Type, -Count), y = Count, fill = SV_Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Distribution of All SV Types in Chordoma (MAF ≤ 5%)", 
       subtitle = "After SVTYPE correction",
       x = "SV Type", y = "Count") +
  theme(legend.position = "none") +
  geom_text(aes(label = Count), vjust = -0.5, size = 4) +
  scale_fill_brewer(palette = "Set2")

ggsave("plots/chordoma_sv_type_distribution_maf5_all_CORRECTED.pdf", p1, width = 8, height = 6)

# 2. SV type distribution - gene-overlapping only
gene_sv_type_df <- as.data.frame(gene_sv_type_dist)
colnames(gene_sv_type_df) <- c("SV_Type", "Count")

p1b <- ggplot(gene_sv_type_df, aes(x = reorder(SV_Type, -Count), y = Count, fill = SV_Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "SV Types Overlapping Protein-Coding Genes (MAF ≤ 5%)", 
       subtitle = "After SVTYPE correction",
       x = "SV Type", y = "Count") +
  theme(legend.position = "none") +
  geom_text(aes(label = Count), vjust = -0.5, size = 4) +
  scale_fill_brewer(palette = "Set2")

ggsave("plots/chordoma_sv_type_distribution_maf5_gene_overlapping_CORRECTED.pdf", p1b, width = 8, height = 6)

# 3. Gene impact distribution
if (nrow(gene_disruptive_svs) > 0) {
  impact_df <- as.data.frame(table(gene_disruptive_svs$gene_impact_type))
  colnames(impact_df) <- c("Impact_Type", "Count")
  
  p2 <- ggplot(impact_df, aes(x = Impact_Type, y = Count, fill = Impact_Type)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Gene Impact Types (Protein-Coding Genes, MAF ≤ 5%)", 
         subtitle = "After SVTYPE correction",
         x = "Impact Type", y = "Count") +
    scale_fill_manual(values = c("pLoF" = "#e41a1c", "CG" = "#377eb8", "IED" = "#4daf4a")) +
    theme(legend.position = "none") +
    geom_text(aes(label = Count), vjust = -0.5, size = 4)
  
  ggsave("plots/chordoma_gene_impact_distribution_maf5_protein_coding_CORRECTED.pdf", p2, width = 6, height = 6)
}

# 4. Comparison plot - all SVs vs gene-overlapping SVs
comparison_data <- data.frame(
  Category = rep(c("All SVs", "Gene-overlapping SVs"), each = length(sv_summary)),
  SV_Type = rep(names(sv_summary), 2),
  Count = c(as.numeric(sv_summary), 
            sapply(names(sv_summary), function(x) ifelse(x %in% names(gene_sv_type_dist), gene_sv_type_dist[x], 0)))
)

p3 <- ggplot(comparison_data, aes(x = SV_Type, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "SV Distribution: All vs Gene-Overlapping", 
       subtitle = "After SVTYPE correction",
       x = "SV Type", y = "Count") +
  scale_fill_manual(values = c("All SVs" = "#999999", "Gene-overlapping SVs" = "#E69F00")) +
  theme(legend.position = "top") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5, size = 3)

ggsave("plots/chordoma_sv_comparison_all_vs_gene_overlapping_CORRECTED.pdf", p3, width = 10, height = 6)

# Step 6: Generate final report -----------------------------------------------
cat("\n=== STEP 6: Generating Final Report ===\n")

# Calculate statistics
all_affected_genes <- unique(unlist(strsplit(
  gene_disruptive_svs$affected_genes[!is.na(gene_disruptive_svs$affected_genes)], ";"
)))
gene_freq <- sort(table(unlist(strsplit(gene_disruptive_svs$affected_genes, ";"))), 
                  decreasing = TRUE)

# Create report
sink("results/chordoma_sv_analysis_summary_maf5_protein_coding_CORRECTED.txt")
cat("CHORDOMA STRUCTURAL VARIANT ANALYSIS REPORT (CORRECTED)\n")
cat("=======================================================\n")
cat("Analysis Date:", Sys.Date(), "\n")
cat("Script: chordoma_sv_analysis_maf5_protein_coding.R (CORRECTED VERSION)\n")
cat("VCF File:", basename(vcf_file), "\n")
cat("Gene Set: Protein-coding genes only\n")
cat("Reference Genome: hg19 (GRCh37)\n\n")

cat("IMPORTANT NOTE\n")
cat("--------------\n")
cat("The simple_SVTYPE_classification function had a bug that misclassified\n")
cat("many SVs (DEL, DUP, INV) as INS. This has been corrected based on the\n")
cat("original VCF INFO field.\n\n")

cat("FILTERING APPLIED\n")
cat("-----------------\n")
cat("MAF filtering: ≤ 5%\n")
cat("Gene type: Protein-coding only\n\n")

cat("OVERALL STATISTICS\n")
cat("------------------\n")
cat("Total SVs after MAF filtering:", nrow(sv_annotated), "\n")
cat("\nSV Type Distribution (All SVs - CORRECTED):\n")
print(table(sv_annotated$SVTYPE))
cat("\nPercentages:\n")
for(type in names(sv_summary)) {
  cat(sprintf("  %s: %d (%.1f%%)\n", type, sv_summary[type], 
              100 * sv_summary[type] / sum(sv_summary)))
}

cat("\n\nGENE ANNOTATION RESULTS\n")
cat("-----------------------\n")
cat("Total protein-coding genes used:", nrow(gene_bed_simple), "\n")
cat("SVs overlapping protein-coding genes:", svs_with_genes, 
    "(", round(100 * svs_with_genes / nrow(sv_annotated), 1), "%)\n")
cat("Potential gene fusions:", nrow(gene_fusions), "\n")

cat("\n\nSV TYPE DISTRIBUTION FOR GENE-OVERLAPPING SVs\n")
cat("----------------------------------------------\n")
print(gene_sv_type_dist)
cat("\nPercentages of gene-overlapping SVs:\n")
for(type in names(gene_sv_type_dist)) {
  cat(sprintf("  %s: %d (%.1f%%)\n", type, gene_sv_type_dist[type], 
              100 * gene_sv_type_dist[type] / sum(gene_sv_type_dist)))
}

cat("\n\nFUNCTIONAL IMPACT ANALYSIS\n")
cat("--------------------------\n")
cat("Gene-disruptive SVs (pLoF/CG/IED):", nrow(gene_disruptive_svs), "\n")
cat("\nImpact Type Distribution:\n")
print(table(gene_disruptive_svs$gene_impact_type))

cat("\nDetailed SV type to impact mapping:\n")
# Create cross-tabulation
impact_by_type <- table(gene_disruptive_svs$SVTYPE, gene_disruptive_svs$gene_impact_type)
print(impact_by_type)

cat("\nTotal unique protein-coding genes affected:", length(all_affected_genes), "\n")

cat("\n\nTOP 20 AFFECTED PROTEIN-CODING GENES\n")
cat("------------------------------------\n")
if (length(gene_freq) >= 20) {
  print(head(gene_freq, 20))
} else {
  print(gene_freq)
}


sink()

# Save workspace
save.image("chordoma_sv_analysis_maf5_protein_coding_CORRECTED.RData")

cat("\n=== ANALYSIS COMPLETE! ===\n")
cat("\nOutput files:\n")
cat("- Results: results/\n")
list.files("results", pattern = "CORRECTED")
cat("\n- Plots: plots/\n")
list.files("plots", pattern = "CORRECTED")
cat("\n- Workspace: chordoma_sv_analysis_maf5_protein_coding_CORRECTED.RData\n")
cat("\nSummary report: results/chordoma_sv_analysis_summary_maf5_protein_coding_CORRECTED.txt\n")

# Print key statistics to console
cat("\n=== KEY STATISTICS (CORRECTED) ===\n")
cat("Total SVs: ", nrow(sv_annotated), "\n")
cat("SVs overlapping genes: ", svs_with_genes, " (", round(100 * svs_with_genes / nrow(sv_annotated), 1), "%)\n")
cat("\nSV types in gene-overlapping set:\n")
for(type in names(gene_sv_type_dist)) {
  cat(sprintf("  %s: %d\n", type, gene_sv_type_dist[type]))
}
cat("\nFunctional impacts:\n")
print(table(gene_disruptive_svs$gene_impact_type))
