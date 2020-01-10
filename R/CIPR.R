#' Cluster Identity Predictor
#'
#' @description CIPR (Cluster Identity PRedictor) is a pipeline that helps annotating
#' unknown single cell clusters in single cell RNA sequencing (scRNAseq experiments).
#' This function scores unknown cluster gene expression signatures against known reference
#' datasets using user-selected analytical approaches to facilitate scRNAseq analysis.
#'
#' @param input_dat Data frame containing normalized log2-transformed gene
#' expression values per cluster OR a table of differentially expressed genes
#' per cluster
#'
#' @param comp_method Method to use for identity score calculations. It accepts
#' one of the following: "logfc_dot_product" (default), "logfc_spearman",
#' "logfc_pearson", "all_genes_spearman", "all_genes_pearson"
#'
#' @param reference Reference data frame containing gene expression data from
#' known cell types. It accepts one of the following: "immgen" (default),
#' "mmrnaseq", "blueprint", "hpca", "dice", "hema", "hsrnaseq", "custom"
#'
#' @param select_ref_subsets The names of cell subsets to be included in the
#' analysis. For using the entire reference dataset use "all", or
#' provide a character vector of cell types of interest. Defaults to "all"
#'
#' @param custom_reference A data frame containing custom reference. There must
#' be a column named 'gene' and other columns contain normalized gene expression
#' data from known samples. Defaults to NULL
#'
#' @param custom_ref_annot A data frame containing custom reference metadata.
#' This is optional to get more informative results from CIPR. The data
#' frame must contain columns named 'short_name' (must match column names in
#' custom reference), 'long_name' (human readable names for reference samples),
#' 'description' (details such as positive and negative sorting markers),
#' 'reference_cell_type' (e.g. T cell, B cell, NK)
#'
#' @param keep_top_var Top n percent of highly variant reference genes to
#' include in the analysis. It accepts a numeric value smaller than or equal
#' to 100 (default). The value of 100 results in keeping all the genes in the
#' reference dataset)
#'
#' @param plot_ind Logical value. Set it to TRUE to plot identity scores for
#' each cluster. Defaults to FALSE
#'
#' @param plot_top Logical value. set it to TRUE to plot top scoring reference
#' cell types for each cluster. Defaults ot TRUE.
#'
#' @param top_num A numeric value determining how many top scoring reference
#' cell types will be plotted for each cluster. Defaults to 5.
#'
#' @param save_png Logical value. Set it to TRUE if you would like to export png
#' images of the results. Defaults to FALSE
#'
#' @param global_plot_obj Logical value. Set it to TRUE if you would like to
#' keep the plots as an object in the global environment. This can be useful
#' for accessing and manipulating the graphs. Defaults to TRUE.
#'
#' @param global_results_obj Logical value. Set it to TRUE if you would like to
#' keep the analysis results as a global object. Defaults to TRUE.
#'
#' @param ... arguments to pass to theme() (for graph manipulation)
#'
#' @return Graphical outputs and/or data frames of identity scores calculated
#' for each cluster in the input data.
#'
#' @examples
#'
#' # Example of using CIPR in conjunction with Seurat
#' library(Seurat)
#' allmarkers <- FindAllMarkers(seurat_object)
#' avgexp <- AverageExpression(seurat_object)
#'
#' # Using built-in immgen as reference and logfc dot product method
#' CIPR(input_dat = allmarkers,
#' comp_method = "logfc_dot_product",
#' reference="immgen",
#' keep_top_var = 100,
#' global_results_obj = T,
#' plot_top = T)
#'
#' #' # Using built-in immgen as reference and all genes spearman method
#'
#' CIPR(input_dat = avgexp,
#' comp_method = "all_genes_spearman",
#' reference="immgen",
#' keep_top_var = 100,
#' global_results_obj = T,
#' plot_top = T)
#'
#'
#' # Using built-in dice reference and logFC spearman method
#' # Variance threshold of top 50%
#'
#' CIPR(input_dat = allmarkers,
#' comp_method = "logfc_spearman",
#' reference="dice",
#' keep_top_var = 50,
#' global_results_obj = T,
#' plot_top = T)
#'
#'
#' # Using a custom reference
#'
#' CIPR(input_dat = allmarkers,
#' comp_method = "logfc_dot_product",
#' reference="custom",
#' custom_ref_dat_path = custom_ref_df,
#' custom_ref_annot_path = custom_annot_df,
#' keep_top_var = 100,
#' global_results_obj = T,
#' plot_top = T)
#'
#'
#' # Using a blueprint-encode reference and limiting the analysis
#' # to "Pericytes", "Skeletal muscle", "Smooth muscle"
#'
#' CIPR(input_dat = allmarkers,
#' comp_method = "logfc_dot_product",
#' reference="blueprint-encode",
#' select_ref_subsets = c("Pericytes", "Skeletal muscle", "Smooth muscle")
#' keep_top_var = 100,
#' global_results_obj = T,
#' plot_top = T)
#'
#' # Using built in example data (logFC signatures per cluster)
#' CIPR(input_dat = example_logfc_data,
#' comp_method = "logfc_dot_product",
#' reference="blueprint-encode",
#' select_ref_subsets = c("Pericytes", "Skeletal muscle", "Smooth muscle")
#' keep_top_var = 100,
#' global_results_obj = T,
#' plot_top = T)
#'
#'
#' # Using built in example data (average expression)
#' CIPR(input_dat = example_avgexp_data,
#' comp_method = "all_spearman",
#' reference="immgen",
#' select_ref_subsets = "all",
#' keep_top_var = 100,
#' global_results_obj = T,
#' plot_top = T)
#'



CIPR <- function(input_dat,
                 comp_method = "logfc_dot_product",
                 reference = NULL,
                 select_ref_subsets = "all",
                 custom_reference = NULL,
                 custom_ref_annot = NULL,
                 keep_top_var = 100,
                 plot_ind = F,
                 plot_top = T,
                 top_num = 5,
                 save_png = F,
                 global_plot_obj = T,
                 global_results_obj = T,
                 update_ref = T,
                 ...
){



  suppressMessages({
    require(ggpubr)
    require(gtools)
    require(tibble)
    require(dplyr)
  })



  ######################### Prepare input_dat   #########################

  message("Preparing input data")

  if(grepl("logfc", comp_method)){

    # Define column names to allow flexibility in case and close matches in column names
    gene_column <- grep("gene", colnames(input_dat), ignore.case = T, value = T)
    logFC_column <- grep("logfc", colnames(input_dat), ignore.case = T, value = T)
    cluster_column <- grep("cluster", colnames(input_dat), ignore.case = T, value = T)


    if(length(c(gene_column, logFC_column, cluster_column)) != 3) stop("Check column names of the input data. Data frame must have columns named as 'gene', 'logfc', and 'cluster'")

    # Convert gene symbols to lower case letters to allow mouse-vs-human comparisons
    input_dat[,gene_column] <- tolower(input_dat[,gene_column])

    # input_dat <- input_dat[!duplicated(input_dat[,gene_column]),]



  } else {

    gene_column <- grep("gene", colnames(input_dat), ignore.case = T, value = T)

    input_dat[,gene_column] <- tolower(input_dat[,gene_column])

    input_dat <- input_dat[!duplicated(input_dat[,gene_column]),]

    if(length(gene_column) != 1) stop("Check column names of the input data. Data frame must have one column named as 'gene'")

  }




  ######################### Prepare ref_dat   #########################

  # if(update_ref == T){  # WORK ON CACHING TO SPEED UP THE PIPELINE

  message("Preparing reference data")

    if(reference == "immgen"){

    message("Reading ImmGen reference data")

    # # Read reference dataset
    # load("data/immgen_expr.rda")
    # ref_dat <- get("immgen_expr")
    # rm(immgen_expr)
    #
    # # Read immgen annotation file for explanations of cell types
    # load("data/immgen_samples.rda")
    # ref_annot <- get("immgen_samples")
    # rm(immgen_samples)

      ref_dat <- immgen_expr
      ref_annot <- immgen_samples


  } else if(reference == "mmrnaseq"){

    message("Reading MmRNAseq reference data")

    ref_dat <- mmrnaseq_expr
    ref_annot <- mmrnaseq_samples

  } else if(reference == "blueprint"){

    message("Reading Blueprint-ENCODE reference data")

    ref_dat <- blueprint_expr
    ref_annot <- blueprint_samples

  } else if(reference == "hpca"){

    message("Reading HCPA reference data")

    ref_dat <- hpca_expr
    ref_annot <- hpca_samples

  } else if(reference == "dice"){

    message("Reading DICE reference data")

    ref_dat <- dice_expr
    ref_annot <- dice_samples

  } else if(reference == "hema"){

    message("Reading hema reference data")

    ref_dat <- hema_expr
    ref_annot <- hema_samples

  } else if(reference == "hsrnaseq"){

    message("Reading HsRNAseq reference data")

    ref_dat <- hsrnaseq_expr
    ref_annot <- hsrnaseq_samples

  } else if(reference == "custom"){

    message("Reading custom reference data")

    # Read reference dataset
    ref_dat <- custom_reference

    # Read immgen annotation file for explanations of cell types
    ref_annot <- custom_ref_annot


  } else {

  stop(message('"reference" argument accepts of the following: "mmrnaseq",  "blueprint", "hpca", "dice", "hema", "hsrnaseq", "custom"'))

  }


  ref_gene_column <- grep("^gene", colnames(ref_dat), ignore.case = T, value = T)

  if(length(ref_gene_column) != 1) stop("Check column names of the input data. Data frame must have one column named as 'gene'")

  ref_dat[, ref_gene_column] <- tolower(ref_dat[, ref_gene_column])





  # Select relevant subsets from the reference

  if(select_ref_subsets != "all"){

    message("Subsetting reference data")

    sel_positions <- which(ref_annot[, "reference_cell_type"] %in% select_ref_subsets)

    if(length(sel_positions) == 0) stop("Selected subset is not present in the reference dataset. Please check spelling.")

    select_ref_subsets <- as.character(ref_annot[sel_positions, "short_name"])

    ref_dat <- ref_dat[, c(ref_gene_column, select_ref_subsets)]

  }






  # Apply quantile filtering

  if(keep_top_var != 100){


  message("Applying variance filtering")

  var_vec <- apply(ref_dat[, ! colnames(ref_dat) %in% ref_gene_column], 1, var, na.rm=T)

  keep_var <- quantile(var_vec, probs = 1-keep_top_var/100, na.rm = T)

  keep_genes <- var_vec >= keep_var

  ref_dat <- ref_dat[keep_genes, ]

  }


  if(grepl("logfc", comp_method)){


    # Calculate row means for each gene (mean expression across the reference cell types)
    gene_avg <- rowMeans(ref_dat[, !colnames(ref_dat) %in% ref_gene_column])

    # Log scale data
    reference_ratio <- sweep(ref_dat[,!colnames(ref_dat) %in% ref_gene_column], 1, FUN="-", gene_avg)

    # Combine gene names and the log fold change in one data frame
    ref_dat <- cbind(tolower(ref_dat[, ref_gene_column]), reference_ratio)

    colnames(ref_dat)[1] <- ref_gene_column

  }

  # } else if(exists("ref_dat") & exists("ref_annot") & update_ref == F){
  #
  #   ref_dat <- get("ref_dat", envir = .GlobalEnv)
  #   ref_annot <- get("ref_annot", envir = .GlobalEnv)
  #
  # }



    ######################## Define clusters ###############################

    if(grepl("logfc", comp_method)){

      clusters <- gtools::mixedsort(
        levels(
          as.factor(
            pull(input_dat, grep("cluster", x = colnames(input_dat),
                                 ignore.case = T, value = T)
            )
          )
        )
      )

    } else {

      clusters <- gtools::mixedsort(
        levels(
          as.factor(
            colnames(input_dat)[!grepl("gene", colnames(input_dat),
                                       ignore.case = T)]
          )
        )
      )}




    ######################### Compare input_dat against ref_dat #############################

    message("Analyzing cluster signatures")


    if(comp_method == "logfc_dot_product"){    ####################################################################################

      # Initiate a master data frame to store the results
      master_df <- data.frame()

      # Iterate over clusters to calculate a distinct identity score for each reference cell type
      for (i in clusters) {

        # Increment the progress bar, and update the detail text.
        # message(paste("Analyzing cluster", i))

        # Subset on the cluster in iteration
        sel_clst <- input_dat %>%
          filter(!!rlang::sym(cluster_column) == i) %>%
          select(c(!!sym(gene_column), !!sym(logFC_column)))


        # Merge SCseq cluster log FC value with immgen log FC for shared genes
        merged <- merge(sel_clst, ref_dat, by.x = gene_column, by.y = ref_gene_column)

        if(dim(merged)[1] < 2) next

        # Calculate a scoring matrix by multiplying log changes of clusters and immgen cells
        reference_scoring <- data.frame(apply(merged[,3:dim(merged)[2]],2,function(x){x*merged[,2]}), check.names = FALSE)

        # Calculate the aggregate score of each immgen cell type by adding
        score_sum <- colSums(reference_scoring)

        # Store identity scores in a data frame
        df <- data.frame(identity_score = score_sum)

        df <- rownames_to_column(df, var="reference_id")


        df <- left_join(df, ref_annot, by=c("reference_id" = "short_name"))


        # Store cluster information in a column
        df$cluster <- i

        # Add confidence-of-prediction calculations here and append to the df
        # Calculate the mean and standard deviation of the aggregate scores per reference cell type
        mean_score_sum <- mean(df$identity_score)
        score_sum_sd <- sd(df$identity_score)

        # Calculate the distance of the identity score from population mean (how many std devs apart?)
        df$z_score <- (df$identity_score - mean_score_sum)/score_sum_sd

        # Calculate the proportion of the genes changing in the same direction between unknown cluster and reference cell type
        df$percent_pos_correlation <- {

          ngenes <- dim(reference_scoring)[1]

          pos_corr_vector <- numeric()

          for(i in 1:dim(reference_scoring)[2]){

            # Calculate number of genes positively correlated (upregulated or downregulated in both unk cluster and reference)
            pos_cor <- ( sum(reference_scoring[, i] > 0) / ngenes ) * 100

            pos_corr_vector <- c(pos_corr_vector, pos_cor)

          } #close for loop

          pos_corr_vector

        } # close expression


        # Add calculation results under the master data frame to have a composite results file
        master_df <- rbind(master_df,df)



      } # close for loop that iterates over clusters1

    } else if(comp_method == "logfc_spearman" | comp_method == "logfc_pearson"){  ########################################################

      # Initiate master data frame to store results
      master_df <- data.frame()


      # Iterate analysis for each cluster. The loop below will calculate a distinct correlation
      # coefficient for each cluster-reference cell pairs
      for (i in clusters) {


        trim_dat <- input_dat %>%
          filter(!!rlang::sym(cluster_column) == i)

        dat_genes <- trim_dat[gene_column] %>% pull() %>% as.character
        ref_genes <- ref_dat[ref_gene_column] %>% pull() %>% as.character

        common_genes <- intersect(dat_genes, ref_genes)


        trim_dat <- trim_dat %>%
          filter(!!rlang::sym(gene_column) %in% common_genes) %>%
          arrange(!!rlang::sym(gene_column)) %>%
          select(- !!rlang::sym(gene_column))


        trim_ref <- ref_dat %>%
          filter(!!rlang::sym(ref_gene_column) %in% common_genes) %>%
          arrange(!!rlang::sym(ref_gene_column)) %>%
          select(- !!rlang::sym(ref_gene_column))


        # Calculate correlation between the the cluster (single column in trimmed input data) and each of the
        # reference cell subsets (columns of the trimmed reference data)
        cor_df <- cor(trim_dat[logFC_column], trim_ref, method = gsub("logfc_", "", comp_method))

        # Store results in a data frame
        df <- data.frame(identity_score = cor_df[1,])

        df <- rownames_to_column(df, var="reference_id")

        # Combine results with reference annotations
        if(reference != "custom"){

          df <- left_join(df, ref_annot, by=c("reference_id" = "short_name"))





        } else if (reference == "custom" & !is.null(custom_ref_annot_path)){

          df <- left_join(df, ref_annot, by=c("reference_id" = "short_name"))


        } else if(reference == "custom" & is.null(custom_ref_annot_path)){

          # Fill in with reminder if annotation file is not updated
          df$reference_cell_type <- rep("Upload annotation file", dim(ref_dat)[2]-1)
          df$short_name <- colnames(ref_dat)[!colnames(ref_dat) %in% ref_gene_column]
          df$long_name <- rep("Upload annotation file", dim(ref_dat)[2]-1)
          df$description <- rep("Upload annotation file", dim(ref_dat)[2]-1)

        }



        # Store cluster information in a column
        df$cluster <- i

        # Add confidence-of-prediction calculations here and append to the df
        # Calculate the mean and standard deviation of the aggregate scores per reference cell type
        mean_cor_coeff <- mean(df$identity_score)
        cor_coeff_sd <- sd(df$identity_score)

        # Calculate the distance of the identity score from population mean (how many std devs apart?)
        df$z_score <- (df$identity_score - mean_cor_coeff)/cor_coeff_sd

        # Add all the results to the master data frame
        master_df <- rbind(master_df, df)


      } # close for loop that iterates over clusters

    } else if(comp_method == "all_genes_spearman" | comp_method == "all_genes_pearson"){  ################################################



      dat_genes <- input_dat[gene_column] %>% pull() %>% as.character
      ref_genes <- ref_dat[ref_gene_column] %>% pull() %>% as.character

      common_genes <- intersect(dat_genes, ref_genes)

      trim_dat <- input_dat %>%
        filter(!!rlang::sym(gene_column) %in% common_genes) %>%
        arrange(!!rlang::sym(gene_column)) %>%
        select_(.dots= paste0("-", gene_column))

      trim_ref <- ref_dat %>%
        filter(!!rlang::sym(ref_gene_column) %in% common_genes) %>%
        arrange(!!rlang::sym(ref_gene_column)) %>%
        select_(.dots=paste0("-", ref_gene_column))

      clusters <- colnames(trim_dat)


      master_df <- data.frame()

      comp_method <- gsub("all_genes_", "", comp_method)

      for (i in clusters) {


        cor_df <- cor(trim_dat[i], trim_ref, method = comp_method)


        df <- data.frame(identity_score = cor_df[1,])

        df <- rownames_to_column(df, var="reference_id")



        if(reference != "custom"){

          df <- left_join(df, ref_annot, by=c("reference_id" = "short_name"))





        } else if (reference == "custom" & !is.null(custom_ref_annot_path)){

          df <- left_join(df, ref_annot, by=c("reference_id" = "short_name"))


        } else if(reference == "custom" & is.null(custom_ref_annot_path)){

          df$reference_cell_type <- rep("Upload annotation file", dim(ref_dat)[2]-1)
          df$short_name <- colnames(ref_dat)[!colnames(ref_dat) %in% ref_gene_column]
          df$long_name <- rep("Upload annotation file", dim(ref_dat)[2]-1)
          df$description <- rep("Upload annotation file", dim(ref_dat)[2]-1)

        }




        df$cluster <- i

        # Add confidence-of-prediction calculations here and append to the df
        # Calculate the mean and standard deviation of the aggregate scores per reference cell type
        mean_cor_coeff <- mean(df$identity_score)
        cor_coeff_sd <- sd(df$identity_score)

        # Calculate the distance of the identity score from population mean (how many std devs apart?)
        df$z_score <- (df$identity_score - mean_cor_coeff)/cor_coeff_sd


        master_df <- rbind(master_df,df)

      }
    }








    if(global_results_obj == T) CIPR_all_results <<- master_df


    #prep individual plots
    if(plot_ind == T){

      ind_clu_plots <- list()

      for (i in clusters) {


        # Extract results calculated for individual clusters
        df_plot <- master_df %>%
          filter(cluster == i)

        # Calculate mean and sd deviation for adding confidence bands to graphs
        score_mean <- mean(df_plot$identity_score)
        score_sd <- sd(df_plot$identity_score)


        plotname <- paste("cluster", i, sep="")

        # Plot identity scores per cluster per reference cell type and add confidence bands
        ind_clu_plots[[plotname]] <- ggdotplot(df_plot, x = "reference_id", y="identity_score",
                                               fill = "reference_cell_type", xlab=F, ylab="Reference identity score",
                                               font.y = c(14, "bold", "black"), size=1, x.text.angle=90,
                                               title = paste("Cluster:",i), font.title = c(15, "bold.italic"),
                                               font.legend = c(15, "plain", "black"))+
          theme(axis.text.x = element_text(size=10, vjust=0.5, hjust=1))+
          geom_hline(yintercept=score_mean)+
          annotate("rect", xmin = 1, xmax = length(df_plot$reference_id),
                   ymin = score_mean-score_sd, ymax = score_mean+score_sd,
                   fill = "gray50", alpha = .1)+
          annotate("rect", xmin = 1, xmax = length(df_plot$reference_id),
                   ymin = score_mean-2*score_sd, ymax = score_mean+2*score_sd,
                   fill = "gray50", alpha = .1)+
          theme(...)


      }

      if(global_plot_obj == T) ind_clu_plots <<-ind_clu_plots


      if(save_png == T) {
        ggexport(filename = "CIPR_individual_clusters.png", plotlist = ind_clu_plots, ncol = 1, width = 1800, height = 360 * length(clusters))
      }
      else {
        print(ggarrange(plotlist = ind_clu_plots, ncol = 1, common.legend = T))
      }

    }

    ################################################################################################################################
    # Prepare top5 summary plots
    # This plot will show the 5 highest scoring reference cell types for each cluster.

    if(plot_top == T){

      # Extract top5 hits from the reuslts
      top_df <- master_df %>%
        group_by(cluster) %>%    #cluster
        top_n(top_num, wt = identity_score) %>%
        arrange(cluster, desc(identity_score))

      # Index variable helps keeping the results for clusters separate and helps ordered outputs
      top_df$index <- 1:nrow(top_df)


      # Order clusters levels for ordered plotting
      ordered_cluster_levels <- gtools::mixedsort(levels(as.factor(top_df$cluster)))


      top_df$cluster <- factor(top_df$cluster, levels = ordered_cluster_levels)



      # Extract relevant columns
      top_df <- select(top_df, cluster,
                       reference_cell_type,
                       reference_id,
                       long_name,
                       description,
                       identity_score,
                       index, everything())

      if(global_results_obj == T) CIPR_top_results <<- top_df

      p <- ggdotplot(top_df, x="index", y="identity_score",
                     fill = "cluster", size=1, x.text.angle=90,
                     font.legend = c(15, "plain", "black")) +
        scale_x_discrete(labels=top_df$reference_id)+
        theme(axis.text.x = element_text(vjust=0.5, hjust=1))+
        theme(...)

      if(global_plot_obj == T) top_plots <<- p

      if(save_png == T) {

        ggexport(p, filename = "CIPR_top_hits.png", ncol = 1, width = 150 * length(clusters), height = 300)

      } else {

        print(p)

      }

    }

  } # close function
