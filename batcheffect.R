library(proBatch)
library(gridExtra)

QCeffectPlot<-function(quantdata,sample_annotation,batch_col, order_col, factor_cols, correctmethod){
if(is.null(quantdata)){return()}
colnames(sample_annotation)<-c("FullRunName",colnames(sample_annotation)[-1])
QCeffectoutput<-list()
QCeffect.ls<-list()
log_transformed_matrix <- log_transform_dm(quantdata,log_base = 2, offset = 1)
color_list <- sample_annotation_to_colors(sample_annotation,
                                          factor_columns = c(batch_col,factor_cols),
                                          numeric_columns = order_col)
QCeffect.ls[["sample mean plot"]]<-plot_sample_mean(log_transformed_matrix, sample_annotation, order_col = order_col,
                 batch_col = batch_col, color_by_batch = TRUE,
                 color_scheme = color_list[[batch_col]])
log_transformed_long <- matrix_to_long(log_transformed_matrix)

QCeffect.ls[["sample box plot"]]<-plot_boxplot(log_transformed_long, sample_annotation,
             batch_col = batch_col, color_scheme = color_list[[batch_col]])

median_normalized_matrix = normalize_data_dm(log_transformed_matrix,
                                             normalize_func = 'median')
par(mar=c(5, 5, 2, 5))
plot_hierarchical_clustering(median_normalized_matrix,
                             sample_annotation = sample_annotation,
                             color_list = color_list,
                             factors_to_plot = factor_cols,
                             distance = 'euclidean', agglomeration = 'complete',
                             label_samples = FALSE)
QCeffectoutput[["clustering"]]<-recordPlot()
QCeffectoutput[["heatmap"]]<-plot_heatmap_diagnostic(median_normalized_matrix, sample_annotation,
                             factors_to_plot = factor_cols,
                             cluster_cols = TRUE,
                             color_list = color_list,
                             show_rownames = FALSE, show_colnames = FALSE)
QCeffect.ls[["PCA"]]<-plot_PCA(median_normalized_matrix, sample_annotation,
                      color_by = batch_col, plot_title = batch_col)
batch_corrected_matrix <-  long_to_matrix(correct_batch_effects_df(df_long = matrix_to_long(median_normalized_matrix),
                                                   sample_annotation = sample_annotation,
                                                   discrete_func = correctmethod,
                                                   continuous_func = 'loess_regression',
                                                   batch_col = batch_col,
                                                   order_col=order_col
                                                   ))

QCeffectoutput[["plot"]]<-QCeffect.ls
QCeffectoutput[["batchcorrectedmatrix"]]<-unlog_dm(batch_corrected_matrix, log_base = 2, offset = 1)
return(QCeffectoutput)
}

