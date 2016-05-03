
################################################################################
########################     Impulse DE package     ############################
################################################################################

### Version:  0.99.0
### Date:     2016-04-06
### Author:   Jil Sander

################################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                             Functions to call                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Annotation preparation    +++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Prepares the annotation table for internal use
annotation_preparation <- function(data_annotation, data_table = NULL,
  colname_time = NULL, colname_condition = NULL, control_timecourse = FALSE,
  control_name = NULL, case_name = NULL){

   if (is.null(data_table)) {
      stop("ERROR: table of data values must be specified")
   }
   if ((nrow(data_annotation) != ncol(data_table)) ||
      (FALSE %in% (sort(rownames(data_annotation)) ==
            sort(colnames(data_table))))) {
                 stop("ERROR: column names of data table must be the same as
                    row names of annotation table!")
   }
   if (is.null(colname_time)) {
      stop("ERROR: name of the column containing the timepoints must be
           specified as character!")
   }
   if (is.null(colname_condition)) {
      stop("ERROR: name of the column containing the condition(s) must be
           specified as character!")
   }
   if (is.numeric(data_annotation[,colname_time]) == FALSE) {
      stop("ERROR: time variable must be numeric and is not allowed to
           contain characters!")
   }
   if (control_timecourse == FALSE &
      length(summary(as.factor(data_annotation[,colname_condition]))) >= 2) {
        stop("ERROR: if you have only one time course do not provide more than
             one condition!")
   }
   if (control_timecourse == TRUE & (length(summary(as.factor(data_annotation[,
        colname_condition]))) > 2) & is.null(case_name)) {
            stop("ERROR: please specify case and control names of interest
                 since you provide more than three conditions!")
   }
   annot <- data_annotation[,c(colname_time,colname_condition)]
   annot <- data_annotation[colnames(data_table),]
   colnames(annot) <- c("Time","Condition")
   annot$Condition <- as.character(annot$Condition)

   if (length(summary(as.factor(data_annotation[,colname_condition]))) > 2) {
      print(paste("Case condition: ", case_name, sep  = ""))
   } else {
      print(paste("Case condition: ",  data_annotation[!(data_annotation[,
        colname_condition] %in% control_name),colname_condition][1], sep = ""))
   }
   if (control_timecourse == TRUE) {
      print(paste("Control condition: ", control_name, sep  = ""))
   }
   if (control_timecourse == TRUE & (length(summary(as.factor(data_annotation[,
        colname_condition]))) > 2)) {
             annot <- annot[annot$Condition %in% c(control_name, case_name),]
   }
   return(annot)
}

################################################################################

#' Impulse model value prediction
#'
#' Calculates impulse model values for given timepoints and predicted
#' impulse parameters.
#' @aliases calc_impulse
#' @param theta numerical vector of impulse parameters with the order
#' beta, h0, h1, h2, t1, t2.
#' @param timepoints numercial vector of time point(s).
#' @return The predicted impulse model values for the given time point(s).
#' @seealso \code{\link{impulse_DE}}, \code{\link{plot_impulse}}.
#' @author Jil Sander
#' @references Chechik, G. and Koller, D. (2009) Timing of Gene Expression
#' Responses to Envi-ronmental Changes. J. Comput. Biol., 16, 279-290.
#' @examples
#' #' theta vector in the order beta, h0, h1, h2, t1, t2
#' theta <- c(9.9, 14.7, 17.0, 16.9, -0.1, 37.0)
#' #' time points
#' timepoints <- c(0, 2, 4, 6, 8, 18, 24, 32, 48, 72)
#' #' calculate impulse values
#' impulse_values <- calc_impulse(theta, timepoints)
#' @export
calc_impulse <- function(theta,timepoints){
  beta1 = theta[1]
  h0 = theta[2]
  h1 = theta[3]
  h2 = theta[4]
  t1 = theta[5]
  t2 = theta[6]
  res = NULL
   for (i in 1:length(timepoints)) {
    res[i] = (1/h1) * (h0 + (h1 - h0) * (1/(1 + exp(-beta1*(timepoints[i] - t1))))) *
        (h2 + (h1 - h2) * (1/(1 + exp(beta1*(timepoints[i] - t2)))))
  }
  res = unlist(res)
  return(res)
}

################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++     Calculation of input for Optimization    ++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### specifies the function of squared errors which is minimized to fit
### parameters
two_impulses <- function(theta,x_vec,y_vec){
  beta1 = theta[1]
  h0 = theta[2]
  h1 = theta[3]
  h2 = theta[4]
  t1 = theta[5]
  t2 = theta[6]
   f = sum((unlist(lapply(x_vec, function(x) {(1/h1) * (h0 + (h1 - h0) *
        (1/(1 + exp(-beta1*(x - t1))))) * (h2 + (h1 - h2) *
        (1/(1 + exp(beta1*(x - t2)))))})) - y_vec)^2)
   return(f)
}

################################################################################

## compile simple functions to make them quicker
calc_impulse_comp <- cmpfun(calc_impulse)
two_impulses_comp <- cmpfun(two_impulses)


################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Clustering    +++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### performs 2-step clustering on expression data or background data
### 1. Step: clustering genes based on correlation to compare the shape
### 2. Step: split the clusters from step 1 into finer clusters based on the
###          eucledian distance to compare the magnitude
cluster_genes_for_impulse <- function(data_table, data_annotation,
    control_timecourse = FALSE, control_name = NULL, plot_clusters = FALSE,
    no_of_clusters = NULL, n_genes = NULL){

  #' @import amap
  corr_thres = 0.85
  eucl_thres = 1.8

  ### assume only one timecourse as default and perform only 1 run
  results = list()
  runs = 1
  label = "case"

  ### perform 3 runs if control timecourse data is present
  if (control_timecourse == TRUE & is.null(control_name) == FALSE) {
     runs = 3
     label = c("combined","case","control")
  }
  ### clustering for different runs
  for (c_runs in 1:runs) {
    if (c_runs == 1) {                # combined or only case if no control
      print(paste("Clustering of ",label[c_runs]," data set", sep = ""))
      dat_pre = data_table
      dat_pre = dat_pre[which(apply(dat_pre,1,function(x){100*
              sd(x)/(mean(x))}) >= 2.5),]
      dat = t(apply(dat_pre,1,function(y){unlist(lapply(as.list(unique(
          as.numeric(as.character(data_annotation$Time)))),function(x){
          mean(y[data_annotation$Time == x])}))}))
      colnames(dat) <- unique(as.numeric(as.character(data_annotation$Time)))
    } else if(c_runs == 2) {         # case
      print(paste("Clustering of ",label[c_runs]," data set", sep = ""))
      dat_pre = data_table[,!(data_annotation$Condition %in% control_name)]
      dat_pre = dat_pre[which(apply(dat_pre,1,function(x){100*
              sd(x)/(mean(x))}) >= 2.5),]
      dat = t(apply(dat_pre,1,function(y){unlist(lapply(as.list(unique(
          as.numeric(as.character(data_annotation[!(data_annotation$Condition
          %in% control_name),"Time"])))),function(x){mean(
          y[data_annotation[!(data_annotation$Condition %in% control_name),
          "Time"] == x])}))}))
      colnames(dat) <- unique(as.numeric(as.character(data_annotation[
          !(data_annotation$Condition %in% control_name),"Time"])))
    } else {                        # control
      print(paste("Clustering of ",label[c_runs]," data set", sep = ""))
      dat_pre = data_table[,data_annotation$Condition %in% control_name]
      dat_pre = dat_pre[which(apply(dat_pre,1,function(x){100*
              sd(x)/(mean(x))}) >= 2.5),]
      dat = t(apply(dat_pre,1,function(y){unlist(lapply(as.list(unique(
          as.numeric(as.character(data_annotation[data_annotation$Condition
          %in% control_name,"Time"])))),function(x){mean(
          y[data_annotation[data_annotation$Condition %in% control_name,
         "Time"] == x])}))}))
      colnames(dat) <- unique(as.numeric(as.character(data_annotation[
          (data_annotation$Condition %in% control_name),"Time"])))
    }
    print(paste("- ", (nrow(data_table) - nrow(dat_pre)),
            " genes were excluded due to very low variation", sep = ""))

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   #~~~~~~~~~~~~~   Step 1: clustering based on correlation   ~~~~~~~~~~~~~~~~~#
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    max_n_pre_clus = 25   # restrict maximum number of pre-clusters to 25
    ### z-transformation of data for clustering
    dat_ztrans <- t(scale(t(dat)))

    ### ---> clustering of expression data
    if (is.null(no_of_clusters)) {
      file_part = "genes"
      n_clus = 0
      for (i in 1:max_n_pre_clus) {

        kmeans_clus <- Kmeans(dat_ztrans,i,200,method = "correlation")

    ### accept clustering if mean correlation between genes of the cluster and
    ### centroid of the cluster is higher than 'corr_thres' (default: 0.75)
        for (j in 1:i) {
          tmp <- dat_ztrans[kmeans_clus$cluster %in% j,]
          indt1 <- NULL
          if (i > 1 & sort(summary(as.factor(kmeans_clus$cluster)))[1]  < 10) {
               indt1 <- "I"
               n_clus = i - 1
               kmeans_clus <- kmeans_clus_pre
               break
          }
          if (is.null(nrow(tmp))) {    # if cluster contains only 1 gene
             average_pearson <- cor(tmp, kmeans_clus$centers[j,])
          } else {     # if cluster contains > 1 genes
             average_pearson <- mean(apply(tmp,1,function(x){cor(x,
                    kmeans_clus$centers[j,])}))
          }
        ### if one cluster does not fufill the correlation condition leave loop,
        ### increase numbers of potential clusters and try again
          if (average_pearson < corr_thres) {
            break
          }
          ### accept clustering and leave loop if all clusters fufill condition
          if (j == i) {
            n_clus = i
            break
          }
        }
        if (is.null(indt1) == FALSE) {
           break
        }
        if (n_clus != 0) {
          break
        ### if maximum number of allowed clusters is reached accept this number
        } else if ((n_clus == 0) & (i == max_n_pre_clus)) {
          n_clus = i
          break
        }
        kmeans_clus_pre = kmeans_clus
      }
      print(paste("--- Number of correlation-based pre-clusters: ",
            n_clus,sep = ""))

    ### clustering of background data based on clusters from expression data
    } else {
      file_part = "bg"
      n_clus <- max(round(no_of_clusters[[c_runs*2 - 1]]*nrow(data_table)/
            n_genes),1)
      kmeans_clus <- Kmeans(dat_ztrans,n_clus,200,method = "correlation")
      n_fine_clust_rand <- round(sample( no_of_clusters[[c_runs*2]] ,n_clus,
            replace = TRUE)*nrow(data_table)/n_genes)
      n_fine_clust_rand[n_fine_clust_rand < 1] = 1
    }
    kmeans_pre_clus_corr <- kmeans_clus
    n_pre_clus <- n_clus
    kmeans_clus_final <- kmeans_pre_clus_corr$cluster
    print(paste("------ Number of genes in pre-clusters: ",
      paste(paste("C",names(summary(as.factor(kmeans_clus_final))),": ",
      summary(as.factor(kmeans_clus_final)), sep = ""),collapse = ", "),sep = ""))


   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   #~~~~~~~~~   Step 2: clustering based on eucledian distance   ~~~~~~~~~~~~~~#
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


    n_fine_clusters <- NULL
    cluster_means <-  NULL
    # restrict maximum number of fine clusters to 5 per pre-cluster
    max_n_fine_clus <- 10
    kmeans_clus <- NULL

     tempo2 <- do.call(c,as.data.frame(dat))
     tempo3 <- (tempo2 - median(tempo2))/mad(tempo2)
     dat_median_ztrans <- matrix(tempo3, nrow(dat), ncol(dat))
     rownames(dat_median_ztrans) = rownames(dat)
     colnames(dat_median_ztrans) = colnames(dat)

    ### cluster each pre-cluster into fine clusters
    for (cl in 1:n_pre_clus) {


     ### ---> clustering of expression data
       if (is.null(no_of_clusters)) {
        if (length(kmeans_pre_clus_corr$cluster[kmeans_pre_clus_corr$cluster %in%
                cl]) >  max_n_fine_clus) {
            n_clus = 0

            for (i in 1:max_n_fine_clus) {

              kmeans_clus <- kmeans(dat_median_ztrans[names(
                  kmeans_pre_clus_corr$cluster[kmeans_pre_clus_corr$cluster
                  %in% cl]),],i,200)
            ### accept clustering if mean eucledian distance between genes of
            ### the cluster and centroid of the cluster is lower than
            ### 'eucl_thres' (default: 5)
              for (j in 1:i) {
              indt2 <- NULL
              tmp <- dat_median_ztrans[names(kmeans_pre_clus_corr$cluster[
                  kmeans_pre_clus_corr$cluster %in% cl]),][kmeans_clus$cluster
                  %in% j,]
               # if smallest cluster would have less than 5 genes take previous
               # number of clusters
              if (i > 1 & sort(summary(as.factor(kmeans_clus$cluster)))[1] < 5) {
                  indt2 = "I"
                  n_clus = i - 1
                  kmeans_clus <- kmeans_clus_pre
                  break
              }

                if (is.null(nrow(tmp))) {    # if cluster contains only 1 gene
                   average_eucledian <- 0
                } else {  # if cluster contains > 1 genes
                   average_eucledian <- sd(apply(tmp,1,function(x){dist(rbind(x,
                        kmeans_clus$centers[j,]))}))
                }

                ### if one cluster does not fufill the eucledian condition leave
                ### loop,increase numbers of potential clusters and try again
                if (average_eucledian > eucl_thres) {
                  break
                }
                ### accept clustering and leave loop if all clusters fulfill
                ### condition
                if (j == i) {
                  n_clus = i
                  break
                }
              }
              if (is.null(indt2) == FALSE) {
                 break
              }
              if (n_clus != 0) {
                break
              ### if maximum number of allowed clusters is reached accept this
              ### number
              } else if ((n_clus == 0) & (i == max_n_fine_clus)) {
                n_clus = i
                break
              }
              kmeans_clus_pre = kmeans_clus
            }
            n_fine_clusters <- c(n_fine_clusters,n_clus)
          } else {
            n_clus = 1
            n_fine_clusters <- c(n_fine_clusters,n_clus)
            kmeans_clus$cluster <- rep(1,length(kmeans_pre_clus_corr$cluster[
                kmeans_pre_clus_corr$cluster %in% cl]))
            names(kmeans_clus$cluster) <- names(kmeans_pre_clus_corr$cluster[
                kmeans_pre_clus_corr$cluster %in% cl])
          }
        ### clustering of background data based on clusters from expression data
        } else {
          if (floor(length(kmeans_pre_clus_corr$cluster[
                kmeans_pre_clus_corr$cluster %in% cl]) / n_fine_clust_rand[cl])
                >= 2) {
            n_clus <- n_fine_clust_rand[cl]
            kmeans_clus <- kmeans(dat_median_ztrans[names(
                  kmeans_pre_clus_corr$cluster[kmeans_pre_clus_corr$cluster
                  %in% cl]),], n_clus,200)
          } else if (n_fine_clust_rand[cl] > 1 & floor(length(
                kmeans_pre_clus_corr$cluster[kmeans_pre_clus_corr$cluster %in%
                cl]) / (n_fine_clust_rand[cl] - 1)) >= 2) {
            n_clus <- n_fine_clust_rand[cl] - 1
            kmeans_clus <- kmeans(dat_median_ztrans[names(
                kmeans_pre_clus_corr$cluster[kmeans_pre_clus_corr$cluster
                %in% cl]),], n_clus,200)
          } else {
            n_clus = 1
            kmeans_clus$cluster <- rep(1,length(
                kmeans_pre_clus_corr$cluster[kmeans_pre_clus_corr$cluster
                %in% cl]))
            names(kmeans_clus$cluster) <- names(kmeans_pre_clus_corr$cluster[
                kmeans_pre_clus_corr$cluster %in% cl])
          }
          n_fine_clusters <- c(n_fine_clusters,n_clus)
        }

#      ### calculate mean values for each final cluster
      for (i in 1:n_clus) {
         tmp <- dat_pre[names(kmeans_clus$cluster)[kmeans_clus$cluster %in% i],]
      if (is.null(dim(tmp))) {       # if cluster contains only 1 gene
          cluster_means <- rbind(cluster_means,tmp)
        } else { # if cluster contains > 1 genes
          cluster_means <- rbind(cluster_means,colMeans(tmp))
        }
      }
      ### overwrite pre-clusters by fine clusters to have cluster information
      ### for each gene
      if (cl == 1) {
        kmeans_clus_final[names(kmeans_pre_clus_corr$cluster[
            kmeans_pre_clus_corr$cluster %in% cl])] <-  kmeans_clus$cluster
      } else {
        kmeans_clus_final[names(kmeans_pre_clus_corr$cluster[
            kmeans_pre_clus_corr$cluster %in% cl])] <-  kmeans_clus$cluster +
            sum(n_fine_clusters[1:(cl - 1)])
      }
    }

    if (file_part != "bg") {
      write.table(kmeans_clus_final,paste("kmeans_clus_final_",label[c_runs],
            "_",file_part,".txt",sep = ""),sep = "\t",quote = FALSE,
            row.names = TRUE, col.names = FALSE)
    }
    rownames(cluster_means)  <- 1:sum(n_fine_clusters)
    n_clus = sum(n_fine_clusters)

    print(paste("--- Final number of clusters: ",n_clus,sep = ""))
    print(paste("------ Number of genes in final clusters: ",
      paste(paste("C",names(summary(as.factor(kmeans_clus_final))),": ",
      summary(as.factor(kmeans_clus_final)), sep = ""),collapse = ", "),sep = ""))


   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   #~~~~~~~~~~~~~~~~   Plot clusters if plot_clusters == TRUE ~~~~~~~~~~~~~~~~~#
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    if (plot_clusters == TRUE & file_part != "bg") {
      pdf(paste("clusters_",label[c_runs],"_",file_part,".pdf",sep = ""),
          height = 10.0,width = 16.0)

      for (i in 1:n_clus) {

        ### 1 or 2 clusters
        if ((i == 1) & (n_clus <= 2)) {
            par(mfrow = c(1,2))
        } else if ((i == 1) & (n_clus <= 4)) {
            par(mfrow = c(2,2))
        ### print all clusters on one page if number of clusters is <= 12
        } else if ((i == 1) & (n_clus <= 12)) {
          par(mfrow = c(3,4))

        ### if number of clusters is > 12 plot always 20 clusters on one page
        } else if (i == 1 || ((i - 1) %% 20 == 0)) {
          par(mfrow = c(4,5))
        }
        ### expression values of genes in cluster
        tmp <- dat_ztrans[kmeans_clus_final %in% i,]
        # if cluster contains only 1 gene
        if (length(kmeans_clus_final[kmeans_clus_final %in% i]) == 1) {
          n_g = 1
        } else {   # if cluster contains > 1 genes
          n_g = nrow(tmp)
        }

        ### order gene expression data according to timepoints
        for (j in 1:n_g) {
          # if cluster contains only 1 gene
          if (length(kmeans_clus_final[kmeans_clus_final %in% i]) == 1) {
            tmp_plot = tmp
          } else {       # if cluster contains > 1 genes
            tmp_plot = tmp[j,]
          }
          ### plot clusters
          if (j == 1) {                     # if cluster contains only 1 gene
            plot(tmp_plot,type = "l",xlab = "Time points", ylab = "Z-score",
                    main = paste("Cluster ",i,sep = ""),ylim = c(min(tmp),max(tmp)),
                    xaxt = "n")
          axis(1, at = 1:length(colnames(tmp)), labels = colnames(tmp), las = 3)
          } else {   # if cluster contains > 1 genes
            points(colnames(tmp),tmp_plot,type = "l")
          }
        }
      }
    dev.off()
    }
    results[[c_runs*4 - 3]] <- kmeans_clus_final
    results[[c_runs*4 - 2]] <- cluster_means
    results[[c_runs*4 - 1]] <- n_pre_clus
    results[[c_runs*4]] <- n_fine_clusters
    names(results)[c((c_runs*4 - 3):(c_runs*4))] <- paste(c("kmeans_clus",
         "cluster_means", "pre_clus","fine_clus"),label[c_runs],sep = "_")
  }

  return(results)
}


################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Impulse model fit    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### fits impulse model to a timecourse dataset
### genes are suppposed to be in rows, samples in columns
### timepoints need to be specified as numbers (integers or float)
### function is splitted into 3 parts
###### 1. fitting impulse model to a single gene
###### 2. fitting impulse model to matrix of genes, which calls 1.
###### 3. prepare data and fit model by calling 2.
### This separation was chosen to be able to use 'apply' to replace for-loops

impulse_fit <- function(data_input, data_annotation, n_iter = 100,
    control_timecourse = FALSE, control_name = NULL, cluster_results = NULL,
    start_values = NULL, fit_backg = FALSE, n_proc = 4, ...){

  ### check dataset for consistency
  if (ncol(data_annotation) != 2 || FALSE %in% (colnames(data_annotation) ==
      c("Time","Condition"))) {
            stop("Please use function 'annotation_preparation' to prepare
                 the annotation table")
  }
  if (control_timecourse == TRUE & is.null(control_name)) {
     stop("Please specify the name of the control in the column 'Condition'")
  }
  if (control_timecourse == TRUE & (is.null(start_values) == FALSE)) {
    start_values <- start_values[c(1,3,5)]
  }
  if (control_timecourse == FALSE & (is.null(control_name) == FALSE)) {
     start_values <- start_values[c(1,3)]
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~  Fit impulse model to gene ~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  ### fits an impulse model to expression values of a single gene

  impulse_fit_gene_wise <- function(vec, timepoints, NPARAM, n_iters = 100,
      fit_to_clusters = FALSE, start_val = NULL, fit_bg = FALSE, ...){

      mat1 = c(0,abs(vec[2:length(vec)] - vec[1:(length(vec) - 1)]) *
            vec[2:length(vec)])
      mn_beta = 0
      mx_beta = 10
      middle_ind = which(timepoints > min(timepoints) & timepoints <
            max(timepoints))
      min_ind = which(timepoints == min(timepoints))[1]
      max_ind = which(timepoints == max(timepoints))[1]
      mx_tm = max(timepoints)
      mn_tm = min(timepoints)
      tmp_ind = which(mat1 == max(mat1))[1]
      peaktime = timepoints[tmp_ind]
      peak = min(which(timepoints == peaktime))
      beta1 = abs(log(abs((vec[peak] - vec[min_ind])/(peaktime -
            timepoints[min_ind]))));

      ### set beta to 1 if calculcated value is infinite
      if (is.finite(beta1) == FALSE || is.na(beta1)) { beta1 = 1 }

      ### define start value theta based on the gene's data
      orig_theta = c(beta1,
         vec[min_ind],                    # h0
         vec[peak],                       # h1
         vec[max_ind],                    # h2
         (peaktime - timepoints[min_ind])/2,            # t1
         peaktime + (timepoints[max_ind] - peaktime)/2) # t2
      names(orig_theta) = c("beta","h0","h1","h2","t1","t2")
      theta = orig_theta

      ### replace theta estimates of 0 by small value, because optimization
      ### function has problems with 0 as input
      rand_num <- runif(length(which(theta == 0)),0,0.0001)
      theta[names(which(theta == 0))] <- rand_num
      names(theta) = c("beta","h0","h1","h2","t1","t2")

      ### ---> fit model to clusters
      if (fit_to_clusters == TRUE & is.null(start_val)) {

          if (fit_bg == FALSE) {
            thetas <- cbind(as.numeric(theta),apply(matrix(rep(runif((n_iters -
                    1)*NPARAM)),NPARAM,(n_iters - 1)),2,
              function(x){c(mn_beta + (mx_beta - mn_beta)*x[1],
                            min(vec) + (max(vec) - min(vec))*x[2:4],
                            mn_tm + (mx_tm - mn_tm)*x[5],
                            (mn_tm + (mx_tm - mn_tm)*x[5]) + (mx_tm - (mn_tm +
                            (mx_tm - mn_tm)*x[5]))*x[6])}))
           } else {
              thetas <- apply(matrix(rep(runif((n_iters)*NPARAM)),NPARAM,
                    (n_iters)),2,
              function(x){c(mn_beta + (mx_beta - mn_beta)*x[1],
                            min(vec) + (max(vec) - min(vec))*x[2:4],
                            mn_tm + (mx_tm - mn_tm)*x[5],
                            (mn_tm + (mx_tm - mn_tm)*x[5]) + (mx_tm - (mn_tm +
                            (mx_tm - mn_tm)*x[5]))*x[6])})
          }
        tmm1 <- system.time({
        fmin_outs <- cbind(
        # fitting via quasi-Newton method (optim, "BFGS")
           apply(thetas[,1:(floor(n_iters/2))],2,function(x){unlist(optim(x,
             two_impulses, x_vec = timepoints,
             y_vec = vec,method = "BFGS")[c("par","value")])}),
        # fitting via PORT routines (nlminb)
           apply(thetas[, c(1,(floor(n_iters/2) + 1):n_iters)],2,function(x){
             unlist(nlminb(x, two_impulses, x_vec = timepoints,
             y_vec = vec)[c("par","objective")])})
           )
        })


        ### return results of 3 best fits, which will serve as start values for
        ### the fits to the genes
        pvec_and_SSE = as.vector(fmin_outs[,order(fmin_outs["value",])][,1:3])

      ### ---> fit model to genes
      } else if (fit_to_clusters == FALSE & is.null(start_val) == FALSE) {
        tmm2 <- system.time({

      ### use first 3 best hits from fit to clusters plus initial guess 'theta'
      ### as start values
        if (fit_bg == TRUE) {
          toll <- matrix(start_val,7,3)[1:NPARAM,]
        } else {
          toll <- cbind(theta,matrix(start_val,7,3)[1:NPARAM,][,1:3])
        }
        fmin_outs <- cbind(
        # fitting via quasi-Newton method (optim, "BFGS")
        apply(toll,2,
          function(x){unlist(optim(x, two_impulses, x_vec = timepoints,
          y_vec = as.numeric(vec),method = "BFGS" )[c("par","value")])}),
        # again fitting via PORT routines (nlminb)
        apply(toll,2,
          function(x){unlist(nlminb(x, two_impulses, x_vec = timepoints,
          y_vec = as.numeric(vec) )[c("par","objective")])})
        )
        })
        if (is.null(dim(fmin_outs[ ,fmin_outs["value",] == min(fmin_outs["value",
            ])])) == TRUE) {
          pvec_and_SSE = fmin_outs[ ,fmin_outs["value",] ==
                min(fmin_outs["value",])]
       } else {
          ### if two or more randomization have the same minimum Sum of
          ### Squared Errors (SSE), choose the first one
          pvec_and_SSE = fmin_outs[ ,fmin_outs["value",] ==
                min(fmin_outs["value",])][,1]
       }
      }
     return(pvec_and_SSE)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~  Fit impulse model to matrix ~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  ### fits an impulse model to all genes of a dataset
  impulse_fit_matrix <- function(data_tab, timepoints, n_it = 100,
        ctrl_tc = FALSE, ctrl = NULL, fit_to_clus = FALSE, start_val = NULL,
        fit_bg = FALSE, n_process = 4, ...){

    NPARAM = 6
    mat1 = NULL
    res <- list()
    for (j in 2:ncol(data_tab)) {
        mat1 = cbind(mat1, abs(data_tab[,j] - data_tab[,j - 1]) * data_tab[,j])
    }

        mc <- min(detectCores() - 1, n_process)
     if (fit_bg == FALSE && nrow(data_tab) > max(2*mc,10)) {

    ## do not split if less than 10 genes
        ind_list = list()
        bord = floor(nrow(data_tab)/mc)
        for (i in 1:mc) {
           ind_list[[i]] <- ((i - 1)*bord + 1):(i*bord)
        }
        if (mc*bord < nrow(data_tab)) {
             ind_list[[mc]] <-  c(ind_list[[mc]],(mc*bord + 1):nrow(data_tab))
        }

        cl <- makeCluster(mc, outfile = "clus_out2.txt")
        my.env <- new.env()
        assign("data_tab", data_tab, envir = my.env)
        assign("calc_impulse_comp", calc_impulse_comp, envir = my.env)
        assign("n_it", calc_impulse_comp, envir = my.env)
        assign("cluster_genes_for_impulse", cluster_genes_for_impulse,
             envir = my.env)
        assign("impulse_fit", impulse_fit, envir = my.env)
        assign("two_impulses", two_impulses, envir = my.env)
        assign("impulse_fit_gene_wise",impulse_fit_gene_wise, envir = my.env)
        assign("impulse_fit_matrix",impulse_fit_matrix, envir = my.env)
        assign("timepoints", timepoints, envir = my.env)
        assign("fit_to_clus", fit_to_clus, envir = my.env)
        assign("start_val", start_val, envir = my.env)
        assign("fit_bg", fit_bg, envir = my.env)
        assign("NPARAM", NPARAM, envir = my.env)

        clusterExport(cl = cl, varlist = c("data_tab",
         "calc_impulse_comp","n_it","cluster_genes_for_impulse","impulse_fit",
         "two_impulses", "impulse_fit_gene_wise",
         "impulse_fit_matrix", "timepoints","fit_to_clus","start_val",
         "fit_bg","NPARAM"), envir = my.env)

        ### fit impulse model to each gene of matrix and get impulse parameters
        resi <- clusterApply(cl, ind_list, function(z){t(apply(data_tab[z,],1,
            function(x){impulse_fit_gene_wise(x,timepoints,
            NPARAM,n_it,fit_to_clus, start_val, fit_bg)}))})
        resmat = do.call(rbind,resi)
        resmat <- resmat[rownames(data_tab),]
        stopCluster(cl)
      } else {
      ### fit impulse model to each gene of matrix and get impulse parameters
        resmat <- t(apply(data_tab,1,function(x){impulse_fit_gene_wise(x,
            timepoints, NPARAM,n_it,fit_to_clus, start_val, fit_bg)}))
      }


    ### use obtained impulse parameters to calculate impulse fit values
    ### ---> if fit to genes
    if (fit_to_clus == FALSE) {
     colnames(resmat) <- c("beta","h0","h1","h2","t1","t2","SSE")
     if (nrow(resmat) == 1) {      # if matrix contains only one gene
       resmat2 <- as.data.frame(t(calc_impulse_comp(resmat[,1:NPARAM],
          unique(sort(timepoints)))),
          row.names = rownames(resmat))
       colnames(resmat2) = unique(sort(timepoints))
     } else {   # if matrix contains > 1 genes
        resmat2 <- t(apply(resmat[,1:NPARAM],1,function(x){calc_impulse_comp(x,
                unique(sort(timepoints)))}))
        colnames(resmat2) = unique(sort(timepoints))
     }

    ### ---> if fit to clusters
    ### more complex because results from 3 best fits need to be saved and later
    ### used for the fit to the genes
    } else {
     colnames(resmat) <- c(paste(c("beta","h0","h1","h2","t1","t2","SSE"),1,
          sep = "_"), paste(c("beta","h0","h1","h2","t1","t2","SSE"),2,sep = "_"),
          paste(c("beta","h0","h1","h2","t1","t2","SSE"),3,sep = "_"))
     resmat2 <- t(apply(resmat,1,function(x){apply(matrix(x,7,3)[1:NPARAM,],2,
          function(y){calc_impulse_comp(y,unique(sort(timepoints)))})}))
     colnames(resmat2) <- c(paste(unique(sort(timepoints)),1, sep = "_"),
          paste(unique(sort(timepoints)),2, sep = "_"),
          paste(unique(sort(timepoints)),3, sep = "_"))
    }
    res[[1]] <- resmat
    res[[2]] <- resmat2
    names(res) <- c("impulse_parameters","impulse_fits")
    return(res)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~  Prepare data for impulse model fit ~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  ### assume only one timecourse as default and perform only 1 run
  g_names = rownames(data_input)
  runs = 1
  label = "case"
  results = list()
  if (control_timecourse == TRUE) {
    runs = 3
    label = c("combined","case","control")
  }

  ### ---> fit to genes
  if (is.null(cluster_results) == FALSE & is.null(start_values) == FALSE) {
    file_ext = "genes"
    c1 = "genes"

    ### perform 3 runs if control timecourse data is present
    if (control_timecourse == TRUE) {
      ### generate a data list containing the combined, case and control data if
      ### control timecourse is present
      dat1 = as.data.frame(data_input)[names(cluster_results[[1]]),]
      dat1_n = as.data.frame(data_input)[!(rownames(data_input) %in%
            rownames(dat1)),]
      dat1_n = dat1_n[,colnames(dat1)]
      dat2 = as.data.frame(data_input)[,!(data_annotation$Condition %in%
            control_name)]
      dat2 = dat2[names(cluster_results[[5]]),]
      dat2_n = as.data.frame(data_input)[!(rownames(data_input) %in%
            rownames(dat2)),]
      dat2_n = dat2_n[,colnames(dat2)]
      dat3 = as.data.frame(data_input)[,(data_annotation$Condition %in%
            control_name)]
      dat3 = dat3[names(cluster_results[[9]]),]
      dat3_n = as.data.frame(data_input)[!(rownames(data_input) %in%
            rownames(dat3)),]
      dat3_n = dat3_n[,colnames(dat3)]
      dat_ind <- c("dat1_n","dat2_n", "dat3_n")
      data_input <- list(dat1, dat2, dat3)

    } else if (control_timecourse == FALSE) {
      dat1 = as.data.frame(data_input)[names(cluster_results[[1]]),]
      dat1_n = as.data.frame(data_input)[!(rownames(data_input) %in%
            rownames(dat1)),]
      dat_ind <- c("dat1_n")
      data_input <- list(dat1)
    }


  ### ---> fit to clusters
  }  else if (is.null(cluster_results) & is.null(start_values)) {
    file_ext = "clusters"
    c1 = "cluster"
    ### input are the means of the gene expression values of the clusters
    ### 2 is combined, 6 is case and 10 is control
    ### if runs == 1 or runs == 2, then 6 and 10 or just 10 are empty
    data_input <- data_input[c(2,6,10)[1:runs]]
  }

  ### fitting for different runs
  for (c_runs in 1:runs) {
    imp_res <- NULL

    ### ---> fit to genes
    if (is.null(cluster_results) == FALSE & is.null(start_values) == FALSE) {

      ### split genes into the clusters and fit impulse model to the
      ### genes of a cluster
      cluster_res_list <- split(data_input[[c_runs]],
            cluster_results[[(c_runs - 1)*4 + 1]])
      ind <- split(1:max(cluster_results[[(c_runs - 1)*4 + 1]]),
            1:max(cluster_results[[(c_runs - 1)*4 + 1]]))

      imp_res_list <- lapply(ind,
          function(x){impulse_fit_matrix(cluster_res_list[[x]],
          as.numeric(as.character(data_annotation[colnames(data_input[[c_runs]])
            ,"Time"])), n_it = n_iter, ctrl_tc = control_timecourse,
            ctrl = control_name, fit_to_clus = FALSE,
            start_val = start_values[[c_runs]][x,],
            fit_bg = fit_backg, n_process = n_proc)})

       if (nrow(get(dat_ind[c_runs])) != 0) {
        tump1  <- do.call(rbind,lapply(1:(length(imp_res_list)),
            function(x) imp_res_list[[x]][[1]]))
        tump2  <- do.call(rbind,lapply(1:(length(imp_res_list)),
            function(x) imp_res_list[[x]][[2]]))

        tmpp1 <- t(apply(get(dat_ind[c_runs]),1,function(x){rep(mean(x),
            ncol(tump2))}))
        colnames(tmpp1) <- colnames(tump2)
        tump2 <- rbind(tump2, tmpp1)
        tmpp2 <-  cbind( matrix(NA,nrow(get(dat_ind[c_runs])),ncol(tump1) - 1),
            apply(cbind(get(dat_ind[c_runs]), tmpp1[,1]),1,
            function(x){sum((x[1:(length(x) - 1)] - x[length(x)])^2)}))
        rownames(tmpp2) <- rownames(get(dat_ind[c_runs]))
        tump1 <- rbind(tump1, tmpp2)

        tempi1 <- tump1
        tempi2 <- tump2
      } else {
         tempi1  <- do.call(rbind,lapply(1:(length(imp_res_list)),
             function(x) imp_res_list[[x]][[1]]))
         tempi2  <- do.call(rbind,lapply(1:(length(imp_res_list)),
             function(x) imp_res_list[[x]][[2]]))
      }
      imp_res$impulse_parameters <- tempi1[g_names,]
      imp_res$impulse_fits <- tempi2[g_names,]

    ### ---> fit to clusters
    } else if (is.null(cluster_results) & is.null(start_values)) {
      ### fit impulse model to the mean of the cluster
      imp_res <-  impulse_fit_matrix(data_input[[c_runs]],
        as.numeric(as.character(data_annotation[colnames(data_input[[c_runs]]),
        "Time"])),
          n_it = n_iter, ctrl_tc = control_timecourse, ctrl = control_name,
                fit_to_clus = TRUE, fit_bg = fit_backg, n_process = n_proc)
    }
    names(imp_res) <- paste(c("impulse_parameters","impulse_fits"),
        label[c_runs], sep = "_")
    results[[names(imp_res)[1]]] <- imp_res[[1]]
    results[[names(imp_res)[2]]] <- imp_res[[2]]
  }

  ### for naming of the output files
  if (fit_backg == TRUE) { c2 = "bg" } else {c2 = "data"}
  return(results)
}

################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Plot impulse fits    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### plots impulse fits to timecourse data and also the control data if present

#' Plot impulse model fits
#'
#' Plots impulse model fits for the specified gene IDs. In the case of two
#' time courses, the fits for the combined, case and control data are plotted.
#' @aliases plot_impulse
#' @param gene_IDs character vector of gene names to be plotted; must be part
#' of the \code{rownames} of \code{data_table}
#' @param data_table numeric matrix of expression values; genes should be in
#' rows, samples in columns. Data should be properly normalized and
#' log2-transformed as well as filtered for present or variable genes.
#' @param data_annotation table providing co-variables for the samples including
#' condition and time points. Time points must be numeric numbers.
#' @param imp_fit_genes list of  fitted impulse model values and parameters as
#' produced by \code{impulse_DE} (list element \code{impulse_fit_results});
#' another possibility is to load the saved fitting object
#' \code{load("impulse_fit_genes.RData")}.
#' @param control_timecourse logical indicating whether a control time
#' timecourse is part of the data set (\code{TRUE}) or not (\code{FALSE}).
#' Default is \code{FALSE}.
#' @param control_name character string specifying the name of the control
#' condition in \code{annotation_table}.
#' @param case_name character string specifying the name of the case
#' condition in \code{annotation_table}. Should be set if more than two
#' conditions are present in \code{data_annotation}.
#' @param file_name_part character string to be used as file extention.
#' @param title_line character string to be used as title for each plot.
#' @param sub_line character string to be used as subtitle for each plot.
#' @return PDF-files containing the plots of the impulse model fits for the
#' specified gene IDs.
#' @seealso \code{\link{impulse_DE}}, \code{\link{calc_impulse}}.
#' @author Jil Sander
#' @references Chechik, G. and Koller, D. (2009) Timing of Gene Expression
#' Responses to Envi-ronmental Changes. J. Comput. Biol., 16, 279-290.
#' @examples
#' #' Install package longitudinal and load it
#' library(longitudinal)
#' #' Attach datasets
#' data(tcell)
#' #' check dimension of data matrix of interest
#' dim(tcell.10)
#' #' generate a proper annotation table
#' annot <- as.data.frame(cbind("Time" =
#'      sort(rep(get.time.repeats(tcell.10)$time,10)),
#'      "Condition" = "activated"), stringsAsFactors = FALSE)
#' #' Time columns must be numeric
#' annot$Time <- as.numeric(annot$Time)
#' #' rownames of annotation table must appear in data table
#' rownames(annot) = rownames(tcell.10)
#' #' since genes must be in rows, transpose data matrix using t()
#' #' consider 6 genes for now only
#' genes <- c("SIVA","CD69","ZNFN1A1","IL4R","MAP2K4","JUND")
#' tcell.10.filtered <- t(tcell.10[,genes])
#' #' generate a list object having the form of the output of impulse_DE
#' #' first the parameter fits and SSEs
#' impulse_parameters_case <- matrix(
#'     c(0.6,	18.6,	17.2,	17.4,	5.1,	40.2,	3.5,
#'       0.3,	-464.9,	18.3,	17.3,	-17.2,	35.3,	17.5,
#'       23.2,	18,	18.8,	18.5,	3,	37,	13.2,
#'       NA,	NA,	NA,	NA,	NA,	NA,	3.1,
#'       NA,	NA,	NA,	NA,	NA,	NA,	9.6,
#'       9.5,	17.5,	18.7,	17.5,	8,	48,	46.7),length(genes),7, byrow = TRUE)
#' rownames(impulse_parameters_case) <- genes
#' colnames(impulse_parameters_case) <- c("beta", "h0", "h1", "h2", "t1", "t2", "SSE")
#' #' then the fitted values for the time points
#' impulse_fits_case <- matrix(c(
#'     18.55,	18.43,	18.15,	17.73,	17.43,	17.24,	17.24,	17.24,	17.38,	17.38,
#'     16.22,	17.18,	17.7,	17.97,	18.12,	18.27,	18.26,	18.03,	17.3,	17.28,
#'     18,	18,	18.82,	18.82,	18.82,	18.82,	18.82,	18.82,	18.48,	18.48,
#'     15.93,	15.93,	15.93,	15.93,	15.93,	15.93,	15.93,	15.93,	15.93,	15.93,
#'     17.62,	17.62,	17.62,	17.62,	17.62,	17.62,	17.62,	17.62,	17.62,	17.62,
#'     17.5,	17.5,	17.5,	17.5,	18.18,	18.67,	18.67,	18.67,	17.98,	17.53)
#'     ,length(genes),length(unique(annot$Time)), byrow = TRUE)
#' rownames(impulse_fits_case) <- genes
#' colnames(impulse_fits_case) <- unique(annot$Time)
#' #' finalize list object
#' impulse_fit_genes <- list("impulse_parameters_case" = impulse_parameters_case,
#'                           "impulse_fits_case" = impulse_fits_case)
#' #' Plot expression values
#' plot_impulse(genes, tcell.10.filtered, annot, impulse_fit_genes)
#' @export
plot_impulse <- function(gene_IDs, data_table, data_annotation,imp_fit_genes,
    control_timecourse = FALSE, control_name = NULL, case_name = NULL,
    file_name_part = "", title_line = "", sub_line = ""){

  print("---Plotting genes")

  if (length(grep("[a-zA-Z]",rownames(data_table))) == 0) {
       rownames(data_table) <- paste(rownames(data_table),"G", sep = "_")
       gene_IDs <- paste(gene_IDs, "G", sep = "_")
  }
  print(gene_IDs)

  ### if control timecourse is present split data into case and control data
  if (control_timecourse == TRUE) {
      fmat_case <- data_table[,!(data_annotation$Condition %in% control_name)]
      fmat_ctrl <- data_table[,data_annotation$Condition %in% control_name]
      timep_case <- as.numeric(as.character(data_annotation[colnames(fmat_case),
            "Time"]))
      timep_ctrl <- as.numeric(as.character(data_annotation[colnames(fmat_ctrl),
            "Time"]))
  }
  timep <- as.numeric(as.character(data_annotation[colnames(data_table),"Time"]))
  pdf(paste("impulse_fit_genes_",file_name_part,".pdf", sep = ""),
      height = 6.0,width = 9.0)
   if (length(gene_IDs) == 1) {
      par(mfrow = c(1,1))
   } else if (length(gene_IDs) <= 4) {
      par(mfrow = c(2,2))
   } else if (length(gene_IDs) <= 6) {
      par(mfrow = c(2,3))
   } else {
    par(mfrow = c(3,3))
   }
   x_vec <- seq(0,max(timep),0.1)
   for (i in 1:length(gene_IDs)) {
     ### if there is no control data plot only case time course data and fit
     if (control_timecourse == FALSE) {
         if (TRUE %in% is.na(imp_fit_genes$impulse_parameters_case[gene_IDs[i],])) {
             calc_case = imp_fit_genes$impulse_fits_case[gene_IDs[i],1]
         } else {
             calc_case = calc_impulse_comp(imp_fit_genes$impulse_parameters_case[gene_IDs[i],1:6],x_vec)
         }
         plot(timep,data_table[gene_IDs[i],],col = "blue",pch = 3,xlim = c(0,max(timep)),
            ylim = c(min(c(as.numeric(data_table[gene_IDs[i],]),as.numeric(calc_case))) - 0.5,
               max(c(as.numeric(data_table[gene_IDs[i],]),as.numeric(calc_case))) + 0.5),
            xlab = "Time", ylab = "Impulse fit or log2 expression value",
            main = paste(gene_IDs[i]," ",title_line, sep = ""),sub = sub_line)
         if (TRUE %in% is.na(imp_fit_genes$impulse_parameters_case[gene_IDs[i],])) {
            abline(h = calc_case , col = "blue")
         } else {
           points(x_vec, calc_case, col = "blue", type = "l")
         }
         legend(x = "bottomright",as.character(data_annotation[1,"Condition"]),
                fill = c("blue"), cex = 0.6)

     ### if there is control data
     } else if (control_timecourse == TRUE) {

         if (TRUE %in% is.na(imp_fit_genes$impulse_parameters_case[gene_IDs[i],])) {
            calc_case = imp_fit_genes$impulse_fits_case[gene_IDs[i],1]
            status_case = FALSE
         } else {
            calc_case = calc_impulse_comp(imp_fit_genes$impulse_parameters_case[gene_IDs[i],1:6],x_vec)
            status_case = TRUE
         }
         if (TRUE %in% is.na(imp_fit_genes$impulse_parameters_control[gene_IDs[i],])) {
            calc_ctrl = imp_fit_genes$impulse_fits_control[gene_IDs[i],1]
            status_ctrl = FALSE
         } else {
            calc_ctrl = calc_impulse_comp(imp_fit_genes$impulse_parameters_control[gene_IDs[i],1:6],x_vec)
            status_ctrl = TRUE
         }
         if (TRUE %in% is.na(imp_fit_genes$impulse_parameters_combined[gene_IDs[i],])) {
            calc_comb = imp_fit_genes$impulse_fits_combined[gene_IDs[i],1]
            status_comb = FALSE
         } else {
            calc_comb = calc_impulse_comp(imp_fit_genes$impulse_parameters_combined[gene_IDs[i],1:6],x_vec)
            status_comb = TRUE
         }

         plot(timep_case,fmat_case[gene_IDs[i],],col = "blue",pch = 3,
              xlim = c(0,max(timep)),
              ylim = c(min(c(as.numeric(data_table[gene_IDs[i],]),
               as.numeric(calc_case), as.numeric(calc_ctrl),
               as.numeric(calc_comb))) - 0.5,
               max(c(as.numeric(data_table[gene_IDs[i],]),as.numeric(calc_case),
               as.numeric(calc_ctrl), as.numeric(calc_comb))) + 0.5),
            xlab = "Time", ylab = "Impulse fit or log2 expression value",
            main = paste(gene_IDs[i]," ",title_line,sep = ""),sub = sub_line)

         points(timep_ctrl,fmat_ctrl[gene_IDs[i],],col = "red",pch = 4)

         if (status_case == FALSE) {
           abline(h = calc_case , col = "blue")
         } else {
           points(x_vec,calc_case, col = "blue", type = "l")
         }
         if (status_ctrl == FALSE) {
            abline(h = calc_ctrl , col = "red")
         } else {
           points(x_vec,calc_ctrl, col = "red", type = "l")
         }
         if (status_comb == FALSE) {
            abline(h = calc_comb , col = "grey")
         } else {
           points(x_vec,calc_comb, col = "grey", type = "l")
         }
         legend(x = "bottomright",
            c(as.character(data_annotation$Condition[data_annotation$Condition !=
            control_name][1]),control_name,"combined"),
            fill = c("blue","red","grey"), cex = 0.6)
     }
   }
   dev.off()
}


#################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Background generation   +++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### generates the background for the DE analysis
### can be parallelized onto 'n_proc' nodes

generate_background <- function(data_table, data_annotation, n_iter = 100,
    imp_fit_genes, control_timecourse = FALSE, control_name = NULL,
    no_of_clus = NULL, n_rands = 50000, n_proc = 4){

  ### get number of pre- and fine clusters as input from the genes
  ### also exclude genes for which the model was not fitted due to almost
  ### no variation
  if (control_timecourse == TRUE) {
    no_of_clus = no_of_clus[c(3,4,7,8,11,12)]
    data_table <- data_table[!(is.na(imp_fit_genes[[1]][,"beta"]) |
        is.na(imp_fit_genes[[3]][,"beta"]) |
        is.na(imp_fit_genes[[5]][,"beta"])),]
  } else if (control_timecourse == FALSE & is.null(control_name)) {
    no_of_clus = no_of_clus[c(3,4)]
    data_table <- data_table[!(is.na(imp_fit_genes[[1]][,"beta"])),]
  }
    imp_fit_genes <- lapply(imp_fit_genes,
        function(x){x <- x[rownames(data_table),]})

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~  Fit background ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  impulse_bg <- function(data_tab, data_annot,  n_it = 100, imp_fit,
        ctrl_tk = FALSE, ctrl_name = NULL,
        n_clust = NULL, n_randoms_tot = 50000, no_proc = 4){

    ### divide randomizations by no. of nodes
    n_rand <- round(n_randoms_tot/no_proc)

    background = rep(0,n_rand)
    timep <- as.numeric(as.character(data_annot[colnames(data_tab),"Time"]))
    NPARAM = 6

    ### take 'n_rand' random entries from data table
    i_s <- sample(1:nrow(data_tab),n_rand, replace = TRUE)
    data_tab_rand <- data_tab[i_s,]
    rownames(data_tab_rand) <- paste("bg_",c(1:length(i_s)))

    ### ---> if control data is present
    if (ctrl_tk == TRUE) {
      data_list <- list(data_tab_rand, data_tab_rand[,
          !(data_annot$Condition %in% ctrl_name)],
          data_tab_rand[,data_annot$Condition %in% ctrl_name])

      ### null model: fits to the combined dataset
      y_i_null_s = (imp_fit[[2]][rownames(data_tab),])[i_s,
            as.character(data_annot[colnames(data_list[[1]]),"Time"])]
      calc_impulse_comp_data_tab_rand_case <- (imp_fit[[4]][rownames(data_tab),
            ])[i_s,as.character(data_annot[colnames(data_list[[2]]),"Time"])]

      ### if real control timecourse is present
        calc_impulse_comp_data_tab_rand_ctrl <- (imp_fit[[6]][rownames(data_tab),
            ])[i_s, as.character(data_annot[colnames(data_list[[3]]),"Time"])]
        RESD_1_s = cbind((data_list[[2]] - calc_impulse_comp_data_tab_rand_case),
          (data_list[[3]] - calc_impulse_comp_data_tab_rand_ctrl))
        RESD_1_s = RESD_1_s[,colnames(data_list[[1]])]

      ### draw random resdiuals and add them to the null model
      Y_star_s = y_i_null_s + t(apply(RESD_1_s,1,
          function(x){x[sample(1:ncol(RESD_1_s),
          ncol(data_tab_rand),replace = TRUE)]}))
      colnames(Y_star_s) <- colnames(data_tab_rand)
      rownames(Y_star_s) <- paste("bg_",c(1:length(i_s)))

      ### cluster Y_star_s data by using the clusters from the genes
      tm_c <- system.time({
        clustering_results_background <- cluster_genes_for_impulse(Y_star_s,
            data_annot, ctrl_tk, ctrl_name, plot_clusters = TRUE,
            no_of_clusters = n_clust, n_genes = nrow(data_tab))
      })
      print(paste("Consumed time bg clustering: ",
            round(tm_c["elapsed"]/60,2)," min", sep = ""))

      ### fit impulse model to the clusters of Y_star_s
      tm_f <- system.time({
        impulse_fit_clusters_bg <- impulse_fit(clustering_results_background,
            data_annot, n_it, control_timecourse = ctrl_tk,
            ctrl_name, fit_backg = TRUE)
      })
      print(paste("Consumed time bg clus fit: ",
            round(tm_f["elapsed"]/60,2)," min", sep = ""))

      ### fit impulse model to each Y_star
      tm_fg <- system.time({
        impulse_fit_genes_bg <- impulse_fit(Y_star_s, data_annot, 1, ctrl_tk,
            ctrl_name, clustering_results_background, impulse_fit_clusters_bg,
            fit_backg = TRUE )
      })

      print(paste("Consumed time bg gene fit: ",
            round(tm_fg["elapsed"]/60,2)," min", sep = ""))

      ### calculate background
        RESD_1_s_b = cbind((Y_star_s[rownames(data_list[[1]]),
            colnames(data_list[[2]])] -
            impulse_fit_genes_bg[[4]][rownames(data_list[[1]]),
            as.character(data_annot[colnames(data_list[[2]]),"Time"])]),
            (Y_star_s[rownames(data_list[[1]]),colnames(data_list[[3]])] -
            impulse_fit_genes_bg[[6]][rownames(data_list[[1]]),
            as.character(data_annot[colnames(data_list[[3]]),"Time"])]))

        colnames(RESD_1_s_b) <- c(colnames(data_list[[2]]),colnames(data_list[[3]]))
        RESD_1_s_b = RESD_1_s_b[,colnames(data_list[[1]])]
        RESD_0_s_b = Y_star_s[rownames(data_list[[1]]),
            colnames(data_list[[1]])] -
            impulse_fit_genes_bg[[2]][rownames(data_list[[1]]),
            as.character(data_annot[colnames(data_list[[1]]),"Time"])]

       SS_1_s_b = apply(RESD_1_s_b,1,function(x){sum(x^2)})
       SS_0_s_b = apply(RESD_0_s_b,1,function(x){sum(x^2)})

    ### ---> if no control data is present
    } else {
      y_i_null_s = rowMeans(data_tab_rand)
      calc_impulse_comp_data_tab_rand <- (imp_fit[[2]][rownames(data_tab),])[i_s,
            as.character(data_annot[,"Time"])]
      RESD_1_s = data_tab_rand - calc_impulse_comp_data_tab_rand

      ### draw random resdiuals and add them to the null model
      Y_star_s = t(y_i_null_s + apply(RESD_1_s,1,function(x){x[sample(1:ncol(RESD_1_s),
          ncol(data_tab_rand),replace = TRUE)]}))
      colnames(Y_star_s) <- colnames(data_tab_rand)

      ### cluster Y_star_s data by using the clusters from the genes
      tm_c <- system.time({
         clustering_results_background <- cluster_genes_for_impulse(Y_star_s,
            data_annot, ctrl_tk, ctrl_name, plot_clusters = TRUE,
            no_of_clusters = n_clust, n_genes = nrow(data_tab))
      })
      print(paste("Consumed time bg clustering: ",
                  round(tm_c["elapsed"]/60,2)," min",sep = ""))

      ### fit impulse model to the clusters of Y_star_s
      tm_f <- system.time({
       impulse_fit_clusters_bg <- impulse_fit(clustering_results_background,
            data_annot, n_it, control_timecourse = ctrl_tk, ctrl_name,
            fit_backg = TRUE)
      })
      print(paste("Consumed time bg clust fit: ",
                  round(tm_f["elapsed"]/60,2)," min",sep = ""))

      ### fit impulse model to each Y_star
      tm_fg <- system.time({
        impulse_fit_genes_bg <- impulse_fit(Y_star_s, data_annot, 1, ctrl_tk, ctrl_name,
            clustering_results_background, impulse_fit_clusters_bg, fit_backg = TRUE )
      })
      print(paste("Consumed time bg gene fit: ",
        round(tm_fg["elapsed"]/60,2)," min", sep = ""))

      ### calculate background
       RESD_1_s_b = Y_star_s - impulse_fit_genes_bg[[2]][rownames(Y_star_s),
            as.character(data_annot[colnames(Y_star_s),"Time"])]
       SS_1_s_b = apply(RESD_1_s_b,1,function(x){sum(x^2)})
       RESD_0_s_b = Y_star_s - rowMeans(Y_star_s)
       SS_0_s_b = apply(RESD_0_s_b,1,function(x){sum(x^2)})
    }
    background = (SS_0_s_b - SS_1_s_b) / SS_1_s_b
    return(background)
  }

#   define a new environment, which is used for parallelization
  my.env <- new.env()
  assign("n_rands", n_rands, envir = my.env)
  assign("n_proc", n_proc, envir = my.env)
  assign("data_table", data_table, envir = my.env)
  assign("data_annotation", data_annotation, envir = my.env)
  assign("n_iter",n_iter,envir = my.env)
  assign("control_timecourse", control_timecourse, envir = my.env)
  assign("control_name", control_name, envir = my.env)
  assign("imp_fit_genes", imp_fit_genes, envir = my.env)
  assign("no_of_clus", no_of_clus, envir = my.env)
  assign("calc_impulse_comp", calc_impulse_comp, envir = my.env)
  assign("cluster_genes_for_impulse", cluster_genes_for_impulse, envir = my.env)
  assign("impulse_fit", impulse_fit, envir = my.env)
  assign("two_impulses", two_impulses, envir = my.env)
  environment(impulse_bg) <- my.env

#' @import parallel
#' @import boot
#
#  ### set number of nodes to maximum - 1
  mc <- min(detectCores() - 1, n_proc)
  assign("mc", mc, envir = my.env)
  cl <- makeCluster(mc, outfile = "cluster_out_random.txt")
  clusterExport(cl = cl, varlist = c("impulse_bg","n_rands","n_proc","data_table",
    "data_annotation","n_iter","control_timecourse","control_name",
    "imp_fit_genes","mc", "no_of_clus","calc_impulse_comp",
    "cluster_genes_for_impulse","impulse_fit", "two_impulses"), envir = my.env)
  junk <- clusterEvalQ(cl,library(boot))
  clusterSetRNGStream(cl,123)
  res <- clusterEvalQ(cl,impulse_bg(data_table, data_annotation, n_iter,
     imp_fit_genes, control_timecourse, control_name, no_of_clus, n_rands, mc))
  environment(res) <- my.env
  res2 <- do.call(c,res)
  stopCluster(cl)
  return(res2)
}

################################################################################


DE_analysis <- function(data_table,data_annotation,impulse_fit_results,
    background, control_timecourse = FALSE, control_name = NULL,
    e_type = "Array", Q = 0.01){

  ### if control timecourse is present split data into case and control data
  if (control_timecourse == TRUE) {
      fmat_case <- as.matrix(data_table[,!(data_annotation$Condition %in%
            control_name)])
      fmat_ctrl <- as.matrix(data_table[,data_annotation$Condition %in%
            control_name])
     timep_case <-  as.numeric(as.character(data_annotation[colnames(fmat_case),
            "Time"]))
      timep_ctrl <- as.numeric(as.character(data_annotation[colnames(fmat_ctrl),
            "Time"]))
  }
  timep <- as.numeric(as.character(data_annotation[colnames(data_table),
        "Time"]))

  ### correct background values
  background[background < 0] = 0
  background = background[background < 100]

  if (control_timecourse == TRUE) {
    y_i_null = impulse_fit_results[[2]][rownames(fmat_case),as.character(timep)]
     RESD_1_s = cbind(fmat_case - impulse_fit_results[[4]][rownames(fmat_case),
        as.character(timep_case)], (fmat_ctrl[rownames(fmat_case),] -
        impulse_fit_results[[6]][rownames(fmat_case),as.character(timep_ctrl)]))
     RESD_1_s = RESD_1_s[,colnames(data_table)]

    ### if there is no control data
  } else if (control_timecourse == FALSE) {
      y_i_null = rowMeans(data_table)
      RESD_1_s = data_table - impulse_fit_results[[2]][rownames(data_table),
            as.character(timep)]
  }
  RESD_0_s = data_table - y_i_null
  SS_0_s = apply(RESD_0_s, 1, function(x){sum(x^2)})
  SS_1_s = apply(RESD_1_s, 1, function(x){sum(x^2)})

  F_stats =  (SS_0_s - SS_1_s) / SS_1_s
  ### avoid cases where SS_1 = 0  & SS_0 = 0--> Quotient will be NaN
  F_stats[ which(SS_0_s == 0 & SS_1_s == 0)] = 0

  ### bootstrap the p-values
  p  =  unlist(lapply(F_stats,function(x){length(which(background >=
        x))/length(background)}))
  p_scaled = p.adjust(p, method = "BH")
  p_scaled_orig = p_scaled


  # calculate FCs

  # FCs timepoint vs earlier
  if (control_timecourse == FALSE) {
       means <- t(apply(data_table, 1,function(x){
       tmp = NULL
       for (i in 1:length(unique(timep))) {
       tmp = c(tmp,mean(2^x[which(data_annotation[colnames(data_table),
            "Time"] == sort(unique(timep))[i])]))
       }
       return(tmp)
       }))
       colnames(means) <- sort(unique(timep))
       Ratios_TPs <- matrix(rep(0,nrow(data_table)*(length(unique(timep)) - 1)),
            nrow(data_table),(length(unique(timep)) - 1))
       for (tt in 1:(length(unique(timep)) - 1)) {
          Ratios_TPs[,tt] <- means[,tt + 1]/means[,tt]
       }

       FC_TP_DEindex <- apply(Ratios_TPs,1,function(x){
          if (length(which(x > 2 || x < 0.5)) >= 2) { "DE" } else {"not_DE"}
       })
  }


  if (control_timecourse == TRUE) {

      # FCs case vs control
      means_case <- t(apply(fmat_case, 1,function(x){
      tmp = NULL
      for (i in 1:length(unique(timep_case))) {
      tmp = c(tmp, mean(2^x[which(data_annotation[colnames(fmat_case),
        "Time"] == sort(unique(timep_case))[i])]))
      }
      return(tmp)
      }))
      colnames(means_case) <- sort(unique(timep_case))

      means_control <- t(apply(fmat_ctrl, 1,function(x){
      tmp = NULL
      for (i in 1:length(unique(timep_ctrl))) {
      tmp = c(tmp, mean(2^x[which(data_annotation[colnames(fmat_ctrl),
            "Time"] == sort(unique(timep_ctrl))[i])]))
      }
      return(tmp)
      }))
      colnames(means_control) <- sort(unique(timep_ctrl))

      common_timepoints = sort(unique(timep_case))[sort(unique(timep_case)) %in%
            unique(timep_ctrl)]
      if (length(common_timepoints) == 0) {
         FC_DEindex = rep("not_detectable",nrow(means_case))
      } else {
       meanRatios <- t(apply(cbind(means_case[,as.character(common_timepoints)],
          means_control[,as.character(common_timepoints)]),1,function(x){
          tmp <- x[1:length(unique(common_timepoints))]/
              (x[(length(unique(common_timepoints)) + 1):(2*
              length(unique(common_timepoints)))])
          return(tmp)
        }))
        FCs <- t(apply(meanRatios,1,function(x){
          y = x
          y[x < 1] = -1/x[x < 1]
          return(y)
        }))
        FC_DEindex <- apply(FCs,1,function(x){
          if (length(which(abs(x) > 2)) >= 2) { "DE" } else {"not_DE"}
        })
      }

      # ANOVA
      anovas <- t(apply(data_table,1,function(x){summary(aov(as.numeric(x) ~
            data_annotation$Time + data_annotation$Condition +
            data_annotation$Time * data_annotation$Condition))[[1]][2:3,
            "Pr(>F)"]}))
      anovas_FDR =  cbind(p.adjust(anovas[,1], method = "BH"),
            p.adjust(anovas[,2], method = "BH"))
      anovas_pre = anovas_FDR[,1] < Q | anovas_FDR[,2] < Q
      anovas_index <- unlist(lapply(anovas_pre, function(x){if (x == TRUE) {"DE"
          } else {"not_DE"}}))
  }

  # SS1 Flag
  SS1_flag = SS_1_s < 20
  error_index = unlist(lapply(SS1_flag,function(x){if (x == FALSE) {"fit might be
      unstable due to large variation"} else {""}}))

  ### exit the function without error if no DE genes are detected
  if (!(TRUE %in% (p_scaled <= Q))) {
    warning("No DE genes were detected. Maybe amount of background genes is
            too low.")
    return(NULL)

  ### if DE genes are detected, finish FDR correction by using the cutoff
  } else {

    ### if control data is present but not as a timecourse, add t-test
    ### p-values to the results
    result =   as.data.frame(cbind("Gene" = row.names(data_table),
            "adj.p" = as.numeric(p_scaled_orig), stringsAsFactors = FALSE))
    result$adj.p <- as.numeric(as.character(result$adj.p))
    result = result[order(result$adj.p),]
    print(paste("Found ",nrow(result[result$adj.p <= Q,])," DE genes",sep = ""))


    if (control_timecourse == TRUE) {
       write.table(as.data.frame(cbind("Gene" = row.names(data_table),
            "adj.p" = p_scaled_orig, "at least 2 TPs for case vs. control" =
            FC_DEindex, "ANOVA condition or time*condition" = anovas_index,
            "prediction error" = error_index)),"pvals_and_flags.txt",
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
    } else {
       if (e_type == "Array") {
         write.table(as.data.frame(cbind("Gene" = row.names(data_table),
            "adj.p" = p_scaled_orig, "at least 2 consecutive TP Ratios" =
            FC_TP_DEindex, "prediction error" = error_index)),
            "pvals_and_flags.txt", quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = NA)
       } else {
         write.table(as.data.frame(cbind("Gene" = row.names(data_table),
            "adj.p" = p_scaled_orig,"prediction error" = error_index)),
            "pvals_and_flags.txt", quote = FALSE, sep = "\t",
            row.names = TRUE, col.names = NA)
       }
    }
    return(result[as.numeric(result$adj.p) <= Q,])
    }
}


################################################################################
################################################################################
################################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                     Final function calling all others
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Differential expression analysis using impulse models
#'
#' Fits an impulse model to time course data and uses this model as a basis
#' to detect differentially expressed (DE) genes. If a single time course
#' data set is given, DE genes are detected over time, whereas if an
#' additional control time course data set is present, DE genes are
#' detected between both datasets.
#' @aliases impulse_DE impulseDE
#' @param expression_table numeric matrix of expression values; genes should
#' be in rows, samples in columns. Data should be properly normalized and
#' log2-transformed as well as filtered for present or variable genes.
#' @param annotation_table table providing co-variables for the samples
#' including condition and time points. Time points must be numeric numbers.
#' @param colname_time character string specifying the column name of the
#' co-variable "\code{Time}" in \code{annotation_table}
#' @param colname_condition character string specifying the column name of
#' the co-variable "\code{Condition}" in \code{annotation_table}
#' @param control_timecourse logical indicating whether a control time
#' timecourse is part of the data set (\code{TRUE}) or not (\code{FALSE}).
#' Default is \code{FALSE}.
#' @param control_name character string specifying the name of the control
#' condition in \code{annotation_table}.
#' @param case_name character string specifying the name of the case
#' condition in \code{annotation_table}. Should be set if more than two
#' conditions are present in \code{annotation_table}.
#' @param expr_type character string with allowed values "\code{Array}" or
#' "\code{Seq}". Default is "\code{Array}".
#' @param plot_clusters logical indicating whether to plot the clusters
#' (\code{TRUE}) or not (\code{FALSE}). Default is \code{TRUE}.
#' @param n_iter numeric value specifying the number of iterations, which are
#' performed to fit the impulse model to the clusters. Default is \code{100}.
#' @param n_randoms numeric value specifying the number of generated randomized
#' background iterations, which are used for differential expression analysis.
#' Default is \code{50000} and this value should not be decreased.
#' @param n_process numeric value indicating the number of processes, which can
#' be used on the machine to run calculations in parallel. Default
#' is \code{4}. The specified value is internally changed to
#' \code{min(detectCores() - 1, n_process)} using the \code{detectCores}
#' function from the package \code{parallel} to avoid overload.
#' @param Q_value numeric value specifying the cutoff to call genes
#' significantly differentially expressed after FDR correction (adjusted
#' p-value). Default is \code{0.01}.
#' @return List containing the following elements:
#' \itemize{
#' \item \code{impulse_fit_results} List containing fitted values and model
#' parameters:
#' \itemize{
#' \item \code{impulse_parameters_combined} Matrix of fitted impulse model
#' parameters and sum of squared fitting errors for the combined dataset.
#' Not existing in the case of a single time course experiment.
#' \item \code{impulse_fits_combined} Matrix of impulse values calculated based
#' on the analyzed time points and the fitted model parameters for the combined
#' dataset. Not existing in the case of a single time course experiment.
#' \item \code{impulse_parameters_case} Matrix of fitted impulse model
#' parameters and sum of squared fitting errors for the case dataset.
#' \item \code{impulse_fits_case}  Matrix of impulse values calculated based
#' on the analyzed time points and the fitted model parameters for the case
#' dataset.
#' \item \code{impulse_parameters_control} Matrix of fitted impulse model
#' parameters and sum of squared fitting errors for the control dataset.
#' Not existing in the case of a single time course experiment.
#' \item \code{impulse_fits_control} Matrix of impulse values calculated based
#' on the analyzed time points and the fitted model parameters for the control
#' dataset. Not existing in the case of a single time course experiment.
#' }
#' \item \code{DE_results} Matrix containing the names of genes being called
#' as differentially expressed according to the specified cutoff \code{Q_value}
#' together with the adjusted p-values.
#' }
#' Additionally, \code{ImpulseDE} saves the following objects and tables into
#' the working directory:
#' \itemize{
#' \item \code{prepared_annotation.RData} Object containing the internally used
#' modified version of \code{annotation_table}
#' \item \code{clus_out2.txt} Text-file saving std out from the multi-threading.
#' Can be ignored.
#' \item \code{kmeans_clus_final_combined_genes.txt} Text-file containing
#' each gene together with its corresponding cluster number for the combined
#' dataset. Not existing in the case of a single time course experiment.
#' \item \code{kmeans_clus_final_case_genes.txt} Text-file containing
#' each gene together with its corresponding cluster number for the case
#' dataset.
#' \item \code{kmeans_clus_final_control_genes.txt} Text-file containing
#' each gene together with its corresponding cluster number for the control
#' dataset. Not existing in the case of a single time course experiment.
#' \item \code{clusters_combined_genes.pdf} PDF-file containing
#' the clusters for the combined dataset. Not existing in the case of a single
#' time course experiment.
#' \item \code{clusters_case_genes.pdf} PDF-file containing
#' the clusters for the case dataset.
#' \item \code{clusters_control_genes.pdf} PDF-file containing
#' the clusters for the control dataset. Not existing in the case of a single
#' time course experiment.
#' \item \code{impulse_fit_clusters.RData} Object containing a list of the
#' fitted impulse model values and paramaters to the cluster means; structure
#' is the same as for the list element \code{impulse_fit_results} of the output
#' value.
#' \item \code{impulse_fit_genes.RData} Object containing a list of the
#' fitted impulse model values and paramaters to the genes; structure is the
#' same as for the list element \code{impulse_fit_results} of the output value.
#' \item \code{cluster_out_random.txt} Text-file saving std out from the
#' multi-threading related to the randomized data. Can be ignored.
#' \item \code{background_results.RData} F-score values generated based on the
#' fits to the randomized data as used for the differential expression analysis.
#' \item \code{impulse_DE_genes.RData} Object containing names of genes being
#' called as differentially expressed according to the specified cutoff
#' \code{Q_value} together with the adjusted p-values; same as the list element
#' \code{DE_results} of the output value.
#' \item \code{pvals_and_flags.txt} Text-file containing all gene names
#' together with the adjusted p-values and flags for differential expression
#' according to additional tests.
#' }
#' @details \code{ImpulseDE} is based on the impulse model proposed by
#' Chechik and Koller, which reflects a two-step behavior of genes within a cell
#' responding to environmental changes (Chechik and Koller, 2009). To detect
#' differentially expressed genes, a five-step workflow is followed:
#' \enumerate{
#' \item The genes are clustered into a limited number of groups using k-means
#' clustering. If \code{plot_clusters} = \code{TRUE}, PDF documents are
#' generated, which contain plots of each cluster. Additionally, a text-file is
#' produced containing each gene together with its corresponding cluster number.
#' \item The impulse model is fitted to the mean expression profiles of the
#' clusters. The best parameter sets are then used for the next step.
#' \item The impulse model is fitted to each gene separately using the parameter
#' sets from step 2 as optimal start point guesses.
#' \item The impulse model is fitted to a randomized dataset (bootstrap), which
#' is essential to detect significantly differentially expressed genes
#' (Storey et al., 2005).
#' \item Detection of differentially expressed genes utilizing the fits to the
#' real and randomized data sets. FDR-correction is performed to obtain adjusted
#' p-values (Benjamini and Hochberg, 1995).
#' }
#' @examples
#' #' Install package longitudinal and load it
#' library(longitudinal)
#' #' Attach datasets
#' data(tcell)
#' #' check dimension of data matrix of interest
#' dim(tcell.10)
#' #' generate a proper annotation table
#' annot <- as.data.frame(cbind("Time" =
#'    sort(rep(get.time.repeats(tcell.10)$time,10)),
#'    "Condition" = "activated"), stringsAsFactors = FALSE)
#' #' Time columns must be numeric
#' annot$Time <- as.numeric(annot$Time)
#' #' rownames of annotation table must appear in data table
#' rownames(annot) = rownames(tcell.10)
#' #' apply ImpulseDE in single time course mode
#' #' since genes must be in rows, transpose data matrix using t()
#' #' For the example, reduce random iterations to 100 and number of
#' #' used processors to 1
#' impulse_results <- impulse_DE(t(tcell.10), annot, "Time", "Condition",
#'    n_randoms = 50, n_process = 1)
#' @seealso \code{\link{plot_impulse}}, \code{\link{calc_impulse}}.
#' @author Jil Sander
#' @references Benjamini, Y. and Hochberg, Y. (1995) Controlling the false
#' discovery rate: a practical and powerful approach to multiple testing.
#' J. R. Stat. Soc. Series B Stat. Methodol., 57, 289-300.
#' @references Storey, J.D. et al. (2005) Significance analysis of time course
#' microarray experiments. Proc. Natl. Acad. Sci. USA, 102, 12837-12841.
#' @references Rangel, C., Angus, J., Ghahramani, Z., Lioumi, M., Sotheran, E.,
#' Gaiba, A., Wild, D.L., Falciani, F. (2004) Modeling T-cell activation using
#' gene expression profiling and state-space models. Bioinformatics, 20(9),
#' 1361-72.
#' @references Chechik, G. and Koller, D. (2009) Timing of Gene Expression
#' Responses to Envi-ronmental Changes. J. Comput. Biol., 16, 279-290.
#' @references Yosef, N. et al. (2013) Dynamic regulatory network controlling
#' TH17 cell differentiation. Nature, 496, 461-468.
#' @export
impulse_DE <- function(expression_table = NULL, annotation_table = NULL,
    colname_time = NULL, colname_condition = NULL, control_timecourse = FALSE,
    control_name = NULL, case_name = NULL, expr_type = "Array",
    plot_clusters = TRUE, n_iter = 100, n_randoms = 50000, n_process = 4,
    Q_value = 0.01){

  tm_tot <- system.time({

    #' @import compiler
    #' @importFrom grDevices dev.off pdf
    #' @importFrom graphics abline axis legend par plot points
    #' @importFrom stats aov cor dist kmeans mad median nlminb optim p.adjust
    #' @importFrom stats runif sd
    #' @importFrom utils head write.table

    ## prepare annotation table for internal usage
    print("START: Prepare annotation table for internal usage")
    print("-------------------------------------------------------------------")

    prepared_annotation <- annotation_preparation(annotation_table,
          expression_table,
          colname_time ,colname_condition, control_timecourse, control_name,
          case_name)
    prepared_annotation <- prepared_annotation[order(
            prepared_annotation$Condition),]
    prepared_annotation <- prepared_annotation[order(prepared_annotation$Time),]
    save(prepared_annotation, file = file.path(getwd(),
            "/prepared_annotation.RData"))
    print("DONE")
    print("###################################################################")

    expression_table <- as.matrix(expression_table)    # to have numeric values
    expression_table <- expression_table[,rownames(prepared_annotation)]

    # exclude genes with missing values(NAs)
    indx <- apply(expression_table,1,function(x){TRUE %in% is.na(x)})
    expression_table <- expression_table[!(indx),]

    # if rownames are just 1,2,3 or if there are no rownames
    if (is.null(rownames(expression_table))) {
       rownames(expression_table) <- paste("G", 1:nrow(expression_table),
            sep = "_")
    } else if (length(grep("[a-zA-Z]",rownames(expression_table))) == 0) {
       rownames(expression_table) <- paste(rownames(expression_table),"G",
            sep = "_")
    }

    ### cluster genes to reduce efforts for fitting the impulse model
    print("START: Clustering genes for Impulse model fit")
    print("-------------------------------------------------------------------")
    load(file.path(getwd(),"/prepared_annotation.RData"))
    tm_clust <- system.time({
      clustering_results <- cluster_genes_for_impulse(expression_table,
          prepared_annotation, control_timecourse, control_name, plot_clusters)
    })
    save(clustering_results, file = file.path(getwd(),"clustering_results.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_clust["elapsed"]/60,2)," min",sep = ""))
    print("##################################################################")

    # fit Impulse model to the clusters
    print("START: Fitting Impulse model to the clusters")
    print("-------------------------------------------------------------------")
    load(file.path(getwd(),"clustering_results.RData"))
    tm_imp_fit_clus <- system.time({
      impulse_fit_clusters <- impulse_fit(clustering_results,prepared_annotation,
          n_iter, control_timecourse, control_name, n_proc = n_process)
    })
    save(impulse_fit_clusters,file = file.path(getwd(),
        "impulse_fit_clusters.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_imp_fit_clus["elapsed"]/60,2)," min",
        sep = ""))
    print("###################################################################")

    ###  fit Impule model to each gene by using the cluster fits as start values
    print("START: Fitting Impulse model to the genes")
    print("-------------------------------------------------------------------")
    load(file.path(getwd(),"impulse_fit_clusters.RData"))
    tm_imp_fit_gen <- system.time({
      impulse_fit_genes <- impulse_fit(expression_table, prepared_annotation,
        n_iter, control_timecourse, control_name, clustering_results,
        impulse_fit_clusters, n_proc = n_process)
    })
    save(impulse_fit_genes,file = file.path(getwd(),"impulse_fit_genes.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_imp_fit_gen["elapsed"]/60,2),
        " min",sep = ""))
    print("###################################################################")

    ### generate background for the DE analysis
    print("START: Generate background")
    print("-------------------------------------------------------------------")
    load(file.path(getwd(),"impulse_fit_genes.RData"))
    tm_bg <- system.time({
      background_results <-  generate_background(expression_table,
        prepared_annotation, n_iter, impulse_fit_genes, control_timecourse,
        control_name, clustering_results, n_randoms, n_process)
    })
    save(background_results,file = file.path(getwd(),"background_results.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_bg["elapsed"]/60,2)," min",sep = ""))
    print("###################################################################")

    ## detect differentially expressed genes
    print("START: DE analysis")
    print("-------------------------------------------------------------------")
    load(file.path(getwd(),"background_results.RData"))
    tm_DE <- system.time({
      impulse_DE_genes <- DE_analysis(expression_table,prepared_annotation,
            impulse_fit_genes, background_results, control_timecourse,
            control_name, expr_type, Q_value)
    })
    ## plot the top DE genes
    if (control_timecourse == TRUE) {
      case_ind =
          as.character(prepared_annotation$Condition[prepared_annotation$Condition
          != control_name][1])
    } else {case_ind = NULL}
    if (is.null(impulse_DE_genes) == FALSE & length(impulse_DE_genes) > 1) {
      if (is.list(impulse_DE_genes) == TRUE &
         is.data.frame(impulse_DE_genes) == FALSE) {
            plot_impulse(as.character(impulse_DE_genes[[1]][,
              "Gene"])[1:(min(nrow(impulse_DE_genes[[1]]),2*18))],
              expression_table, prepared_annotation,impulse_fit_genes,
              control_timecourse, control_name, case_ind, file_name_part = "DE",
              title_line = "", sub_line = "")
        } else {
            plot_impulse(as.character(impulse_DE_genes[,
              "Gene"])[1:(min(nrow(impulse_DE_genes),2*18))],
              expression_table, prepared_annotation, impulse_fit_genes,
              control_timecourse, control_name, case_ind, file_name_part = "DE",
              title_line = "", sub_line = "")
      }
    }
    save(impulse_DE_genes,file = file.path(getwd(),"impulse_DE_genes.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_DE["elapsed"]/60,2)," min",sep = ""))
    print("##################################################################")
  })
  print(paste("TOTAL consumed time: ",round(tm_tot["elapsed"]/60,2),
        " min",sep = ""))
  return(list("impulse_fit_results" = impulse_fit_genes,
        "DE_results" = impulse_DE_genes))
}
