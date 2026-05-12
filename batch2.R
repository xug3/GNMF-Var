library(R.matlab)
library(FAST,lib.loc="/gpfs/gibbs/project/wang_zuoheng/gx28/R_lib")
library(openxlsx)


results=read.xlsx("Freyja_cov_24_25_final_v3.xlsx",detectDates = TRUE )
#results[,4]=as.Date(results[,4], format = "%m/%d/%Y")
results=results[order(results[,4]),]
colnames(results)
#results=results[,-c(23,32,33,30,19,37,40,25)]
dim(results)
results[results$GE.1>0.15,c("Date","GE.1")]
results[1:10,]



new_sample=read.csv("mutation_frequency_matrix_processed_redup_correctdatev2.csv")
dim(new_sample)
new_sample[1:10,1:10]
setdiff(new_sample[,1],results[,1])
which(results[,1]=="Feb24a_SARS64")
results[436,]
str(new_sample)
# check out data structure
countmatrixt = t(new_sample[,-c(1:4)])
#countmatrixt = t(cross_dsads[[1]])
dim(countmatrixt)
colsumsall=rowSums(countmatrixt)
rowcounts=apply(countmatrixt,1,function(x) length(which(x>0)))
table(rowcounts)
mutationname=rownames(countmatrixt)
sample_name=new_sample[,1]

date_sample=as.Date(new_sample[,4])
table(date_sample)

adjacent_matrix=diag(dim(new_sample)[1])
for(i in 1:(dim(new_sample)[1]-1)){
  locationi=new_sample[i,3]
  datei=as.Date(date_sample[i])
  for(j in (i+1):dim(new_sample)[1]){
    locationj=new_sample[j,3]
    datej=as.Date(date_sample[j])
    date_diff <- datej - datei
    # Check if the difference is greater than one week (7 days)
    if((date_diff <=7)&(locationi==locationj)){
      adjacent_matrix[i,j]=adjacent_matrix[j,i]=1
    } 
  }
}
adjacent_matrix[1:20,1:20]


# ============================================================
# Batch 2 GNMF runtime loop
# Put this after countmatrixt, date_sample, adjacent_matrix are created
# ============================================================

# -----------------------------
# Quarter information for batch 2
# -----------------------------
quarter_info <- data.frame(
  quarter = c("2023-Q4", "2024-Q1", "2024-Q2", "2024-Q3", "2024-Q4"),
  start_date = as.Date(c("2023-10-20", "2024-01-01", "2024-04-01",
                         "2024-07-01", "2024-10-01")),
  end_date = as.Date(c("2023-12-31", "2024-03-31", "2024-06-30",
                       "2024-09-30", "2025-01-15")),
  K = c(20, 20, 20, 20, 20)
)

# -----------------------------
# Output folder
# -----------------------------
out_dir <- "GNMF_runtime_batch2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

runtime_results <- data.frame()

# -----------------------------
# Loop through quarters
# -----------------------------
for (q in seq_len(nrow(quarter_info))) {
  
  qname <- quarter_info$quarter[q]
  start_q <- quarter_info$start_date[q]
  end_q <- quarter_info$end_date[q]
  K_q <- quarter_info$K[q]
  
  message("Running GNMF for ", qname, " ...")
  
  # Select samples in this quarter
  idx <- which(date_sample >= start_q & date_sample <= end_q)
  
  if (length(idx) == 0) {
    warning("No samples found for ", qname)
    next
  }
  
  # Subset mutation matrix and adjacency matrix
  subweek_countmatrixt <- countmatrixt[, idx, drop = FALSE]
  sub_adjacent_matrix <- adjacent_matrix[idx, idx, drop = FALSE]
  
  # Filter mutations present in at least 3 samples
  rowcounts <- rowSums(subweek_countmatrixt > 0)
  subweek_countmatrixt <- subweek_countmatrixt[rowcounts > 2, , drop = FALSE]
  
  # Set row/column names
  rownames(subweek_countmatrixt) <- rownames(subweek_countmatrixt)
  colnames(subweek_countmatrixt) <- paste0(new_sample[idx, 3], "_", new_sample[idx, 4])
  
  # GNMF config
  config <- dmain_config
  config$r <- K_q
  config$lambda_1 <- 0.1   # GNMF graph regularization
  config$lambda_2 <- 0.1
  config$n_iter <- 500000
  
  # Run and time GNMF
  set.seed(1)
  gc()
  
  t_run <- system.time({
    res <- dmain(
      subweek_countmatrixt,
      sub_adjacent_matrix,
      config
    )
  })
  
  saveRDS(res,file.path(out_dir,paste0(qname,".rds")))
  # Save runtime information
  runtime_i <- data.frame(
    quarter = qname,
    start_date = start_q,
    end_date = end_q,
    sample_size = ncol(subweek_countmatrixt),
    mutations_after_filter = nrow(subweek_countmatrixt),
    K = K_q,
    lambda_1 = config$lambda_1,
    lambda_2 = config$lambda_2,
    n_iter = config$n_iter,
    elapsed_sec = as.numeric(t_run["elapsed"]),
    elapsed_min = as.numeric(t_run["elapsed"]) / 60,
    input_matrix_MB = as.numeric(object.size(subweek_countmatrixt)) / 1024^2,
    adjacency_matrix_MB = as.numeric(object.size(sub_adjacent_matrix)) / 1024^2
  )
  
  runtime_results <- rbind(runtime_results, runtime_i)
  
  # Save partial runtime table after each quarter
  write.csv(
    runtime_results,
    file = file.path(out_dir, "GNMF_batch2_runtime_results_partial.csv"),
    row.names = FALSE
  )
  
  
}

# -----------------------------
# Save final runtime table
# -----------------------------
write.csv(
  runtime_results,
  file = file.path(out_dir, "GNMF_batch2_runtime_results.csv"),
  row.names = FALSE
)


runtime_results
