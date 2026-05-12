library(R.matlab)
library(openxlsx)
library(dplyr)
library(tidyr)
library(FAST,lib.loc="/gpfs/gibbs/project/wang_zuoheng/gx28/R_lib")


results=read.xlsx("XZ_run_v2_rename.xlsx")
results <- results %>%
  mutate(across(where(is.logical), as.numeric)) %>%  # 全NA列从 logical 转成 numeric
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
x_chr <- sprintf("%06d", as.integer(results[[3]]))   # 变字符并左侧补零到6位
results[,3]=as.Date(x_chr, format = "%m%d%y")
results=results[order(results[,3]),]
v_res=results[,1]
results[1:10,]

results[which(results$XBB.1.16>0.15)[1],]

first_date_table=NULL
for(i in 1:18) {
  variantsi=colnames(results)[i+5]
  datei=results[results[,i+5]>0.15,3][1]
  first_date_table=rbind(first_date_table,data.frame(variantsi,datei))
}
#order_first_date_table=data.frame(first_date_table[order(first_date_table[,2]),])  
colnames(first_date_table)=c("variants","date_freyja")
first_date_table


cross_dsads <- readMat("Yf_all_ivar_variant_sep21_nov23_50xcoverage_gt_80_vRefine_vNoDup_github.mat")

# check out data structure
countmatrixt = t(cross_dsads[[1]])
dim(countmatrixt)
colsumsall=rowSums(countmatrixt)
rowcounts=apply(countmatrixt,1,function(x) length(which(x>0)))
#countmatrixt =countmatrixt[-which(rowcounts<6),]
dim(countmatrixt)
samplenames=cross_dsads[[2]]


mutationname=cross_dsads[[3]]
mutationname=c(unlist(mutationname))

str(samplenames)
samplenames[1,2][[1]][[1]][1,1]
samplenames[1,1]

adjacent_matrix=diag(rep(1,1397))
for(i in 1:1396){
  if(!is.null(samplenames[i,2][[1]])){
    spliti=strsplit(samplenames[i,2][[1]][[1]][1,1],split="-")
    locationi=spliti[[1]][1]
    datei=as.Date(spliti[[1]][length(spliti[[1]])], format = "%m%d%y")
    
    for(j in (i+1):1397){
      if(!is.null(samplenames[j,2][[1]])){
        locationj=strsplit(samplenames[j,2][[1]][[1]][1,1],split="-")[[1]][1]
        splitj=strsplit(samplenames[j,2][[1]][[1]][1,1],split="-")
        locationj=splitj[[1]][1]
        datej=as.Date(splitj[[1]][length(splitj[[1]])], format = "%m%d%y")
        date_diff <- datej - datei
        # Check if the difference is greater than one week (7 days)
        if((date_diff <=7)&(locationi==locationj)){
          adjacent_matrix[i,j]=adjacent_matrix[j,i]=1
        } 
      }
    }
  }
}
adjacent_matrix[1:20,1:20]


date_sample=rep(0,1397)
loc_date=rep(0,1397)
sampleid=rep(0,1397)
for(k in 1:1397) date_sample[k]=samplenames[k,][[3]][[1]][1,1]
for(k in 1:1397) sampleid[k]=samplenames[k,][[1]][[1]][1,1]
date_location=rep(0,1397)
for(k in 1:1397) {
  if(!is.null(samplenames[k,][[2]][[1]][1,1])){
    spliti=strsplit(samplenames[k,2][[1]][[1]][1,1],split="-")
    date_location[k]=spliti[[1]][1]
    date_veci=as.Date(spliti[[1]][length(spliti[[1]])], format = "%m%d%y")
    date_veci2=spliti[[1]][length(spliti[[1]])]
    loc_date[k]=paste0(date_location[k],"-",date_veci2)
  }
}
date_location[1:20]
date_sample[1:20]
length(loc_date)
dim(countmatrixt)


matchindx=match(v_res,sampleid)
colnames(countmatrixt)=loc_date
rownames(countmatrixt)=mutationname
countmatrixt=countmatrixt[,matchindx]
sampleidsub=sampleid[matchindx]
dim(countmatrixt)
adjacent_matrix=adjacent_matrix[matchindx,matchindx]
samplenames=samplenames[matchindx,]
date_sample=date_sample[matchindx]
date_sample[date_sample=="07312323"]="073123"



# -----------------------------
# Quarter information for batch 1
# -----------------------------
quarter_info <- data.frame(
  quarter = c("2021-Q4", "2022-Q1", "2022-Q2", "2022-Q3",
              "2022-Q4", "2023-Q1", "2023-Q2", "2023-Q3"),
  start_date = as.Date(c("2021-09-06", "2022-01-01", "2022-04-01", "2022-07-01",
                         "2022-10-01", "2023-01-01", "2023-04-01", "2023-07-01")),
  end_date = as.Date(c("2021-12-31", "2022-03-31", "2022-06-30", "2022-09-30",
                       "2022-12-31", "2023-03-31", "2023-06-30", "2023-10-30")),
  K = c(4, 4, 20, 20, 20, 20, 20, 20)
)

# -----------------------------
# Extract sample date from column names, e.g. Facility5-071122
# -----------------------------
get_sample_date <- function(x) {
  parts <- strsplit(x, "-")[[1]]
  date_chr <- parts[length(parts)]
  
  # Fix occasional duplicated year strings, e.g. 07312323 -> 073123
  if (nchar(date_chr) > 6) {
    date_chr <- substr(date_chr, 1, 6)
  }
  
  as.Date(date_chr, format = "%m%d%y")
}

sample_dates <- as.Date(sapply(colnames(countmatrixt), get_sample_date))

# Optional check
print(summary(sample_dates))
print(table(is.na(sample_dates)))

# -----------------------------
# Runtime results table
# -----------------------------
runtime_results <- data.frame()

# -----------------------------
# Output folder
# -----------------------------
out_dir <- "GNMF_runtime_batch1"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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
  idx <- which(sample_dates >= start_q & sample_dates <= end_q)
  
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
  # Collect runtime information
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
  
  # Save partial results after each quarter
  write.csv(
    runtime_results,
    file = file.path(out_dir, "GNMF_batch1_runtime_results_partial.csv"),
    row.names = FALSE
  )
  
}

# -----------------------------
# Save final runtime table
# -----------------------------
write.csv(
  runtime_results,
  file = file.path(out_dir, "GNMF_batch1_runtime_results.csv"),
  row.names = FALSE
)


runtime_results