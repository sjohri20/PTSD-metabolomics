#Executing disjoint metabolomics 
setwd("C:/Users/shj622/Desktop/Shreya/Metabolomics_Clinical_Data/")
options(stringsAsFactors=FALSE)

library(limma)
#metabolite data
metab.T1 = read.delim("FtCampbellT1-Metabolomics-Plasma-USACEHR-Metabolon-20171207-raw.txt", sep = "\t")
msamples.T1 = colnames(metab.T1)
colnames(metab.T1) = sapply(1:length(msamples.T1), function(x)substr(msamples.T1[x],1,9))
# metab.T1 = log10(metab.T1)
# metab.T1 = normalizeMedianAbsValues(metab.T1)
# metab.T1 = as.data.frame(metab.T1)
#metab.T1 = log2(metab.T1)

metab.T2 = read.delim("FtCampbellT3-Metabolomics-Plasma-USACEHR-Metabolon-20171207-raw.txt", sep = "\t")
msamples.T2 = colnames(metab.T2)
colnames(metab.T2) = sapply(1:length(msamples.T2), function(x)substr(msamples.T2[x],1,9))
# metab.T2 = log10(metab.T2)
# metab.T2 = normalizeMedianAbsValues(metab.T2)
# metab.T2 = as.data.frame(metab.T2)
#metab.T2 = log2(metab.T2)

#clinical data
clin.T1 = read.delim("FtCampbellT1-Clinical-All-All-ScoredData-20160620-raw.txt", sep = "\t")
csamples.T1 = colnames(clin.T1)
colnames(clin.T1) = sapply(1:length(csamples.T1), function(x)substr(csamples.T1[x],1,9))

clin.T2 = read.delim("FtCampbellT3-Clinical-All-All-ScoredData-20181004-raw.txt", sep = "\t")
csamples.T2 = colnames(clin.T2)
colnames(clin.T2) = sapply(1:length(csamples.T2), function(x)substr(csamples.T2[x],1,9))

#com_metab_names = intersect(colnames(metab.T1), colnames(metab.T2))
com_clin_metab_names_1 = intersect(colnames(clin.T1), colnames(metab.T1))
com_clin_metab_names_2 = intersect(colnames(clin.T2), colnames(metab.T2))

common_metab_1 = metab.T1[,com_clin_metab_names_1]
common_metab_2 = metab.T2[,com_clin_metab_names_2]

common_clin_1 = clin.T1[,com_clin_metab_names_1]
common_clin_2 = clin.T2[,com_clin_metab_names_2]

#to get PCL scores as row 4 in each dataset
common_row = intersect(rownames(common_clin_1), rownames(common_clin_2))
common_clin_1 = common_clin_1[common_row,]
common_clin_2 = common_clin_2[common_row,]

#TBI patients to be removed
a = colnames(common_clin_1[which(common_clin_1[7,]=="mild TBI")])
a_temp = common_clin_1[-which(is.na(common_clin_1[9,])==TRUE)]
a = intersect(a,colnames(a_temp[which((as.numeric(a_temp[9,])>=8)==TRUE)]))
b = colnames(common_clin_1[which(common_clin_1[7,]=="moderate TBI")])
c = colnames(common_clin_1[which(common_clin_1[7,]=="severe TBI")])
d = colnames(common_clin_2[which(common_clin_2[7,]=="mild TBI")])
d_temp = common_clin_2[-which(is.na(common_clin_2[9,])==TRUE)]
d = intersect(d,colnames(d_temp[which((as.numeric(d_temp[9,])>=8)==TRUE)]))
e = colnames(common_clin_2[which(common_clin_2[7,]=="moderate TBI")])
f = colnames(common_clin_2[which(common_clin_2[7,]=="severe TBI")])

#removing TBI patients
common_clin_1 = common_clin_1[-(which(colnames(common_clin_1) %in% intersect(a,colnames(common_clin_1))))]
common_clin_1 = common_clin_1[-(which(colnames(common_clin_1) %in% intersect(b,colnames(common_clin_1))))]
common_clin_1 = common_clin_1[-(which(colnames(common_clin_1) %in% intersect(c,colnames(common_clin_1))))]

common_clin_2 = common_clin_2[-(which(colnames(common_clin_2) %in% intersect(d,colnames(common_clin_2))))]
common_clin_2 = common_clin_2[-(which(colnames(common_clin_2) %in% intersect(e,colnames(common_clin_2))))]
#common_clin_2 = common_clin_2[-(which(colnames(common_clin_2) %in% intersect(f,colnames(common_clin_2))))]

common_metab_1 = common_metab_1[-(which(colnames(common_metab_1) %in% intersect(a,colnames(common_metab_1))))]
common_metab_1 = common_metab_1[-(which(colnames(common_metab_1) %in% intersect(b,colnames(common_metab_1))))]
common_metab_1 = common_metab_1[-(which(colnames(common_metab_1) %in% intersect(c,colnames(common_metab_1))))]
#normalising data
common_metab_1 = log2(common_metab_1)
common_metab_1 = normalizeMedianAbsValues(common_metab_1)
common_metab_1 = as.data.frame(common_metab_1)

common_metab_2 = common_metab_2[-(which(colnames(common_metab_2) %in% intersect(d,colnames(common_metab_2))))]
common_metab_2 = common_metab_2[-(which(colnames(common_metab_2) %in% intersect(e,colnames(common_metab_2))))]
# common_metab_2 = common_metab_2[-(which(colnames(common_metab_2) %in% intersect(f,colnames(common_metab_2))))]
common_metab_2 = log2(common_metab_2)
common_metab_2 = normalizeMedianAbsValues(common_metab_2)
common_metab_2 = as.data.frame(common_metab_2)

#Removing PCL-scores = NA in clinical data
trial_clin1 = common_clin_1[,-which(is.na(as.integer(common_clin_1[4,])) == TRUE)]
trial_clin2 = common_clin_2[,-which(is.na(as.integer(common_clin_2[4,])) == TRUE)]  

#segregating cases and controls in clinical data
#for T1
#getting controls
temp_control_wna=trial_clin1
count1=0
for (i in 1:ncol(trial_clin1))
{
  if (as.integer(trial_clin1[4,i])!=0)
  {
    temp_control_wna=temp_control_wna[,-(i-count1)]
    count1=count1+1
  }
}

#write.csv(temp_control_wna,"clin1_PTSD_controls_intersect_without_na.csv")

#getting cases
temp_cases_wna = trial_clin1
count2=0
for (i in 1:ncol(trial_clin1))
{
  if (as.integer(trial_clin1[4,i])<31)
  {
    temp_cases_wna=temp_cases_wna[,-(i-count2)]
    count2=count2+1
  }
}
#write.csv(temp_cases_wna,"clin1_PTSD_cases_intersect_without_na.csv")

#for T3
#getting controls
temp_control_wna_2=trial_clin2
count3=0
for (i in 1:ncol(trial_clin2))
{
  if (as.integer(trial_clin2[4,i])!=0)
  {
    temp_control_wna_2=temp_control_wna_2[,-(i-count3)]
    count3=count3+1
  }
}

#write.csv(temp_control_wna_2,"clin2_PTSD_controls_intersect_without_na.csv")

#getting cases
temp_cases_wna_2 = trial_clin2
count4=0
for (i in 1:ncol(trial_clin2))
{
  if (as.integer(trial_clin2[4,i])<31)
  {
    temp_cases_wna_2=temp_cases_wna_2[,-(i-count4)]
    count4=count4+1
  }
}
#write.csv(temp_cases_wna_2,"clin2_PTSD_cases_intersect_without_na.csv")


#analysing metabolites for T1
mean_T1 = list()
sd_T1 = list()    #indices
cv_T1 = list()    #indices
metab_remove_cv = list()
metab_remove_iqr = list()
iqr_T1 = list()
count1=0
count2=0
count3=0
for (i in 1:nrow(common_metab_1))
{
  mean_T1[i] = mean(as.numeric(common_metab_1[i,2:ncol(common_metab_1)]))
  sd_T1[i] = sd(as.numeric(common_metab_1[i,2:ncol(common_metab_1)]))
  if (as.numeric(mean_T1[i])==0)
  {
    count1=count1+1
    metab_remove_cv[count1] = i
    cv_T1[i] = 99
  }
  # }else
  # {
  #   cv_T1[i] = as.numeric(sd_T1[i])/as.numeric(mean_T1[i])
  #   if (as.numeric(cv_T1[i])<=0.05)
  #   {
  #     count1=count1+1
  #     metab_remove_cv[count1] = i
  #   }
  # }
  # iqr_T1[i] = IQR(as.numeric(common_metab_1[i,2:ncol(common_metab_1)]), na.rm=FALSE)
  # if (as.numeric(iqr_T1[i])<=0.1)
  # {
  #   count2=count2+1
  #   metab_remove_iqr[count2] = i
  # }
  # if (as.numeric(sd_T1[i])==0)
  # {
  #   count3=count3+1
  # }
}

#Analysing distribution of cv_T1
h1 = hist(as.numeric(cv_T1), main="Histogram for distribution of CV", xlab="CV", ylim = c(0,1500))
text(h1$mids,h1$counts,labels=h1$counts, adj=c(0.5, -0.5))

h2 = hist(as.numeric(iqr_T1), main="Histogram for distribution of IQR", xlab="IQR", ylim = c(0,1000))
text(h2$mids,h2$counts,labels=h2$counts, adj=c(0.5, -0.5))

remove_metab = unique(c(c(as.numeric(metab_remove_cv)), c(as.numeric(metab_remove_iqr))))

#t_tests between cases and controls for each metabolite at T1
new_metab_list = common_metab_1[-remove_metab,]

new_metab_cases = new_metab_list[,colnames(temp_cases_wna)]
new_metab_controls = new_metab_list[,colnames(temp_control_wna)]

tr_new_metab_cases.T1 = t(new_metab_cases)
tr_new_metab_controls.T1 = t(new_metab_controls)
# t_test.T1 = data.frame()
# cname = colnames(tr_new_metab_cases.T1)
# count_pval_1=0
# p_val.T1 = list()
# count1=0
# up.T1 = list()
# count2=0
# down.T1 = list()
# for (i in 1:ncol(tr_new_metab_cases.T1))
# {
#   ttest = t.test(tr_new_metab_cases.T1[,i], tr_new_metab_controls.T1[,i])
#   #write(ttest$p.value, file = "t_tests_T1.txt", append=TRUE)
#   t_test.T1[1,i] = as.numeric(ttest$p.value)
#   t_test.T1[2,i] = as.numeric(ttest$statistic)
#   if (as.numeric(ttest$p.value)<=0.05)
#   {
#     count_pval_1=count_pval_1+1
#     p_val.T1[count_pval_1] = cname[i]
#   }
#   if (as.numeric(ttest$statistic)>0)
#   {
#     count1=count1+1
#     up.T1[count1]= cname[i]
#   } else if (as.numeric(ttest$statistic)<0)
#   {
#     count2=count2+1
#     down.T1[count2]= cname[i]
#   }
# }
# colnames(t_test.T1) = cname
# rownames(t_test.T1) = c("p-value","t-statistic")
# up_regulated_metabs.T1 = t_test.T1[,intersect(as.character(up.T1), as.character(p_val.T1))]
# down_regulated_metabs.T1 = t_test.T1[,intersect(as.character(down.T1), as.character(p_val.T1))]
# #write.csv(t_test.T1,"T-test-T1-only-0.csv")
# #write.csv(up_regulated_metabs.T1,"upregulated_T1.csv")
# #write.csv(down_regulated_metabs.T1,"downregulated_T1.csv")


#analysing metabolites for T3
mean_T2 = list()
sd_T2 = list()    #indices
cv_T2 = list()    #indices
metab_remove_cv_2 = list()
metab_remove_iqr_2 = list()
iqr_T2 = list()
count1=0
count2=0
count3=0
for (i in 1:nrow(common_metab_2))
{
  mean_T2[i] = mean(as.numeric(common_metab_2[i,2:ncol(common_metab_2)]))
  sd_T2[i] = sd(as.numeric(common_metab_2[i,2:ncol(common_metab_2)]))
  cv_T2[i] = as.numeric(sd_T2[i])/as.numeric(mean_T2[i])
  iqr_T2[i] = IQR(as.numeric(common_metab_2[i,2:ncol(common_metab_2)]), na.rm=FALSE)
  if (mean_T2[i]==0)
  {
    count1=count1+1
    metab_remove_cv_2[count1] = i
    cv_T2[i] = 99
  }
  # }else
  # {
  #   cv_T2[i] = as.numeric(sd_T2[i])/as.numeric(mean_T2[i])
  #   if (as.numeric(cv_T2[i])<=0.1)
  #   {
  #     count1=count1+1
  #     metab_remove_cv[count1] = i
  #   }
  # }
  # if (as.numeric(iqr_T2[i])<=0.1)
  # {
  #   count2=count2+1
  #   metab_remove_iqr_2[count2] = i
  # }
  # if (as.numeric(sd_T2[i])==0)
  # {
  #   count3=count3+1
  # }
}

#Analysisng distribution of cv_T2
h1 = hist(as.numeric(cv_T2), main="Histogram for distribution of CV", xlab="CV", ylim = c(0,1500))
text(h1$mids,h1$counts,labels=h1$counts, adj=c(0.5, -0.5))

h2 = hist(as.numeric(iqr_T2), main="Histogram for distribution of IQR", xlab="IQR", ylim = c(0,1500))
text(h2$mids,h2$counts,labels=h2$counts, adj=c(0.5, -0.5))

remove_metab_2 = unique(c(c(as.numeric(metab_remove_cv_2)), c(as.numeric(metab_remove_iqr_2))))

#t_tests between cases and controls for each metabolite at T3
new_metab_list_2 = common_metab_2[-remove_metab_2,]

new_metab_cases_2 = new_metab_list_2[,colnames(temp_cases_wna_2)]
new_metab_controls_2 = new_metab_list_2[,colnames(temp_control_wna_2)]

tr_new_metab_cases.T2 = t(new_metab_cases_2)
tr_new_metab_controls.T2 = t(new_metab_controls_2)
# t_test.T2 = data.frame()
# cname_2 = colnames(tr_new_metab_cases.T2)
# count_pval_2=0
# p_val.T2 = list()
# count1=0
# up.T2 = list()
# count2=0
# down.T2 = list()
# for (i in 1:ncol(tr_new_metab_cases.T2))
# {
#   ttest = t.test(tr_new_metab_cases.T2[,i], tr_new_metab_controls.T2[,i])
#   #write(ttest$p.value, file = "t_tests_T3.txt", append=TRUE)
#   t_test.T2[1,i] = as.numeric(ttest$p.value)
#   t_test.T2[2,i] = as.numeric(ttest$statistic)
#   if (as.numeric(ttest$p.value)<=0.05)
#   {
#     count_pval_2=count_pval_2+1
#     p_val.T2[count_pval_2] = cname_2[i]
#   }
#   if (as.numeric(ttest$statistic)>0)
#   {
#     count1=count1+1
#     up.T2[count1]= cname_2[i]
#   } else if (as.numeric(ttest$statistic)<0)
#   {
#     count2=count2+1
#     down.T2[count2]= cname_2[i]
#   }
# }
# colnames(t_test.T2) = cname_2
# rownames(t_test.T2) = c("p-value","t-statistic")
# up_regulated_metabs.T2 = t_test.T2[,intersect(as.character(up.T2), as.character(p_val.T2))]
# down_regulated_metabs.T2 = t_test.T2[,intersect(as.character(down.T2), as.character(p_val.T2))]
# #write.csv(t_test.T2,"T-test-T3-only-0.csv")
# #write.csv(up_regulated_metabs.T2,"upregulated_T3.csv")
# #write.csv(down_regulated_metabs.T2,"downregulated_T3.csv")
# 
# 
# #calculating mean for cases and controls for each metab in T1
# frame1 = data.frame()
# 
# up_list_mean_cases.T1 = list()
# up_list_sd_cases.T1 = list()
# up_list_mean_controls.T1 = list()
# up_list_sd_controls.T1 = list()
# for (i in 1:ncol(up_regulated_metabs.T1))
# {
#   up_list_mean_cases.T1[i] = mean(as.numeric(new_metab_cases[colnames(up_regulated_metabs.T1)[i],]))
#   up_list_sd_cases.T1[i] = sd(as.numeric(new_metab_cases[colnames(up_regulated_metabs.T1)[i],]))
#   up_list_mean_controls.T1[i] = mean(as.numeric(new_metab_controls[colnames(up_regulated_metabs.T1)[i],]))
#   up_list_sd_controls.T1[i] = sd(as.numeric(new_metab_controls[colnames(up_regulated_metabs.T1)[i],]))
#   frame1[i,1] = as.numeric(up_list_mean_cases.T1[i])
#   frame1[i,2] = as.numeric(up_list_sd_cases.T1[i])
#   frame1[i,3] = as.numeric(up_list_mean_controls.T1[i])
#   frame1[i,4] = as.numeric(up_list_sd_controls.T1[i])
#   frame1[i,5] = as.numeric(up_regulated_metabs.T1[1,i])
#   frame1[i,6] = as.numeric(up_regulated_metabs.T1[2,i])
#   rownames(frame1)[i] = as.character(intersect(as.character(up.T1), as.character(p_val.T1)))[i]
# }
# down_list_mean_cases.T1 = list()
# down_list_sd_cases.T1 = list()
# down_list_mean_controls.T1 = list()
# down_list_sd_controls.T1 = list()
# for (j in 1:ncol(down_regulated_metabs.T1))
# {
#   down_list_mean_cases.T1[j] = mean(as.numeric(new_metab_cases[colnames(down_regulated_metabs.T1)[j],]))
#   down_list_sd_cases.T1[j] = sd(as.numeric(new_metab_cases[colnames(down_regulated_metabs.T1)[j],]))
#   down_list_mean_controls.T1[j] = mean(as.numeric(new_metab_controls[colnames(down_regulated_metabs.T1)[j],]))
#   down_list_sd_controls.T1[j] = sd(as.numeric(new_metab_controls[colnames(down_regulated_metabs.T1)[j],]))
#   frame1[i+j,1] = as.numeric(down_list_mean_cases.T1[j])
#   frame1[i+j,2] = as.numeric(down_list_sd_cases.T1[j])
#   frame1[i+j,3] = as.numeric(down_list_mean_controls.T1[j])
#   frame1[i+j,4] = as.numeric(down_list_sd_controls.T1[j])
#   frame1[i+j,5] = as.numeric(down_regulated_metabs.T1[1,j])
#   frame1[i+j,6] = as.numeric(down_regulated_metabs.T1[2,j])
#   rownames(frame1)[i+j] = as.character(intersect(as.character(down.T1), as.character(p_val.T1)))[j]
# }
# 
# #calculating mean for cases and controls for each metab in T3
# up_list_mean_cases.T2 = list()
# up_list_sd_cases.T2 = list()
# up_list_mean_controls.T2 = list()
# up_list_sd_controls.T2 = list()
# for (k in 1:ncol(up_regulated_metabs.T2))
# {
#   up_list_mean_cases.T2[k] = mean(as.numeric(new_metab_cases_2[colnames(up_regulated_metabs.T2)[k],]))
#   up_list_sd_cases.T2[k] = sd(as.numeric(new_metab_cases_2[colnames(up_regulated_metabs.T2)[k],]))
#   up_list_mean_controls.T2[k] = mean(as.numeric(new_metab_controls_2[colnames(up_regulated_metabs.T2)[k],]))
#   up_list_sd_controls.T2[k] = sd(as.numeric(new_metab_controls_2[colnames(up_regulated_metabs.T2)[k],]))
#   frame1[i+j+k,1] = as.numeric(up_list_mean_cases.T2[k])
#   frame1[i+j+k,2] = as.numeric(up_list_sd_cases.T2[k])
#   frame1[i+j+k,3] = as.numeric(up_list_mean_controls.T2[k])
#   frame1[i+j+k,4] = as.numeric(up_list_sd_controls.T2[k])
#   frame1[i+j+k,5] = as.numeric(up_regulated_metabs.T2[1,k])
#   frame1[i+j+k,6] = as.numeric(up_regulated_metabs.T2[2,k])
#   rownames(frame1)[i+j+k] = as.character(intersect(as.character(up.T2), as.character(p_val.T2)))[k]
# }
# down_list_mean_cases.T2 = list()
# down_list_sd_cases.T2 = list()
# down_list_mean_controls.T2 = list()
# down_list_sd_controls.T2 = list()
# for (m in 1:ncol(down_regulated_metabs.T2))
# {
#   down_list_mean_cases.T2[m] = mean(as.numeric(new_metab_cases[colnames(down_regulated_metabs.T2)[m],]))
#   down_list_sd_cases.T2[m] = sd(as.numeric(new_metab_cases[colnames(down_regulated_metabs.T2)[m],]))
#   down_list_mean_controls.T2[m] = mean(as.numeric(new_metab_controls[colnames(down_regulated_metabs.T2)[m],]))
#   down_list_sd_controls.T2[m] = sd(as.numeric(new_metab_controls[colnames(down_regulated_metabs.T2)[m],]))
#   frame1[i+j+k+m,1] = as.numeric(down_list_mean_cases.T2[m])
#   frame1[i+j+k+m,2] = as.numeric(down_list_sd_cases.T2[m])
#   frame1[i+j+k+m,3] = as.numeric(down_list_mean_controls.T2[m])
#   frame1[i+j+k+m,4] = as.numeric(down_list_sd_controls.T2[m])
#   frame1[i+j+k+m,5] = as.numeric(down_regulated_metabs.T2[1,m])
#   frame1[i+j+k+m,6] = as.numeric(down_regulated_metabs.T2[2,m])
#   rownames(frame1)[i+j+k+m] = as.character(intersect(as.character(down.T2), as.character(p_val.T2)))[m]
# }
# colnames(frame1) = c("Mean_Cases", "SD_Cases", "Mean_Controls", "SD_Controls", "p-value", "t-statistic")
# #write.csv(frame1,"table.csv")
# 
# #incorporating HMDB and KEGG IDs into the table
# library(readxl)
# frame_info = read_excel("table.xlsx", sheet=3)
# frame_id = read_excel("IDs.xlsx", sheet=1)
# merged_table = merge(frame_info, frame_id, by="BIOCHEMICAL")
# #write.csv(merged_table,"table_with_IDs_T3.csv")

#Non parametric tests
#Mann whitney U test

w_test.T1 = data.frame()
cname = colnames(tr_new_metab_cases.T1)
count_wpval_1=0
wp_val.T1 = list()
count1=0
w_up.T1 = list()
count2=0
w_down.T1 = list()
for (i in 1:ncol(tr_new_metab_cases.T1))
{
  wtest = wilcox.test(tr_new_metab_cases.T1[,i], tr_new_metab_controls.T1[,i])
  #write(ttest$p.value, file = "t_tests_T1.txt", append=TRUE)
  median_case = median(tr_new_metab_cases.T1[,i])
  median_control = median(tr_new_metab_controls.T1[,i])
  w_test.T1[1,i] = as.numeric(wtest$p.value)
  w_test.T1[2,i] = as.numeric(wtest$statistic)
  if (as.numeric(wtest$p.value)<=0.05)
  {
    count_wpval_1=count_wpval_1+1
    wp_val.T1[count_wpval_1] = cname[i]
    if (median_control<median_case)
    {
      count1=count1+1
      w_up.T1[count1]= cname[i]
    } else if (median_control>median_case)
    {
      count2=count2+1
      w_down.T1[count2]= cname[i]
    }
  }
  
}
colnames(w_test.T1) = cname
rownames(w_test.T1) = c("p-value","w-statistic")
w_up_regulated_metabs.T1 = w_test.T1[,intersect(as.character(w_up.T1), as.character(wp_val.T1))]
w_down_regulated_metabs.T1 = w_test.T1[,intersect(as.character(w_down.T1), as.character(wp_val.T1))]
# write.csv(w_up_regulated_metabs.T1,"w_upregulated_T1.csv")
# write.csv(w_down_regulated_metabs.T1,"w_downregulated_T1.csv")


#T3
w_test.T2 = data.frame()
cname_2 = colnames(tr_new_metab_cases.T2)
count_wpval_2=0
wp_val.T2 = list()
count3=0
w_up.T2 = list()
count4=0
w_down.T2 = list()
for (i in 1:ncol(tr_new_metab_cases.T2))
{
  wtest = wilcox.test(tr_new_metab_cases.T2[,i], tr_new_metab_controls.T2[,i])
  #write(ttest$p.value, file = "t_tests_T1.txt", append=TRUE)
  median_case = median(tr_new_metab_cases.T2[,i])
  median_control = median(tr_new_metab_controls.T2[,i])
  w_test.T2[1,i] = as.numeric(wtest$p.value)
  w_test.T2[2,i] = as.numeric(wtest$statistic)
  if (as.numeric(wtest$p.value)<=0.05)
  {
    count_wpval_2=count_wpval_2+1
    wp_val.T2[count_wpval_2] = cname_2[i]
    if (median_control<median_case)
    {
      count3=count3+1
      w_up.T2[count3]= cname_2[i]
    } else if (median_control>median_case)
    {
      count4=count4+1
      w_down.T2[count4]= cname_2[i]
    }
  }
}
colnames(w_test.T2) = cname_2
rownames(w_test.T2) = c("p-value","w-statistic")
w_up_regulated_metabs.T2 = w_test.T2[,intersect(as.character(w_up.T2), as.character(wp_val.T2))]
w_down_regulated_metabs.T2 = w_test.T2[,intersect(as.character(w_down.T2), as.character(wp_val.T2))]
# write.csv(w_up_regulated_metabs.T2,"w_upregulated_T3.csv")
# write.csv(w_down_regulated_metabs.T2,"w_downregulated_T3.csv")


# #Q values
# library(qvalue)
# qobj_param.T1 = qvalue(t_test.T1[1,])
# qobj_nparam.T1 = qvalue(w_test.T1[1,])
# qobj.T1 = p.adjust(t_test.T1[1,], method = "fdr")
# 
# qobj_param.T2 = qvalue(t_test.T2[1,])
# qobj_nparam.T2 = qvalue(w_test.T2[1,])
# #qobj.T2 = p.adjust(t_test.T2[1,], method = "fdr")
# 
# qobj_log.T1 = qvalue(log_t_test.T1[1,])
# # qobj_log.T1 = p.adjust(log_t_test.T1[1,], method = "fdr")
# 
# qobj_log.T2 = qvalue(log_t_test.T2[1,])
# # qobj_log.T2 = p.adjust(log_t_test.T2[1,], method = "fdr")
# 

#drawing comparisons
x = unique(c(colnames(w_down_regulated_metabs.T1), colnames(w_down_regulated_metabs.T2),
             colnames(w_up_regulated_metabs.T1), colnames(w_up_regulated_metabs.T2)))
x = x[-235]
ladder_frame = data.frame()
transition_patients = read.csv("transition patients.csv")
for (i in 1:length(x))
{
  ladder_frame[i,1] = mean(as.numeric(metab.T1[as.character(x[i]),transition_patients$Pateints]))
  ladder_frame[i,2] = mean(as.numeric(metab.T2[as.character(x[i]),transition_patients$Pateints]))
  if (as.numeric(ladder_frame[i,1])>as.numeric(ladder_frame[i,2]))
  {
    ladder_frame[i,3] = "Going DOWN"
  } else
  {
    ladder_frame[i,3] = "Going UP"
  }
}
colnames(ladder_frame) = c("T1", "T3")
rownames(ladder_frame) = x

id_names = as.data.frame(read_excel("IDs.xlsx"))
rownames(id_names) = id_names$BIOCHEMICAL
ids_matched = id_names[as.character(x),]

#plotting graph of cases vs controls diff in metab levels at T1 and T3
median_cases.T1 = data.frame()
for (i in 1:ncol(tr_new_metab_cases.T1))
{
  median_cases.T1[1,i] = median(tr_new_metab_cases.T1[,i])
}
colnames(median_cases.T1) = colnames(tr_new_metab_cases.T1)

median_cases.T2 = data.frame()
for (i in 1:ncol(tr_new_metab_cases.T2))
{
  median_cases.T2[1,i] = median(tr_new_metab_cases.T2[,i])
}
colnames(median_cases.T2) = colnames(tr_new_metab_cases.T2)

median_controls.T1 = data.frame()
for (i in 1:ncol(tr_new_metab_controls.T1))
{
  median_controls.T1[1,i] = median(tr_new_metab_controls.T1[,i])
}
colnames(median_controls.T1) = colnames(tr_new_metab_controls.T1)

median_controls.T2 = data.frame()
for (i in 1:ncol(tr_new_metab_controls.T2))
{
  median_controls.T2[1,i] = median(tr_new_metab_controls.T2[,i])
}
colnames(median_controls.T2) = colnames(tr_new_metab_controls.T2)

delta.T1 = rbind(median_cases.T1,median_controls.T1)
rownames(delta.T1) = c("Cases", "Controls")

#cases-controls
for (i in 1:ncol(delta.T1))
{
  delta.T1[1,i] = delta.T1[1,i] - delta.T1[2,i]
}
delta.T1 = delta.T1[-2,]
rownames(delta.T1) = "Cases-Controls@T1"

delta.T2 = rbind(median_cases.T2,median_controls.T2)
rownames(delta.T2) = c("Cases", "Controls")
for (i in 1:ncol(delta.T2))
{
  delta.T2[1,i] = delta.T2[1,i] - delta.T2[2,i]
}
delta.T2 = delta.T2[-2,]
rownames(delta.T2) = "Cases-Controls@T3"

# delta.T1 = delta.T1[intersect(colnames(delta.T1),colnames(delta.T2))]
# delta.T2 = delta.T2[intersect(colnames(delta.T1),colnames(delta.T2))]

delta.T1 = delta.T1[as.character(x)]
delta.T2 = delta.T2[as.character(x)]

graph_frame = cbind(t(delta.T1),t(delta.T2))
plot(graph_frame, main = "Mean differences" )
abline(h=0.0,col="red")
abline(v=0.0,col="red")
# raw_file_1 = read.delim("hi.csv", sep = ",")
# t1_data = raw_file_1
# count1=0
# t2_data = raw_file_1
# count2=0
# t3_data = raw_file_1
# count3=0
# for (i in 1:ncol(raw_file_1))
# {
#   if (as.integer(raw_file_1[1,i])==1)
#   {
#     t2_data = t2_data[,-as.integer(i-count2)]
#     t3_data = t3_data[,-as.integer(i-count3)]
#     count2=count2+1
#     count3=count3+1
#   } else if(as.integer(raw_file_1[1,i])==2) {
#     t1_data = t1_data[,-as.integer(i-count1)]
#     t3_data = t3_data[,-as.integer(i-count3)]
#     count1=count1+1
#     count3=count3+1
#   } else if(as.integer(raw_file_1[1,i])==3) {
#     t2_data = t2_data[,-as.integer(i-count2)]
#     t1_data = t1_data[,-as.integer(i-count1)]
#     count2=count2+1
#     count1=count1+1
#   }
# }
# t1_data = t1_data[-1,]
# t2_data = t2_data[-1,]
# t3_data = t3_data[-1,]
# colnames(t1_data)=colnames(t1_data)[1,1:4]
# # write.csv(t1_data[-1,],"T1_new.csv")
# # write.csv(t2_data[-1,],"T2_new.csv")
# # write.csv(t3_data[-1,],"T3_new.csv")
# colnames(t1_data) = colnames(metab.T1)
# colnames(t3_data) = colnames(metab.T2)

#sapply(1:length(t1_data), function(x)substr(t1_data[x],1,4)+str_pad(substr(t1_data,5,),5,side="left", pad = "0"))


#library(glmnet)

#kolmogorov-smirnov test 
#T1: cases vs controls
ks_frame.T1 = data.frame()
ks_sign.T1 = data.frame()
count_ks_pval.T1 = 0
for (i in 1:nrow(new_metab_cases))
{
  ks = ks.test(as.numeric(new_metab_cases[as.character(rownames(new_metab_cases)[i]),]), as.numeric(new_metab_controls[as.character(rownames(new_metab_controls)[i]),]))
  ks_frame.T1[i,1] = ks$p.value
  ks_frame.T1[i,2] = ks$statistic
  ks_frame.T1[i,3] = ks$method
  if (ks$p.value<=0.05)
  {
    count_ks_pval.T1= count_ks_pval.T1+1 
    ks_sign.T1[count_ks_pval.T1,1] = rownames(new_metab_cases)[i]
    ks_sign.T1[count_ks_pval.T1,2] = ks$p.value
    ks_sign.T1[count_ks_pval.T1,3] = ks$statistic
    ks_sign.T1[count_ks_pval.T1,4] = id_names[rownames(new_metab_cases)[i],2]
    ks_sign.T1[count_ks_pval.T1,5] = id_names[rownames(new_metab_cases)[i],3]
  }
}
rownames(ks_frame.T1) = rownames(new_metab_cases)
colnames(ks_frame.T1) = c("P-value", "Statistic", "Method")

#T3: cases vs controls
ks_frame.T2 = data.frame()
ks_sign.T2 = data.frame()
count_ks_pval.T2 = 0
for (i in 1:nrow(new_metab_cases_2))
{
  ks = ks.test(as.numeric(new_metab_cases_2[as.character(rownames(new_metab_cases_2)[i]),]), as.numeric(new_metab_controls_2[as.character(rownames(new_metab_controls_2)[i]),]))
  ks_frame.T2[i,1] = ks$p.value
  ks_frame.T2[i,2] = ks$statistic
  ks_frame.T2[i,3] = ks$method
  if (ks$p.value<=0.05)
  {
    count_ks_pval.T2= count_ks_pval.T2+1 
    ks_sign.T2[count_ks_pval.T2,1] = rownames(new_metab_cases_2)[i]
    ks_sign.T2[count_ks_pval.T2,2] = ks$p.value
    ks_sign.T2[count_ks_pval.T2,3] = ks$statistic
    ks_sign.T2[count_ks_pval.T2,4] = id_names[rownames(new_metab_cases_2)[i],2]
    ks_sign.T2[count_ks_pval.T2,5] = id_names[rownames(new_metab_cases_2)[i],3]
  }
}
rownames(ks_frame.T2) = rownames(new_metab_cases_2)
colnames(ks_frame.T2) = c("P-value", "Statistic", "Method")

#EMD test
#T1: cases vs controls

#T3: cases vs controls

#library(EMDomics)
library(corrplot)

correlations_down.T2 <- cor(t(common_metab_2[as.character(colnames(w_down_regulated_metabs.T2)),colnames(new_metab_cases_2)]))
corrplot(correlations_down.T2, method="circle")

correlations_up.T2 <- cor(t(common_metab_2[as.character(colnames(w_up_regulated_metabs.T2)),colnames(new_metab_cases_2)]))
corrplot(correlations_up.T2, method="circle")

#differential expression variance
deg_tr_new_metab_cases.T1 = tr_new_metab_cases.T1
for (i in 1:ncol(tr_new_metab_cases.T1))
{
  mean_col = mean(tr_new_metab_cases.T1[,i])
  for ( j in 1:nrow(deg_tr_new_metab_cases.T1))
  {
    deg_tr_new_metab_cases.T1[j,i] = abs(deg_tr_new_metab_cases.T1[j,i] - mean_col)
  }
}

deg_tr_new_metab_controls.T1 = tr_new_metab_controls.T1
for (i in 1:ncol(tr_new_metab_controls.T1))
{
  mean_col = mean(tr_new_metab_controls.T1[,i])
  for (j in 1:nrow(deg_tr_new_metab_controls.T1))
  {
    deg_tr_new_metab_controls.T1[j,i] = abs(deg_tr_new_metab_controls.T1[j,i] - mean_col)
  }
}

deg_tr_new_metab_cases.T2 = tr_new_metab_cases.T2
for (i in 1:ncol(tr_new_metab_cases.T2))
{
  mean_col = mean(tr_new_metab_cases.T2[,i])
  for ( j in 1:nrow(deg_tr_new_metab_cases.T2))
  {
    deg_tr_new_metab_cases.T2[j,i] = abs(deg_tr_new_metab_cases.T2[j,i] - mean_col)
  }
}

deg_tr_new_metab_controls.T2 = tr_new_metab_controls.T2
for (i in 1:ncol(tr_new_metab_controls.T2))
{
  mean_col = mean(tr_new_metab_controls.T2[,i])
  for (j in 1:nrow(deg_tr_new_metab_controls.T2))
  {
    deg_tr_new_metab_controls.T2[j,i] = abs(deg_tr_new_metab_controls.T2[j,i] - mean_col)
  }
}

t_test.T1 = data.frame()
cname = colnames(deg_tr_new_metab_cases.T1)
count_pval_1=0
p_val.T1 = list()
count1=0
up.T1 = list()
count2=0
down.T1 = list()
for (i in 1:ncol(deg_tr_new_metab_cases.T1))
{
  ttest = t.test(deg_tr_new_metab_cases.T1[,i], deg_tr_new_metab_controls.T1[,i])
  #write(ttest$p.value, file = "t_tests_T1.txt", append=TRUE)
  t_test.T1[1,i] = as.numeric(ttest$p.value)
  t_test.T1[2,i] = as.numeric(ttest$statistic)
  if (as.numeric(ttest$p.value)<=0.05)
  {
    count_pval_1=count_pval_1+1
    p_val.T1[count_pval_1] = cname[i]
  }
  if (as.numeric(ttest$statistic)>0)
  {
    count1=count1+1
    up.T1[count1]= cname[i]
  } else if (as.numeric(ttest$statistic)<0)
  {
    count2=count2+1
    down.T1[count2]= cname[i]
  }
}
colnames(t_test.T1) = cname
rownames(t_test.T1) = c("p-value","t-statistic")
up_regulated_metabs.T1 = t_test.T1[,intersect(as.character(up.T1), as.character(p_val.T1))]
down_regulated_metabs.T1 = t_test.T1[,intersect(as.character(down.T1), as.character(p_val.T1))]

t_test.T2 = data.frame()
cname_2 = colnames(deg_tr_new_metab_cases.T2)
count_pval_2=0
p_val.T2 = list()
count1=0
up.T2 = list()
count2=0
down.T2 = list()
for (i in 1:ncol(deg_tr_new_metab_cases.T2))
{
  ttest = t.test(deg_tr_new_metab_cases.T2[,i], deg_tr_new_metab_controls.T2[,i])
  #write(ttest$p.value, file = "t_tests_T3.txt", append=TRUE)
  t_test.T2[1,i] = as.numeric(ttest$p.value)
  t_test.T2[2,i] = as.numeric(ttest$statistic)
  if (as.numeric(ttest$p.value)<=0.05)
  {
    count_pval_2=count_pval_2+1
    p_val.T2[count_pval_2] = cname_2[i]
  }
  if (as.numeric(ttest$statistic)>0)
  {
    count1=count1+1
    up.T2[count1]= cname_2[i]
  } else if (as.numeric(ttest$statistic)<0)
  {
    count2=count2+1
    down.T2[count2]= cname_2[i]
  }
}
colnames(t_test.T2) = cname_2
rownames(t_test.T2) = c("p-value","t-statistic")
up_regulated_metabs.T2 = t_test.T2[,intersect(as.character(up.T2), as.character(p_val.T2))]
down_regulated_metabs.T2 = t_test.T2[,intersect(as.character(down.T2), as.character(p_val.T2))]

library(neuralnet)
library(InformationValue)
####Framework####
for (j in 1:6)
{
  frame_data = as.data.frame(read_excel("this.xlsx",sheet=10))
  rownames(frame_data) = frame_data$...1
  frame_data = frame_data[-1]
  
  train = frame_data[1:180,]
  test = frame_data[181:218,]
  # names_frame = c(as.character(colnames(train)))
  # a = as.formula(paste('Output ~ ', paste(names_frame,collapse = '+')))
  #write(j,file = "summary.txt", append=TRUE)
  for (i in 2:10)
  {
    nn = neuralnet(Classifier ~., data=train, hidden = c(81,33,2), act.fct = "logistic",linear.output = FALSE)
    predict = compute(nn,test)
    #predict$net.result
    prob = predict$net.result
    pred = ifelse(prob>0.5,1,0)
    results = cbind(test["Classifier"],pred)
    cm = confusionMatrix(pred,test["Classifier"])
    cm
    plotROC(test["Classifier"],pred)
    # print("problem  in write")
    # write(cm$table,file = "summary.txt", append = TRUE)
    # write(cm$overall, file="summary.txt", append = TRUE)
  }
  # write("----------------------------------------------------", file = "summary.txt", append = TRUE)
}
plot(nn)

####Converting to Modules ####
module_frame_cases.T2=data.frame()

for (i in 1:nrow(tr_new_metab_cases.T2))
{
  for (j in 2:6)
  {
    sheet = read_excel("grouping metabs.xlsx", sheet=j)
    name_list = as.character(sheet$Names)
    module_frame_cases.T2[i,j-1] = 0
    for (k in 1:length(name_list))
    {
      module_frame_cases.T2[i,j-1] = module_frame_cases.T2[i,j-1] + tr_new_metab_cases.T2[i,name_list[k]]
    }
    module_frame_cases.T2[i,j-1] = module_frame_cases.T2[i,j-1]/length(name_list)
  }
}
rownames(module_frame_cases.T2) = rownames(tr_new_metab_cases.T2)


module_frame_controls.T2=data.frame()

for (i in 1:nrow(tr_new_metab_controls.T2))
{
  for (j in 2:6)
  {
    sheet = read_excel("grouping metabs.xlsx", sheet=j)
    name_list = as.character(sheet$Names)
    module_frame_controls.T2[i,j-1] = 0
    for (k in 1:length(name_list))
    {
      module_frame_controls.T2[i,j-1] = module_frame_controls.T2[i,j-1] + tr_new_metab_controls.T2[i,name_list[k]]
    }
    module_frame_controls.T2[i,j-1] = module_frame_controls.T2[i,j-1]/length(name_list)
  }
}
rownames(module_frame_controls.T2) = rownames(tr_new_metab_controls.T2)

for (i in 1:5)
{
  tt = t.test(module_frame_cases.T2[,i], module_frame_controls.T2[,i])
  print(tt$p.value)
}

# training SVM
library(e1071)
library(caret)
library(readxl)
svm_data = as.data.frame(read_excel("this.xlsx", sheet = 8))
rownames(svm_data)=svm_data$...1
svm_data=svm_data[-1]
train = svm_data[1:180,]
test = svm_data[181:218,]
svm_model = svm(Status~.,data=train, type = "C-classification", kernel = "linear")
y_pred = predict(svm_model,test)
results = cbind(test["Status"],y_pred)
cm =confusionMatrix(test["Status"], y_pred)


#training logistic regression
logistic_model = glm(Classifier~.,data=train, family = "binomial")
# summary(logistic_model)
pred = predict(logistic_model,test,type = "response")
library(InformationValue)
opt_cutoff = optimalCutoff(test["Status"],pred)
y_pred = ifelse(pred>opt_cutoff,1,0)
results = cbind(test["Classifier"],y_pred)
#library(caret)
confusionMatrix(y_pred,test["Classifier"], threshold = opt_cutoff) 
plotROC(test["Classifier"], y_pred)
Concordance(test["Classifier"],y_pred)
specificity(test["Classifier"],y_pred, threshold = opt_cutoff)
sensitivity(test["Classifier"],y_pred, threshold = opt_cutoff)

#PCA
princomp = prcomp(svm_data, center=TRUE)
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
ggbiplot(princomp)
