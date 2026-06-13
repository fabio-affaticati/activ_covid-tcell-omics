### ----------------------------
### Helper function to preprocess ANY (CD4/CD8) dataset
### ----------------------------

process_cytometry <- function(df) {
  
  df_processed <- df %>%
    mutate(
      RunNumber = str_match(Name, "^RUN(\\d+)")[,2],
      Middle = str_match(Name, "^RUN\\d+_([^_].*?)_neg|^RUN\\d+_([^_].*?)$")[,2]
    ) %>%
    separate(Middle, into = c("InternalID", "ID2", "ID3", "Timepoint"), 
             sep = "_", remove = FALSE, fill = "right") %>%
    mutate(
      Timepoint = case_when(
        str_detect(ID3, "^T\\d") ~ ID3,
        str_detect(Timepoint, "^T\\d") ~ Timepoint,
        TRUE ~ "T1"
      ),
      Extract = ifelse(!is.na(ID3) & str_detect(ID3, "^T\\d"),
                       paste(ID2),
                       paste(ID2, ID3, sep = "_"))
    ) %>%
    mutate(
      ID2 = ifelse(is.na(ID2), "T1", ID2),
      ID3 = ifelse(is.na(ID3), "T1", ID3),
      Timepoint = ifelse(is.na(Timepoint), "T1", Timepoint),
      Extract = ifelse(is.na(Extract), paste(ID2, ID3, sep = "_"), Extract)
    ) %>%
    select(-Middle)
  
  df_final <- df_processed %>%
    mutate(
      sample_id = case_when(
        is.na(ID3) ~ ID2,
        ID3 %in% c("T1","T2","T3","T4") ~ ID2,
        TRUE ~ paste0(ID2, "_", ID3)
      )
    ) %>%
    select(sample_id, Timepoint, Metacluster, percentage, RunNumber) %>%
    mutate(sample_id = str_replace(sample_id, "_", "")) %>%
    drop_na() %>%
    mutate(
      sample_id = if_else(str_starts(sample_id, "VZV"), paste0("ERC", sample_id), sample_id),
      sample_id = if_else(str_starts(sample_id, "VPK"), paste0("ERC", sample_id), sample_id),
      sample_id = if_else(str_starts(sample_id, "DCP"), paste0("ERC", sample_id), sample_id),
      sample_id = if_else(str_starts(sample_id, "PT"), paste0("Covid", sample_id), sample_id),
      sample_id = if_else(str_starts(sample_id, "NICU"), paste0("ERC", sample_id), sample_id)
    )
  
  # -----------------------------
  # 🔥 NEW RULES
  # -----------------------------
  
  df_final <- df_final %>%
    # drop VZV/TYF samples at T2 or T3
    filter(!(str_detect(sample_id, "VZV|TYF") & Timepoint %in% c("T2", "T3"))) %>%
    
    # change T4 → T2 for VZV/TYF
    mutate(
      Timepoint = if_else(
        str_detect(sample_id, "VZV|TYF") & Timepoint == "T4",
        "T2",
        Timepoint
      )
    )
  
  df_final
}



### ----------------------------
### Helper: build heatmap matrix
### ----------------------------

make_matrix <- function(df_tp) {
  df_tp %>%
    select(sample_id, Metacluster, percentage) %>%
    pivot_wider(names_from = sample_id, values_from = percentage) %>%
    arrange(Metacluster) %>%
    column_to_rownames("Metacluster") %>%
    as.matrix()
}


### ----------------------------
### Helper: Remove entire samples if ANY row has percentage > 50
### ----------------------------

remove_samples_percentage <- function(df_tp) {
  
  bad_samples <- df_tp %>%
    group_by(sample_id) %>%
    summarize(has_big = any(percentage > 50)) %>%
    filter(has_big) %>%
    pull(sample_id)
  
  df_tp %>% filter(!sample_id %in% bad_samples)
}


### ----------------------------
### Helper: Add missing metacluster rows
### ----------------------------

library(dplyr)
library(tidyr)

enforce_55_metaclusters <- function(df, total_mc = 55) {
  # ensure Metacluster is integer (or numeric)
  df <- df %>% mutate(Metacluster = as.integer(Metacluster))
  
  # columns that must not be auto-filled with cluster numbers
  # (we'll treat percentage specially and fill other metadata with the group's first non-NA)
  meta_cols <- setdiff(names(df), c("sample_id", "Timepoint", "Metacluster", "percentage"))
  
  df_fixed <- df %>%
    # group per sample_id x Timepoint
    group_by(sample_id, Timepoint) %>%
    # ensure a row for every Metacluster 1:total_mc
    complete(Metacluster = seq_len(total_mc)) %>%
    # percentage: replace NA (newly created rows) with 0
    mutate(percentage = replace_na(percentage, 0)) %>%
    # for all other metadata columns, fill missing values with the first non-NA value in the group
    # if a whole group had NA for a meta column (unlikely), it stays NA
    mutate(across(all_of(meta_cols), ~ {
      # find first non-NA in this group
      first_val <- { tmp <- .; tmp_non_na <- tmp[!is.na(tmp)]; if (length(tmp_non_na)) tmp_non_na[1] else NA }
      # replace NAs with that first_val
      replace_na(., first_val)
    })) %>%
    ungroup() %>%
    arrange(sample_id, Timepoint, Metacluster)
  
  return(df_fixed)
}

