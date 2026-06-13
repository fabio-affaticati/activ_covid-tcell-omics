### ----------------------------
### Helper function to preprocess CD4/CD8 cytometry tables
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
  
  # Validation samples use different visit labels; keep comparable baseline/follow-up visits.
  df_final <- df_final %>%
    filter(!(str_detect(sample_id, "VZV|TYF") & Timepoint %in% c("T2", "T3"))) %>%
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
### Helper: remove samples with one metacluster above 50 percent
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
### Helper: complete the expected FlowSOM metacluster grid
### ----------------------------

library(dplyr)
library(tidyr)

enforce_55_metaclusters <- function(df, total_mc = 55) {
  # Metacluster IDs must be numeric before completing the expected range.
  df <- df %>% mutate(Metacluster = as.integer(Metacluster))
  
  # Preserve per-sample metadata while adding zero rows for absent FlowSOM metaclusters.
  meta_cols <- setdiff(names(df), c("sample_id", "Timepoint", "Metacluster", "percentage"))
  
  df_fixed <- df %>%
    group_by(sample_id, Timepoint) %>%
    complete(Metacluster = seq_len(total_mc)) %>%
    mutate(percentage = replace_na(percentage, 0)) %>%
    mutate(across(all_of(meta_cols), ~ {
      first_val <- { tmp <- .; tmp_non_na <- tmp[!is.na(tmp)]; if (length(tmp_non_na)) tmp_non_na[1] else NA }
      replace_na(., first_val)
    })) %>%
    ungroup() %>%
    arrange(sample_id, Timepoint, Metacluster)
  
  return(df_fixed)
}
