
#' Dupilumab codes
get_dupilumab_codes <- function(drug_name='dupilumab',
                                exclude_drug='',
                                drug_codes=c('1876376','1876399','2375326','1876400'),
                                concept_ancestor=vocabulary_tbl('concept_ancestor'),
                                drug_tbl=cdm_tbl('drug_exposure'),
                                concept_tbl=vocabulary_tbl('concept')
) {
  
  rxnorm_codes_first <- concept_tbl %>%
    filter(concept_code %in% drug_codes) %>%
    filter(vocabulary_id %in% c('RxNorm','RxNorm Extension')) %>%
    distinct(concept_id)
  
  rxnorm_codes <- rxnorm_codes_first %>%
    left_join(concept_ancestor, by=c('concept_id'='ancestor_concept_id')) %>%
    distinct(descendant_concept_id) %>%
    rename(concept_id=descendant_concept_id) %>%
    dplyr::union(rxnorm_codes_first) %>%
    distinct(concept_id) %>%
    left_join(vocabulary_tbl('concept'), by='concept_id') %>%
    mutate(found_via_rxnav=1L) %>%
    compute_new()
  
  all_codes <- concept_tbl %>%
    anti_join(rxnorm_codes, by='concept_id') %>%
    filter(str_detect(lower(concept_name),drug_name)) %>%
    mutate(found_via_rxnav=0L) %>%
    dplyr::union(rxnorm_codes) %>%
    filter(!is.na(concept_id)) %>%
    compute_new()
  
  if (nchar(exclude_drug)>0L) {
    all_codes <- all_codes %>%
      filter(!str_detect(lower(concept_name),exclude_drug))
    
    rxnorm_codes <- rxnorm_codes %>%
      filter(!str_detect(lower(concept_name),exclude_drug))
  }
  
  # all_drug <- drug_tbl %>%
  #   inner_join(rxnorm_codes %>% distinct(concept_id, found_via_rxnav), by=c('drug_concept_id'='concept_id')) %>%
  #   compute_new()
  
  full <- list()
  full$rxnorm_codes <- rxnorm_codes
  #full$all_drug <- all_drug
  
  return(full)
  
}





#' @param cohort_tbl tbl with person_id and index_date
#' @param dx_cond tbl with all rows fitting a condition in condition_occurrence
#' @param diff_bw_codes requirement for number of days between earliest and latest condition code
#' 
get_pats_w_dx <- function(cohort_tbl=rslt$dup_dates,
                          dx_cond=results_tbl('persistent_asthma_dx'),
                          diff_bw_codes=180L) {
  
  has_2codes_name <- paste0('has2codes_',as.character(diff_bw_codes),'diff')
  
  combined <- dx_cond %>%
    inner_join(cohort_tbl, by=c('person_id')) %>%
    filter(as.Date(condition_start_date)<=as.Date(index_date)) %>%
    distinct(person_id, index_date, condition_start_date) %>%
    group_by(person_id, index_date) %>%
    summarise(min_dx_date=min(condition_start_date),
              max_dx_date=max(condition_start_date),
              n_dx=n()) %>%
    ungroup() %>%
    mutate(dx_date_diff=max_dx_date-min_dx_date) %>%
    mutate(!!has_2codes_name:=case_when(dx_date_diff>=diff_bw_codes ~ 1L,
                                        TRUE ~ 0L)) %>%
    compute_new(indexes=c('person_id'))
  
}



get_prior_followup_criteria <- function(cohort_tbl=results_tbl('asthma_dx_dup'),
                                        visit_tbl=cdm_tbl('visit_occurrence'),
                                        lookback_max=365L,
                                        lookback_min=7,
                                        lookforward_min=60L,
                                        lookforward_max=365L,
                                        visit_types=c(9201L, 2000001532L, 42898160L, 44814710L,
                                                      9202L, 2000000469L,
                                                      9203L, 2000000048L,
                                                      2000000088L,
                                                      581399,
                                                      44814711L, 2000000104L, 44814653L, 44814649L, 44814650L
                                        )) {
  
  #' Lookback criteria
  pre_criteria <- cohort_tbl %>%
    select(person_id, index_date) %>%
    left_join(visit_tbl, by='person_id') %>%
    filter(visit_concept_id %in% visit_types) %>%
    filter(as.Date(visit_start_date) >= (as.Date(index_date) - lookback_max) &&
             as.Date(visit_start_date) <= (as.Date(index_date) - lookback_min)) %>%
    distinct(person_id) %>%
    mutate(lookback_pre=1L) %>%
    compute_new(indexes=c('person_id'))
  
  post_criteria <- cohort_tbl %>%
    select(person_id, index_date) %>%
    left_join(visit_tbl, by='person_id') %>%
    filter(visit_concept_id %in% visit_types) %>%
    filter(as.Date(visit_start_date) >= (as.Date(index_date) + lookforward_min) &&
             as.Date(visit_start_date) <= (as.Date(index_date) + lookforward_max)) %>%
    distinct(person_id) %>%
    mutate(lookback_post=1L) %>%
    compute_new(indexes=c('person_id'))
  
  pre_post_criteria <- cohort_tbl %>%
    select(person_id, index_date) %>%
    left_join(pre_criteria, by='person_id') %>%
    left_join(post_criteria, by='person_id') %>%
    mutate(lookback_pre_post=case_when((lookback_pre==1L & lookback_post==1L) ~ 1L,
                                       TRUE ~ 0L)) %>%
    compute_new(indexes=c('person_id'))
}







#' Get attrition for dupilumab study
get_attrition <- function(person_tbl=cdm_tbl('person'),
                          asthma_tbl=results_tbl('persistent_asthma_dx'),
                          aderm_tbl=results_tbl('atopic_dermatitis_dx'),
                          asthma_dup=results_tbl('asthma_dx_dup'),
                          aderm_dup=results_tbl('aderm_dx_dup'),
                          study_start_date=as.Date('2018-10-01'),
                          study_end_date=as.Date('2022-06-01'),
                          asthma_dup_lblf=results_tbl('asthma_dx_dup_lblf'),
                          aderm_dup_lblf=results_tbl('aderm_dx_dup_lblf'),
                          age_limit=18L
) {
  
  total_ct <- person_tbl %>% filter(((study_start_date-as.Date(birth_date))/365.25) < age_limit) %>%
    distinct(person_id) %>% count() %>% pull()
  
  #' asthma patients
  asthma <- asthma_tbl %>%
    filter(as.Date(condition_start_date)<=as.Date(study_end_date)) %>%
    group_by(person_id) %>%
    summarise(min_cond_date=min(condition_start_date)) %>%
    ungroup() %>%
    left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
    mutate(age = (min_cond_date - birth_date) / 365.25) %>%
    filter(age < age_limit)
  
  asthma_ct <- asthma %>% count() %>% pull()
  
  #' atopic dermatitis patients
  aderm <- aderm_tbl %>%
    filter(as.Date(condition_start_date)<=as.Date(study_end_date)) %>%
    group_by(person_id) %>%
    summarise(min_cond_date=min(condition_start_date)) %>%
    ungroup() %>%
    left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
    mutate(age = (min_cond_date - birth_date) / 365.25) %>%
    filter(age < age_limit)
  
  aderm_ct <- aderm %>% count() %>% pull()
  
  #' asthma or atopic dermatitis patients
  asthma_aderm_ct <- asthma %>%
    distinct(person_id) %>%
    dplyr::union(aderm %>% distinct(person_id)) %>%
    distinct(person_id) %>%
    count() %>% pull()
  
  #' asthma with dupilumab
  asthma_dup_tot <- asthma_dup %>%
    filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
    mutate(age = (index_date - birth_date) / 365.25) %>%
    filter(age < age_limit)
  
  asthma_dup_ct <- asthma_dup_tot %>% count() %>% pull()
  
  #' atopic dermatitis with dupilumab
  aderm_dup_tot <- aderm_dup %>%
    filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
    mutate(age = (index_date - birth_date) / 365.25) %>%
    filter(age < age_limit)
  
  aderm_dup_ct <- aderm_dup_tot %>% count() %>% pull()
  
  #' asthma or atopic dermatitis with dupilumab
  asthma_aderm_dup_ct <- asthma_dup %>%
    dplyr::union(aderm_dup) %>%
    group_by(person_id) %>%
    summarise(index_date=min(index_date)) %>%
    filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
    mutate(age = (index_date - birth_date) / 365.25) %>%
    filter(age < age_limit) %>%
    distinct(person_id) %>%
    count() %>% pull()
  
  #' asthma with dupilumab and lblf
  asthma_dup_lblf_tot <- asthma_dup_lblf %>%
    filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
    mutate(age = (index_date - birth_date) / 365.25) %>%
    filter(age < age_limit)
  
  asthma_dup_lblf_ct <- asthma_dup_lblf_tot %>% count() %>% pull()
  
  #' atopic dermatitis with dupilumab and lblf
  aderm_dup_lblf_tot <- aderm_dup_lblf %>%
    filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
    mutate(age = (index_date - birth_date) / 365.25) %>%
    filter(age < age_limit)
  
  aderm_dup_lblf_ct <- aderm_dup_lblf_tot %>% count() %>% pull()
  
  #' asthma or atopic dermatitis with dupilumab and lblf
  asthma_aderm_dup_lblf_ct <- asthma_dup_lblf %>%
    dplyr::union(aderm_dup_lblf) %>%
    group_by(person_id) %>%
    summarise(index_date=min(index_date)) %>%
    filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
    mutate(age = (index_date - birth_date) / 365.25) %>%
    filter(age < age_limit) %>%
    distinct(person_id) %>%
    count() %>% pull()
  
  desc <- c('01. Patients aged <18 on 10/01/2018 with available EHR data in network',
            '02. Patients aged <18 diagnosed with persistent asthma',
            '03. Patients aged <18 diagnosed with atopic dermatitis',
            '04. Patients aged <18 diagnosed with persistent asthma or atopic dermatitis',
            '05. Patients aged <18 diagnosed with persistent asthma on or before Dupilumab rx',
            '06. Patients aged <18 diagnosed with atopic dermatitis on or before Dupilumab rx',
            '07. Patients aged <18 diagnosed with persistent asthma or atopic dermatitis on or before Dupilumab rx',
            '09. Patients aged <18 diagnosed with persistent asthma on or before Dupilumab rx, with lb/lf critera',
            '09. Patients aged <18 diagnosed with atopic dermatitis on or before Dupilumab rx, with lb/lf criteria',
            '10. Patients aged <18 diagnosed with persistent asthma or atopic dermatitis on or before Dupilumab rx, with lb/lf criteria'
  )
  
  cts_all <- c(total_ct, asthma_ct, aderm_ct, asthma_aderm_ct,
               asthma_dup_ct, aderm_dup_ct, asthma_aderm_dup_ct,
               asthma_dup_lblf_ct, aderm_dup_lblf_ct, asthma_aderm_dup_lblf_ct)
  
  final_att <- cbind(desc, cts_all) %>% as_tibble()
  
}


#' #' Get attrition for dupilumab study, only allowing 2 diagnostic codes at least 6 months apart
#' get_attrition_2codes <- function(person_tbl=cdm_tbl('person'),
#'                           asthma_tbl=results_tbl('persistent_asthma_dx'),
#'                           aderm_tbl=results_tbl('atopic_dermatitis_dx'),
#'                           asthma_dup=results_tbl('asthma_dx_dup'),
#'                           aderm_dup=results_tbl('aderm_dx_dup'),
#'                           study_start_date=as.Date('2018-10-01'),
#'                           study_end_date=as.Date('2022-06-01'),
#'                           asthma_dup_lblf=results_tbl('asthma_dx_dup_lblf'),
#'                           aderm_dup_lblf=results_tbl('aderm_dx_dup_lblf'),
#'                           age_limit=18L
#' ) {
#'   
#'   total_ct <- person_tbl %>% filter(((study_start_date-as.Date(birth_date))/365.25) < age_limit) %>%
#'     distinct(person_id) %>% count() %>% pull()
#'   
#'   #' asthma patients
#'   asthma <- asthma_tbl %>%
#'     filter(as.Date(condition_start_date)<=as.Date(study_end_date)) %>%
#'     group_by(person_id) %>%
#'     summarise(min_cond_date=min(condition_start_date),
#'               max_cond_date=max(condition_start_date)) %>%
#'     ungroup() %>%
#'     mutate(date_diff=max_cond_date-min_cond_date) %>%
#'     filter(date_diff>=180L) %>%
#'     left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
#'     mutate(age = (min_cond_date - birth_date) / 365.25) %>%
#'     filter(age < age_limit)
#'   
#'   asthma_ct <- asthma %>% count() %>% pull()
#'   
#'   #' atopic dermatitis patients
#'   aderm <- aderm_tbl %>%
#'     filter(as.Date(condition_start_date)<=as.Date(study_end_date)) %>%
#'     group_by(person_id) %>%
#'     summarise(min_cond_date=min(condition_start_date),
#'               max_cond_date=max(condition_start_date)) %>%
#'     ungroup() %>%
#'     mutate(date_diff=max_cond_date-min_cond_date) %>%
#'     filter(date_diff>=180L) %>%
#'     left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
#'     mutate(age = (min_cond_date - birth_date) / 365.25) %>%
#'     filter(age < age_limit)
#'   
#'   aderm_ct <- aderm %>% count() %>% pull()
#'   
#'   #' asthma or atopic dermatitis patients
#'   asthma_aderm_ct <- asthma %>%
#'     dplyr::union(aderm) %>%
#'     group_by(person_id) %>%
#'     summarise(index_date=min(min_cond_date)) %>%
#'     left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
#'     mutate(age = (index_date - birth_date) / 365.25) %>%
#'     filter(age < age_limit) %>%
#'     distinct(person_id) %>%
#'     count() %>% pull()
#'   
#'   #' asthma with dupilumab
#'   asthma_dup_tot <- asthma_dup %>%
#'     filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
#'     filter(has2codes_180diff==1L) %>%
#'     left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
#'     mutate(age = (index_date - birth_date) / 365.25) %>%
#'     filter(age < age_limit)
#'   
#'   asthma_dup_ct <- asthma_dup_tot %>% count() %>% pull()
#'   
#'   #' atopic dermatitis with dupilumab
#'   aderm_dup_tot <- aderm_dup %>%
#'     filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
#'     filter(has2codes_180diff==1L) %>%
#'     left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
#'     mutate(age = (index_date - birth_date) / 365.25) %>%
#'     filter(age < age_limit)
#'   
#'   aderm_dup_ct <- aderm_dup_tot %>% count() %>% pull()
#'   
#'   #' asthma or atopic dermatitis with dupilumab
#'   asthma_aderm_dup_ct <- asthma_dup_tot %>%
#'     distinct(person_id) %>%
#'     dplyr::union(aderm_dup_tot %>% distinct(person_id)) %>%
#'     distinct(person_id) %>%
#'     count() %>% pull()
#'   
#'   #' asthma with dupilumab and lblf
#'   asthma_dup_lblf_tot <- asthma_dup_lblf %>%
#'     inner_join(asthma_dup_tot %>% distinct(person_id), by='person_id') %>%
#'     filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
#'     left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
#'     mutate(age = (index_date - birth_date) / 365.25) %>%
#'     filter(age < age_limit)
#'   
#'   asthma_dup_lblf_ct <- asthma_dup_lblf_tot %>% count() %>% pull()
#'   
#'   #' atopic dermatitis with dupilumab and lblf
#'   aderm_dup_lblf_tot <- aderm_dup_lblf %>%
#'     inner_join(aderm_dup_tot %>% distinct(person_id), by='person_id') %>%
#'     filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
#'     left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
#'     mutate(age = (index_date - birth_date) / 365.25) %>%
#'     filter(age < age_limit)
#'   
#'   aderm_dup_lblf_ct <- aderm_dup_lblf_tot %>% count() %>% pull()
#'   
#'   #' asthma or atopic dermatitis with dupilumab and lblf
#'   asthma_aderm_dup_lblf_ct <- asthma_dup_lblf %>%
#'     left_join(asthma_dup %>% distinct(person_id, has2codes_180diff), by='person_id') %>%
#'     filter(has2codes_180diff==1L) %>%
#'     dplyr::union(aderm_dup_lblf %>%
#'                    left_join(aderm_dup %>% distinct(person_id, has2codes_180diff), by='person_id') %>%
#'                    filter(has2codes_180diff==1L)) %>%
#'     group_by(person_id) %>%
#'     summarise(index_date=min(index_date)) %>%
#'     filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
#'     left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
#'     mutate(age = (index_date - birth_date) / 365.25) %>%
#'     filter(age < age_limit) %>%
#'     distinct(person_id) %>%
#'     count() %>% pull()
#'   
#'   desc <- c('01. Patients aged <18 on 10/01/2018 with available EHR data in network',
#'             '02. Patients aged <18 diagnosed with persistent asthma',
#'             '03. Patients aged <18 diagnosed with atopic dermatitis',
#'             '04. Patients aged <18 diagnosed with persistent asthma or atopic dermatitis',
#'             '05. Patients aged <18 diagnosed with persistent asthma on or before Dupilumab rx',
#'             '06. Patients aged <18 diagnosed with atopic dermatitis on or before Dupilumab rx',
#'             '07. Patients aged <18 diagnosed with persistent asthma or atopic dermatitis on or before Dupilumab rx',
#'             '09. Patients aged <18 diagnosed with persistent asthma on or before Dupilumab rx, with lb/lf critera',
#'             '09. Patients aged <18 diagnosed with atopic dermatitis on or before Dupilumab rx, with lb/lf criteria',
#'             '10. Patients aged <18 diagnosed with persistent asthma or atopic dermatitis on or before Dupilumab rx, with lb/lf criteria'
#'   )
#'   
#'   cts_all <- c(total_ct, asthma_ct, aderm_ct, asthma_aderm_ct,
#'                asthma_dup_ct, aderm_dup_ct, asthma_aderm_dup_ct,
#'                asthma_dup_lblf_ct, aderm_dup_lblf_ct, asthma_aderm_dup_lblf_ct)
#'   
#'   final_att <- cbind(desc, cts_all) %>% as_tibble()
#'   
#' }


#' Get attrition for dupilumab study (using 2 codes for computable phenotype definition)
get_attrition_2codes <- function(asthma_tbl=results_tbl('persistent_asthma_dx'),
                                 aderm_tbl=results_tbl('atopic_dermatitis_dx'),
                                 asthma_or_aderm=results_tbl('pasthma_aderm_2codes'),
                                 
                                 asthma_w_dup=results_tbl('asthma_2dx_dup_lblf'),
                                 aderm_w_dup=results_tbl('aderm_2dx_dup_lblf'),
                                 
                                 atopic_ctrl_lblf = results_tbl('pasthma_aderm_2dx_lblf'),
                                 
                                 non_atopic=results_tbl('no_disease'),
                                 non_atopic_lblf=results_tbl('no_diasease_lblf'),
                                 
                                 final_cht=results_tbl('pasthma_aderm_2dx_demog'),
                                 
                                 study_start_date=as.Date('2018-10-01'),
                                 study_end_date=as.Date('2022-06-01'),
                                 age_limit=18L,
                                 person_tbl=cdm_tbl('person')
) {
  
  #' total count
  total_ct <- person_tbl %>% filter(((study_start_date-as.Date(birth_date))/365.25) < age_limit) %>%
    distinct(person_id) %>% count() %>% pull()
  
  #' asthma patients
  asthma <- asthma_tbl %>%
    filter(as.Date(condition_start_date)<=as.Date(study_end_date)) %>%
    group_by(person_id) %>%
    summarise(min_cond_date=min(condition_start_date)) %>%
    ungroup() %>%
    left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
    mutate(age = (min_cond_date - birth_date) / 365.25) %>%
    filter(age < age_limit)
  
  asthma_ct <- asthma %>% count() %>% pull()
  
  #' atopic dermatitis patients
  aderm <- aderm_tbl %>%
    filter(as.Date(condition_start_date)<=as.Date(study_end_date)) %>%
    group_by(person_id) %>%
    summarise(min_cond_date=min(condition_start_date)) %>%
    ungroup() %>%
    left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
    mutate(age = (min_cond_date - birth_date) / 365.25) %>%
    filter(age < age_limit)
  
  aderm_ct <- aderm %>% count() %>% pull()
  
  #' asthma or atopic dermatitis patients
  asthma_aderm_ct <- asthma %>%
    dplyr::union(aderm) %>%
    distinct(person_id) %>% count() %>% pull()
  
  #' asthma or atopic dermatitis patients with 2 codes (computable phenotype)
  asthma_aderm_cp_ct <- asthma_or_aderm %>%
    distinct(person_id) %>%
    count() %>% pull()
  
  ############################################################################
  #' asthma or atopic dermatitis with dupilumab
  asthma_aderm_dup_ct <- asthma_w_dup %>%
    dplyr::union(aderm_w_dup) %>%
    group_by(person_id) %>%
    summarise(index_date=min(index_date)) %>%
    ungroup() %>%
    filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
    mutate(age = (index_date - birth_date) / 365.25) %>%
    filter(age < age_limit) %>%
    distinct(person_id) %>%
    count() %>% pull()
  
  #' asthma or atopic dermatitis with dupilumab and lb
  # asthma_aderm_dup_lb_ct <- asthma_w_dup %>%
  #   filter(lookback_pre==1L) %>%
  #   dplyr::union(aderm_w_dup %>%
  #                  filter(lookback_pre==1L)) %>%
  #   group_by(person_id) %>%
  #   summarise(index_date=min(index_date)) %>%
  #   filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
  #   left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
  #   mutate(age = (index_date - birth_date) / 365.25) %>%
  #   filter(age < age_limit) %>%
  #   distinct(person_id) %>%
  #   count() %>% pull()
  
  asthma_aderm_dup_lb_ct <- final_cht %>%
    filter(cht=='1_atopic_dup') %>%
    distinct(person_id) %>%
    count() %>% pull()
  
  #' asthma or atopic dermatitis with dupilumab and lblf
  asthma_aderm_dup_lblf_ct <- asthma_w_dup %>%
    filter(lookback_pre_post==1L) %>%
    dplyr::union(aderm_w_dup %>%
                   filter(lookback_pre_post==1L)) %>%
    group_by(person_id) %>%
    summarise(index_date=min(index_date)) %>%
    filter(as.Date(index_date)<=as.Date(study_end_date)) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by='person_id') %>%
    mutate(age = (index_date - birth_date) / 365.25) %>%
    filter(age < age_limit) %>%
    distinct(person_id) %>%
    inner_join(final_cht %>% filter(cht=='1_atopic_dup') %>% distinct(person_id), by=c('person_id')) %>%
    count() %>% pull()
  
  ############################################################################
  #' asthma or atopic dermatitis without dupilumab
  atopic_ctrl_ct <- atopic_ctrl_lblf %>%
    distinct(person_id) %>%
    count() %>% pull()
  
  #' asthma or atopic dermatitis without dupilumab and lb
  # atopic_ctrl_lb_ct <- atopic_ctrl_lblf %>%
  #   anti_join(asthma_dup_lblf, by=c('person_id')) %>%
  #   anti_join(aderm_dup_lblf, by=c('person_id')) %>%
  #   filter(lookback_pre==1L) %>%
  #   distinct(person_id) %>%
  #   count() %>% pull()
  
  atopic_ctrl_lb_ct <- final_cht %>%
    filter(cht=='2_atopic_controls') %>%
    distinct(person_id) %>%
    count() %>% pull()
  
  #' asthma or atopic dermatitis without dupilumab and lblf
  atopic_ctrl_lblf_ct <- atopic_ctrl_lblf %>%
    anti_join(asthma_dup_lblf, by=c('person_id')) %>%
    anti_join(aderm_dup_lblf, by=c('person_id')) %>%
    inner_join(final_cht %>% filter(cht=='2_atopic_controls') %>% distinct(person_id), by=c('person_id')) %>%
    filter(lookback_pre_post==1L) %>%
    distinct(person_id) %>%
    count() %>% pull()
  
  ############################################################################
  #' non-atopic controls
  non_atopic_ct <- non_atopic %>%
    anti_join(asthma_tbl, by=c('person_id')) %>%
    anti_join(aderm_tbl, by=c('person_id')) %>%
    distinct(person_id) %>% count() %>% pull()
  
  #' non-atopic controls with lb
  # non_atopic_lb_ct <- non_atopic_lblf %>%
  #   filter(lookback_pre==1L) %>%
  #   anti_join(asthma_tbl, by=c('person_id')) %>%
  #   anti_join(aderm_tbl, by=c('person_id')) %>%
  #   distinct(person_id) %>% count() %>% pull()
  
  non_atopic_lb_ct <- final_cht %>%
    filter(cht=='3_non_atopic_controls') %>%
    distinct(person_id) %>%
    count() %>% pull()
  
  #' non-atopic controls with lblf
  non_atopic_lblf_ct <- non_atopic_lblf %>%
    filter(lookback_pre_post==1L) %>%
    anti_join(asthma_tbl, by=c('person_id')) %>%
    anti_join(aderm_tbl, by=c('person_id')) %>%
    inner_join(final_cht %>% filter(cht=='3_non_atopic_controls') %>% distinct(person_id), by=c('person_id')) %>%
    distinct(person_id) %>% count() %>% pull()
  
  desc <- c('01. Patients aged <18 years at the study start date with a visit within the study period',
            '02.A (From Step 01) At least 1 diagnosis code for persistent asthma',
            '02.B (From Step 01) At least 1 diagnosis code for atopic dermatitis',
            '02. (From Step 01) At least 1 diagnosis code for atopic disease',
            '03. (From Step 02) Meets atopic disease computable phenotype (2 dx at least 6 months apart)',
            
            '04. (From Step 03) Any Dupilumab exposure during the study period',
            '05. (From Step 04) Has prior visit requirement (7 to 365 days before cohort entry date)',
            '06. (From Step 05) Has post visit requirement (60 to 365 days after cohort entry date)',
            
            '07. (From Step 3) No Dupilumab exposure during the study period',
            '08. (From Step 7) Has prior visit requirement (7 to 365 days before cohort entry date)',
            '09. (From Step 8) Has post visit requirement (60 to 365 days after cohort entry date)',
            
            '10. (From Step 1) No atopic px or atopic dx on or before cohort entry, or dupilumab on or before end of follow-up',
            '11. (From Step 10) Has prior visit requirement (7 to 365 days before cohort entry date)',
            '12. (From Step 11) Has post visit requirement (60 to 365 days after cohort entry date)')
  
  cts_all <- c(total_ct, asthma_ct, aderm_ct, asthma_aderm_ct, asthma_aderm_cp_ct,
               asthma_aderm_dup_ct, asthma_aderm_dup_lb_ct, asthma_aderm_dup_lblf_ct,
               atopic_ctrl_ct, atopic_ctrl_lb_ct, atopic_ctrl_lblf_ct,
               non_atopic_ct, non_atopic_lb_ct, non_atopic_lblf_ct)
  
  final_att <- cbind(desc, cts_all) %>% as_tibble()
  
}



#' Get a random visit_occurrence_id/visit date within the study period
get_random_visit <- function(cohort=rslt$pasthma_aderm,
                             visit_tbl=cdm_tbl('visit_occurrence'),
                             study_start_date=as.Date('2018-10-01'),
                             study_end_date=as.Date('2022-06-01')) {
  
  set.seed(18)
  
  rslt$visits <- visit_tbl %>%
    inner_join(cohort, by=c('person_id')) %>%
    filter(visit_start_date >= study_start_date) %>%
    filter(visit_start_date <= study_end_date) %>%
    group_by(person_id) %>%
    slice_sample(visit_occurrence_id, n=1) %>%
    ungroup() %>%
    compute_new(indexes=c('person_id'))
  
  rslt$final <- cohort %>%
    left_join(rslt$visits %>% select(person_id, visit_occurrence_id, visit_start_date), by=c('person_id')) %>%
    mutate(has_visit_in_pd=case_when(!is.na(visit_start_date) ~ 1L,
                                     TRUE ~ 0L)) %>%
    compute_new(indexes=c('person_id'))
  
}


#' Get condition codes (commented out: and only retain records occurring at the earliest possible date?)
get_min_cond <- function(cond_tbl=cdm_tbl('condition_occurrence'),
                         codeset=rslt$eczema_codeset) {
  
  conds <- cond_tbl %>%
    inner_join(codeset %>% distinct(concept_id), by=c('condition_concept_id'='concept_id')) %>%
    mutate(concept='condition_concept_id') %>%
    compute_new(indexes=c('condition_occurrence_id'))
  
  cond_src <- cond_tbl %>%
    anti_join(conds, by=c('condition_occurrence_id')) %>%
    inner_join(codeset %>% distinct(concept_id), by=c('condition_source_concept_id'='concept_id')) %>%
    mutate(concept='condition_source_concept_id') %>%
    dplyr::union(conds) %>%
    group_by(person_id) %>%
    summarise(min_date=min(condition_start_date)) %>%
    ungroup() %>%
    compute_new(indexes=c('person_id'))
  
}



#' Get drug exposure codes (commented out: and only retain records occurring at the earliest possible date?)
get_min_drug <- function(drug_tbl=cdm_tbl('drug_exposure'),
                         codeset=rslt$asthma_aderm_px_codes) {
  
  drugs <- drug_tbl %>%
    inner_join(codeset %>% distinct(concept_id), by=c('drug_concept_id'='concept_id')) %>%
    mutate(concept='drug_concept_id') %>%
    compute_new(indexes=c('drug_exposure_id'))
  
  drug_src <- drug_tbl %>%
    anti_join(drugs, by=c('drug_exposure_id')) %>%
    inner_join(codeset %>% distinct(concept_id), by=c('drug_source_concept_id'='concept_id')) %>%
    mutate(concept='drug_source_concept_id') %>%
    dplyr::union(drugs) %>%
    group_by(person_id) %>%
    summarise(min_date=min(drug_exposure_start_date)) %>%
    ungroup() %>%
    compute_new(indexes=c('person_id'))
  
}


find_condition_codes <- function(codeset=codeset,
                                 cond_tbl=cdm_tbl('condition_occurrence')) {
  
  cond_info <- cond_tbl %>%
    inner_join(codeset %>% select(concept_id), by=c('condition_concept_id'='concept_id')) %>%
    compute_new(indexes=c('condition_occurrence_id'))
  
  
  cond_info_src <- cond_tbl %>%
    anti_join(cond_info, by=c('condition_occurrence_id')) %>%
    inner_join(codeset %>% select(concept_id), by=c('condition_source_concept_id'='concept_id')) %>%
    dplyr::union(cond_info) %>%
    compute_new(indexes=c('condition_occurrence_id'))
  
}




condition_outcomes <- function(outcome_tbl=outcome,
                               study_start=as.Date('2018-01-01'),
                               age_max=18L,
                               person_tbl=cdm_tbl('person')) {
  
  summary <- outcome_tbl %>%
    distinct(person_id, condition_start_date) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    mutate(age=(as.Date(condition_start_date)-as.Date(birth_date))/365.25) %>%
    filter(age<age_max) %>%
    group_by(condition_start_date) %>%
    summarise(n=n()) %>%
    filter(condition_start_date > study_start) %>%
    collect()
  
  # p <- ggplot(summary, aes(x=condition_start_date, y=n)) +
  #   geom_line() +
  #   scale_x_date(date_labels = "%Y %b", date_breaks = "3 months") +
  #   theme(axis.text.x=element_text(angle=60, hjust=1))
  
  spc(data = summary,
      x = condition_start_date,
      y = n,
      chart = "c",
      period="month",
      x.format = '%d-%m-%Y')
  #period = "month")
}



#' Get records of all conditions occurring follow_up_min days after initial dupilumab date
#' And limit to condition concept_ids occurring in at least 10 people
get_dup_cond_records <- function(dupilumab_tbl=results_tbl('dupilumab'),
                                 cond_tbl=cdm_tbl('condition_occurrence'),
                                 concept_tbl=vocabulary_tbl('concept'),
                                 follow_up_min=60L,
                                 min_cond_ct=10L) {
  
  earliest_dup <- dupilumab_tbl %>%
    group_by(person_id) %>%
    summarise(min_dup_date=min(drug_exposure_start_date)) %>%
    ungroup()
  
  conds <- cond_tbl %>%
    inner_join(earliest_dup, by=c('person_id')) %>%
    filter(condition_start_date >= min_dup_date + days(follow_up_min)) %>%
    mutate(cond_from='condition_concept_id') %>%
    compute_new(indexes=c('condition_occurrence_id'))
  
  conds_src <- cond_tbl %>%
    inner_join(earliest_dup, by=c('person_id')) %>%
    anti_join(conds, by=c('condition_occurrence_id')) %>%
    filter(condition_start_date >= min_dup_date + days(follow_up_min)) %>%
    mutate(cond_from='condition_source_concept_id') %>%
    dplyr::union(conds) %>%
    mutate(concept_id=case_when(cond_from=='condition_source_concept_id' ~ condition_source_concept_id,
                                cond_from=='condition_concept_id' ~ condition_concept_id)) %>%
    compute_new(indexes=c('person_id'))
  
  #' Find number of conditions indicated per concept_id
  outcome_conds_summary <- conds_src %>%
    group_by(concept_id) %>%
    summarise(n_cond_records=n()) %>%
    ungroup() %>%
    filter(n_cond_records > min_cond_ct) %>%
    compute_new()
  
  #' Find number of patients with conditions indicated per concept_id
  outcome_conds_summary_psn <- conds_src %>%
    distinct(person_id, concept_id) %>%
    group_by(concept_id) %>%
    summarise(n_cond_people=n()) %>%
    ungroup() %>%
    filter(n_cond_people > min_cond_ct) %>%
    select(concept_id, n_cond_people) %>%
    compute_new()
  
  #' Summarise findings
  outcome_conds_summary_fin <- outcome_conds_summary %>%
    inner_join(outcome_conds_summary_psn, by=c('concept_id')) %>%
    left_join(concept_tbl %>%
                copy_to_new(dest=config('db_src'), name='concept_tbl'), by=c('concept_id')) %>%
    select(concept_id, concept_name, vocabulary_id, n_cond_people, n_cond_records) %>%
    as_data_frame() %>%
    collect() %>%
    arrange(desc(n_cond_people))
  
}




#' Get a summary of most frequent conditions occurring at least follow_up_min days after initial dupilumab exposure
get_dupilumab_outcomes <- function(dupilumab_tbl=results_tbl('dupilumab'),
                                   outcome_codes=load_codeset('outcome_codes/all_outcome_codes_icd_snomed') %>%
                                     mutate(concept_id=as.numeric(concept_id)),
                                   cond_tbl=cdm_tbl('condition_occurrence'),
                                   concept_tbl=vocabulary_tbl('concept'),
                                   follow_up_min=60L,
                                   min_cond_ct=10L
) {
  
  earliest_dup <- dupilumab_tbl %>%
    group_by(person_id) %>%
    summarise(min_dup_date=min(drug_exposure_start_date)) %>%
    ungroup()
  
  #' #' String search for ICD10 concept codes within condition_source_value
  #' icd_outcome_codes <- outcome_codes %>%
  #'   filter(vocabulary_id %in% c('ICD10CM','ICD10')) %>%
  #'   distinct(concept_code) %>%
  #'   mutate(concept_code=lower(concept_code)) %>%
  #'   as_data_frame() %>%
  #'   collect()
  #' 
  #' icd_outcome_list <- icd_outcome_codes[[1]] %>% as.list()
  #' icd_outcomes_all <- paste(icd_outcome_list, collapse = '|')
  #' 
  #' outcome_conds_icd_subset <- cond_tbl %>%
  #'   inner_join(earliest_dup, by=c('person_id')) %>%
  #'   filter(condition_start_date >= min_dup_date + days(follow_up_min)) %>%
  #'   filter(str_detect(lower(condition_source_value),icd_outcomes_all)) %>%
  #'   compute_new(indexes=c('condition_occurrence_id'))
  
  #' Search through condition concept IDs
  outcome_conds <- cond_tbl %>%
    inner_join(earliest_dup, by=c('person_id')) %>%
    inner_join(outcome_codes, by=c('condition_concept_id'='concept_id')) %>%
    filter(condition_start_date >= min_dup_date + days(follow_up_min)) %>%
    compute_new(indexes=c('person_id'))
  
  #' Search through condition source concept IDs and concatenate with condition concept ID search
  outcome_cond_src <- cond_tbl %>%
    inner_join(earliest_dup, by=c('person_id')) %>%
    inner_join(outcome_codes, by=c('condition_source_concept_id'='concept_id')) %>%
    anti_join(outcome_conds, by=c('condition_occurrence_id')) %>%
    filter(condition_start_date >= min_dup_date + days(follow_up_min)) %>%
    mutate(cond_from='condition_source_concept_id') %>%
    dplyr::union(outcome_conds %>% mutate(cond_from='condition_concept_id')) %>%
    mutate(concept_id=case_when(cond_from=='condition_source_concept_id' ~ condition_source_concept_id,
                                cond_from=='condition_concept_id' ~ condition_concept_id)) %>%
    compute_new(indexes=c('person_id'))
  
  #' Find number of conditions indicated per concept_id
  outcome_conds_summary <- outcome_cond_src %>%
    group_by(concept_id, outcome_set_id, outcome_set_name) %>%
    summarise(n_cond_records=n()) %>%
    ungroup() %>%
    filter(n_cond_records > min_cond_ct) %>%
    compute_new()
  
  #' Find number of patients with conditions indicated per concept_id
  outcome_conds_summary_psn <- outcome_cond_src %>%
    distinct(person_id, concept_id, outcome_set_id, outcome_set_name) %>%
    group_by(concept_id) %>%
    summarise(n_cond_people=n()) %>%
    ungroup() %>%
    filter(n_cond_people > min_cond_ct) %>%
    select(concept_id, n_cond_people) %>%
    compute_new()
  
  #' Summarise findings
  outcome_conds_summary_fin <- outcome_conds_summary %>%
    inner_join(outcome_conds_summary_psn, by=c('concept_id')) %>%
    select(concept_id, outcome_set_id, outcome_set_name, n_cond_people, n_cond_records) %>%
    left_join(concept_tbl %>%
                copy_to_new(dest=config('db_src'),name='concept_tbl'),
              by=c('concept_id'='concept_id')) %>%
    as_data_frame() %>%
    collect() %>%
    arrange(outcome_set_id, desc(n_cond_people), n_cond_records,
            concept_id, concept_name, vocabulary_id)
  
}





#' Wrapper to compute all PMCA functions needed--
#' allows separate computation for each group in the cohort
#' @param cohort_tbl cohort table with person_id and observation_date
#' @param description description to be appended onto output pmca tables
#' @param pmca_xwalk_codes codeset with PMCA ICD10 codes to be input into the pmca_lookup fxn
pmca_functions <- function(cohort_tbl=rslt$covid_neg_oriac,
                           description='covnegoriac',
                           pmca_xwalk_codes=load_codeset('pmca_icd10_no_atopic')
                           
) {
  
  pmca_lookup_tbl_name <- paste0('pmca_lookup_',description)
  pmca_summary_tbl_name <- paste0('pmca_summary_',description)
  pmca_category_tbl_name <- paste0('pmca_category_',description)
  pmca_all_flags_tbl_name <- paste0('pmca_all_flags_',description)
  
  #' PMCA lookup tbl
  pmca_lookup <- produce_pmca_lookup(cohort=cohort_tbl,
                                     pmca_xwalk=pmca_xwalk_codes)
  
  output_tbl(pmca_lookup, pmca_lookup_tbl_name,
             indexes=c('person_id'))
  
  #' PMCA summary tbl
  pmca_summary <- compute_pmca_summary(pmca_lookup_tbl = results_tbl(pmca_lookup_tbl_name)) %>%
    compute_new()
  
  output_tbl(pmca_summary,
             pmca_summary_tbl_name,
             indexes=c('person_id'))
  
  #' PMCA categories tbl (most conservative algorithm)
  pmca_category <- compute_pmca_cats_cons(pmca_summary_tbl = pmca_summary,
                                          cohort_tbl = cohort_tbl %>% add_site())
  
  output_tbl(pmca_category,
             pmca_category_tbl_name,
             indexes=c('person_id'))
  
  #' PMCA all flags tbl
  final_pmca_all_flags <- pmca_all_flags(pmca_summary_tbl = pmca_summary,
                                         cohort_tbl = cohort_tbl,
                                         pmca_category = pmca_category) %>%
    left_join(cohort_tbl %>% select(person_id, observation_date), by='person_id')
  
  output_tbl(final_pmca_all_flags,
             pmca_all_flags_tbl_name,
             indexes=list('person_id'))
}


#' Function from pasc25_hazard_ratio / code / cohort_pmca
#' produce PMCA table with 3 year lookback
#'
#' @param cohort cohort of patients with patients and observation_date
#' @param pmca_xwalk codeset that has flags for body systems and whether or not progressive
#' @param condition_tbl formatting of the condition occurrence table
#' @param visit_tbl formatting of the visit occurrence table
#'
#' @return table that has conditions and flags for body systems
#' columns:
#' person_id | observation_date | result_derivation | condition_concept_id | condition_concept_name |
#' condition_source_concept_id | condition_source_value | visit_occurrence_id | condition_start_state |
#' description | body_system | progressive | visit_concept_id | site
#'

produce_pmca_lookup <- function (cohort,
                                 pmca_xwalk=load_codeset('pmca_icd10'),
                                 condition_tbl=cdm_tbl('condition_occurrence'),
                                 visit_tbl=cdm_tbl('visit_occurrence')) {
  only_pmca_conds <-
    cohort %>%
    inner_join(
      select(condition_tbl,
             person_id,
             condition_concept_id,
             condition_concept_name,
             condition_source_concept_id,
             condition_source_concept_name,
             visit_occurrence_id,
             condition_start_date),
      by=c('person_id')
    ) %>% filter(observation_date - condition_start_date <= 1096L &
                   observation_date - condition_start_date >= 0L) %>%
    inner_join(
      pmca_xwalk,
      by=c('condition_source_concept_id'='concept_id')
    ) %>%
    inner_join(
      select(
        visit_tbl,
        visit_occurrence_id,
        visit_concept_id
      )
    ) %>% filter(visit_concept_id %in% c(9202L,9201L,9203L,
                                         2000000048L,2000000088L,581399L)) %>%
    distinct()
  
}



#' Function from pasc25_hazard_ratio / code / cohort_pmca
#' compute patient, body system with visit number, and flags for malignancy and progressive
#'
#' @param pmca_lookup_tbl pmca table output from `produce_pmca_lookup`;
#' must contain `observation_date` and `body_system` and `condition_start_date` and flagged conditions
#'
#' @return computes information to be able to apply algorithms; groups by body system and counts visits,
#' with flags for progressive or malignancy for patients
#'
#' person_id | observation_date | result_derivation | body_system |
#' yr_1 | yr_2 | yr_3 | total_visits | progressive | malignancy
#'

compute_pmca_summary <- function(pmca_lookup_tbl) {
  
  add_year <-
    pmca_lookup_tbl %>%
    filter(
      ! body_system == 'malignancy'
    ) %>% mutate(
      flag_yr = case_when(
        condition_start_date < observation_date &
          condition_start_date >= sql("(observation_date - interval '1 year')::date") ~ 'yr_1',
        condition_start_date < sql("(observation_date - interval '1 year')::date") &
          condition_start_date >= sql("(observation_date - interval '2 years')::date") ~ 'yr_2',
        condition_start_date < sql("(observation_date - interval '2 years')::date") &
          condition_start_date >= sql("(observation_date - interval '3 years')::date") ~ 'yr_3',
        TRUE ~ 'no_yr'
      )) %>%
    group_by(
      person_id,
      observation_date,
      #result_derivation,
      body_system,
      flag_yr
    ) %>% summarise(
      visit_yr_ct=as.integer(n_distinct(condition_start_date))
    ) %>% pivot_wider(names_from = flag_yr,
                      values_from = visit_yr_ct,
                      values_fill = 0L) %>%
    #collect() %>%
    #mutate(total_visits=sum(c_across(starts_with("yr_")), na.rm=T)) %>%
    mutate(total_visits = yr_1 + yr_2 + yr_3) %>%
    ungroup() %>% compute_new(temporary=TRUE,
                              indexes=list('person_id'))
  
  progressive_malignant_pts <-
    pmca_lookup_tbl %>%
    filter(
      progressive == 'yes' |
        body_system == 'malignancy'
    ) %>% mutate(malignancy =
                   case_when(
                     body_system == 'malignancy' ~ 'yes',
                     TRUE ~ 'no'
                   )) %>%
    filter(malignancy=='yes' | progressive == 'yes') %>%
    select(person_id,
           progressive,
           malignancy) %>% distinct %>%
    compute_new(temporary=TRUE,
                indexes=list('person_id'))
  
  all_pts <-
    dplyr::union(
      add_year %>% ungroup() %>% select(person_id),
      progressive_malignant_pts %>% select(person_id)
    )
  
  all_pts %>%
    left_join(add_year) %>%
    left_join(progressive_malignant_pts) %>%
    mutate(
      progressive=case_when(
        is.na(progressive) ~ 'no',
        TRUE ~ progressive
      ),
      malignancy=case_when(
        is.na(malignancy) ~ 'no',
        TRUE ~ malignancy
      )
    ) %>% select(
      person_id,observation_date,#result_derivation,
      body_system,yr_1,yr_2,yr_3,total_visits,progressive,malignancy
    )
  
}



#' compute classification algorithm for *most conservative*: 
#' *complex chronic* is defined as patients with at least 1 visit for two body systems for all three years OR progressive OR malignant
#' *chronic* is defined as having at least 1 visit all three years for just one body system
#' 
#' @param pmca_lookup_tbl output from `produce_pmca_lookup`
#' @param cohort_tbl cohort table with person_id, observation_date and site
#' 
#' @return table that has patients in the *most conservative* category with the following columns:
#' person_id | observation_date | body_system | yr_1 | yr_2 | yr_3 | total_visits |
#' progressive | malignancy | complex_chronic | chronic | non_complex_chronic
#'

compute_pmca_cats_cons <- function(pmca_summary_tbl,
                                   cohort_tbl) {
  
  gt_two_bs <- 
    pmca_summary_tbl %>%
    filter(
      yr_1 > 0 &
        yr_2 > 0 &
        yr_3 > 0
    ) %>%
    group_by(person_id,
             observation_date) %>%
    summarise(body_system_ct=n_distinct(body_system)) %>%
    filter(
      body_system_ct > 1
    ) %>% ungroup()
  
  prog_or_malig <-
    pmca_summary_tbl %>% 
    filter(progressive == 'yes' | malignancy == 'yes') 
  
  complex_pts <- 
    dplyr::union(
      gt_two_bs %>% select(person_id,observation_date),
      prog_or_malig %>% select(person_id,observation_date)
    ) %>% mutate(complex_chronic = 1L) %>% 
    compute_new(temporary=TRUE,
                indexes=list('person_id'))
  
  chronic_pts <- 
    pmca_summary_tbl %>%
    filter(
      yr_1 > 0 &
        yr_2 > 0 &
        yr_3 > 0
    ) %>% anti_join(complex_pts,
                    by=c('person_id','observation_date')) %>%
    distinct(person_id, observation_date) %>% mutate(chronic = 1L) %>%
    compute_new(temporary=TRUE,
                indexes=list('person_id'))
  
  pmca_summary_tbl %>%
    right_join(cohort_tbl,
               by=c('person_id','observation_date')) %>% 
    left_join(complex_pts) %>%
    left_join(chronic_pts) %>%
    mutate(complex_chronic=case_when(is.na(complex_chronic) ~ 0L,
                                     TRUE ~ complex_chronic),
           chronic=case_when(is.na(chronic) ~ 0L,
                             TRUE ~ chronic)) %>%
    mutate(non_complex_chronic = 
             case_when(complex_chronic == 1 | chronic == 1 ~ 0L,
                       TRUE ~ 1L)) %>% select(person_id,
                                              observation_date,
                                              site,
                                              complex_chronic,
                                              chronic,
                                              non_complex_chronic) %>% distinct()
}





#' @param pmca_summary pmca summary table output from compute_pmca_summary. Counts are unique by person_id, body_system, and malignancy/progressive.
#' @param cohort_tbl full cohort tbl
#' @param pmca_category pmca_category table output from compute_pmca_cat_cons function
#' @return table with person-level flags for progressive/malignant and number of body systems affected
#' 
pmca_all_flags <- function(pmca_summary_tbl,
                           cohort_tbl,
                           pmca_category) {
  
  progressive_malignancy_flag <- pmca_summary_tbl %>%
    mutate(progressive_num=case_when(progressive=='yes' ~ 1L,
                                     progressive=='no' ~ 0L),
           malignancy_num=case_when(malignancy=='yes' ~ 1L,
                                    malignancy=='no' ~ 0L)) %>%
    group_by(person_id, observation_date) %>%
    summarise(progressive_ct = sum(as.numeric(progressive_num)),
              malignancy_ct=sum(as.numeric(malignancy_num)))
  
  body_systems_count <- pmca_summary_tbl %>%
    distinct(person_id, observation_date, body_system) %>%
    mutate(body_system_flag=1L) %>%
    ungroup() %>%
    group_by(person_id, observation_date) %>%
    summarise(n_body_systems = count(body_system_flag)) %>%
    ungroup()
  
  
  final_flags <- pmca_summary_tbl %>%
    distinct(person_id, observation_date) %>%
    left_join(body_systems_count) %>%
    left_join(progressive_malignancy_flag) %>%
    right_join(select(cohort_tbl, person_id)) %>%
    right_join(select(pmca_category, person_id, site, complex_chronic, chronic, non_complex_chronic)) %>%
    mutate(
      progressive_ct = case_when(progressive_ct > 0L ~ progressive_ct,
                                 TRUE ~ 0L),
      malignancy_ct = case_when(malignancy_ct > 0L ~ malignancy_ct,
                                TRUE ~ 0L),
      n_body_systems = case_when(n_body_systems > 0L ~ n_body_systems,
                                 TRUE ~ 0L)) %>%
    select(person_id,
           site,
           progressive_ct,
           malignancy_ct,
           complex_chronic,
           chronic,
           non_complex_chronic,
           n_body_systems)
  
}



#' @param cohort cohort: must include person_id
#' @return two figures in a list: one with visit types and one with total visits
get_prior_visit_count <- function(cohort=results_tbl('scd_pat_cohort_long') %>%
                                    mutate(date_of_entry=observation_date) %>%
                                    select(person_id, date_of_entry, group),
                                  visit_tbl=cdm_tbl('visit_occurrence'),
                                  days_prior_start=1095L,
                                  days_prior_end=1L) {
  
  prior_visits <- cohort %>%
    left_join(visit_tbl, by='person_id') %>%
    filter(as.Date(visit_start_date) >= as.Date(date_of_entry) - days(days_prior_start)) %>%
    filter(as.Date(visit_start_date) <= as.Date(date_of_entry) - days(days_prior_end)) %>%
    #filter(as.Date(visit_start_date) < as.Date(date_of_entry)) %>%
    #filter(as.Date(visit_start_date) > as.Date(date_of_entry) - days(days_prior)) %>%
    mutate(visit_type=case_when(visit_concept_id %in% c(9201L, 2000000048L, 2000000088L) ~ 'Inpatient',
                                visit_concept_id %in% c(9202L, 581399L) ~ 'Outpatient',
                                visit_concept_id %in% c(9203L) ~ 'ED',
                                TRUE ~ 'Other')) %>%
    distinct(person_id, date_of_entry, visit_start_date, visit_type, group) %>%
    group_by(person_id, date_of_entry, group, visit_type) %>%
    summarise(n_visits=n()) %>%
    ungroup() %>%
    mutate(n_visit_grp=case_when(n_visits>=1L && n_visits <6L ~ '01 to 05 visits',
                                 n_visits>=6L && n_visits <11L ~ '06 to 10 visits',
                                 n_visits>=11L && n_visits <25L ~ '11 to 24 visits',
                                 n_visits>=25L && n_visits <50L ~ '25 to 49 visits',
                                 n_visits>=50L && n_visits <100L ~ '50 to 99 visits',
                                 n_visits >=100L ~ '100+ visits',
                                 TRUE ~ '00 visits')) %>%
    pivot_wider(names_from=visit_type, values_from=c(n_visits, n_visit_grp)) %>%
    compute_new(indexes=c('person_id'))
  
  prior_visits2 <- prior_visits %>%
    full_join(cohort %>% select(person_id, date_of_entry, group), by=c('person_id','date_of_entry','group')) %>%
    mutate(ED = case_when(!is.na(n_visits_ED) ~ n_visits_ED,
                          TRUE ~ 0L),
           Inpatient = case_when(!is.na(n_visits_Inpatient) ~ n_visits_Inpatient,
                                 TRUE ~ 0L),
           Outpatient = case_when(!is.na(n_visits_Outpatient) ~ n_visits_Outpatient,
                                  TRUE ~ 0L),
           Other = case_when(!is.na(n_visits_Other) ~ n_visits_Other,
                             TRUE ~ 0L)) %>%
    select(-n_visits_ED, -n_visits_Inpatient, -n_visits_Outpatient, -n_visits_Other) %>%
    mutate(Total = ED + Inpatient + Outpatient + Other) %>%
    mutate(total_pre = case_when(Total>=1L && Total <6L ~ '01 to 05 visits',
                                 Total>=6L && Total <11L ~ '06 to 10 visits',
                                 Total>=11L && Total <25L ~ '11 to 24 visits',
                                 Total>=25L && Total <50L ~ '25 to 49 visits',
                                 Total>=50L && Total <100L ~ '50 to 99 visits',
                                 Total >=100L ~ '100+ visits',
                                 TRUE ~ '00 visits')) %>%
    mutate(n_visit_grp_ED = case_when(!is.na(n_visit_grp_ED) ~ n_visit_grp_ED,
                                      TRUE ~ '0 visits'),
           n_visit_grp_Inpatient = case_when(!is.na(n_visit_grp_Inpatient) ~ n_visit_grp_Inpatient,
                                             TRUE ~ '0 visits'),
           n_visit_grp_Outpatient = case_when(!is.na(n_visit_grp_Outpatient) ~ n_visit_grp_Outpatient,
                                              TRUE ~ '0 visits'),
           n_visit_grp_Other = case_when(!is.na(n_visit_grp_Other) ~ n_visit_grp_Other,
                                         TRUE ~ '0 visits')) %>%
    rename(visits_pre_ED=n_visit_grp_ED,
           visits_pre_IP=n_visit_grp_Inpatient,
           visits_pre_OP=n_visit_grp_Outpatient) %>%
    rename(total_visits=Total) %>%
    select(person_id,
           date_of_entry,
           group,
           total_pre,
           total_visits,
           visits_pre_ED,
           visits_pre_IP,
           visits_pre_OP) %>%
    compute_new(indexes=c('person_id')) %>%
    return()
}



#' Use within larger get_insurance_status function to categorize insurance
categorize_insurance <- function(visit_info=visits,
                                 visit_payer=cdm_tbl('visit_payer'),
                                 cohort_tbl=cohort_tbl) {
  
  visit_payer_reduced <- visit_payer %>%
    inner_join(visit_info %>% select(visit_occurrence_id), by='visit_occurrence_id') %>%
    mutate(plan_class_desc=case_when(str_detect(lower(plan_class),'medicare') ~ 'Medicare_Medicaid',
                                     str_detect(lower(plan_class),'medicaid') ~ 'Medicare_Medicaid',
                                     str_detect(lower(plan_class),'public') ~ 'Other_Public',
                                     str_detect(lower(plan_class),'private') ~ 'Private',
                                     str_detect(lower(plan_class),'self') ~ 'Private',
                                     TRUE ~ 'Other_Unknown')) %>%
    filter(!plan_class_desc=='Other_Unknown') %>%
    compute_new(indexes=c('visit_occurrence_id'))
  
  insurance <- visit_info %>%
    left_join(visit_payer_reduced, by='visit_occurrence_id') %>%
    select(person_id, index_date, cohort, visit_occurrence_id, plan_class_desc) %>%
    mutate(plan_class_desc=case_when(is.na(plan_class_desc) ~ 'Other_Unknown',
                                     TRUE ~ plan_class_desc)) %>%
    group_by(person_id, index_date, cohort, plan_class_desc) %>%
    summarise(n=n()) %>%
    ungroup() %>%
    pivot_wider(names_from=plan_class_desc, values_from=n, values_fill=0) %>%
    mutate(plan_class=case_when(Medicare_Medicaid >= 1L ~ 'Medicare_Medicaid', #Prioritizing Public over Private
                                Other_Public >= 1L ~ 'Other_Public',
                                Private >=1L ~ 'Private',
                                Other_Unknown >=1L ~ 'Other_Unknown')) %>%
    full_join(cohort_tbl, by=c('person_id','index_date','cohort')) %>%
    mutate(plan_class=case_when(!is.na(plan_class) ~ plan_class,
                                TRUE ~ 'Other_Unknown')) %>%
    compute_new(indexes=c('person_id'))
  
}

#' Look for insurance status from visits occurring within the study period closest in time prior to the index event
#' Searching through the closest 3 dates within 365 days
get_insurance_status <- function(study_start_date=as.Date('2020-04-01'),
                                 study_end_date=as.Date('2020-12-31'),
                                 visit_tbl=cdm_tbl('visit_occurrence'),
                                 visit_payer=cdm_tbl('visit_payer'),
                                 cohort_tbl=results_tbl('pasthma_aderm_2dx_demog'),
                                 max_time_cap=365L
) {
  
  visits <- visit_tbl %>%
    inner_join(cohort_tbl %>%
                 select(person_id, index_date, cohort), by='person_id') %>%
    filter(as.Date(visit_start_date) >= as.Date(index_date) - max_time_cap) %>%
    filter(as.Date(visit_start_date) <= as.Date(index_date)) %>%
    group_by(person_id, index_date, cohort) %>%
    mutate(closest_date=max(visit_start_date)) %>%
    ungroup() %>%
    filter(visit_start_date==closest_date) %>%
    compute_new(indexes=c('person_id'))
  
  insurance <- categorize_insurance(visit_info=visits,
                                    cohort_tbl=cohort_tbl)
  
  ##############################################################################
  #' Populate with second closest date
  no_plan_class <- insurance %>%
    filter(plan_class=='Other_Unknown') %>%
    distinct(person_id, index_date, cohort) %>%
    compute_new(indexes=c('person_id'))
  
  visits2 <- visit_tbl %>%
    inner_join(no_plan_class, by='person_id') %>%
    anti_join(visits, by='visit_occurrence_id') %>%
    filter(as.Date(visit_start_date) >= as.Date(index_date) - max_time_cap) %>%
    filter(as.Date(visit_start_date) <= as.Date(index_date)) %>%
    group_by(person_id, index_date, cohort) %>%
    mutate(closest_date=max(visit_start_date)) %>%
    ungroup() %>%
    filter(visit_start_date==closest_date) %>%
    compute_new(indexes=c('person_id'))
  
  insurance2 <- categorize_insurance(visit_info=visits2,
                                     cohort_tbl=no_plan_class)
  
  insurance_combo2 <- insurance %>%
    select(person_id, index_date, cohort, plan_class) %>%
    anti_join(insurance2, by=c('person_id','index_date','cohort')) %>%
    dplyr::union(insurance2 %>%
                   select(person_id, index_date, cohort, plan_class)) %>%
    compute_new(indexes=c('person_id'))
  
  ##############################################################################
  #' Populate with third closest date
  no_plan_class2 <- insurance_combo2 %>%
    filter(plan_class=='Other_Unknown') %>%
    distinct(person_id, index_date, cohort) %>%
    compute_new(indexes=c('person_id'))
  
  visits3 <- visit_tbl %>%
    inner_join(no_plan_class2, by='person_id') %>%
    anti_join(visits, by='visit_occurrence_id') %>%
    anti_join(visits2, by='visit_occurrence_id') %>%
    filter(as.Date(visit_start_date) >= as.Date(index_date) - max_time_cap) %>%
    filter(as.Date(visit_start_date) <= as.Date(index_date)) %>%
    group_by(person_id, index_date, cohort) %>%
    mutate(closest_date=max(visit_start_date)) %>%
    ungroup() %>%
    filter(visit_start_date==closest_date) %>%
    compute_new(indexes=c('person_id'))
  
  insurance3 <- categorize_insurance(visit_info=visits3,
                                     cohort_tbl=no_plan_class2)
  
  insurance_combo3 <- insurance_combo2 %>%
    anti_join(insurance2, by=c('person_id','index_date','cohort')) %>%
    dplyr::union(insurance2 %>%
                   select(person_id, index_date, cohort, plan_class)) %>%
    compute_new(indexes=c('person_id'))
  
  return(insurance_combo3)
  
}



#' Get visit type at index event
get_visit_type_at_index <- function(cohort_tbl=results_tbl('full_cohort'),
                                    visit_tbl=cdm_tbl('visit_occurrence')) {
  
  visit_type <- cohort_tbl %>% select(person_id, index_date, cohort) %>%
    inner_join(visit_tbl, by=c('person_id','index_date'='visit_start_date')) %>%
    group_by(person_id) %>%
    slice_min(visit_occurrence_id, n=1) %>%
    ungroup() %>%
    mutate(visit_type=case_when(visit_concept_id %in% c(9201L, 2000000048L, 2000000088L, 2000001532L) ~ 'Inpatient',
                                visit_concept_id %in% c(9202L, 581399L, 2000000469L) ~ 'Outpatient_and_telemedicine',
                                visit_concept_id %in% c(9203L) ~ 'ED',
                                visit_concept_id %in% c(44814711L,
                                                        42898160L,
                                                        44814710L,
                                                        2000000104L) ~ 'Other',
                                
                                TRUE ~ 'Other')) %>%
    select(person_id, index_date, cohort, visit_type, visit_concept_id) %>%
    compute_new()
  
}



#' Get list of when t1 conditions occur prior to a patient's index date
get_prior_t1_outcomes_by_psn <- function(cohort_tbl=results_tbl('pasthma_aderm_2dx_demog'),
                                         outcome_codes=load_codeset('outcome_codes/type1_disease_codes_reviewed_020524') %>%
                                           mutate(concept_id=as.numeric(concept_id)),
                                         cond_tbl=cdm_tbl('condition_occurrence')
) {
  
  prior_t1_outcomes <- list()
  
  for (a in (1:26)) {
    
    outcome_concepts <- outcome_codes %>%
      filter(outcome_set_id==a)
    
    #' Search through condition concept IDs
    prior_outcome_conds <- cond_tbl %>%
      inner_join(cohort_tbl, by=c('person_id')) %>%
      inner_join(outcome_concepts, by=c('condition_concept_id'='concept_id')) %>%
      filter(condition_start_date <= index_date) %>%
      compute_new(indexes=c('person_id'))
    
    #' Search through condition source concept IDs and concatenate with condition concept ID search
    prior_t1_outcomes[[a]] <- cond_tbl %>%
      inner_join(cohort_tbl, by=c('person_id')) %>%
      inner_join(outcome_concepts, by=c('condition_source_concept_id'='concept_id')) %>%
      anti_join(prior_outcome_conds, by=c('condition_occurrence_id')) %>%
      filter(condition_start_date <= index_date) %>%
      dplyr::union(prior_outcome_conds) %>%
      group_by(person_id, index_date, dx_group, outcome_set_id, outcome_set_name) %>%
      summarise(prior_cond_start_date=min(condition_start_date)) %>%
      ungroup() %>%
      compute_new(indexes=c('person_id'))
    
  }
  
  final_prior_t1_outcomes <- reduce(.x=prior_t1_outcomes,
                                    .f=dplyr::union) %>%
    compute_new(indexes=c('person_id'))
  
  print(a)
  
  return(final_prior_t1_outcomes)
}


#' Get list of when t1 conditions occur after a patient's index date (within a particular time-frame)
get_t1_outcomes_by_psn <- function(cohort_tbl=results_tbl('pasthma_aderm_2dx_demog'),
                                   outcome_codes=load_codeset('outcome_codes/type1_disease_codes_reviewed_020524') %>%
                                     mutate(concept_id=as.numeric(concept_id)),
                                   cond_tbl=cdm_tbl('condition_occurrence'),
                                   follow_up_min=60L,
                                   follow_up_max=1461L
) {
  
  t1_outcomes <- list()
  
  for (a in (1:26)) {
    
    outcome_concepts <- outcome_codes %>%
      filter(outcome_set_id==a)
    
    #' Search through condition concept IDs
    outcome_conds <- cond_tbl %>%
      inner_join(cohort_tbl, by=c('person_id')) %>%
      inner_join(outcome_concepts, by=c('condition_concept_id'='concept_id')) %>%
      filter(condition_start_date >= index_date + days(follow_up_min)) %>%
      filter(condition_start_date <= index_date + days(follow_up_max)) %>%
      compute_new(indexes=c('person_id'))
    
    #' Search through condition source concept IDs and concatenate with condition concept ID search
    t1_outcomes[[a]] <- cond_tbl %>%
      inner_join(cohort_tbl, by=c('person_id')) %>%
      inner_join(outcome_concepts, by=c('condition_source_concept_id'='concept_id')) %>%
      anti_join(outcome_conds, by=c('condition_occurrence_id')) %>%
      filter(condition_start_date >= index_date + days(follow_up_min)) %>%
      filter(condition_start_date <= index_date + days(follow_up_max)) %>%
      dplyr::union(outcome_conds) %>%
      group_by(person_id, index_date, dx_group, outcome_set_id, outcome_set_name) %>%
      summarise(min_cond_start_date=min(condition_start_date)) %>%
      ungroup() %>%
      compute_new(indexes=c('person_id'))
    
  }
  
  final_t1_outcomes <- reduce(.x=t1_outcomes,
                              .f=dplyr::union) %>%
    compute_new(indexes=c('person_id'))
  
  print(a)
  
  return(final_t1_outcomes)
}



#' Make wide table indicating dates and outcomes for each t1 disease of interest
make_outcomes_list <- function(t1_outcomes=results_tbl('t1_outcomes'),
                               prior_t1_outcomes=results_tbl('prior_t1_outcomes'),
                               cohort_tbl=results_tbl('pasthma_aderm_2dx_demog')) {
  
  any_prior_t1_outcomes <- prior_t1_outcomes %>%
    distinct(person_id, index_date, dx_group) %>%
    mutate(any_prior_t1_outcome=1L) %>%
    compute_new(indexes=c('person_id'))
  
  cht_without_prior_t1 <- cohort_tbl %>%
    anti_join(any_prior_t1_outcomes, by=c('person_id','index_date','dx_group')) %>%
    select(person_id, index_date, dx_group, cht) %>%
    compute_new(indexes=c('person_id'))
  
  for (a in 1:26) {
    
    t1_outcome <- t1_outcomes %>%
      filter(outcome_set_id==a)
    
    cond_start_date_name <- paste0('min_cond_start_date_',as.character(a))
    outcome_name <- paste0('t1_outcome_',as.character(a))
    
    outcome_info <- cht_without_prior_t1 %>%
      left_join(t1_outcome %>%
                  select(person_id, index_date, dx_group, min_cond_start_date), by=c('person_id','index_date','dx_group')) %>%
      mutate(!!outcome_name:=case_when(!is.na(min_cond_start_date) ~ 1L,
                                       TRUE ~ 0L)) %>%
      rename(!!cond_start_date_name:=min_cond_start_date)
    
    #' Iteratively add info for each t1 disease outcome
    cht_without_prior_t1 <- outcome_info
    
  }
  
  final_cht <- cht_without_prior_t1 %>%
    compute_new(indexes=c('person_id'))
  
  return(final_cht)
  
}



#' Estimate Type 1 disease prevalence by site in patients with asthma and atopic dermatitis
#' NOTE: patient ages are calculated as of July 1 of the indicated year, NOT showing age at diagnosis of type 1 disease or atopic disease
#' However, patients had to be aged <18 at the time of their type 1 diagnosis or atopic disease.
#' Get counts of t1 diseases over counts of patients with persistent asthma/atopic dermatitis by site, year, and age group
get_t1_pct_by_yr <- function(pasthma = results_tbl('persistent_asthma_dx'),
                             aderm = results_tbl('atopic_dermatitis_dx'),
                             outcome_codes=load_codeset('outcome_codes/type1_disease_codes_reviewed_020524') %>%
                               mutate(concept_id=as.numeric(concept_id)),
                             years=c('2017', '2018', '2019', '2020', '2021', '2022', '2023'),
                             person_tbl = cdm_tbl('person'),
                             cond_tbl = cdm_tbl('condition_occurrence'),
                             sites = c('cchmc','chop','colorado','lurie','national','nationwide',
                                       'nemours','seattle','stanford','texas'),
                             ages_start = c(0, 6, 12),
                             ages_end = c(6, 12, 18),
                             max_age = 18L
) {
  
  all_counts <- list()
  
  pasthma_age <- pasthma %>%
    left_join(person_tbl %>% select(person_id, birth_date)) %>%
    mutate(age = (condition_start_date - birth_date) / 325.25) %>%
    filter(age < max_age) %>%
    compute_new(indexes=c('person_id'))
  
  aderm_age <- aderm %>%
    left_join(person_tbl %>% select(person_id, birth_date)) %>%
    mutate(age = (condition_start_date - birth_date) / 325.25) %>%
    filter(age < max_age) %>%
    compute_new(indexes=c('person_id'))
  
  #' Search through condition concept IDs
  outcome_conds <- cond_tbl %>%
    inner_join(outcome_codes %>% distinct(concept_id), by=c('condition_concept_id'='concept_id')) %>%
    compute_new(indexes=c('person_id'))
  
  #' Search through condition source concept IDs and concatenate with condition concept ID search
  t1_age <- cond_tbl %>%
    inner_join(outcome_codes %>% distinct(concept_id), by=c('condition_source_concept_id'='concept_id')) %>%
    anti_join(outcome_conds, by=c('condition_occurrence_id')) %>%
    dplyr::union(outcome_conds) %>%
    group_by(person_id) %>%
    summarise(t1_start_date=min(condition_start_date)) %>%
    ungroup() %>%
    left_join(person_tbl %>% select(person_id, birth_date)) %>%
    mutate(age = (t1_start_date - birth_date) / 325.25) %>%
    filter(age < max_age) %>%
    compute_new(indexes=c('person_id'))
  
  
  global_count <- 0L
  
  for (a in 1:length(years)) {
    
    for (b in 1:length(sites)) {
      
      for (c in 1:length(ages_start)) {
        
        global_count <- global_count + 1L
        
        all_counts_num <- (a-1) + b*c
        
        year <- as.Date(paste0(as.character(years[[a]]),'-07-01'))
        curr_site <- sites[[b]]
        age_start <- ages_start[[c]]
        age_end <- ages_end[[c]]
        age_group <- paste0(as.character(age_start),' to <',as.character(age_end))
        
        #' Persistent asthma
        pasthma_age_yr <- pasthma_age %>%
          filter(age < age_end) %>%
          mutate(current_age=(year - birth_date)/365.25) %>%
          filter(site==curr_site) %>%
          filter(condition_start_date <= year) %>%
          group_by(person_id, current_age) %>%
          summarise(min_cond_start_date=min(condition_start_date),
                    max_cond_start_date=max(condition_start_date)) %>%
          ungroup() %>%
          mutate(cond_date_diff = max_cond_start_date - min_cond_start_date) %>%
          mutate(diff_6m_plus = case_when(cond_date_diff > 180L ~ 1L,
                                          TRUE ~ 0L))
        
        pasthma_1dx_ct <- pasthma_age_yr %>%
          filter(current_age >= age_start) %>%
          filter(current_age < age_end) %>%
          count() %>% pull()
        
        pasthma_2dx_ct <- pasthma_age_yr %>%
          filter(diff_6m_plus==1L) %>%
          filter(current_age >= age_start) %>%
          filter(current_age < age_end) %>%
          count() %>% pull()
        
        #' Atopic dermatitis
        aderm_age_yr <- aderm_age %>%
          mutate(current_age=(year - birth_date)/365.25) %>%
          filter(age < age_end) %>%
          filter(site==curr_site) %>%
          filter(condition_start_date <= year) %>%
          group_by(person_id, current_age) %>%
          summarise(min_cond_start_date=min(condition_start_date),
                    max_cond_start_date=max(condition_start_date)) %>%
          ungroup() %>%
          mutate(cond_date_diff = max_cond_start_date - min_cond_start_date) %>%
          mutate(diff_6m_plus = case_when(cond_date_diff > 180L ~ 1L,
                                          TRUE ~ 0L))
        
        aderm_1dx_ct <- aderm_age_yr %>%
          filter(current_age >= age_start) %>%
          filter(current_age < age_end) %>%
          count() %>% pull()
        
        aderm_2dx_ct <- aderm_age_yr %>%
          filter(diff_6m_plus==1L) %>%
          filter(current_age >= age_start) %>%
          filter(current_age < age_end) %>%
          count() %>% pull()
        
        #' Persistent asthma or atopic dermatitis
        pasthma_aderm_1dx <- pasthma_age_yr %>%
          distinct(person_id) %>%
          full_join(aderm_age_yr %>% distinct(person_id))
        
        pasthma_aderm_1dx_ct <- pasthma_aderm_1dx %>%
          count() %>% pull()
        
        pasthma_aderm_2dx <- pasthma_age_yr %>%
          filter(diff_6m_plus==1L) %>%
          distinct(person_id) %>%
          full_join(aderm_age_yr %>% filter(diff_6m_plus==1L) %>% distinct(person_id))
        
        pasthma_aderm_2dx_ct <- pasthma_aderm_2dx %>%
          count() %>% pull()
        
        #' T1 disease, of patients with 2 dx for persistent asthma or atopic dermatitis
        t1_age_yr <- t1_age %>%
          add_site() %>%
          filter(age < age_end) %>%
          inner_join(pasthma_aderm_2dx, by=c('person_id')) %>%
          mutate(current_age=(year - birth_date)/365.25) %>%
          filter(current_age >= age_start) %>%
          filter(current_age < age_end) %>%
          filter(site==curr_site) %>%
          filter(t1_start_date <= year) %>%
          group_by(person_id) %>%
          summarise(n_annual_t1=n()) %>%
          ungroup()
        
        t1_psn_ct <- t1_age_yr %>%
          count() %>% pull()
        
        #' All counts
        count_description <- c('Persistent asthma patients, 2 codes',
                               'Persistent asthma patients, 1 code',
                               'Atopic dermatitis patients, 2 codes',
                               'Atopic dermatitis patients, 1 code',
                               'Persistent asthma or atopic dermatitis patients, 2 codes',
                               'Persistent asthma or atopic dermatitis patients, 1 code'
        )
        
        counts <- c(pasthma_2dx_ct,
                    pasthma_1dx_ct,
                    aderm_2dx_ct,
                    aderm_1dx_ct,
                    pasthma_aderm_2dx_ct,
                    pasthma_aderm_1dx_ct
        )
        
        all_counts[[global_count]] <- cbind(count_description, counts) %>%
          as_tibble() %>%
          mutate(year=years[[a]]) %>%
          mutate(site=curr_site) %>%
          mutate(age_grp=age_group) %>%
          mutate(t1_pats=as.numeric(t1_psn_ct)) %>%
          mutate(counts=as.numeric(counts)) %>%
          mutate(t1_ct_pct=100*t1_pats/counts)
        
        print(paste0(global_count,': ',years[[a]],' : ',curr_site,' : ',age_group))
        
      }
      
    }
    
  }
  
  final_counts <- reduce(.x=all_counts,
                         .f=dplyr::union)
  
  return(final_counts)
  
}





#' Estimate Type 1 disease prevalence by site in patients IN COHORT with asthma, atopic dermatitis, or no asthma or atopic dermatitis
#' NOTE: patient ages are calculated as of DECEMBER 31ST (previously July 1st) of the indicated year, NOT showing age at diagnosis of type 1 disease or atopic disease
#' However, patients had to be aged <18 at the time of their type 1 diagnosis or atopic disease.
#' Get counts of t1 diseases over counts of patients with persistent asthma/atopic dermatitis by site, year, and age group
get_all_t1_pct_by_yr <- function(dup_cht=results_tbl('pasthma_aderm_dup'),
                                 nodup_cht=results_tbl('pasthma_aderm_2dx_lblf') %>% filter(lookback_pre==1L),
                                 nonatopic_cht=results_tbl('no_diasease_lblf') %>% filter(lookback_pre==1L),
                                 outcome_codes=load_codeset('outcome_codes/type1_disease_codes_reviewed_020524') %>%
                                   mutate(concept_id=as.numeric(concept_id)),
                                 years=c('2017', '2018', '2019', '2020', '2021', '2022', '2023'),
                                 person_tbl = cdm_tbl('person'),
                                 cond_tbl = cdm_tbl('condition_occurrence'),
                                 sites = c('cchmc','chop','colorado','lurie','national','nationwide',
                                           'nemours','seattle','stanford','texas'),
                                 ages_start = c(0, 6, 12),
                                 ages_end = c(6, 12, 18),
                                 max_age = 18L
) {
  
  all_counts <- list()
  
  atopic_pats <- dup_cht %>%
    distinct(person_id) %>%
    dplyr::union(nodup_cht %>%
                   distinct(person_id))
  
  nonatopic_pats <- nonatopic_cht %>%
    distinct(person_id)
  
  #' Search through condition concept IDs
  outcome_conds <- cond_tbl %>%
    inner_join(outcome_codes %>% distinct(concept_id), by=c('condition_concept_id'='concept_id')) %>%
    compute_new(indexes=c('person_id'))
  
  #' Search through condition source concept IDs and concatenate with condition concept ID search
  t1_age <- cond_tbl %>%
    inner_join(outcome_codes %>% distinct(concept_id), by=c('condition_source_concept_id'='concept_id')) %>%
    anti_join(outcome_conds, by=c('condition_occurrence_id')) %>%
    dplyr::union(outcome_conds) %>%
    group_by(person_id) %>%
    summarise(t1_start_date=min(condition_start_date)) %>%
    ungroup() %>%
    left_join(person_tbl %>% select(person_id, birth_date)) %>%
    mutate(age = (t1_start_date - birth_date) / 325.25) %>%
    filter(age < max_age) %>%
    compute_new(indexes=c('person_id'))
  
  
  global_count <- 0L
  
  for (a in 1:length(years)) {
    
    for (b in 1:length(sites)) {
      
      for (c in 1:length(ages_start)) {
        
        global_count <- global_count + 1L
        
        all_counts_num <- (a-1) + b*c
        
        #year <- as.Date(paste0(as.character(years[[a]]),'-07-01'))
        year <- as.Date(paste0(as.character(years[[a]]),'-12-31')) #' Previously set to July 1 of the year; updated to December 31
        curr_site <- sites[[b]]
        age_start <- ages_start[[c]]
        age_end <- ages_end[[c]]
        age_group <- paste0(as.character(age_start),' to <',as.character(age_end))
        
        #' Patients with atopic disease in age category
        atopic_age_yr <- atopic_pats %>%
          add_site() %>%
          left_join(person_tbl %>% select(person_id, birth_date)) %>%
          mutate(current_age = (year - birth_date) / 325.25) %>%
          filter(current_age < max_age) %>%
          filter(current_age >= age_start) %>%
          filter(current_age < age_end) %>%
          filter(site==curr_site) %>%
          distinct(person_id)
        
        atopic_ct <- atopic_age_yr %>% count() %>% pull()
        
        #' Patients with nonatopic disease in age category
        nonatopic_age_yr <- nonatopic_pats %>%
          add_site() %>%
          left_join(person_tbl %>% select(person_id, birth_date)) %>%
          mutate(current_age = (year - birth_date) / 325.25) %>%
          filter(current_age < max_age) %>%
          filter(current_age >= age_start) %>%
          filter(current_age < age_end) %>%
          filter(site==curr_site) %>%
          distinct(person_id)
        
        nonatopic_ct <- nonatopic_age_yr %>% count() %>% pull()
        
        #' T1 disease, of patients with atopic disease
        t1_age_yr_atopic <- t1_age %>%
          add_site() %>%
          filter(age < age_end) %>%
          inner_join(atopic_age_yr, by=c('person_id')) %>%
          mutate(current_age=(year - birth_date)/365.25) %>%
          filter(current_age >= age_start) %>%
          filter(current_age < age_end) %>%
          filter(site==curr_site) %>%
          filter(t1_start_date <= year) %>%
          group_by(person_id) %>%
          summarise(n_annual_t1=n()) %>%
          ungroup()
        
        #' T1 disease, of patients with nonatopic disease
        t1_age_yr_nonatopic <- t1_age %>%
          add_site() %>%
          filter(age < age_end) %>%
          inner_join(nonatopic_age_yr, by=c('person_id')) %>%
          mutate(current_age=(year - birth_date)/365.25) %>%
          filter(current_age >= age_start) %>%
          filter(current_age < age_end) %>%
          filter(site==curr_site) %>%
          filter(t1_start_date <= year) %>%
          group_by(person_id) %>%
          summarise(n_annual_t1=n()) %>%
          ungroup()
        
        t1_atopic_psn_ct <- t1_age_yr_atopic %>%
          count() %>% pull()
        
        t1_nonatopic_psn_ct <- t1_age_yr_nonatopic %>%
          count() %>% pull()
        
        #' All counts
        count_description <- c('Atopic disease patients in cohort',
                               'Nonatopic disease patients in cohort')
        
        cohort_counts <- c(atopic_ct,
                           nonatopic_ct)
        
        t1_counts <- c(t1_atopic_psn_ct,
                       t1_nonatopic_psn_ct)
        
        all_counts[[global_count]] <- cbind(count_description, cohort_counts, t1_counts) %>%
          as_tibble() %>%
          mutate(year=years[[a]]) %>%
          mutate(site=curr_site) %>%
          mutate(age_grp=age_group) %>%
          mutate(cohort_counts=as.numeric(cohort_counts)) %>%
          mutate(t1_counts=as.numeric(t1_counts)) %>%
          mutate(t1_ct_pct=100*t1_counts/cohort_counts)
        
        print(paste0(global_count,': ',years[[a]],' : ',curr_site,' : ',age_group))
        
      }
      
    }
    
  }
  
  final_counts <- reduce(.x=all_counts,
                         .f=dplyr::union)
  
  return(final_counts)
  
}


#' Remove site stratification and perform only for ages 6 to 17 (age as of December 31 of the indicated year)
#' Estimate Type 1 disease prevalence by site in patients IN COHORT with asthma, atopic dermatitis, or no asthma or atopic dermatitis
#' NOTE: patient ages are calculated as of DECEMBER 31ST of the indicated year, NOT showing age at diagnosis of type 1 disease or atopic disease
#' Get counts of t1 diseases over counts of patients with persistent asthma/atopic dermatitis by year
get_all_t1_pct_by_yr_simplified <- function(dup_cht=results_tbl('pasthma_aderm_dup'),
                                            nodup_cht=results_tbl('pasthma_aderm_2dx_lblf') %>% filter(lookback_pre==1L),
                                            nonatopic_cht=results_tbl('no_diasease_lblf') %>% filter(lookback_pre==1L),
                                            outcome_codes=load_codeset('outcome_codes/type1_disease_codes_reviewed_020524') %>%
                                              mutate(concept_id=as.numeric(concept_id)),
                                            years=c('2017', '2018', '2019', '2020', '2021', '2022', '2023'),
                                            person_tbl = cdm_tbl('person'),
                                            cond_tbl = cdm_tbl('condition_occurrence'),
                                            sites = c('cchmc','chop','colorado','lurie','national','nationwide',
                                                      'nemours','seattle','stanford','texas'),
                                            min_age = 6L,
                                            max_age = 18L
) {
  
  all_counts <- list()
  
  atopic_pats <- dup_cht %>%
    distinct(person_id, index_date) %>%
    add_site() %>%
    filter(site %in% sites) %>%
    dplyr::union(nodup_cht %>%
                   distinct(person_id, index_date)) %>%
    distinct(person_id, index_date)
  
  nonatopic_pats <- nonatopic_cht %>%
    distinct(person_id) %>%
    add_site() %>%
    filter(site %in% sites) %>%
    distinct(person_id)
  
  #' Search through condition concept IDs
  outcome_conds <- cond_tbl %>%
    filter(site %in% sites) %>%
    inner_join(outcome_codes %>% distinct(concept_id), by=c('condition_concept_id'='concept_id')) %>%
    compute_new(indexes=c('person_id'))
  
  #' Search through condition source concept IDs and concatenate with condition concept ID search
  t1_age <- cond_tbl %>%
    filter(site %in% sites) %>%
    inner_join(outcome_codes %>% distinct(concept_id), by=c('condition_source_concept_id'='concept_id')) %>%
    anti_join(outcome_conds, by=c('condition_occurrence_id')) %>%
    dplyr::union(outcome_conds) %>%
    group_by(person_id) %>%
    summarise(t1_start_date=min(condition_start_date)) %>%
    ungroup() %>%
    left_join(person_tbl %>% select(person_id, birth_date)) %>%
    mutate(age = (t1_start_date - birth_date) / 325.25) %>%
    filter(age >= min_age) %>%
    filter(age < max_age) %>%
    compute_new(indexes=c('person_id'))
  
  for (a in 1:length(years)) {
    
    year <- as.Date(paste0(as.character(years[[a]]),'-12-31'))
    age_group <- paste0(as.character(min_age),' to <',as.character(max_age))
    
    #' Patients with atopic disease in age category
    atopic_age_yr <- atopic_pats %>%
      left_join(person_tbl %>% select(person_id, birth_date)) %>%
      mutate(current_age = (year - birth_date) / 325.25) %>%
      filter(current_age >= min_age) %>%
      filter(current_age < max_age) %>%
      distinct(person_id)
    
    atopic_ct <- atopic_age_yr %>% count() %>% pull()
    
    #' Patients with nonatopic disease in age category
    nonatopic_age_yr <- nonatopic_pats %>%
      left_join(person_tbl %>% select(person_id, birth_date)) %>%
      mutate(current_age = (year - birth_date) / 325.25) %>%
      filter(current_age >= min_age) %>%
      filter(current_age < max_age) %>%
      distinct(person_id)
    
    nonatopic_ct <- nonatopic_age_yr %>% count() %>% pull()
    
    #' T1 disease, of patients with atopic disease
    t1_age_yr_atopic <- t1_age %>%
      inner_join(atopic_age_yr, by=c('person_id')) %>%
      filter(t1_start_date <= year) %>%
      group_by(person_id) %>%
      summarise(n_annual_t1=n()) %>%
      ungroup()
    
    #' T1 disease, of patients with nonatopic disease
    t1_age_yr_nonatopic <- t1_age %>%
      inner_join(nonatopic_age_yr, by=c('person_id')) %>%
      filter(t1_start_date <= year) %>%
      group_by(person_id) %>%
      summarise(n_annual_t1=n()) %>%
      ungroup()
    
    t1_atopic_psn_ct <- t1_age_yr_atopic %>%
      count() %>% pull()
    
    t1_nonatopic_psn_ct <- t1_age_yr_nonatopic %>%
      count() %>% pull()
    
    #' All counts
    count_description <- c('Atopic disease patients in cohort',
                           'Nonatopic disease patients in cohort')
    
    cohort_counts <- c(atopic_ct,
                       nonatopic_ct)
    
    t1_counts <- c(t1_atopic_psn_ct,
                   t1_nonatopic_psn_ct)
    
    all_counts[[a]] <- cbind(count_description, cohort_counts, t1_counts) %>%
      as_tibble() %>%
      mutate(year=years[[a]]) %>%
      mutate(cohort_counts=as.numeric(cohort_counts)) %>%
      mutate(t1_counts=as.numeric(t1_counts)) %>%
      mutate(t1_ct_pct=100*t1_counts/cohort_counts)
    
    print(paste0(a,': ',years[[a]]))
    
  }
  
  final_counts <- reduce(.x=all_counts,
                         .f=dplyr::union)
  
  return(final_counts)
  
}



#' Get counts of dupilumab over counts of patients with persistent asthma/atopic dermatitis by site, year, and age group
get_dup_pct_by_yr <- function(pasthma = results_tbl('persistent_asthma_dx'),
                              aderm = results_tbl('atopic_dermatitis_dx'),
                              dupilumab = results_tbl('dupilumab'),
                              years=c('2017', '2018', '2019', '2020', '2021', '2022', '2023'),
                              person_tbl = cdm_tbl('person'),
                              sites = c('cchmc','chop','colorado','lurie','national','nationwide',
                                        'nemours','seattle','stanford','texas'),
                              ages_start = c(0, 5, 11),
                              ages_end = c(5, 12, 18)
) {
  
  all_counts <- list()
  
  pasthma_age <- pasthma %>%
    left_join(person_tbl %>% select(person_id, birth_date)) %>%
    mutate(age = (condition_start_date - birth_date) / 325.25) %>%
    compute_new(indexes=c('person_id'))
  
  aderm_age <- aderm %>%
    left_join(person_tbl %>% select(person_id, birth_date)) %>%
    mutate(age = (condition_start_date - birth_date) / 325.25) %>%
    compute_new(indexes=c('person_id'))
  
  dup_age <- dupilumab %>%
    left_join(person_tbl %>% select(person_id, birth_date)) %>%
    mutate(dup_age = (drug_exposure_start_date - birth_date) / 325.25) %>%
    compute_new(indexes=c('person_id'))
  
  global_count <- 0L
  
  for (a in 1:length(years)) {
    
    for (b in 1:length(sites)) {
      
      for (c in 1:length(ages_start)) {
        
        global_count <- global_count + 1L
        
        all_counts_num <- (a-1) + b*c
        
        year <- as.Date(paste0(as.character(years[[a]]),'-12-31'))
        year_start <- as.Date(paste0(as.character(years[[a]]),'-01-01'))
        curr_site <- sites[[b]]
        age_start <- ages_start[[c]]
        age_end <- ages_end[[c]]
        age_group <- paste0(as.character(age_start),' to <',as.character(age_end))
        
        #' Persistent asthma
        pasthma_age_yr <- pasthma_age %>%
          filter(age >= age_start) %>%
          filter(age < age_end) %>%
          filter(site==curr_site) %>%
          filter(condition_start_date <= year) %>%
          group_by(person_id) %>%
          summarise(min_cond_start_date=min(condition_start_date),
                    max_cond_start_date=max(condition_start_date)) %>%
          ungroup() %>%
          mutate(cond_date_diff = max_cond_start_date - min_cond_start_date) %>%
          mutate(diff_6m_plus = case_when(cond_date_diff > 180L ~ 1L,
                                          TRUE ~ 0L))
        
        pasthma_1dx_ct <- pasthma_age_yr %>%
          count() %>% pull()
        
        pasthma_2dx_ct <- pasthma_age_yr %>%
          filter(diff_6m_plus==1L) %>% count() %>% pull()
        
        #' Atopic dermatitis
        aderm_age_yr <- aderm_age %>%
          filter(age >= age_start) %>%
          filter(age < age_end) %>%
          filter(site==curr_site) %>%
          filter(condition_start_date <= year) %>%
          group_by(person_id) %>%
          summarise(min_cond_start_date=min(condition_start_date),
                    max_cond_start_date=max(condition_start_date)) %>%
          ungroup() %>%
          mutate(cond_date_diff = max_cond_start_date - min_cond_start_date) %>%
          mutate(diff_6m_plus = case_when(cond_date_diff > 180L ~ 1L,
                                          TRUE ~ 0L))
        
        aderm_1dx_ct <- aderm_age_yr %>%
          count() %>% pull()
        
        aderm_2dx_ct <- aderm_age_yr %>%
          filter(diff_6m_plus==1L) %>% count() %>% pull()
        
        #' Persistent asthma or atopic dermatitis
        pasthma_aderm_1dx_ct <- pasthma_age_yr %>%
          distinct(person_id) %>%
          full_join(aderm_age_yr %>% distinct(person_id)) %>%
          count() %>% pull()
        
        pasthma_aderm_2dx_ct <- pasthma_age_yr %>%
          filter(diff_6m_plus==1L) %>%
          distinct(person_id) %>%
          full_join(aderm_age_yr %>% filter(diff_6m_plus==1L) %>% distinct(person_id)) %>%
          count() %>% pull()
        
        #' Dupilumab
        dup_age_yr <- dup_age %>%
          filter(dup_age >= age_start) %>%
          filter(dup_age < age_end) %>%
          filter(site==curr_site) %>%
          filter(drug_exposure_start_date <= year) %>%
          filter(drug_exposure_start_date >= year_start) %>%
          group_by(person_id) %>%
          summarise(n_annual_dup=n()) %>%
          ungroup()
        
        dup_psn_ct <- dup_age_yr %>%
          count() %>% pull()
        
        #' All counts
        count_description <- c('Persistent asthma patients, 2 codes',
                               'Persistent asthma patients, 1 code',
                               'Atopic dermatitis patients, 2 codes',
                               'Atopic dermatitis patients, 1 code',
                               'Persistent asthma or atopic dermatitis patients, 2 codes',
                               'Persistent asthma or atopic dermatitis patients, 1 code'
        )
        
        counts <- c(pasthma_2dx_ct,
                    pasthma_1dx_ct,
                    aderm_2dx_ct,
                    aderm_1dx_ct,
                    pasthma_aderm_2dx_ct,
                    pasthma_aderm_1dx_ct
        )
        
        all_counts[[global_count]] <- cbind(count_description, counts) %>%
          as_tibble() %>%
          mutate(year=years[[a]]) %>%
          mutate(site=curr_site) %>%
          mutate(age_grp=age_group) %>%
          mutate(dup_pats=as.numeric(dup_psn_ct)) %>%
          mutate(counts=as.numeric(counts)) %>%
          mutate(dup_ct_pct=100*dup_pats/counts)
        
        print(paste0(global_count,': ',years[[a]],' : ',curr_site,' : ',age_group))
        
      }
      
    }
    
  }
  
  final_counts <- reduce(.x=all_counts,
                         .f=dplyr::union)
  
  return(final_counts)
  
}



#' Get attrition table for atopic disesase without dupilumab      
get_attrition_atopic_nodup <- function(asthma_or_aderm=results_tbl('pasthma_aderm_2codes'),
                                       
                                       pasthma_2dx=results_tbl('preindex_pasthma_2dx'),
                                       aderm_2dx=results_tbl('preindex_aderm_2dx'),
                                       
                                       pasthma_aderm_nodup=results_tbl('pasthma_aderm_2dx_visit'),
                                       pasthma_aderm_nodup_lb=results_tbl('pasthma_aderm_2dx_lblf'),
                                       
                                       t1_outcomes=results_tbl('cht_w_t1_outcomes'),
                                       final_cht=results_tbl('pasthma_aderm_2dx_demog_wcovs'),
                                       study_start_date=as.Date('2018-10-01'),
                                       study_end_date=as.Date('2022-06-01'),
                                       age_limit=18L,
                                       person_tbl=cdm_tbl('person'),
                                       visit_tbl=cdm_tbl('visit_occurrence')
) {
  
  # Patients aged <18 years by the study end date with a visit within the study period
  total_ct <- person_tbl %>% filter(((study_end_date-as.Date(birth_date))/365.25) < age_limit) %>%
    distinct(person_id) %>%
    inner_join(visit_tbl, by=c('person_id')) %>%
    filter(visit_start_date >= study_start_date,
           visit_start_date <= study_end_date) %>%
    distinct(person_id) %>% count() %>% pull()
  
  # Persistent asthma: 2 codes at least 6 months apart on or prior to cohort entry date
  pasthma_ct <- pasthma_2dx %>% distinct(person_id) %>% count() %>% pull()
  
  # Atopic dermatitis: 2 codes at least 6 months apart on or prior to cohort entry date
  aderm_ct <- aderm_2dx %>% distinct(person_id) %>% count() %>% pull()
  
  # Atopic disease: persistent asthma or atopic dermatitis on or prior to cohort entry date
  pasthma_aderm_ct <- pasthma_2dx %>% distinct(person_id) %>%
    full_join(aderm_2dx %>% distinct(person_id)) %>%
    count() %>% pull()
  
  #Not receiving Dupilumab within 4 years following earliest atopic disease diagnosis
  atopic_nodup_ct <- pasthma_aderm_nodup %>% distinct(person_id) %>% count() %>% pull()
  
  #' Visit 7 to 365 days before cohort entry date
  lb_ct <- pasthma_aderm_nodup_lb %>% filter(lookback_pre==1L) %>% distinct(person_id) %>% count() %>% pull()
  
  #' No Type 1 disease prior to cohort entry date
  end_ct <- t1_outcomes %>% select(person_id, index_date, cht) %>%
    mutate(cohort=cht) %>%
    left_join(final_cht, by=c('person_id','index_date','cohort')) %>%
    filter(cohort=='2_atopic_controls') %>% distinct(person_id) %>% count() %>% pull()
  
  
  attrition_cts <- c(total_ct, pasthma_ct, aderm_ct, pasthma_aderm_ct,
                     atopic_nodup_ct, lb_ct, end_ct) %>% as.numeric()
  
  attrition_pcts <- c(100, 100*pasthma_ct/total_ct, 100*aderm_ct/total_ct, 100*pasthma_aderm_ct/total_ct,
                      100*atopic_nodup_ct/pasthma_aderm_ct, 100*lb_ct/atopic_nodup_ct, 100*end_ct/lb_ct)
  
  attrition_steps <- c('1. Patients aged <18 years by the study end date with a visit within the study period',
                       '2. Persistent asthma: 2 codes at least 6 months apart on or prior to cohort entry date',
                       '3. Atopic dermatitis: 2 codes at least 6 months apart on or prior to cohort entry date',
                       '4. Atopic disease: persistent asthma or atopic dermatitis on or prior to cohort entry date',
                       '5. Not receiving Dupilumab within 4 years following earliest atopic disease diagnosis',
                       '6. Has visit 7 to 365 days before cohort entry date',
                       '7. Has no Type 1 disease prior to cohort entry date')
  
  attrition <- cbind(attrition_steps, attrition_cts, attrition_pcts) %>%
    as_tibble() %>%
    mutate(attrition_cts=as.numeric(attrition_cts)) %>%
    mutate(attrition_pcts=as.numeric(attrition_pcts))
  
  return(attrition)
  
}


#' Get attrition table for atopic disesase without dupilumab      
get_attrition_atopic_nodup_6to17 <- function(asthma_or_aderm=results_tbl('pasthma_aderm_2codes'),
                                             
                                             pasthma_2dx=results_tbl('preindex_pasthma_2dx'),
                                             aderm_2dx=results_tbl('preindex_aderm_2dx'),
                                             
                                             pasthma_aderm_nodup=results_tbl('pasthma_aderm_2dx_visit'),
                                             pasthma_aderm_nodup_lb=results_tbl('pasthma_aderm_2dx_lblf'),
                                             
                                             t1_outcomes=results_tbl('cht_w_t1_outcomes_rev'),
                                             final_cht=results_tbl('pasthma_aderm_2dx_demog_wcovs'),
                                             study_start_date=as.Date('2018-10-01'),
                                             study_end_date=as.Date('2022-06-01'),
                                             age_lower_limit=6L,
                                             age_upper_limit=18L,
                                             person_tbl=cdm_tbl('person'),
                                             visit_tbl=cdm_tbl('visit_occurrence')
) {
  
  # Patients aged 6 to <18 years by the study end date with a visit within the study period
  total_ct <- person_tbl %>% filter(((study_end_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((study_end_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    distinct(person_id) %>%
    inner_join(visit_tbl, by=c('person_id')) %>%
    filter(visit_start_date >= study_start_date,
           visit_start_date <= study_end_date) %>%
    distinct(person_id) %>% count() %>% pull()
  
  # Persistent asthma: 2 codes at least 6 months apart on or prior to cohort entry date
  pasthma_ct <- pasthma_2dx %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    distinct(person_id) %>% count() %>% pull()
  
  # Atopic dermatitis: 2 codes at least 6 months apart on or prior to cohort entry date
  aderm_ct <- aderm_2dx %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    distinct(person_id) %>% count() %>% pull()
  
  # Atopic disease: persistent asthma or atopic dermatitis on or prior to cohort entry date
  pasthma_aderm_ct <- pasthma_2dx %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    distinct(person_id) %>%
    full_join(aderm_2dx %>%
                left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
                filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
                filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
                distinct(person_id)) %>%
    count() %>% pull()
  
  #Not receiving Dupilumab within 4 years following earliest atopic disease diagnosis
  atopic_nodup_ct <- pasthma_aderm_nodup %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    distinct(person_id) %>% count() %>% pull()
  
  #' Visit 7 to 365 days before cohort entry date
  lb_ct <- pasthma_aderm_nodup_lb %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    filter(lookback_pre==1L) %>% distinct(person_id) %>% count() %>% pull()
  
  #' No Type 1 disease prior to cohort entry date
  end_ct <- t1_outcomes %>% select(person_id, index_date, cht) %>%
    mutate(cohort=cht) %>%
    left_join(final_cht, by=c('person_id','index_date','cohort')) %>%
    select(person_id, index_date, cohort) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    filter(cohort=='2_atopic_controls') %>% distinct(person_id) %>% count() %>% pull()
  
  
  attrition_cts <- c(total_ct, pasthma_ct, aderm_ct, pasthma_aderm_ct,
                     atopic_nodup_ct, lb_ct, end_ct) %>% as.numeric()
  
  attrition_pcts <- c(100, 100*pasthma_ct/total_ct, 100*aderm_ct/total_ct, 100*pasthma_aderm_ct/total_ct,
                      100*atopic_nodup_ct/pasthma_aderm_ct, 100*lb_ct/atopic_nodup_ct, 100*end_ct/lb_ct)
  
  attrition_steps <- c('1. Patients aged <18 years by the study end date with a visit within the study period',
                       '2. Persistent asthma: 2 codes at least 6 months apart on or prior to cohort entry date',
                       '3. Atopic dermatitis: 2 codes at least 6 months apart on or prior to cohort entry date',
                       '4. Atopic disease: persistent asthma or atopic dermatitis on or prior to cohort entry date',
                       '5. Not receiving Dupilumab within 4 years following earliest atopic disease diagnosis',
                       '6. Has visit 7 to 365 days before cohort entry date',
                       '7. Has no Type 1 disease prior to cohort entry date')
  
  attrition <- cbind(attrition_steps, attrition_cts, attrition_pcts) %>%
    as_tibble() %>%
    mutate(attrition_cts=as.numeric(attrition_cts)) %>%
    mutate(attrition_pcts=as.numeric(attrition_pcts))
  
  return(attrition)
  
}


#' Get attrition table for atopic disesase with dupilumab      
get_attrition_atopic_dup <- function(asthma_or_aderm=results_tbl('pasthma_aderm_2codes'),
                                     pasthma_2dx=results_tbl('preindex_pasthma_2dx'),
                                     aderm_2dx=results_tbl('preindex_aderm_2dx'),
                                     
                                     pasthma_dup=results_tbl('asthma_2dx_dup_lblf'),
                                     aderm_dup=results_tbl('aderm_2dx_dup_lblf'),
                                     
                                     pasthma_aderm_lb=results_tbl('pasthma_aderm_dup'),
                                     dup_tbl=results_tbl('dupilumab2'),
                                     
                                     t1_outcomes=results_tbl('cht_w_t1_outcomes_rev'),
                                     final_cht=results_tbl('pasthma_aderm_2dx_demog_wcovs'),
                                     study_start_date=as.Date('2018-10-01'),
                                     study_end_date=as.Date('2022-06-01'),
                                     age_limit=18L,
                                     person_tbl=cdm_tbl('person'),
                                     visit_tbl=cdm_tbl('visit_occurrence')
) {
  
  # Patients aged <18 years by the study end date with a visit within the study period
  total_ct <- person_tbl %>% filter(((study_end_date-as.Date(birth_date))/365.25) < age_limit) %>%
    distinct(person_id) %>%
    inner_join(visit_tbl, by=c('person_id')) %>%
    filter(visit_start_date >= study_start_date,
           visit_start_date <= study_end_date) %>%
    distinct(person_id) %>% count() %>% pull()
  
  # Persistent asthma: 2 codes at least 6 months apart on or prior to cohort entry date
  pasthma_ct <- pasthma_2dx %>% distinct(person_id) %>% count() %>% pull()
  
  # Atopic dermatitis: 2 codes at least 6 months apart on or prior to cohort entry date
  aderm_ct <- aderm_2dx %>% distinct(person_id) %>% count() %>% pull()
  
  # Atopic disease: persistent asthma or atopic dermatitis on or prior to cohort entry date
  pasthma_aderm_ct <- pasthma_2dx %>% distinct(person_id) %>%
    full_join(aderm_2dx %>% distinct(person_id)) %>%
    count() %>% pull()
  
  #Receiving Dupilumab within 4 years following earliest atopic disease diagnosis
  atopic_dup_ct <- pasthma_dup %>%
    left_join(person_tbl %>% select(person_id, birth_date)) %>%
    mutate(age = (index_date-birth_date) / 365.25) %>%
    filter(age < age_limit) %>%
    filter(index_date >= study_start_date) %>%
    filter(index_date <= study_end_date) %>%
    distinct(person_id) %>%
    full_join(aderm_dup %>% left_join(person_tbl %>% select(person_id, birth_date)) %>%
                mutate(age = (index_date-birth_date) / 365.25) %>%
                filter(age < age_limit) %>%
                filter(index_date >= study_start_date) %>%
                filter(index_date <= study_end_date) %>% distinct(person_id)) %>%
    count() %>% pull()
  
  #' Visit 7 to 365 days before cohort entry date
  lb_ct <- pasthma_aderm_lb %>% distinct(person_id) %>% count() %>% pull()
  
  #' No Type 1 disease prior to cohort entry date
  end_ct <- t1_outcomes %>% select(person_id, index_date, cht) %>%
    mutate(cohort=cht) %>%
    left_join(final_cht, by=c('person_id','index_date','cohort')) %>%
    filter(cohort=='1_atopic_dup') %>% distinct(person_id) %>% count() %>% pull()
  
  # Aged <18 with any dupilumab exposure during the study period
  init_dup_ct <- dup_tbl %>%
    group_by(person_id) %>%
    summarise(min_dup_date=min(drug_exposure_start_date)) %>%
    ungroup() %>%
    filter(min_dup_date >= study_start_date,
           min_dup_date <= study_end_date) %>%
    left_join(person_tbl %>% select(person_id, birth_date)) %>%
    mutate(age = (birth_date - min_dup_date) / 365.25) %>%
    filter(age < age_limit) %>% distinct(person_id) %>% count() %>% pull()
  
  attrition_cts <- c(total_ct, pasthma_ct, aderm_ct, pasthma_aderm_ct,
                     atopic_dup_ct, lb_ct, end_ct, init_dup_ct) %>% as.numeric()
  
  attrition_pcts <- c(100, 100*pasthma_ct/total_ct, 100*aderm_ct/total_ct, 100*pasthma_aderm_ct/total_ct,
                      100*atopic_dup_ct/pasthma_aderm_ct, 100*lb_ct/atopic_dup_ct, 100*end_ct/lb_ct, NA)
  
  attrition_steps <- c('1. Patients aged <18 years by the study end date with a visit within the study period',
                       '2. Persistent asthma: 2 codes at least 6 months apart on or prior to cohort entry date',
                       '3. Atopic dermatitis: 2 codes at least 6 months apart on or prior to cohort entry date',
                       '4. Atopic disease: persistent asthma or atopic dermatitis on or prior to cohort entry date',
                       '5. Receiving Dupilumab within 4 years following earliest atopic disease diagnosis',
                       '6. Has visit 7 to 365 days before cohort entry date',
                       '7. Has no Type 1 disease prior to cohort entry date',
                       'Extra note: Patients aged <18 with any Dupilumab exposure during the study period')
  
  attrition <- cbind(attrition_steps, attrition_cts, attrition_pcts) %>%
    as_tibble() %>%
    mutate(attrition_cts=as.numeric(attrition_cts)) %>%
    mutate(attrition_pcts=as.numeric(attrition_pcts))
  
  return(attrition)
  
}



#' Get attrition table for atopic disesase with dupilumab      
get_attrition_atopic_dup_6to17 <- function(asthma_or_aderm=results_tbl('pasthma_aderm_2codes'),
                                           pasthma_2dx=results_tbl('preindex_pasthma_2dx'),
                                           aderm_2dx=results_tbl('preindex_aderm_2dx'),
                                           
                                           pasthma_dup=results_tbl('asthma_2dx_dup_lblf'),
                                           aderm_dup=results_tbl('aderm_2dx_dup_lblf'),
                                           
                                           pasthma_aderm_lb=results_tbl('pasthma_aderm_dup'),
                                           dup_tbl=results_tbl('dupilumab2'),
                                           
                                           t1_outcomes=results_tbl('cht_w_t1_outcomes_rev'),
                                           final_cht=results_tbl('pasthma_aderm_2dx_demog_wcovs'),
                                           study_start_date=as.Date('2018-10-01'),
                                           study_end_date=as.Date('2022-06-01'),
                                           age_lower_limit=6L,
                                           age_upper_limit=18L,
                                           person_tbl=cdm_tbl('person'),
                                           visit_tbl=cdm_tbl('visit_occurrence')
) {
  
  # Patients aged <18 years by the study end date with a visit within the study period
  total_ct <- person_tbl %>% filter(((study_end_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((study_end_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    distinct(person_id) %>%
    inner_join(visit_tbl, by=c('person_id')) %>%
    filter(visit_start_date >= study_start_date,
           visit_start_date <= study_end_date) %>%
    distinct(person_id) %>% count() %>% pull()
  
  # Persistent asthma: 2 codes at least 6 months apart on or prior to cohort entry date
  pasthma_ct <- pasthma_2dx %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    distinct(person_id) %>% count() %>% pull()
  
  # Atopic dermatitis: 2 codes at least 6 months apart on or prior to cohort entry date
  aderm_ct <- aderm_2dx %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    distinct(person_id) %>% count() %>% pull()
  
  # Atopic disease: persistent asthma or atopic dermatitis on or prior to cohort entry date
  pasthma_aderm_ct <- pasthma_2dx %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    distinct(person_id) %>%
    full_join(aderm_2dx %>%
                left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
                filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
                filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
                distinct(person_id)) %>%
    count() %>% pull()
  
  #Receiving Dupilumab within 4 years following earliest atopic disease diagnosis
  atopic_dup_ct <- pasthma_dup %>%
    left_join(person_tbl %>% select(person_id, birth_date)) %>%
    mutate(age = (index_date-birth_date) / 365.25) %>%
    filter(age < age_upper_limit) %>%
    filter(age >= age_lower_limit) %>%
    filter(index_date >= study_start_date) %>%
    filter(index_date <= study_end_date) %>%
    distinct(person_id) %>%
    full_join(aderm_dup %>% left_join(person_tbl %>% select(person_id, birth_date)) %>%
                mutate(age = (index_date-birth_date) / 365.25) %>%
                filter(age < age_upper_limit) %>%
                filter(age >= age_lower_limit) %>%
                filter(index_date >= study_start_date) %>%
                filter(index_date <= study_end_date) %>% distinct(person_id)) %>%
    count() %>% pull()
  
  #' Visit 7 to 365 days before cohort entry date
  lb_ct <- pasthma_aderm_lb %>%
    select(person_id, index_date) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    distinct(person_id) %>% count() %>% pull()
  
  #' No Type 1 disease prior to cohort entry date
  end_ct <- t1_outcomes %>% select(person_id, index_date, cht) %>%
    mutate(cohort=cht) %>%
    left_join(final_cht, by=c('person_id','index_date','cohort')) %>%
    select(person_id, index_date, cohort) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    filter(cohort=='1_atopic_dup') %>% distinct(person_id) %>% count() %>% pull()
  
  # Aged <18 with any dupilumab exposure during the study period
  init_dup_ct <- dup_tbl %>%
    group_by(person_id) %>%
    summarise(min_dup_date=min(drug_exposure_start_date)) %>%
    ungroup() %>%
    filter(min_dup_date >= study_start_date,
           min_dup_date <= study_end_date) %>%
    left_join(person_tbl %>% select(person_id, birth_date)) %>%
    mutate(age = (min_dup_date - birth_date) / 365.25) %>%
    filter(age < age_upper_limit) %>%
    filter(age >= age_lower_limit) %>%
    distinct(person_id) %>% count() %>% pull()
  
  attrition_cts <- c(total_ct, pasthma_ct, aderm_ct, pasthma_aderm_ct,
                     atopic_dup_ct, lb_ct, end_ct, init_dup_ct) %>% as.numeric()
  
  attrition_pcts <- c(100, 100*pasthma_ct/total_ct, 100*aderm_ct/total_ct, 100*pasthma_aderm_ct/total_ct,
                      100*atopic_dup_ct/pasthma_aderm_ct, 100*lb_ct/atopic_dup_ct, 100*end_ct/lb_ct, NA)
  
  attrition_steps <- c('1. Patients aged 6 to <18 years by the study end date with a visit within the study period',
                       '2. Persistent asthma: 2 codes at least 6 months apart on or prior to cohort entry date',
                       '3. Atopic dermatitis: 2 codes at least 6 months apart on or prior to cohort entry date',
                       '4. Atopic disease: persistent asthma or atopic dermatitis on or prior to cohort entry date',
                       '5. Receiving Dupilumab within 4 years following earliest atopic disease diagnosis',
                       '6. Has visit 7 to 365 days before cohort entry date',
                       '7. Has no Type 1 disease prior to cohort entry date',
                       'Extra note: Patients aged 6 to <18 with any Dupilumab exposure during the study period')
  
  attrition <- cbind(attrition_steps, attrition_cts, attrition_pcts) %>%
    as_tibble() %>%
    mutate(attrition_cts=as.numeric(attrition_cts)) %>%
    mutate(attrition_pcts=as.numeric(attrition_pcts))
  
  return(attrition)
  
}


#' Get attrition table for atopic disesase without dupilumab      
get_attrition_nonatopic_nodup <- function(nonatopic=results_tbl('no_disease'),
                                          nonatopic_lb=results_tbl('no_diasease_lblf'),
                                          t1_outcomes=results_tbl('cht_w_t1_outcomes'),
                                          final_cht=results_tbl('pasthma_aderm_2dx_demog_wcovs'),
                                          study_start_date=as.Date('2018-10-01'),
                                          study_end_date=as.Date('2022-06-01'),
                                          age_limit=18L,
                                          person_tbl=cdm_tbl('person'),
                                          visit_tbl=cdm_tbl('visit_occurrence')
) {
  
  # Patients aged <18 years by the study end date with a visit within the study period
  total_ct <- person_tbl %>% filter(((study_end_date-as.Date(birth_date))/365.25) < age_limit) %>%
    distinct(person_id) %>%
    inner_join(visit_tbl, by=c('person_id')) %>%
    filter(visit_start_date >= study_start_date,
           visit_start_date <= study_end_date) %>%
    distinct(person_id) %>% count() %>% pull()
  
  # Non-atopic disease or prescriptions and no dupilumab
  nonatopic_ct <- nonatopic %>% distinct(person_id) %>% count() %>% pull()
  
  #' Visit 7 to 365 days before cohort entry date
  lb_ct <- nonatopic_lb %>% filter(lookback_pre==1L) %>% distinct(person_id) %>% count() %>% pull()
  
  #' No Type 1 disease prior to cohort entry date
  end_ct <- t1_outcomes %>% select(person_id, index_date, cht) %>%
    mutate(cohort=cht) %>%
    left_join(final_cht, by=c('person_id','index_date','cohort')) %>%
    filter(cohort=='3_non_atopic_controls') %>% distinct(person_id) %>% count() %>% pull()
  
  
  attrition_cts <- c(total_ct, nonatopic_ct, lb_ct, end_ct) %>% as.numeric()
  
  attrition_pcts <- c(100, 100*nonatopic_ct/total_ct, 100*lb_ct/nonatopic_ct, 100*end_ct/lb_ct)
  
  attrition_steps <- c('1. Patients aged <18 years by the study end date with a visit within the study period',
                       '2. No atopic px or dx on or before cohort entry, and no dupilumab on or before the end of follow-up',
                       '3. Has in-person visit 7-365 days before cohort entry date',
                       '4. Has no Type 1 disease prior to cohort entry date')
  
  attrition <- cbind(attrition_steps, attrition_cts, attrition_pcts) %>%
    as_tibble() %>%
    mutate(attrition_cts=as.numeric(attrition_cts)) %>%
    mutate(attrition_pcts=as.numeric(attrition_pcts))
  
  return(attrition)
  
}


#' Get attrition table for atopic disesase without dupilumab      
get_attrition_nonatopic_nodup_6to17 <- function(nonatopic=results_tbl('no_disease'),
                                                nonatopic_lb=results_tbl('no_disease_lblf'),
                                                t1_outcomes=results_tbl('cht_w_t1_outcomes_rev'),
                                                final_cht=results_tbl('pasthma_aderm_2dx_demog_wcovs'),
                                                study_start_date=as.Date('2018-10-01'),
                                                study_end_date=as.Date('2022-06-01'),
                                                age_lower_limit=6L,
                                                age_upper_limit=18L,
                                                person_tbl=cdm_tbl('person'),
                                                visit_tbl=cdm_tbl('visit_occurrence')
) {
  
  # Patients aged 6 to <18 years by the study end date with a visit within the study period
  total_ct <- person_tbl %>% filter(((study_end_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((study_end_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    distinct(person_id) %>%
    inner_join(visit_tbl, by=c('person_id')) %>%
    filter(visit_start_date >= study_start_date,
           visit_start_date <= study_end_date) %>%
    distinct(person_id) %>% count() %>% pull()
  
  # Non-atopic disease or prescriptions and no dupilumab
  nonatopic_ct <- nonatopic %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((visit_start_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((visit_start_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    distinct(person_id) %>% count() %>% pull()
  
  #' Visit 7 to 365 days before cohort entry date
  lb_ct <- nonatopic_lb %>% filter(lookback_pre==1L) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    distinct(person_id) %>% count() %>% pull()
  
  #' No Type 1 disease prior to cohort entry date
  end_ct <- t1_outcomes %>% select(person_id, index_date, cht) %>%
    mutate(cohort=cht) %>%
    left_join(final_cht, by=c('person_id','index_date','cohort')) %>%
    select(person_id, index_date, cohort) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    filter(((index_date-as.Date(birth_date))/365.25) < age_upper_limit) %>%
    filter(((index_date-as.Date(birth_date))/365.25) >= age_lower_limit) %>%
    filter(cohort=='3_non_atopic_controls') %>% distinct(person_id) %>% count() %>% pull()
  
  
  attrition_cts <- c(total_ct, nonatopic_ct, lb_ct, end_ct) %>% as.numeric()
  
  attrition_pcts <- c(100, 100*nonatopic_ct/total_ct, 100*lb_ct/nonatopic_ct, 100*end_ct/lb_ct)
  
  attrition_steps <- c('1. Patients aged 6 to <18 years by the study end date with a visit within the study period',
                       '2. No atopic px or dx on or before cohort entry, and no dupilumab on or before the end of follow-up',
                       '3. Has in-person visit 7-365 days before cohort entry date',
                       '4. Has no Type 1 disease prior to cohort entry date')
  
  attrition <- cbind(attrition_steps, attrition_cts, attrition_pcts) %>%
    as_tibble() %>%
    mutate(attrition_cts=as.numeric(attrition_cts)) %>%
    mutate(attrition_pcts=as.numeric(attrition_pcts))
  
  return(attrition)
  
}



#' Get persistent asthma: at least 1 prescription of ICS plus LABA on at least 2 occasions separated by
#' at least 6 months
get_persistent_asthma_from_drug <- function(asthma_meds=load_codeset('all_asthma_codes_combined'),
                                            drug_tbl=cdm_tbl('drug_exposure'),
                                            person_tbl=cdm_tbl('person'),
                                            days_apart=180,
                                            age_max=18L
) {
  
  #' Note: only checking drug_concept_id, not drug_source_concept_id currently
  asthma_drugs <- drug_tbl %>%
    inner_join(asthma_meds %>%
                 distinct(concept_id, category), by=c('drug_concept_id'='concept_id')) %>%
    compute_new(indexes=list('drug_exposure_id'))
  
  #' ICS, LABA and combination
  asthma_drug_ics <- asthma_drugs %>%
    filter(category=='ICS') %>%
    distinct(person_id, drug_exposure_start_date) %>%
    mutate(ics=1L) %>%
    compute_new(indexes=list('person_id'))
  
  asthma_drug_laba <- asthma_drugs %>%
    filter(category=='LABA') %>%
    distinct(person_id, drug_exposure_start_date) %>%
    mutate(laba=1L) %>%
    compute_new(indexes=list('person_id'))
  
  asthma_drug_ics_laba <- asthma_drug_ics %>%
    full_join(asthma_drug_laba, by=c('person_id','drug_exposure_start_date')) %>%
    filter(ics==1L && laba==1L) %>%
    group_by(person_id) %>%
    summarise(min_drug_date = min(drug_exposure_start_date),
              max_drug_date = max(drug_exposure_start_date)) %>%
    ungroup() %>%
    filter(max_drug_date - min_drug_date >= days_apart) %>%
    compute_new(indexes=list('person_id'))
  
  asthma_drug_combo <- asthma_drugs %>%
    filter(category=='combo') %>%
    group_by(person_id) %>%
    summarise(min_drug_date = min(drug_exposure_start_date),
              max_drug_date = max(drug_exposure_start_date)) %>%
    ungroup() %>%
    filter(max_drug_date - min_drug_date >= days_apart) %>%
    compute_new(indexes=list('person_id'))
  
  #' Get counts of patients with all types of categories and classify
  asthma_drug_pats <- asthma_drug_ics_laba %>%
    dplyr::union(asthma_drug_combo) %>%
    group_by(person_id) %>%
    summarise(min_drug_date_final=min(min_drug_date),
              max_drug_date_final=max(max_drug_date)) %>%
    ungroup() %>%
    rename(min_drug_date=min_drug_date_final,
           max_drug_date=max_drug_date_final) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    mutate(age = (min_drug_date - birth_date) / 365.25) %>%
    filter(age < age_max) %>%
    compute_new(indexes=list('person_id'))
  
  return(asthma_drug_pats)
  
}


#' Get moderate/severe atopic dermatitis
get_moderate_severe_aderm_from_drug <- function(aderm_meds=load_codeset('eczema_topical_codes_050324'),
                                                drug_tbl=cdm_tbl('drug_exposure'),
                                                person_tbl=cdm_tbl('person'),
                                                days_apart=180,
                                                age_max=18L
) {
  
  #' Note: only checking drug_concept_id, not drug_source_concept_id currently
  aderm_drugs <- drug_tbl %>%
    inner_join(aderm_meds %>%
                 distinct(concept_id, drug, severity), by=c('drug_concept_id'='concept_id')) %>%
    compute_new(indexes=list('drug_exposure_id'))
  
  #' Moderate and severe drugs
  aderm_drug <- aderm_drugs %>%
    filter(severity %in% c('moderate','severe')) %>%
    distinct(person_id, drug_exposure_start_date) %>%
    group_by(person_id) %>%
    summarise(min_drug_date = min(drug_exposure_start_date),
              max_drug_date = max(drug_exposure_start_date)) %>%
    ungroup() %>%
    filter(max_drug_date - min_drug_date >= days_apart) %>%
    compute_new(indexes=list('person_id'))
  
  #' Get counts of patients with all types of categories and classify
  aderm_drug_pats <- aderm_drug %>%
    group_by(person_id) %>%
    summarise(min_drug_date_final=min(min_drug_date),
              max_drug_date_final=max(max_drug_date)) %>%
    ungroup() %>%
    rename(min_drug_date=min_drug_date_final,
           max_drug_date=max_drug_date_final) %>%
    left_join(person_tbl %>% select(person_id, birth_date), by=c('person_id')) %>%
    mutate(age = (min_drug_date - birth_date) / 365.25) %>%
    filter(age < age_max) %>%
    compute_new(indexes=list('person_id'))
  
  return(aderm_drug_pats)
  
}


#' Get asthma severity only from what's available in drug_exposure
#' Asthma severity categories: mild, moderate/severe, other (does not fit in mild or moderate/severe categories), or no information
#' 5 year time limit
get_asthma_severity_from_drug <- function(cohort_tbl=results_tbl('pasthma_aderm_2dx_demog_wcovs'),
                                          asthma_meds=load_codeset('all_asthma_codes_combined'),
                                          drug_tbl=cdm_tbl('drug_exposure'),
                                          days_prior=1826L
) {
  
  #' Note: only checking drug_concept_id, not drug_source_concept_id currently
  asthma_drugs <- drug_tbl %>%
    inner_join(asthma_meds %>%
                 distinct(concept_id, category), by=c('drug_concept_id'='concept_id')) %>%
    compute_new(indexes=list('drug_exposure_id'))
  
  #' Table of all asthma medications within time frame for patients in the cohort
  asthma_drug_cht <- cohort_tbl %>%
    select(person_id, index_date, cohort) %>%
    inner_join(asthma_drugs, by=c('person_id')) %>%
    filter(drug_exposure_start_date <= index_date) %>%
    filter(drug_exposure_start_date >= index_date-days_prior) %>%
    compute_new(indexes=list('drug_exposure_id'))
  
  #' Get counts of patients with all types of categories and classify
  asthma_drug_categories <- asthma_drug_cht %>%
    distinct(person_id, index_date, cohort, category) %>%
    mutate(n=1L) %>%
    pivot_wider(names_from='category', values_from='n', values_fill='0') %>%
    mutate(moderate_severe=case_when(injection==1L ~ 1L,
                                     combo==1L ~ 1L,
                                     ICS==1L && (LABA+LTRA+theophylline+zileuton > 0L) ~ 1L,
                                     TRUE ~ 0L)) %>%
    mutate(mild_strict=case_when(ICS==1L && (LABA+LTRA+Cromolyn+theophylline+zileuton+injection+combo+SABA==0L) ~ 1L, #' this rule says ONLY ONE
                                 Cromolyn==1L && (LABA+LTRA+ICS+theophylline+zileuton+injection+combo+SABA==0L) ~ 1L,
                                 LTRA==1L && (LABA+ICS+Cromolyn+theophylline+zileuton+injection+combo+SABA==0L) ~ 1L,
                                 theophylline==1L && (LABA+LTRA+Cromolyn+ICS+zileuton+injection+combo+SABA==0L) ~ 1L,
                                 SABA==1L && (LABA+LTRA+Cromolyn+theophylline+zileuton+injection+combo+ICS==0L) ~ 1L,
                                 TRUE ~ 0L
    )) %>%
    mutate(mild_loose=case_when((ICS+Cromolyn+LTRA+theophylline+SABA > 0L) && (injection+combo+LABA+zileuton==0L) ~ 1L,
                                TRUE ~ 0L)) %>%
    mutate(asthma_severity=case_when(moderate_severe==1L && mild_strict==1L ~ 'moderate_severe',
                                     moderate_severe==1L && mild_strict==0L ~ 'moderate_severe',
                                     moderate_severe==0L && mild_strict==1L ~ 'mild_strict',
                                     moderate_severe==0L && mild_strict==0L ~ 'other',
                                     TRUE ~ 'other'
    )) %>%
    mutate(asthma_severity_loose=case_when(moderate_severe==1L && mild_loose==1L ~ 'moderate_severe',
                                           moderate_severe==1L && mild_loose==0L ~ 'moderate_severe',
                                           moderate_severe==0L && mild_loose==1L ~ 'mild_loose',
                                           moderate_severe==0L && mild_loose==0L ~ 'other',
                                           TRUE ~ 'other'
    )) %>%
    full_join(cohort_tbl %>% select(person_id, index_date, cohort), by=c('person_id','index_date','cohort')) %>%
    mutate(asthma_severity2=case_when(is.na(asthma_severity) ~ 'no_drug_information',
                                      TRUE ~ asthma_severity)) %>%
    select(-asthma_severity) %>%
    mutate(asthma_severity_loose2=case_when(is.na(asthma_severity_loose) ~ 'no_drug_information',
                                            TRUE ~ asthma_severity_loose)) %>%
    select(-asthma_severity_loose) %>%
    rename(asthma_severity=asthma_severity2) %>%
    rename(asthma_severity_loose=asthma_severity_loose2) %>%
    compute_new(indexes=list('person_id'))
  
  return(asthma_drug_categories)
  
}



#' Get asthma severity only from what's available in drug_exposure
#' Asthma severity categories: mild, moderate/severe, other (does not fit in mild or moderate/severe categories), or no information
#' 5 year time limit
get_aderm_severity_from_drug <- function(cohort_tbl=results_tbl('pasthma_aderm_2dx_demog_wcovs'),
                                         aderm_meds=load_codeset('eczema_topical_codes_050324'),
                                         drug_tbl=cdm_tbl('drug_exposure'),
                                         days_prior=1826L
) {
  
  aderm_drugs <- drug_tbl %>%
    inner_join(aderm_meds %>%
                 distinct(concept_id, drug, severity), by=c('drug_concept_id'='concept_id')) %>%
    mutate(concept_id_type='drug_concept_id') %>%
    compute_new(indexes=list('drug_exposure_id'))
  
  aderm_drugs_src <- drug_tbl %>%
    anti_join(aderm_drugs, by=c('drug_exposure_id')) %>%
    inner_join(aderm_meds %>%
                 distinct(concept_id, drug, severity), by=c('drug_source_concept_id'='concept_id')) %>%
    mutate(concept_id_type='drug_source_concept_id') %>%
    compute_new(indexes=list('drug_exposure_id'))
  
  aderm_drugs_all <- aderm_drugs %>%
    dplyr::union(aderm_drugs_src) %>%
    compute_new(indexes=list('drug_exposure_id'))
  
  #' Table of all aderm medications within time frame for patients in the cohort
  aderm_drug_cht <- cohort_tbl %>%
    select(person_id, index_date, cohort) %>%
    inner_join(aderm_drugs_all, by=c('person_id')) %>%
    filter(drug_exposure_start_date <= index_date) %>%
    filter(drug_exposure_start_date >= index_date-days_prior) %>%
    compute_new(indexes=list('drug_exposure_id'))
  
  #' Get counts of patients with all types of categories and classify
  aderm_drug_categories <- aderm_drug_cht %>%
    distinct(person_id, index_date, cohort, severity) %>%
    mutate(n=1L) %>%
    pivot_wider(names_from='severity', values_from='n') %>%
    mutate(mild=case_when(!is.na(mild) ~ mild,
                          TRUE ~ 0L)) %>%
    mutate(moderate=case_when(!is.na(moderate) ~ moderate,
                              TRUE ~ 0L)) %>%
    mutate(severe=case_when(!is.na(severe) ~ severe,
                            TRUE ~ 0L)) %>%
    mutate(aderm_severity=case_when(severe > 0L ~ 'severe',
                                    moderate > 0L ~ 'moderate',
                                    mild > 0L ~ 'mild')) %>%
    rename(aderm_severe=severe,
           aderm_moderate=moderate,
           aderm_mild=mild) %>%
    compute_new(indexes=list('person_id'))
  
  return(aderm_drug_categories)
  
}






