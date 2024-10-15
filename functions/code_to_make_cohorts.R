require(ggplot2)
require(rocqi)
require(tidyr)

# By convention, accumulate execution results in a list rather than as
# independent variables, in order to make returning the entire set easier
rslt <- list()

################################################################################
#' Checking dupilumab drug exposure outcomes

rslt$dup_codes2 <- load_codeset('dupilumab_codes2')

rslt$dupilumab_drugs2 <- cdm_tbl('drug_exposure') %>%
  inner_join(rslt$dup_codes2 %>% select(concept_id), by=c('drug_concept_id'='concept_id')) %>%
  compute_new(indexes=c('drug_exposure_id'))

rslt$dupilumab_drugs_src2 <- cdm_tbl('drug_exposure') %>%
  anti_join(rslt$dupilumab_drugs2, by=c('drug_exposure_id')) %>%
  inner_join(rslt$dup_codes2 %>% select(concept_id), by=c('drug_source_concept_id'='concept_id')) %>%
  dplyr::union(rslt$dupilumab_drugs2) %>%
  compute_new(indexes=c('drug_exposure_id'))

rslt$dupilumab_src_value2 <- cdm_tbl('drug_exposure') %>%
  anti_join(rslt$dupilumab_drugs_src2, by=c('drug_exposure_id')) %>%
  filter(str_detect(lower(drug_source_value),'dupilumab')) %>%
  dplyr::union(rslt$dupilumab_drugs_src2) %>%
  compute_new(indexes=c('drug_exposure_id'))

output_tbl(rslt$dupilumab_src_value2, 'dupilumab2', indexes=c('drug_exposure_id'))

##############################################################################
#' Identify patients with persistent asthma diagnosis codes, based only on the diagnosis code (not on asthma phenotype)

rslt$pas_codeset <- load_codeset('persistent_asthma_codeset')

rslt$pas_conds <- cdm_tbl('condition_occurrence') %>%
  inner_join(rslt$pas_codeset %>% distinct(concept_id), by=c('condition_concept_id'='concept_id')) %>%
  mutate(concept='condition_concept_id') %>%
  compute_new(indexes=c('condition_occurrence_id'))

rslt$pas_cond_src <- cdm_tbl('condition_occurrence') %>%
  anti_join(rslt$pas_conds, by=c('condition_occurrence_id')) %>%
  inner_join(rslt$pas_codeset %>% distinct(concept_id), by=c('condition_source_concept_id'='concept_id')) %>%
  mutate(concept='condition_source_concept_id') %>%
  dplyr::union(rslt$pas_conds) %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$pas_cond_src, 'persistent_asthma_dx', indexes=c('person_id'))

#' Add severity information
rslt$pas_cond <- results_tbl('persistent_asthma_dx') %>%
  filter(concept=='condition_concept_id') %>%
  left_join(vocabulary_tbl('concept') %>% select(concept_id, concept_name), by=c('condition_concept_id'='concept_id')) %>%
  mutate(pas_dx_severity=case_when(str_detect(lower(concept_name),'mild persistent') ~ 'mild',
                                   str_detect(lower(concept_name),'moderate persistent') ~ 'moderate',
                                   str_detect(lower(concept_name),'severe persistent') ~ 'severe',
                                   str_detect(lower(concept_name),'severe controlled persistent') ~ 'severe',
                                   str_detect(lower(concept_name),'severe uncontrolled persistent') ~ 'severe',
                                   TRUE ~ 'unknown'
  )) %>%
  compute_new(indexes=list(c('person_id')))

rslt$pas_cond_src <- results_tbl('persistent_asthma_dx') %>%
  filter(concept=='condition_source_concept_id') %>%
  left_join(vocabulary_tbl('concept') %>% select(concept_id, concept_name), by=c('condition_source_concept_id'='concept_id')) %>%
  mutate(pas_dx_severity=case_when(str_detect(lower(concept_name),'mild persistent') ~ 'mild',
                                   str_detect(lower(concept_name),'moderate persistent') ~ 'moderate',
                                   str_detect(lower(concept_name),'severe persistent') ~ 'severe',
                                   str_detect(lower(concept_name),'severe controlled persistent') ~ 'severe',
                                   str_detect(lower(concept_name),'severe uncontrolled persistent') ~ 'severe',
                                   TRUE ~ 'unknown'
  )) %>%
  compute_new(indexes=list(c('person_id')))

rslt$pas_cond_all <- rslt$pas_cond %>%
  dplyr::union(rslt$pas_cond_src) %>%
  output_tbl('persistent_asthma_dx_severity', indexes=c('person_id'))

##############################################################################
#' Identify patients with any asthma diagnosis codes, based only on the diagnosis code (not on asthma phenotype)

rslt$asth_codeset <- load_codeset('asthma_snomed_icd_codes')

rslt$asth_conds <- cdm_tbl('condition_occurrence') %>%
  inner_join(rslt$asth_codeset %>% distinct(concept_id), by=c('condition_concept_id'='concept_id')) %>%
  mutate(concept='condition_concept_id') %>%
  compute_new(indexes=c('condition_occurrence_id'))

rslt$asth_cond_src <- cdm_tbl('condition_occurrence') %>%
  anti_join(rslt$asth_conds, by=c('condition_occurrence_id')) %>%
  inner_join(rslt$asth_codeset %>% distinct(concept_id), by=c('condition_source_concept_id'='concept_id')) %>%
  mutate(concept='condition_source_concept_id') %>%
  dplyr::union(rslt$asth_conds) %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$asth_cond_src, 'any_asthma_dx', indexes=c('person_id'))

##############################################################################
#' Identify patients with atopic dermatitis diagnosis codes

rslt$aderm_codeset <- load_codeset('atopic_dermatitis_codeset3')

rslt$aderm_conds <- cdm_tbl('condition_occurrence') %>%
  inner_join(rslt$aderm_codeset %>% distinct(concept_id), by=c('condition_concept_id'='concept_id')) %>%
  mutate(concept='condition_concept_id') %>%
  compute_new(indexes=c('condition_occurrence_id'))

rslt$aderm_cond_src <- cdm_tbl('condition_occurrence') %>%
  anti_join(rslt$aderm_conds, by=c('condition_occurrence_id')) %>%
  inner_join(rslt$aderm_codeset %>% distinct(concept_id), by=c('condition_source_concept_id'='concept_id')) %>%
  mutate(concept='condition_source_concept_id') %>%
  dplyr::union(rslt$aderm_conds) %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$aderm_cond_src,'atopic_dermatitis_dx',indexes=c('person_id'))

##############################################################################
##############################################################################
##############################################################################
#' Group 1: Non-Atopic Disease, No Dupilumab
################################################################################

#' Identify patients with other atopic disease (eczema, IgE-mediated food allergy,
#' allergic rhinitis, eosinophilic esophagitis)

rslt$eczema_codeset <- load_codeset('eczema_concepts')
rslt$food_allergy_codeset <- load_codeset('food_allergy_icd_snomed') #' not igG-mediated
rslt$allergic_rhinitis_codeset <- load_codeset('allergic_rhinitis_icd_snomed')
rslt$eosinophilic_esophagitis_codeset <- load_codeset('esophagitis_icd_snomed') #' not eosinophilic

rslt$eczema_conds <- get_min_cond(codeset=rslt$eczema_codeset) %>%
  mutate(type='eczema') %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$eczema_conds,'eczema', indexes=c('person_id'))

rslt$food_allergy_conds <- get_min_cond(codeset=rslt$food_allergy_codeset) %>%
  mutate(type='food_allergy') %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$food_allergy_conds,'food_allergy', indexes=c('person_id'))

rslt$allergic_rhinitis_conds <- get_min_cond(codeset=rslt$allergic_rhinitis_codeset) %>%
  mutate(type='allergic_rhinitis') %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$allergic_rhinitis_conds,'allergic_rhinitis', indexes=c('person_id'))

rslt$eosinophilic_esophagitis_conds <- get_min_cond(codeset=rslt$eosinophilic_esophagitis_codeset) %>%
  mutate(type='esophagitis') %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$eosinophilic_esophagitis_conds,'esophagitis', indexes=c('person_id'))

##############################################################################
#' Identify patients with other asthma or atopic dermatitis prescriptions

rslt$asthma_aderm_px_codes <- load_codeset('atopic_px')

rslt$asthma_aderm_px <- get_min_drug(codeset=rslt$asthma_aderm_px_codes) %>%
  mutate(type='asthma_aderm_px')

output_tbl(rslt$asthma_aderm_px, 'asthma_aderm_px', indexes=c('person_id'))

rslt$other_atopic_or_px <- results_tbl('eczema') %>%
  dplyr::union(results_tbl('food_allergy')) %>%
  dplyr::union(results_tbl('allergic_rhinitis')) %>%
  dplyr::union(results_tbl('esophagitis')) %>%
  dplyr::union(results_tbl('asthma_aderm_px')) %>%
  group_by(person_id) %>%
  summarise(min_atopic_or_px_date=min(min_date)) %>%
  ungroup() %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$other_atopic_or_px, 'other_atopic_or_px',
           indexes=c('person_id'))

#' Get random visit date for all patients in between study start and end dates
rslt$rand_any <- get_random_visit(cohort=cdm_tbl('person'),
                                  study_start_date=as.Date('2018-10-01'),
                                  study_end_date=as.Date('2022-06-01'))

#' Limit to patients who are aged <18 as of random visit date
rslt$rand_under18 <- rslt$rand_any %>%
  left_join(cdm_tbl('person') %>% select(person_id, birth_date)) %>%
  filter(((visit_start_date-as.Date(birth_date))/365.25) < 18L) %>%
  compute_new(indexes=c('person_id'))

#' Remove patients with other atopic/asthma/aderm prescription or diagnosis on or before the CED
rslt$earliest_asthma <- results_tbl('any_asthma_dx') %>%
  group_by(person_id) %>%
  summarise(min_asthma_date=min(condition_start_date)) %>%
  ungroup() %>%
  compute_new(indexes=c('person_id'))

rslt$earliest_aderm <- results_tbl('atopic_dermatitis_dx') %>%
  group_by(person_id) %>%
  summarise(min_aderm_date=min(condition_start_date)) %>%
  ungroup() %>%
  compute_new(indexes=c('person_id'))

rslt$earliest_dup <- results_tbl('dupilumab2') %>%
  group_by(person_id) %>%
  summarise(min_dupilumab_date=min(drug_exposure_start_date)) %>%
  ungroup() %>%
  compute_new(indexes=c('person_id'))

rslt$no_disease <- rslt$rand_under18 %>%
  left_join(results_tbl('other_atopic_or_px'), by=c('person_id')) %>%
  left_join(rslt$earliest_asthma, by=c('person_id')) %>%
  left_join(rslt$earliest_aderm, by=c('person_id')) %>%
  left_join(rslt$earliest_dup, by=c('person_id')) %>%
  mutate(other_atopic_or_px_dx=case_when(min_atopic_or_px_date <= visit_start_date ~ 1L,
                                         TRUE ~ 0L)) %>%
  mutate(asthma=case_when(min_asthma_date <= visit_start_date ~ 1L,
                          TRUE ~ 0L)) %>%
  mutate(aderm=case_when(min_aderm_date <= visit_start_date ~ 1L,
                         TRUE ~ 0L)) %>%
  mutate(dup=case_when(min_dupilumab_date <= visit_start_date + days(1461) ~ 1L,
                       TRUE ~ 0L)) %>%
  filter(other_atopic_or_px_dx == 0L) %>%
  filter(asthma == 0L) %>%
  filter(aderm == 0L) %>%
  filter(dup == 0L) %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$no_disease %>%
             select(person_id, visit_start_date), 'no_disease', indexes=c('person_id'))

#' Apply lookback/lookforward criteria for in-person visits
rslt$no_disease_lblf <- get_prior_followup_criteria(cohort_tbl=results_tbl('no_disease') %>%
                                                      rename(index_date=visit_start_date),
                                                    visit_tbl=cdm_tbl('visit_occurrence'),
                                                    lookback_max=365L,
                                                    lookback_min=7L,
                                                    lookforward_min=60L,
                                                    lookforward_max=365L)

rslt$no_disease_lblf_out <- rslt$no_disease_lblf %>%
  output_tbl('no_disease_lblf')

rslt$no_disease_lblf_fmt <- results_tbl('no_disease_lblf') %>%
  filter(lookback_pre==1L) %>%
  mutate(persistent_asthma=0L) %>%
  mutate(atopic_dermatitis=0L) %>%
  mutate(dup_duration=0L) %>%
  mutate(min_dup_date=as.Date(NA)) %>%
  mutate(max_dup_date=as.Date(NA)) %>%
  mutate(dup_ct=0L) %>%
  mutate(daily_dup_dur=0L) %>%
  mutate(monthly_dup_ct=0L) %>%
  mutate(yearly_dup_ct=0L) %>%
  left_join(cdm_tbl('person') %>% select(person_id, birth_date, gender_concept_id, race_concept_id, ethnicity_concept_id), by='person_id') %>%
  mutate(age=(as.Date(index_date)-as.Date(birth_date))/365.25) %>%
  mutate(dx_group='non_atopic_controls') %>%
  select(person_id, index_date, atopic_dermatitis, persistent_asthma, dx_group, birth_date, gender_concept_id,
         race_concept_id, ethnicity_concept_id, min_dup_date, max_dup_date, dup_ct, daily_dup_dur, monthly_dup_ct,
         yearly_dup_ct, age) %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$no_disease_lblf_fmt, 'no_disease_fmt', indexes=c('person_id'))

################################################################################
################################################################################
################################################################################
#' Group 2: Atopic Disease, No Dupilumab: Using 2 condition codes
################################################################################

#' Identify patients with 2 persistent asthma or 2 atopic dermatitis codes at least 6 months apart
#' who were aged <18 at time of earliest dx

rslt$earliest_pasthma_2codes <- results_tbl('persistent_asthma_dx') %>%
  group_by(person_id) %>%
  summarise(min_pasthma_date=min(condition_start_date),
            max_pasthma_date=max(condition_start_date)) %>%
  ungroup() %>%
  mutate(pasthma_date_diff=max_pasthma_date-min_pasthma_date) %>%
  filter(pasthma_date_diff >= 180L) %>%
  left_join(cdm_tbl('person') %>% select(person_id, birth_date), by=c('person_id')) %>%
  mutate(age=(as.Date(min_pasthma_date)-as.Date(birth_date))/365.25) %>%
  filter(age<18L) %>%
  mutate(pasthma_2dx=1L) %>%
  compute_new(indexes=c('person_id'))

rslt$earliest_aderm_2codes <- results_tbl('atopic_dermatitis_dx') %>%
  group_by(person_id) %>%
  summarise(min_aderm_date=min(condition_start_date),
            max_aderm_date=max(condition_start_date)) %>%
  ungroup() %>%
  mutate(aderm_date_diff=max_aderm_date-min_aderm_date) %>%
  filter(aderm_date_diff >= 180L) %>%
  left_join(cdm_tbl('person') %>% select(person_id, birth_date), by=c('person_id')) %>%
  mutate(age=(as.Date(min_aderm_date)-as.Date(birth_date))/365.25) %>%
  filter(age<18L) %>%
  mutate(aderm_2dx=1L) %>%
  compute_new(indexes=c('person_id'))

rslt$pasthma_aderm_2codes <- rslt$earliest_pasthma_2codes %>%
  select(person_id, min_pasthma_date, pasthma_2dx) %>%
  full_join(rslt$earliest_aderm_2codes %>%
              select(person_id, min_aderm_date, aderm_2dx), by=c('person_id')) %>%
  mutate(pasthma_2dx=case_when(pasthma_2dx==1L ~ 1L,
                               TRUE ~ 0L)) %>%
  mutate(aderm_2dx=case_when(aderm_2dx==1L ~ 1L,
                             TRUE ~ 0L)) %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$pasthma_aderm_2codes, 'pasthma_aderm_2codes', indexes=c('person_id'))

################################################################################
#' Find patients with at least 1 prescription of inhaled corticosteroid plus LABA prescription
#' on at least 2 occasions separated by at least 180 days:
#' Looking for a combo prescription OR an ICS + LABA prescription on the same day to count as 1 occasion
#' NOT IMPLEMENTED: Exclude patients with a mild intermittent asthma code occurring after indication of persistent asthma

rslt$persistent_asthma_drug_pats <- get_persistent_asthma_from_drug(asthma_meds=load_codeset('all_asthma_codes_combined'),
                                                                    days_apart=180,
                                                                    age_max=18L)

output_tbl(rslt$persistent_asthma_drug_pats, 'persistent_asthma_drug_pats', indexes=c('person_id'))

################################################################################
#' Find patients with at least 1 prescription of class V or highter potency topical steroids
#' on at least 2 occasions separated by at least 180 days
#' NOT IMPLEMENTED: Exclude patients with a mild eczema code occurring after the indication of moderate-to-severe atopic dermatitis

rslt$moderate_severe_aderm_drug_pats <- get_moderate_severe_aderm_from_drug(aderm_meds=load_codeset('eczema_topical_codes_050324'),
                                                                            days_apart=180,
                                                                            age_max=18L)

output_tbl(rslt$moderate_severe_aderm_drug_pats, 'moderate_severe_aderm_drug_pats', indexes=c('person_id'))

################################################################################
#' Add drug info

rslt$pasthma_2codes_wdrugs <- results_tbl('pasthma_aderm_2codes') %>%
  filter(pasthma_2dx==1L) %>%
  select(person_id, min_pasthma_date) %>%
  dplyr::union(results_tbl('persistent_asthma_drug_pats') %>% select(person_id, min_drug_date) %>%
                 rename(min_pasthma_date=min_drug_date)) %>%
  group_by(person_id) %>%
  summarise(min_pasthma_date=min(min_pasthma_date)) %>%
  ungroup() %>%
  left_join(results_tbl('pasthma_aderm_2codes') %>%
              filter(pasthma_2dx==1L) %>% distinct(person_id, pasthma_2dx), by=c('person_id')) %>%
  left_join(results_tbl('persistent_asthma_drug_pats') %>%
              mutate(pasthma_2px=1L) %>% distinct(person_id, pasthma_2px), by=c('person_id')) %>%
  compute_new()

rslt$aderm_2codes_wdrugs <- results_tbl('pasthma_aderm_2codes') %>%
  filter(aderm_2dx==1L) %>%
  select(person_id, min_aderm_date) %>%
  dplyr::union(results_tbl('moderate_severe_aderm_drug_pats') %>% select(person_id, min_drug_date) %>%
                 rename(min_aderm_date=min_drug_date)) %>%
  group_by(person_id) %>%
  summarise(min_aderm_date=min(min_aderm_date)) %>%
  ungroup() %>%
  left_join(results_tbl('pasthma_aderm_2codes') %>%
              filter(aderm_2dx==1L) %>% distinct(person_id, aderm_2dx), by=c('person_id')) %>%
  left_join(results_tbl('moderate_severe_aderm_drug_pats') %>%
              mutate(aderm_2px=1L) %>% distinct(person_id, aderm_2px), by=c('person_id')) %>%
  compute_new()

rslt$pasthma_aderm_2codes_wdrugs <- rslt$pasthma_2codes_wdrugs %>%
  full_join(rslt$aderm_2codes_wdrugs, by=c('person_id')) %>%
  mutate(pasthma_2dx = case_when(pasthma_2dx == 1L ~ 1L,
                                 TRUE ~ 0L)) %>%
  mutate(pasthma_2px = case_when(pasthma_2px == 1L ~ 1L,
                                 TRUE ~ 0L)) %>%
  mutate(aderm_2dx = case_when(aderm_2dx == 1L ~ 1L,
                               TRUE ~ 0L)) %>%
  mutate(aderm_2px = case_when(aderm_2px == 1L ~ 1L,
                               TRUE ~ 0L)) %>%
  compute_new()

output_tbl(rslt$pasthma_aderm_2codes_wdrugs, 'pasthma_aderm_2codes_wdrugs', indexes=c('person_id'))

################################################################################
#' Get random visit date in study period for patients with any persistent asthma or atopic dermatitis dx

rslt$rand_2dx <- get_random_visit(cohort=results_tbl('pasthma_aderm_2codes_wdrugs'),
                                  study_start_date=as.Date('2018-10-01'),
                                  study_end_date=as.Date('2022-06-01'))

#' Narrow to patients who have a visit in the study period and are <18 at random visit date
rslt$pasthma_aderm_2dx_visit <- rslt$rand_2dx %>%
  filter(has_visit_in_pd==1L) %>%
  rename(index_date=visit_start_date) %>%
  left_join(cdm_tbl('person') %>% select(person_id, birth_date), by=c('person_id')) %>%
  mutate(age=(as.Date(index_date)-as.Date(birth_date))/365.25) %>%
  filter(age<18L) %>%
  compute_new(indexes=c('person_id'))

#' Identify patients whose earliest atopic dermatitis or persistent asthma dates occurred on or before their index date
rslt$preindex_aderm_2dx <- rslt$pasthma_aderm_2dx_visit %>%
  filter(aderm_2dx==1L | aderm_2px==1L) %>%
  select(person_id, min_aderm_date, index_date, aderm_2dx, aderm_2px) %>%
  filter(min_aderm_date <= index_date) %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$preindex_aderm_2dx, 'preindex_aderm_2dx', indexes=c('person_id'))

rslt$preindex_pasthma_2dx <- rslt$pasthma_aderm_2dx_visit %>%
  filter(pasthma_2dx==1L | pasthma_2px==1L) %>%
  select(person_id, min_pasthma_date, index_date, pasthma_2dx, pasthma_2px) %>%
  filter(min_pasthma_date <= index_date) %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$preindex_pasthma_2dx, 'preindex_pasthma_2dx', indexes=c('person_id'))

rslt$pasthma_aderm_2dx_vfin <- rslt$pasthma_aderm_2dx_visit %>%
  select(person_id, visit_occurrence_id, index_date, has_visit_in_pd) %>%
  left_join(rslt$preindex_aderm_2dx, by=c('person_id','index_date')) %>%
  left_join(rslt$preindex_pasthma_2dx, by=c('person_id','index_date')) %>%
  filter(!is.na(min_pasthma_date) || !is.na(min_aderm_date)) %>%
  mutate(pasthma_2dx=case_when(pasthma_2dx==1L ~ 1L,
                               TRUE ~ 0L)) %>%
  mutate(aderm_2dx=case_when(aderm_2dx==1L ~ 1L,
                             TRUE ~ 0L)) %>%
  mutate(pasthma_2px=case_when(pasthma_2px==1L ~ 1L,
                               TRUE ~ 0L)) %>%
  mutate(aderm_2px=case_when(aderm_2px==1L ~ 1L,
                             TRUE ~ 0L)) %>%
  compute_new(indexes=c('person_id'))

#' Remove patients who had a dupilumab prescription within 4 years after their index date or earlier
rslt$earliest_dup <- results_tbl('dupilumab2') %>%
  group_by(person_id) %>%
  summarise(min_dupilumab_date=min(drug_exposure_start_date)) %>%
  ungroup() %>%
  compute_new(indexes=c('person_id'))

rslt$pasthma_aderm_2dx_vfin_nodup <- rslt$pasthma_aderm_2dx_vfin %>%
  left_join(rslt$earliest_dup, by=c('person_id')) %>%
  mutate(dup_within_study = case_when(min_dupilumab_date <= (index_date + days(1461)) ~ 1L,
                                      TRUE ~ 0L)) %>%
  filter(dup_within_study==0L) %>%
  select(-dup_within_study) %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$pasthma_aderm_2dx_vfin_nodup, 'pasthma_aderm_2dx_visit', indexes=c('person_id'))

#' Apply lookback/lookforward criteria
rslt$pasthma_aderm_2dx_visit_lblf <- get_prior_followup_criteria(cohort_tbl=rslt$pasthma_aderm_2dx_vfin_nodup,
                                                                 lookback_max=365L,
                                                                 lookback_min=7L,
                                                                 lookforward_min=60L,
                                                                 lookforward_max=365L)

rslt$pasthma_aderm_2dx_cht <- rslt$pasthma_aderm_2dx_vfin %>%
  left_join(rslt$pasthma_aderm_2dx_visit_lblf) %>%
  #filter(lookback_pre_post==1L) %>%
  output_tbl('pasthma_aderm_2dx_lblf', indexes=c('person_id'))



################################################################################
################################################################################
################################################################################
#' Group 3: Atopic Disease, Dupilumab-Treated: Using 2 condition codes
################################################################################

#' Dupilumab patients with persistent asthma dx or atopic dermatitis dx on or before dupilumab date

#' Get earliest dupilumab dates per person aged <18
rslt$dup_dates <- results_tbl('dupilumab2') %>%
  select(person_id, drug_exposure_start_date) %>%
  left_join(cdm_tbl('person') %>% select(person_id, birth_date), by='person_id') %>%
  mutate(age=(as.Date(drug_exposure_start_date)-as.Date(birth_date))/365.25) %>%
  filter(age<18L) %>%
  group_by(person_id) %>%
  summarise(index_date=min(drug_exposure_start_date)) %>%
  ungroup()

#' Get patients with persistent asthma (2 dx or px at least 6 months apart) and dupilumab exposure on or after diagnosis
#' EARLIEST condition or px code needs to occur before the dupilumab exposure

rslt$asthma_pats_w_dup <- results_tbl('pasthma_aderm_2codes_wdrugs') %>%
  filter(pasthma_2dx==1L | pasthma_2px==1L) %>%
  inner_join(rslt$dup_dates, by=c('person_id')) %>%
  filter(min_pasthma_date <= index_date) %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$asthma_pats_w_dup, 'asthma_2dx_dup', indexes=c('person_id'))

rslt$aderm_pats_w_dup <- results_tbl('pasthma_aderm_2codes_wdrugs') %>%
  filter(aderm_2dx==1L | aderm_2px==1L) %>%
  inner_join(rslt$dup_dates, by=c('person_id')) %>%
  filter(min_aderm_date <= index_date) %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$aderm_pats_w_dup, 'aderm_2dx_dup', indexes=c('person_id'))

##############################################################################
#' Asthma and atopic dermatitis with dupilumab: IN-PERSON visit required: lookback/lookforward

rslt$asthma_dup_lblf <- get_prior_followup_criteria(cohort_tbl=results_tbl('asthma_2dx_dup'),
                                                    visit_tbl=cdm_tbl('visit_occurrence'),
                                                    lookback_max=365L,
                                                    lookback_min=7L,
                                                    lookforward_min=60L,
                                                    lookforward_max=365L)

rslt$asthma_dup_lblf_out <- rslt$asthma_dup_lblf %>%
  output_tbl('asthma_2dx_dup_lblf')


rslt$aderm_dup_lblf <- get_prior_followup_criteria(cohort_tbl=results_tbl('aderm_2dx_dup'),
                                                   visit_tbl=cdm_tbl('visit_occurrence'),
                                                   lookback_max=365L,
                                                   lookback_min=7L,
                                                   lookforward_min=60L,
                                                   lookforward_max=365L)

rslt$aderm_dup_lblf_out <- rslt$aderm_dup_lblf %>%
  output_tbl('aderm_2dx_dup_lblf')

################################################################################
################################################################################
################################################################################
#' Combine patients from different cohorts into one group (atopic with dupilumab,
#' atopic without dupilumab, non-atopic)
################################################################################

rslt$pasthma_aderm_dup <- results_tbl('asthma_2dx_dup_lblf') %>%
  #filter(lookback_pre_post==1L) %>%
  filter(lookback_pre==1L) %>%
  distinct(person_id, index_date) %>%
  mutate(persistent_asthma=1L) %>%
  full_join(results_tbl('aderm_2dx_dup_lblf') %>%
              filter(lookback_pre==1L) %>%
              #filter(lookback_pre_post==1L) %>%
              distinct(person_id, index_date) %>%
              mutate(atopic_dermatitis=1L)) %>%
  group_by(person_id) %>%
  summarise(index_date=min(index_date),
            atopic_dermatitis=sum(atopic_dermatitis),
            persistent_asthma=sum(persistent_asthma)) %>%
  mutate(persistent_asthma=case_when(persistent_asthma >= 1L ~ 1L,
                                     TRUE ~ 0L)) %>%
  mutate(atopic_dermatitis=case_when(atopic_dermatitis >= 1L ~ 1L,
                                     TRUE ~ 0L)) %>%
  mutate(dx_group=case_when(atopic_dermatitis==1L && persistent_asthma==1L ~ 'aderm_and_pasthma',
                            atopic_dermatitis==1L && persistent_asthma==0L ~ 'aderm_no_pasthma',
                            atopic_dermatitis==0L && persistent_asthma==1L ~ 'pasthma_no_aderm',
                            atopic_dermatitis==0L && persistent_asthma==0L ~ 'no_aderm_or_pasthma')) %>%
  filter(as.Date(index_date)<=as.Date('2022-06-01')) %>%
  left_join(cdm_tbl('person') %>% select(person_id, birth_date, gender_concept_id, race_concept_id, ethnicity_concept_id), by='person_id') %>%
  left_join(results_tbl('dup_duration2'), by=c('person_id')) %>%
  mutate(age = (index_date - birth_date) / 365.25) %>%
  filter(age < 18L) %>%
  output_tbl('pasthma_aderm_dup', indexes=c('person_id'))

################################################################################
#' Remove patients who had dupilumab exposure from pasthma/aderm overall list
#' Add back in dupilumab exposure patients

#' Cohort 1: dupilumab-treated
rslt$pasthma_aderm_dup_2dx_cht <- results_tbl('pasthma_aderm_dup') %>%
  select(person_id, index_date, atopic_dermatitis, persistent_asthma, dx_group,
         min_dup_date, max_dup_date, dup_ct, daily_dup_dur, monthly_dup_ct, yearly_dup_ct) %>%
  mutate(visit_occurrence_id=as.numeric(NA)) %>%
  mutate(trt='dupilumab') %>%
  mutate(cht='1_atopic_dup') %>%
  compute_new(indexes=c('person_id'))

#' Cohort 2: atopic controls
rslt$pasthma_aderm_2dx_cht <- results_tbl('pasthma_aderm_2dx_lblf') %>%
  anti_join(rslt$pasthma_aderm_dup_2dx_cht, by=c('person_id')) %>%
  filter(lookback_pre==1L) %>%
  mutate(persistent_asthma = case_when(pasthma_2px + pasthma_2dx >= 1L ~ 1L,
                                       TRUE ~ 0L)) %>%
  mutate(atopic_dermatitis = case_when(aderm_2dx + aderm_2px >= 1L ~ 1L,
                                       TRUE ~ 0L)) %>%
  select(person_id, visit_occurrence_id, index_date, persistent_asthma, atopic_dermatitis) %>%
  mutate(dx_group=case_when(atopic_dermatitis==1L && persistent_asthma==1L ~ 'aderm_and_pasthma',
                            atopic_dermatitis==1L && persistent_asthma==0L ~ 'aderm_no_pasthma',
                            atopic_dermatitis==0L && persistent_asthma==1L ~ 'pasthma_no_aderm',
                            atopic_dermatitis==0L && persistent_asthma==0L ~ 'no_aderm_or_pasthma')) %>%
  mutate(trt='no_dupilumab') %>%
  mutate(min_dup_date=as.Date(NA),
         max_dup_date=as.Date(NA),
         dup_ct=0L,
         daily_dup_dur=0L,
         monthly_dup_ct=0L,
         yearly_dup_ct=0L) %>%
  mutate(cht='2_atopic_controls') %>%
  compute_new(indexes=c('person_id'))

#' Cohort 3: non-atopic controls
rslt$no_disease_cht <- results_tbl('no_disease_fmt') %>%
  anti_join(rslt$pasthma_aderm_dup_2dx_cht, by=c('person_id')) %>%
  anti_join(rslt$pasthma_aderm_2dx_cht, by=c('person_id')) %>%
  select(person_id, index_date, atopic_dermatitis, persistent_asthma, dx_group,
         min_dup_date, max_dup_date, dup_ct, daily_dup_dur, monthly_dup_ct, yearly_dup_ct) %>%
  mutate(visit_occurrence_id=as.numeric(NA)) %>%
  mutate(trt='no_dupilumab') %>%
  mutate(cht='3_non_atopic_controls') %>%
  compute_new(indexes=c('person_id'))

#' Combine cohorts 1, 2 and 3
rslt$cohort_w_demog <- rslt$pasthma_aderm_dup_2dx_cht %>%
  dplyr::union(rslt$pasthma_aderm_2dx_cht) %>%
  dplyr::union(rslt$no_disease_cht) %>%
  filter(index_date >= as.Date('2018-10-01')) %>%
  filter(index_date <= as.Date('2022-06-01')) %>%
  left_join(cdm_tbl('person') %>% select(person_id, birth_date, gender_concept_id, race_concept_id, ethnicity_concept_id), by='person_id') %>%
  mutate(age = (index_date - birth_date) / 365.25) %>%
  filter(age < 18L) %>%
  mutate(
    sex_cat=case_when(#gender_concept_id==8507L ~ 'Male',
      gender_concept_id==8532L ~ 'Female',
      TRUE ~ 'Male or Other/unknown/ambiguous')) %>%
  mutate(
    raceth_cat=case_when(ethnicity_concept_id == 38003563L ~ 'Hispanic',
                         ethnicity_concept_id == 38003564L & race_concept_id == 8527L ~ 'NH_White',
                         ethnicity_concept_id == 38003564L & race_concept_id == 8516L ~ 'NH_Black/AA',
                         ethnicity_concept_id == 38003564L & race_concept_id == 8515L ~ 'NH_Asian',
                         ethnicity_concept_id == 38003564L & race_concept_id == 44814659L ~ 'NH_Other_or_Multiple_Race', #NH multiple race
                         ethnicity_concept_id == 38003564L & race_concept_id %in% c(8567L, 8557L) ~ 'NH_Other_or_Multiple_Race', #NH other
                         TRUE ~ 'Other/Unknown')) %>%
  mutate(
    age_group = case_when(#age < 1 ~ '<1',
      age < 5 ~ '00 to 04',
      age < 10 ~ '05 to 09',
      age < 15 ~ '10 to 14',
      age < 18 ~ '15 to 17',
      TRUE ~ 'Old or Unknown')) %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$cohort_w_demog, 'pasthma_aderm_2dx_demog', indexes=c('person_id'))

##############################################################################
#' Get outcomes
##############################################################################

outcome_codes_rev <- load_codeset('outcome_codes/type1_disease_codes_05_17_2024') %>%
  distinct(concept_id, outcome_set_id, outcome_set_name) %>%
  left_join(vocabulary_tbl('concept') %>% select(concept_id, concept_name, vocabulary_id),
            by=c('concept_id')) %>%
  mutate(concept_id=as.numeric(concept_id))

rslt$prior_t1_outcomes_rev <- get_prior_t1_outcomes_by_psn(cohort_tbl=results_tbl('pasthma_aderm_2dx_demog'),
                                                           outcome_codes=outcome_codes_rev)

output_tbl(rslt$prior_t1_outcomes_rev, 'prior_t1_outcomes_rev', indexes=c('person_id'))

rslt$t1_outcomes_rev <- get_t1_outcomes_by_psn(cohort_tbl=results_tbl('pasthma_aderm_2dx_demog'),
                                               outcome_codes=outcome_codes_rev,
                                               follow_up_min=60L,
                                               follow_up_max=1461L
)

output_tbl(rslt$t1_outcomes_rev, 't1_outcomes_rev', indexes=c('person_id'))

rslt$outcomes_list <- make_outcomes_list(t1_outcomes=results_tbl('t1_outcomes_rev'),
                                         prior_t1_outcomes=results_tbl('prior_t1_outcomes_rev'),
                                         cohort_tbl=results_tbl('pasthma_aderm_2dx_demog'))

output_tbl(rslt$outcomes_list, 'cht_w_t1_outcomes_rev', indexes=c('person_id'))

##############################################################################
#' Add covariates
#############################################################################

#' I'm currently setting a precedent with public at the top (medicare, medicaid, other public, commercial, self pay, other/unknown)
rslt$insurance <- get_insurance_status(study_start_date=as.Date('2018-10-01'),
                                       study_end_date=as.Date('2022-06-01'),
                                       cohort_tbl=results_tbl('pasthma_aderm_2dx_demog') %>%
                                         rename(cohort=cht),
                                       max_time_cap=365L)

rslt$insurance_output <- rslt$insurance %>%
  output_tbl('insurance_status', indexes=c('person_id'))

#############################################################################
message('Visit type at index event')

rslt$visit_type <- get_visit_type_at_index(cohort_tbl=results_tbl('pasthma_aderm_2dx_demog') %>%
                                             rename(cohort=cht)) %>%
  rename(incident_visit_concept_id=visit_concept_id) %>%
  compute_new(indexes=c('person_id'))

output_tbl(rslt$visit_type, 'visit_types', indexes=c('person_id'))

#############################################################################
message('Prior visit count')

rslt$prior_visits <- get_prior_visit_count(cohort=results_tbl('pasthma_aderm_2dx_demog') %>%
                                             rename(group=cht) %>%
                                             rename(date_of_entry=index_date),
                                           visit_tbl=cdm_tbl('visit_occurrence'),
                                           days_prior_start=365L,
                                           days_prior_end=1L)

output_tbl(rslt$prior_visits, 'prior_visits', indexes=c('person_id'))

##############################################################################
#' Estimate asthma severity based on asthma drugs

rslt$asthma_severity_from_drug <- get_asthma_severity_from_drug(cohort_tbl=results_tbl('pasthma_aderm_2dx_demog') %>%
                                                                  mutate(cohort=cht),
                                                                asthma_meds=load_codeset('all_asthma_codes_combined'),
                                                                days_prior=1826L)

output_tbl(rslt$asthma_severity_from_drug, 'asthma_severity_from_drug',
           indexes=list(c('person_id')))

#' Estimate asthma severity based on most recent asthma diagnosis (looking within 5 years prior to index date)
rslt$pas_severity <- results_tbl('persistent_asthma_dx_severity') %>%
  inner_join(results_tbl('pasthma_aderm_2dx_demog') %>% select(person_id, index_date), by=c('person_id')) %>%
  filter(condition_start_date <= index_date) %>%
  filter(condition_start_date >= (index_date - days(1826))) %>%
  group_by(person_id, index_date) %>%
  mutate(max_pasthma_date=max(condition_start_date)) %>%
  ungroup() %>%
  filter(condition_start_date==max_pasthma_date) %>%
  group_by(person_id, index_date, pas_dx_severity) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  pivot_wider(names_from=c('pas_dx_severity'), values_from=n, values_fill=0) %>%
  rename(asthma_mild=mild,
         asthma_moderate=moderate,
         asthma_severe=severe) %>%
  mutate(asthma_dx_severity=case_when(asthma_severe>0L ~ 'severe',
                                      asthma_moderate>0L && asthma_severe==0L ~ 'moderate',
                                      asthma_mild>0L && asthma_moderate+asthma_severe==0L ~ 'mild',
                                      TRUE ~ 'unknown'
  )) %>%
  select(person_id, index_date, asthma_dx_severity, asthma_severe, asthma_moderate, asthma_mild) %>%
  compute_new(indexes=list(c('person_id')))

output_tbl(rslt$pas_severity, 'asthma_dx_severity',
           indexes=list(c('person_id')))

##############################################################################
#' Estimate atopic dermatitis severity based on atopic dermatitis drugs 
rslt$aderm_severity_from_drug <- get_aderm_severity_from_drug(cohort_tbl=results_tbl('pasthma_aderm_2dx_demog') %>%
                                                                rename(cohort=cht),
                                                              aderm_meds=load_codeset('eczema_topical_codes_050324'),
                                                              days_prior=1826L)

output_tbl(rslt$aderm_severity_from_drug, 'aderm_severity_from_drug',
           indexes=list('person_id'))

#############################################################################
##############################################################################
message('Steps 1 and 2: Calculate PMCA lookup table and summary table for all sites in cohort')

#' Remove all references to asthma, atopic dermatitis and atopic disease from the pmca ICD10 codeset
rslt$all_atopic_codes <- load_codeset('asthma_snomed_icd_codes') %>% select(concept_id) %>%
  dplyr::union(load_codeset('persistent_asthma_codeset') %>% select(concept_id)) %>%
  dplyr::union(load_codeset('atopic_dermatitis_codeset3') %>% select(concept_id)) %>%
  dplyr::union(load_codeset('eczema_concepts') %>% select(concept_id)) %>%
  dplyr::union(load_codeset('food_allergy_icd_snomed') %>% select(concept_id)) %>%
  dplyr::union(load_codeset('allergic_rhinitis_icd_snomed') %>% select(concept_id)) %>%
  dplyr::union(load_codeset('esophagitis_icd_snomed') %>% select(concept_id))

rslt$pmca_icd10_no_atopic <- load_codeset('pmca_icd10') %>%
  anti_join(rslt$all_atopic_codes, by=c('concept_id')) %>%
  as_data_frame() %>%
  write_csv('specs/pmca_icd10_no_atopic.csv')

#' The specified cohort should have person_id and observation_date; a 3-year lookback from
#' the observation date is applied
rslt$pmca_cht1 <- results_tbl('pasthma_aderm_2dx_demog') %>%
  filter(cht=='1_atopic_dup') %>%
  rename(observation_date=index_date) %>%
  select(person_id, observation_date)

rslt$pmca_cht2 <- results_tbl('pasthma_aderm_2dx_demog') %>%
  filter(cht=='2_atopic_controls') %>%
  rename(observation_date=index_date) %>%
  select(person_id, observation_date)

rslt$pmca_cht3 <- results_tbl('pasthma_aderm_2dx_demog') %>%
  filter(cht=='3_non_atopic_controls') %>%
  rename(observation_date=index_date) %>%
  select(person_id, observation_date)

rslt$pmca_functions1 <- pmca_functions(cohort_tbl=rslt$pmca_cht1,
                                       description='1atopdup_noatopcodes',
                                       pmca_xwalk_codes=load_codeset('pmca_icd10_no_atopic'))

rslt$pmca_functions2 <- pmca_functions(cohort_tbl=rslt$pmca_cht2,
                                       description='2atopctrl_noatopcodes',
                                       pmca_xwalk_codes=load_codeset('pmca_icd10_no_atopic'))

rslt$pmca_functions3 <- pmca_functions(cohort_tbl=rslt$pmca_cht3,
                                       description='3natopctrl_noatopcodes',
                                       pmca_xwalk_codes=load_codeset('pmca_icd10_no_atopic'))

rslt$pmca_all_flags_final <- results_tbl('pmca_all_flags_1atopdup_noatopcodes') %>%
  mutate(cohort='1_atopic_dup') %>%
  dplyr::union(results_tbl('pmca_all_flags_2atopctrl_noatopcodes') %>%
                 mutate(cohort='2_atopic_controls')) %>%
  dplyr::union(results_tbl('pmca_all_flags_3natopctrl_noatopcodes') %>%
                 mutate(cohort='3_non_atopic_controls')) %>%
  select(-site) %>%
  rename(index_date=observation_date) %>%
  compute_new(indexes=c('person_id'))

##############################################################################
#' Find patients who died

rslt$death <- cdm_tbl('death') %>%
  group_by(person_id) %>%
  summarise(min_death_date=min(death_date)) %>%
  ungroup() %>%
  rename(death_date=min_death_date) %>%
  mutate(died=1L) %>%
  inner_join(results_tbl('pasthma_aderm_2dx_demog') %>% select(person_id), by=c('person_id')) %>%
  compute_new(indexes=c('person_id'))

##############################################################################
#' Find patients who had incidence of cancer within 60 days to 1461 days after the index date

rslt$cancer_codes <- load_codeset('dx_2022_12_Cancer_V1') %>%
  distinct(concept_id)

rslt$cancer_concept <- results_tbl('cht_w_t1_outcomes_rev') %>%
  select(person_id, index_date, cht) %>%
  left_join(cdm_tbl('condition_occurrence'), by=c('person_id')) %>%
  inner_join(rslt$cancer_codes, by=c('condition_concept_id'='concept_id')) %>%
  compute_new(indexes=list(c('condition_occurrence_id')))

rslt$cancer_concept_src <- results_tbl('cht_w_t1_outcomes_rev') %>%
  select(person_id, index_date, cht) %>%
  left_join(cdm_tbl('condition_occurrence'), by=c('person_id')) %>%
  filter(condition_concept_id==0L) %>%
  inner_join(rslt$cancer_codes, by=c('condition_source_concept_id'='concept_id')) %>%
  dplyr::union(rslt$cancer_concept) %>%
  group_by(person_id, index_date, cht) %>%
  summarise(min_cancer_date=min(condition_start_date)) %>%
  ungroup() %>%
  compute_new(indexes=list(c('person_id')))

rslt$incident_cancer <- rslt$cancer_concept_src %>%
  mutate(incident_cancer = case_when(!is.na(min_cancer_date) && (index_date - min_cancer_date) >= 60L &&
                                       (index_date - min_cancer_date) >= 1461L ~ 1L,
                                     TRUE ~ 0L)) %>%
  full_join(results_tbl('cht_w_t1_outcomes_rev') %>%
              select(person_id, index_date, cht), by=c('person_id','index_date','cht')) %>%
  mutate(incident_cancer = case_when(!is.na(incident_cancer) ~ incident_cancer,
                                     TRUE ~ 0L)) %>%
  select(person_id, index_date, incident_cancer, min_cancer_date) %>%
  output_tbl('incident_cancer',
             indexes=list(c('person_id')))

##############################################################################
#' Make table including all covariate information

rslt$cht_with_covariates <- results_tbl('pasthma_aderm_2dx_demog') %>%
  add_site() %>%
  rename(cohort=cht) %>%
  right_join(results_tbl('cht_w_t1_outcomes_rev'), by=c('person_id','index_date','dx_group')) %>%
  left_join(rslt$pmca_all_flags_final, by=c('person_id','index_date','cohort')) %>%
  left_join(rslt$death, by=c('person_id')) %>%
  left_join(results_tbl('incident_cancer'), by=c('person_id','index_date')) %>%
  left_join(results_tbl('insurance_status'), by=c('person_id','index_date','cohort')) %>%
  left_join(results_tbl('visit_types'), by=c('person_id','index_date','cohort')) %>%
  left_join(results_tbl('prior_visits') %>%
              rename(cohort=group) %>%
              rename(index_date=date_of_entry), by=c('person_id','index_date','cohort')) %>%
  left_join(results_tbl('asthma_severity_from_drug') %>%
              select(person_id, index_date, cohort, asthma_severity_loose), by=c('person_id','index_date','cohort')) %>%
  left_join(results_tbl('asthma_dx_severity'), by=c('person_id','index_date')) %>%
  left_join(results_tbl('aderm_severity_from_drug') %>%
              select(person_id, index_date, cohort, aderm_severity), by=c('person_id','index_date','cohort')) %>%
  mutate(cohort=case_when(cohort=='1_atopic_dup' ~ 'atopic_dx_w_dupilumab',
                          cohort=='2_atopic_controls' ~ 'atopic_dx_no_dupilumab',
                          cohort=='3_non_atopic_controls' ~ 'non_atopic_dx_no_dupilumab',
                          TRUE ~ NA)) %>%
  mutate(asthma_severity_loose=case_when(!is.na(asthma_severity_loose) ~ asthma_severity_loose,
                                         TRUE ~ 'no_drug_information')) %>%
  mutate(asthma_dx_severity=case_when(!is.na(asthma_dx_severity) ~ asthma_dx_severity,
                                      TRUE ~ 'unknown')) %>%
  mutate(asthma_severity_dx_drug=case_when(asthma_severity_loose=='moderate_severe' | asthma_dx_severity=='severe' | asthma_dx_severity=='moderate' ~ 'moderate_severe',
                                           asthma_severity_loose=='mild_loose' | asthma_dx_severity=='mild' ~ 'mild',
                                           TRUE ~ 'other_unknown'
  )) %>%
  mutate(asthma_severity_in_asthma=case_when(dx_group %in% c('aderm_and_pasthma','pasthma_no_aderm') ~ asthma_severity_dx_drug,
                                             TRUE ~ 'no_asthma')) %>%
  mutate(aderm_severity=case_when(!is.na(aderm_severity) ~ aderm_severity,
                                  TRUE ~ 'other_unknown')) %>%
  mutate(aderm_severity_in_aderm=case_when(dx_group %in% c('aderm_and_pasthma','aderm_no_pasthma') ~ aderm_severity,
                                           TRUE ~ 'no_atopic_dermatitis')) %>%
  mutate(asthma_aderm_severity=case_when(dx_group %in% c('aderm_no_pasthma','aderm_and_pasthma','pasthma_no_aderm') && aderm_severity=='severe' | aderm_severity=='moderate' | asthma_severity_dx_drug=='moderate_severe' ~ 'moderate_severe',
                                         dx_group %in% c('aderm_no_pasthma','aderm_and_pasthma','pasthma_no_aderm') && aderm_severity=='mild' && asthma_severity_dx_drug=='mild' ~ 'mild',
                                         dx_group %in% c('aderm_no_pasthma','aderm_and_pasthma','pasthma_no_aderm') && aderm_severity=='mild' && asthma_severity_dx_drug=='other_unknown' ~ 'mild',
                                         dx_group %in% c('aderm_no_pasthma','aderm_and_pasthma','pasthma_no_aderm') && aderm_severity=='other_unknown' && asthma_severity_dx_drug=='mild' ~ 'mild',
                                         TRUE ~ 'other_unknown'
  )) %>%
  mutate(complex_chronic_flag = case_when(complex_chronic == 1L ~ 'complex_chronic',
                                          non_complex_chronic == 1L ~ 'not_chronic_or_complex_chronic',
                                          chronic == 1L ~ 'chronic',
                                          TRUE ~ 'not_chronic_or_complex_chronic')) %>%
  mutate(progressive_flag = case_when(progressive_ct > 0L ~ 'yes',
                                      TRUE ~ 'no'),
         malignancy_flag = case_when(malignancy_ct > 0L ~ 'yes',
                                     TRUE ~ 'no')) %>%
  mutate(n_body_systems=case_when(is.na(n_body_systems) ~ 0L,
                                  TRUE ~ n_body_systems)) %>%
  mutate(nbs=case_when(
    n_body_systems>=5L~"05-17 body systems",
    n_body_systems>=3L~"03-04 body systems",
    n_body_systems==2L~"02 body systems",
    n_body_systems==1L~"01 body system",
    n_body_systems==0L~"No body systems")) %>%
  mutate(index_period=case_when(index_date >= as.Date('2018-10-01') && index_date < as.Date('2019-04-01') ~ '01_Oct_Mar_2019',
                                index_date >= as.Date('2019-04-01') && index_date < as.Date('2019-10-01') ~ '02_Apr_Sep_2019',
                                index_date >= as.Date('2019-10-01') && index_date < as.Date('2020-04-01') ~ '03_Oct_Mar_2020',
                                index_date >= as.Date('2020-04-01') && index_date < as.Date('2020-10-01') ~ '04_Apr_Sep_2020',
                                index_date >= as.Date('2010-10-01') && index_date < as.Date('2021-04-01') ~ '05_Oct_Mar_2021',
                                index_date >= as.Date('2021-04-01') && index_date < as.Date('2021-10-01') ~ '06_Apr_Sep_2021',
                                index_date >= as.Date('2021-10-01') && index_date < as.Date('2022-04-01') ~ '07_Oct_Mar_2022',
                                index_date >= as.Date('2022-04-01') && index_date <= as.Date('2022-06-01') ~ '08_Apr_Jun1_2022'
  )) %>%
  mutate(total_pre_1y=case_when(total_visits <= 4 ~ '01_to_04_visits',
                                total_visits > 4 && total_visits <= 9 ~ '05_to_09_visits',
                                total_visits > 9 && total_visits <= 16 ~ '10_to_16_visits',
                                total_visits > 16 && total_visits <= 30 ~ '17_to_30_visits',
                                total_visits > 30 ~ '31+ visits')) %>%
  mutate(remove_asthma_aderm=case_when(cohort=='non_atopic_dx_no_dupilumab' &&
                                         asthma_severity_loose %in% c('mild_loose','moderate_severe','other') ~ 1L,
                                       cohort=='non_atopic_dx_no_dupilumab' &&
                                         aderm_severity %in% c('mild','moderate','severe') ~ 1L,
                                       TRUE ~ 0L)) %>%
  mutate(
    sex_cat=case_when(#gender_concept_id==8507L ~ 'Male',
      gender_concept_id==8532L ~ 'Female',
      TRUE ~ 'Male or Other/unknown/ambiguous')) %>%
  mutate(
    raceth_cat=case_when(ethnicity_concept_id == 38003563L ~ 'Hispanic',
                         ethnicity_concept_id == 38003564L & race_concept_id == 8527L ~ 'NH_White',
                         ethnicity_concept_id == 38003564L & race_concept_id == 8516L ~ 'NH_Black/AA',
                         ethnicity_concept_id == 38003564L & race_concept_id == 8515L ~ 'NH_Asian',
                         ethnicity_concept_id == 38003564L & race_concept_id == 44814659L ~ 'NH_Other_or_Multiple_Race', #NH multiple race
                         ethnicity_concept_id == 38003564L & race_concept_id %in% c(8567L, 8557L) ~ 'NH_Other_or_Multiple_Race', #NH other
                         TRUE ~ 'Other/Unknown')) %>%
  mutate(
    age_group = case_when(#age < 1 ~ '<1',
      age < 6 ~ '00 to 05',
      age < 10 ~ '06 to 09',
      age < 14 ~ '10 to 13',
      age < 18 ~ '14 to 17',
      TRUE ~ 'Old or Unknown')) %>%
  #' site_mask information is removed
  mutate(t1_outcome_any=case_when(t1_outcome_1+t1_outcome_2+t1_outcome_3+t1_outcome_4+
                                    t1_outcome_5+t1_outcome_6+t1_outcome_7+t1_outcome_8+
                                    t1_outcome_9+t1_outcome_10+t1_outcome_11+t1_outcome_12+
                                    t1_outcome_13+t1_outcome_14+t1_outcome_15+t1_outcome_16+
                                    t1_outcome_17+t1_outcome_18+t1_outcome_19+t1_outcome_20+
                                    t1_outcome_21+t1_outcome_22+t1_outcome_23+t1_outcome_24+
                                    t1_outcome_25+t1_outcome_26 > 0L ~ 1L,
                                  TRUE ~ 0L)) %>%
  mutate(t1_wo_skin=case_when(t1_outcome_1+t1_outcome_2+t1_outcome_4+
                                t1_outcome_5+t1_outcome_6+t1_outcome_7+t1_outcome_8+
                                t1_outcome_9+t1_outcome_10+t1_outcome_11+
                                t1_outcome_14+t1_outcome_15+t1_outcome_16+
                                t1_outcome_17+t1_outcome_18+t1_outcome_19+t1_outcome_20+
                                t1_outcome_21+t1_outcome_22+t1_outcome_23+t1_outcome_24+
                                t1_outcome_25+t1_outcome_26 > 0L ~ 1L,
                              TRUE ~ 0L)) %>%
  mutate(t1_skin=case_when(t1_outcome_3+t1_outcome_13+t1_outcome_12>0L ~ 1L,
                           TRUE ~ 0L)) %>%
  mutate(t1_endocrine=case_when(t1_outcome_6+t1_outcome_5+t1_outcome_10+t1_outcome_25>0L ~ 1L,
                                TRUE ~ 0L)) %>%
  mutate(t1_rheumatologic=case_when(t1_outcome_15+t1_outcome_11+t1_outcome_17+t1_outcome_4+t1_outcome_16+t1_outcome_21+t1_outcome_18+t1_outcome_23+t1_outcome_20>0L ~ 1L,
                                    TRUE ~ 0L)) %>%
  mutate(t1_gi=case_when(t1_outcome_1+t1_outcome_2+t1_outcome_19>0L ~ 1L,
                         TRUE ~ 0L)) %>%
  mutate(t1_hematologic=case_when(t1_outcome_14+t1_outcome_26>0L ~ 1L,
                                  TRUE ~ 0L)) %>%
  mutate(t1_neuro=case_when(t1_outcome_9+t1_outcome_7+t1_outcome_8+t1_outcome_22+t1_outcome_24>0L ~ 1L,
                            TRUE ~ 0L)) %>%
  filter(remove_asthma_aderm==0L) %>%
  compute_new(indexes=c('person_id'))


output_tbl(rslt$cht_with_covariates, 'pasthma_aderm_2dx_demog_wcovs',
           indexes=c('person_id'))

##############################################################################
##############################################################################
#' Get attrition
##############################################################################

#' Attrition for ages 6 to <18
rslt$attrition_atopic_dup_6to17 <- get_attrition_atopic_dup_6to17(asthma_or_aderm=results_tbl('pasthma_aderm_2codes'),
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
                                                                  age_upper_limit=18L)

output_tbl(rslt$attrition_atopic_dup_6to17, 'attrition_atopic_dup_6to17')

#' Attrition for ages 6 to <18
rslt$attrition_atopic_nodup_6to17 <- get_attrition_atopic_nodup_6to17(asthma_or_aderm=results_tbl('pasthma_aderm_2codes'),
                                                                      
                                                                      pasthma_2dx=results_tbl('preindex_pasthma_2dx'),
                                                                      aderm_2dx=results_tbl('preindex_aderm_2dx'),
                                                                      
                                                                      pasthma_aderm_nodup=results_tbl('pasthma_aderm_2dx_visit'),
                                                                      pasthma_aderm_nodup_lb=results_tbl('pasthma_aderm_2dx_lblf'),
                                                                      
                                                                      t1_outcomes=results_tbl('cht_w_t1_outcomes_rev'),
                                                                      final_cht=results_tbl('pasthma_aderm_2dx_demog_wcovs'),
                                                                      study_start_date=as.Date('2018-10-01'),
                                                                      study_end_date=as.Date('2022-06-01'),
                                                                      age_lower_limit=6L,
                                                                      age_upper_limit=18L)

output_tbl(rslt$attrition_atopic_nodup_6to17, 'attrition_atopic_nodup_6to17')

#' Attrition for ages 6 to <18
rslt$attrition_nonatopic_nodup_6to17 <- get_attrition_nonatopic_nodup_6to17(nonatopic=results_tbl('no_disease'),
                                                                            nonatopic_lb=results_tbl('no_disease_lblf'),
                                                                            t1_outcomes=results_tbl('cht_w_t1_outcomes_rev'),
                                                                            final_cht=results_tbl('pasthma_aderm_2dx_demog_wcovs'),
                                                                            study_start_date=as.Date('2018-10-01'),
                                                                            study_end_date=as.Date('2022-06-01'),
                                                                            age_lower_limit=6L,
                                                                            age_upper_limit=18L,
                                                                            person_tbl=cdm_tbl('person'),
                                                                            visit_tbl=cdm_tbl('visit_occurrence'))

output_tbl(rslt$attrition_nonatopic_nodup_6to17, 'attrition_nonatopic_nodup_6to17')
