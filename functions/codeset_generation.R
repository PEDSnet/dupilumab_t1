##############################################################################
#' Create Omalizumab/Xolair codeset

rslt$omalizumab <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('302379', '1159566')) %>%
  distinct(concept_id)

rslt$omalizumab_str <- vocabulary_tbl('concept') %>%
  filter(str_detect(lower(concept_name),'omalizumab|xolair')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  distinct(concept_id)

rslt$omalizumab_desc <- rslt$omalizumab %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  dplyr::union(rslt$omalizumab) %>%
  dplyr::union(rslt$omalizumab_str) %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  compute_new()

write_csv(rslt$omalizumab_desc %>% as_data_frame(), 'specs/omalizumab_codes.csv')

##############################################################################
#' Create atopic dermatitis codeset

rslt$aderm <- load_codeset('allergies') %>%
  filter(sub_cluster=='Atopic_dermatitis')

rslt$aderm_codes_other <- vocabulary_tbl('concept') %>%
  filter(str_detect(lower(concept_name),'atopic')) %>%
  filter(str_detect(lower(concept_name),'dermatitis')) %>%
  filter(vocabulary_id %in% c('ICD10CM','IC10','SNOMED'))

rslt$aderm_codes_combo <- rslt$aderm_codes %>%
  dplyr::union(rslt$aderm_codes_other) %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  compute_new()

rslt$aderm_codes_mapped <- rslt$aderm_codes_combo %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept_relationship'), by=c('concept_id'='concept_id_1')) %>%
  filter(relationship_id=='Maps to') %>%
  distinct(concept_id_2) %>%
  rename(concept_id=concept_id_2) %>%
  dplyr::union(rslt$aderm_codes_combo %>% distinct(concept_id)) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('ICD10CM','IC10','SNOMED')) %>%
  compute_new()

rslt$aderm_codes_desc <- rslt$pas_codes_mapped %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept_ancestor'),by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  dplyr::union(rslt$aderm_codes_mapped %>% distinct(concept_id)) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('ICD10CM','IC10','SNOMED')) %>%
  filter(!is.na(concept_id)) %>%
  compute_new()

rslt$aderm_codes_desc %>% as.data.frame() %>% write_csv('specs/atopic_dermatitis_codeset.csv')

################################################################################
################################################################################
#' Checking creation of Dupilumab codeset
################################################################################  
################################################################################
#' Dupilumab codeset creation, uncrestricted by RxNorm

rslt$dup2 <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('1876376','1876399','2375326','1876400')) %>%
  distinct(concept_id)

rslt$dup_str <- vocabulary_tbl('concept') %>%
  filter(str_detect(lower(concept_name),'dupilumab')) %>%
  distinct(concept_id)

rslt$dup2_desc <- rslt$dup2 %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  select(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  dplyr::union(rslt$dup2) %>%
  dplyr::union(rslt$dup_str) %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  compute_new()

write_csv(rslt$dup2_desc %>% as.data.frame(), 'specs/dupilumab_codes2.csv')


##############################################################################
#' Create codeset of patients with persistent asthma

rslt$pas_codes <- vocabulary_tbl('concept') %>%
  filter(str_detect(lower(concept_code),'j45.3|j45.4|j45.5')) %>%
  filter(vocabulary_id %in% c('ICD10CM','IC10')) %>%
  compute_new()

rslt$pas_codes_other <- vocabulary_tbl('concept') %>%
  filter(str_detect(lower(concept_name),'persistent')) %>%
  filter(str_detect(lower(concept_name),'asthma')) %>%
  filter(!str_detect(lower(concept_name),'intermittent')) %>%
  filter(vocabulary_id %in% c('ICD10CM','IC10','SNOMED')) %>%
  compute_new()

rslt$pas_codes_combo <- rslt$pas_codes %>%
  dplyr::union(rslt$pas_codes_other) %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  compute_new()

rslt$pas_codes_mapped <- rslt$pas_codes_combo %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept_relationship'), by=c('concept_id'='concept_id_1')) %>%
  filter(relationship_id=='Maps to') %>%
  distinct(concept_id_2) %>%
  rename(concept_id=concept_id_2) %>%
  dplyr::union(rslt$pas_codes_combo %>% distinct(concept_id)) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('ICD10CM','IC10','SNOMED')) %>%
  compute_new()

rslt$pas_codes_desc <- rslt$pas_codes_mapped %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept_ancestor'),by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  dplyr::union(rslt$pas_codes_mapped %>% distinct(concept_id)) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(!str_detect(lower(concept_name),'intermittent')) %>%
  filter(vocabulary_id %in% c('ICD10CM','IC10','SNOMED')) %>%
  filter(!is.na(concept_id)) %>%
  compute_new()

rslt$pas_codes_desc %>% as.data.frame() %>% write_csv('specs/persistent_asthma_codeset.csv')

################################################################################

#' IgE-mediated food allergy, eosinophilic esophagitis, allergic rhinitis and eczema
rslt$food_allergy <- read_csv('specs/food_allergy_concepts.csv') %>% copy_to_new(dest=config('db_src'), name='food_allergy')
rslt$esophagitis <- read_csv('specs/esophagitis.csv') %>% copy_to_new(dest=config('db_src'), name='esophagitis')
rslt$allergic_rhinitis <- read_csv('specs/allergic_rhinitis_codes.csv') %>% copy_to_new(dest=config('db_src'), name='allergic_rhinitis')
rslt$eczema <- read_csv('specs/eczema_concepts.csv') %>% copy_to_new(dest=config('db_src'), name='eczema')

rslt$esophagitis_codes <- rslt$esophagitis %>%
  distinct(concept_id) %>%
  #left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  left_join(vocabulary_tbl('concept_relationship'), by=c('concept_id'='concept_id_1')) %>%
  filter(relationship_id=='Maps to') %>%
  distinct(concept_id_2) %>%
  rename(concept_id=concept_id_2) %>%
  dplyr::union(rslt$esophagitis %>% distinct(concept_id)) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('ICD10CM','IC10','SNOMED')) %>%
  as_data_frame()

write_csv(rslt$esophagitis_codes, 'specs/esophagitis_icd_snomed.csv')

rslt$allergic_rhinitis_codes <- rslt$allergic_rhinitis %>%
  distinct(concept_id) %>%
  #left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  left_join(vocabulary_tbl('concept_relationship'), by=c('concept_id'='concept_id_1')) %>%
  filter(relationship_id %in% c('Maps to','Mapped from')) %>%
  distinct(concept_id_2) %>%
  rename(concept_id=concept_id_2) %>%
  dplyr::union(rslt$allergic_rhinitis %>% distinct(concept_id)) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('ICD10CM','IC10','SNOMED')) %>%
  as_data_frame()

write_csv(rslt$allergic_rhinitis_codes, 'specs/allergic_rhinitis_icd_snomed.csv')

rslt$food_allergy_codes <- rslt$food_allergy %>%
  distinct(concept_id) %>%
  #left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  left_join(vocabulary_tbl('concept_relationship'), by=c('concept_id'='concept_id_1')) %>%
  filter(relationship_id %in% c('Maps to','Mapped from')) %>%
  distinct(concept_id_2) %>%
  rename(concept_id=concept_id_2) %>%
  dplyr::union(rslt$food_allergy %>% distinct(concept_id)) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('ICD10CM','IC10','SNOMED')) %>%
  as_data_frame()

write_csv(rslt$food_allergy_codes, 'specs/food_allergy_icd_snomed.csv')

rslt$eczema_codes <- rslt$eczema %>%
  distinct(concept_id) %>%
  #left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  left_join(vocabulary_tbl('concept_relationship'), by=c('concept_id'='concept_id_1')) %>%
  filter(relationship_id %in% c('Maps to','Mapped from')) %>%
  distinct(concept_id_2) %>%
  rename(concept_id=concept_id_2) %>%
  dplyr::union(rslt$eczema %>% distinct(concept_id)) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('ICD10CM','IC10','SNOMED')) %>%
  as_data_frame()

write_csv(rslt$eczema_codes, 'specs/eczema_icd_snomed.csv')

##############################################################################

rslt$food_allergy_codeset <- load_codeset('food_allergy_icd_snomed')
rslt$allergic_rhinitis_codeset <- load_codeset('allergic_rhinitis_icd_snomed')
rslt$esophagitis_codeset <- load_codeset('esophagitis_icd_snomed')

rslt$food_allergy <- get_cond_records(codeset=rslt$food_allergy_codeset)
rslt$allergic_rhinitis <- get_cond_records(codeset=rslt$allergic_rhinitis_codeset)
rslt$esophagitis <- get_cond_records(codeset=rslt$esophagitis_codeset)


output_tbl(rslt$food_allergy_cond_src, 'food_allergy_dx', indexes=c('person_id'))

##############################################################################

celiac <- load_codeset('outcome_codes/dx_celiac') %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=1L) %>%
  mutate(outcome_set_name='Celiac Disease') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

crohns <- load_codeset('outcome_codes/crohns_dx') %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=2L) %>%
  mutate(outcome_set_name='Crohns Disease') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

psoriasis <- load_codeset('outcome_codes/dx_psoriasis') %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=3L) %>%
  mutate(outcome_set_name='Psoriasis') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

lupus <- load_codeset('outcome_codes/combined_lupus') %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=4L) %>%
  mutate(outcome_set_name='Lupus') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

t1d <- load_codeset('outcome_codes/T1DM') %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=5L) %>%
  mutate(outcome_set_name='Type 1 Diabetes') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

hashimotos <- load_codeset('outcome_codes/hashimoto_prelim') %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=6L) %>%
  mutate(outcome_set_name='Hashimotos') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

ms <- load_codeset('outcome_codes/ms') %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=7L) %>%
  mutate(outcome_set_name='Multiple Sclerosis') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

nmo <- load_codeset('outcome_codes/neuromyelitis_optica') %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=8L) %>%
  mutate(outcome_set_name='Neuromyelitis Optica') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

ad_enc <- load_codeset('outcome_codes/dx_encephalitis') %>%
  filter(str_detect(lower(concept_name),'acute')) %>%
  filter(str_detect(lower(concept_name),'disseminated')) %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=9L) %>%
  mutate(outcome_set_name='Acute Disseminated Encephaloymelitis') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

addison <- load_codeset('outcome_codes/dx_addison') %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=10L) %>%
  mutate(outcome_set_name='Addisons Disease') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

sjogren <- load_codeset('outcome_codes/dx_sjogren') %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=11L) %>%
  mutate(outcome_set_name='Sjorens Disease') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

alopecia <- load_codeset('outcome_codes/dx_alopecia_areata') %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=12L) %>%
  mutate(outcome_set_name='Alopecia') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

vitiligo <- load_codeset('outcome_codes/dx_vitiligo') %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=13L) %>%
  mutate(outcome_set_name='Vitiligo') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

itp <- load_codeset('outcome_codes/dx_itp') %>%
  select(concept_id, concept_code, concept_name, vocabulary_id) %>%
  mutate(outcome_set_id=14L) %>%
  mutate(outcome_set_name='Immune thrombocytopenia') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

arthritis <- load_codeset('outcome_codes/arthritis') %>%
  mutate(type=case_when(str_detect(lower(concept_name),'rheumatoid') ~ 'rheumatoid',
                        str_detect(lower(concept_name),'psoria') ~ 'psoriatic',
                        TRUE ~ 'other_unknown'
  ))

rheumatoid_arthritis <- arthritis %>%
  filter(type=='rheumatoid') %>%
  select(concept_id, concept_name, concept_code, vocabulary_id) %>%
  mutate(outcome_set_id=15L) %>%
  mutate(outcome_set_name='Rheumatoid Arthritis') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

psoriatic_arthritis <- arthritis %>%
  filter(type=='psoriatic') %>%
  select(concept_id, concept_name, concept_code, vocabulary_id) %>%
  mutate(outcome_set_id=16L) %>%
  mutate(outcome_set_name='Psoriatic Arthritis') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))

inflammatory_arthritis <- vocabulary_tbl('concept') %>%
  filter(vocabulary_id %in% c('ICD10','ICD10CM','SNOMED')) %>%
  filter(str_detect(lower(concept_name),'inflammatory arthritis')) %>%
  filter(!str_detect(lower(concept_name),'without')) %>%
  filter(!str_detect(lower(concept_name),'suspected')) %>%
  select(concept_id, concept_name, concept_code, vocabulary_id) %>%
  mutate(outcome_set_id=17L) %>%
  mutate(outcome_set_name='Inflammatory Arthritis') %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))



ad_concepts <- read_csv('specs/outcome_codes/ad_concepts_102423.csv') %>%
  filter(!is.na(Name)) %>%
  rename(concept_set_id=`Concept Set ID`,
         concept_code=`Concept Code`,
         concept_name=`Concept Name`,
         vocabulary_id=Vocabulary,
         concept_set_name=Name) %>%
  select(concept_set_id, concept_set_name, concept_code,
         concept_name, vocabulary_id) %>%
  select(-concept_name) %>%
  copy_to_new(dest=config('db_src'), name='ad_concepts') %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_code','vocabulary_id')) %>%
  mutate(concept_id=as.numeric(concept_id),
         concept_code=as.character(concept_code),
         concept_name=as.character(concept_name),
         vocabulary_id=as.character(vocabulary_id))


ankylosing_spondylitis <- ad_concepts %>%
  filter(concept_set_id==1) %>%
  select(concept_id, concept_name, concept_code, vocabulary_id) %>%
  mutate(outcome_set_id=18L) %>%
  mutate(outcome_set_name='Ankylosing Spondylitis')

ulcerative_colitis <- ad_concepts %>%
  filter(str_detect(lower(concept_name),'ulcer')) %>%
  filter(str_detect(lower(concept_name),'colitis')) %>%
  select(concept_id, concept_name, concept_code, vocabulary_id) %>%
  mutate(outcome_set_id=19L) %>%
  mutate(outcome_set_name='Ulcerative Colitis')

dermatomyositis <- ad_concepts %>%
  filter(str_detect(lower(concept_name),'dermatomyositis')) %>%
  select(concept_id, concept_name, concept_code, vocabulary_id) %>%
  mutate(outcome_set_id=20L) %>%
  mutate(outcome_set_name='Dermatomyositis')

polymyositis <- ad_concepts %>%
  filter(str_detect(lower(concept_name),'polymyositis')) %>%
  select(concept_id, concept_name, concept_code, vocabulary_id) %>%
  mutate(outcome_set_id=21L) %>%
  mutate(outcome_set_name='Polymyositis')

tm <- ad_concepts %>%
  filter(str_detect(lower(concept_name),'transverse')) %>%
  filter(str_detect(lower(concept_name),'myelitis')) %>%
  select(concept_id, concept_name, concept_code, vocabulary_id) %>%
  mutate(outcome_set_id=22L) %>%
  mutate(outcome_set_name='Transverse Myelitis')

behcet <- ad_concepts %>%
  filter(concept_set_id==44) %>%
  select(concept_id, concept_name, concept_code, vocabulary_id) %>%
  mutate(outcome_set_id=23L) %>%
  mutate(outcome_set_name='Behcets Disease')

mg <- ad_concepts %>%
  filter(concept_set_id==31) %>%
  select(concept_id, concept_name, concept_code, vocabulary_id) %>%
  mutate(outcome_set_id=24L) %>%
  mutate(outcome_set_name='Myasthenia Gravis')

graves <- ad_concepts %>%
  filter(concept_set_id==13) %>%
  select(concept_id, concept_name, concept_code, vocabulary_id) %>%
  mutate(outcome_set_id=25L) %>%
  mutate(outcome_set_name='Graves Disease')

anemia <- ad_concepts %>%
  filter(concept_set_id==22) %>%
  select(concept_id, concept_name, concept_code, vocabulary_id) %>%
  mutate(outcome_set_id=26L) %>%
  mutate(outcome_set_name='Anemia')



all_codes <- celiac %>%
  dplyr::union(crohns) %>%
  dplyr::union(psoriasis) %>%
  dplyr::union(lupus) %>%
  dplyr::union(t1d) %>%
  dplyr::union(hashimotos) %>%
  dplyr::union(ms) %>%
  dplyr::union(nmo) %>%
  dplyr::union(ad_enc) %>%
  dplyr::union(addison) %>%
  dplyr::union(sjogren) %>%
  dplyr::union(alopecia) %>%
  dplyr::union(vitiligo) %>%
  dplyr::union(itp) %>%
  dplyr::union(rheumatoid_arthritis) %>%
  dplyr::union(psoriatic_arthritis) %>%
  dplyr::union(inflammatory_arthritis) %>%
  dplyr::union(ankylosing_spondylitis) %>%
  dplyr::union(ulcerative_colitis) %>%
  dplyr::union(dermatomyositis) %>%
  dplyr::union(polymyositis) %>%
  dplyr::union(tm) %>%
  dplyr::union(behcet) %>%
  dplyr::union(mg) %>%
  dplyr::union(graves) %>%
  dplyr::union(anemia) %>%
  collect() %>%
  arrange(outcome_set_id, concept_id) %>%
  as_data_frame()

write_csv(all_codes, 'all_outcome_codes.csv')

all_outcome_codes <- read_csv('specs/outcome_codes/all_outcome_codes.csv') %>%
  copy_to_new(dest=config('db_src'))

map_to_snomed_icd <- all_outcome_codes %>%
  left_join(vocabulary_tbl('concept_relationship'), by=c('concept_id'='concept_id_1')) %>%
  filter(relationship_id %in% c('Maps to','Mapped from')) %>%
  compute_new(indexes=c('concept_id'))

mapped_concepts <- map_to_snomed_icd %>%
  select(concept_id_2, outcome_set_name, outcome_set_id) %>%
  rename(concept_id=concept_id_2) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  select(concept_id, concept_code, concept_name, vocabulary_id, outcome_set_id, outcome_set_name) %>%
  filter(vocabulary_id %in% c('ICD10','ICD10CM','SNOMED')) %>%
  compute_new(indexes=c('concept_id'))

all_codes_mapped <- all_outcome_codes %>%
  dplyr::union(mapped_concepts) %>%
  distinct(concept_id, concept_code, concept_name, vocabulary_id, outcome_set_id, outcome_set_name) %>%
  filter(!is.na(concept_code)) %>%
  filter(!is.na(concept_id)) %>%
  compute_new(indexes=c('concept_id'))

all_codes_mapped_collect <- all_codes_mapped %>%
  collect() %>%
  arrange(outcome_set_id, concept_id) %>%
  as_data_frame()
# 
# write_csv(all_codes_mapped_collect, 'all_outcome_codes_icd_snomed.csv')
#' Note: these codes were subsequently manually edited.


##############################################################################
#' Identify patients with prescriptions for albuterol, inhaled corticosteroid,
#' eczema ointments, epinephrine

#' Albuterol
rslt$albuterol <- load_codeset('rx_albuterol') %>%
  distinct(concept_id) %>%
  mutate(drug_type='albuterol')

#' Inhaled corticosteroids
rslt$asthma_meds <- load_codeset('asthma_meds_severity')

rslt$asthma_ics <- rslt$asthma_meds %>%
  filter(drug_type=='ics') %>%
  distinct(concept_id) %>%
  mutate(drug_type='ics')

#' eczema ointment: https://nationaleczema.org/eczema/treatment/topicals/
#' topical corticosteroids
#' topical PDE4 inhibitors: crisaborole
#' topical JAK inhibitors: opzelura
#' topical calcineurin inhibitors

rslt$topical_corticosteroids <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('2660029', '2659863', '1160817', '1161449', '2659417',
                             '1157140', '1154232', '2659280', '1160295',
                             '2659633', '2659280', '2661133', '2661051', '2662286', '2662650', '2661871','2661113',
                             '2663064', '1154230', '2661133', '1164013', '2661113')) %>%
  filter(vocabulary_id %in% c('RxNorm','RxNorm Extension'))

rslt$topical_corticosteroids_string <- vocabulary_tbl('concept') %>%
  filter(str_detect(lower(concept_name),'betamethasone diproprionate|clobetasol propionate|fluocinonide|flurandrenolide|halobetasol propionate|
                      amcinonide|desoximetasone|diflorasone diacetate|fluocinonide|halcinonide|
                      betamethasone valerate|fluocinolone acetonide|flurandrenolide|fluticasone propionate|hydrocortisone butyrate|hydrocortisone valerate|
                      mometasone furoate|triamcinolone acetonide|alclometasone dipropionate|desonide|fluocinolone acetonide|hydrocortisone')) %>%
  filter(str_detect(lower(concept_name),'topical|cream|gel|ointment|lotion')) %>% 
  filter(domain_id=='Drug') %>%
  filter(vocabulary_id %in% c('RxNorm','RxNorm Extension'))

rslt$topical_corticosteroids_all <- rslt$topical_corticosteroids %>% distinct(concept_id) %>%
  dplyr::union(rslt$topical_corticosteroids_string %>% distinct(concept_id)) %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  distinct(concept_id) %>%
  dplyr::union(rslt$topical_corticosteroids %>% distinct(concept_id)) %>%
  dplyr::union(rslt$topical_corticosteroids_string %>% distinct(concept_id)) %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(str_detect(lower(concept_name),'topical|cream|gel|ointment|lotion')) %>% 
  filter(domain_id=='Drug') %>%
  filter(vocabulary_id %in% c('RxNorm','RxNorm Extension')) %>%
  distinct(concept_id) %>%
  mutate(drug_type='topical_corticosteroids') %>%
  compute_new(indexes=c('concept_id'))

rslt$eczema_ointment_pde4_jak <- vocabulary_tbl('concept') %>%
  filter(str_detect(lower(concept_name),'crisaborole|opzelura')) %>%
  filter(vocabulary_id %in% c('RxNorm','RxNorm Extension'))

rslt$eczema_ointment_pde4_jak_all <- rslt$eczema_ointment_pde4_jak %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  distinct(concept_id) %>%
  dplyr::union(rslt$eczema_ointment_pde4_jak %>% distinct(concept_id)) %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept'),by=c('concept_id')) %>%
  distinct(concept_id) %>%
  compute_new(indexes=c('concept_id'))

rslt$eczema_ointment_cnis <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('1158739', '1164204'))

rslt$eczema_ointment_cnis_all <- rslt$eczema_ointment_cnis %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  distinct(concept_id) %>%
  dplyr::union(rslt$eczema_ointment_cnis %>% distinct(concept_id)) %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept'),by=c('concept_id'))

rslt$eczema_ointment <- rslt$topical_corticosteroids_all %>%
  dplyr::union(rslt$eczema_ointment_pde4_jak_all) %>%
  dplyr::union(rslt$eczema_ointment_cnis_all) %>%
  distinct(concept_id) %>%
  mutate(drug_type='eczema_topical') %>%
  compute_new(indexes=c('concept_id'))


#' Epinephrine
rslt$epinephrine <- vocabulary_tbl('concept') %>%
  filter(concept_id %in% c(44814535L, 1344014L, 40038600L,
                           19076867L, 19076899L))

rslt$epinephrine2 <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('3992', '1163886', '1163887', '1163888',
                             '1163889', '1163891', '1163892')) %>%
  filter(vocabulary_id %in% c('RxNorm','RxNorm Extension'))

rslt$epinephrine_all <- rslt$epinephrine %>% distinct(concept_id) %>%
  dplyr::union(rslt$epinephrine2 %>% distinct(concept_id)) %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(str_detect(lower(concept_name),'epinephrine')) %>%
  distinct(concept_id) %>%
  dplyr::union(rslt$epinephrine %>% distinct(concept_id)) %>%
  dplyr::union(rslt$epinephrine2 %>% distinct(concept_id)) %>%
  distinct(concept_id) %>%
  mutate(drug_type='epinephrine')

rslt$asthma_aderm_px <- rslt$albuterol %>%
  dplyr::union(rslt$asthma_ics) %>%
  dplyr::union(rslt$eczema_ointment) %>%
  dplyr::union(rslt$epinephrine_all) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id'))

write_csv(rslt$asthma_aderm_px %>% as_data_frame(),'specs/atopic_px.csv')

##############################################################################
#' Try getting a limited set of all eczema ointment concepts,
#' only to aid with ease of codeset review

rslt$topical_corticosteroids <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('2660029', '2659863', '1160817', '1161449', '2659417',
                             '1157140', '1154232', '2659280', '1160295',
                             '2659633', '2659280', '2661133', '2661051', '2662286', '2662650', '2661871','2661113',
                             '2663064', '1154230', '2661133', '1164013', '2661113')) %>%
  filter(vocabulary_id %in% c('RxNorm','RxNorm Extension'))

rslt$topical_corticosteroids_string <- vocabulary_tbl('concept') %>%
  filter(str_detect(lower(concept_name),'betamethasone diproprionate|clobetasol propionate|fluocinonide|flurandrenolide|halobetasol propionate|
                      amcinonide|desoximetasone|diflorasone diacetate|fluocinonide|halcinonide|
                      betamethasone valerate|fluocinolone acetonide|flurandrenolide|fluticasone propionate|hydrocortisone butyrate|hydrocortisone valerate|
                      mometasone furoate|triamcinolone acetonide|alclometasone dipropionate|desonide|fluocinolone acetonide|hydrocortisone')) %>%
  filter(str_detect(lower(concept_name),'topical|cream|gel|ointment|lotion')) %>% 
  filter(domain_id=='Drug') %>%
  filter(vocabulary_id %in% c('RxNorm','RxNorm Extension')) %>%
  compute_new(indexes=c('concept_id'))

rslt$topical_corticosteroids_string_nodesc <- rslt$topical_corticosteroids_string %>%
  filter(!str_detect(lower(concept_name),'mg/mg')) %>% 
  filter(!str_detect(lower(concept_name),'mg/ml')) %>% 
  compute_new(indexes=c('concept_id'))

rslt$topical_corticosteroids_all <- rslt$topical_corticosteroids %>% distinct(concept_id) %>%
  dplyr::union(rslt$topical_corticosteroids_string_nodesc %>% distinct(concept_id)) %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  mutate(category='Topical Corticosteroid') %>%
  compute_new(indexes=c('concept_id'))

rslt$eczema_ointment_pde4 <- vocabulary_tbl('concept') %>%
  filter(str_detect(lower(concept_name),'crisaborole')) %>%
  filter(vocabulary_id %in% c('RxNorm','RxNorm Extension'))

rslt$eczema_ointment_pde4_all <- rslt$eczema_ointment_pde4 %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept'),by=c('concept_id')) %>%
  mutate(category='Topical PDE4 Inhibitor') %>%
  compute_new(indexes=c('concept_id'))

rslt$eczema_ointment_jak <- vocabulary_tbl('concept') %>%
  filter(str_detect(lower(concept_name),'opzelura')) %>%
  filter(vocabulary_id %in% c('RxNorm','RxNorm Extension'))

rslt$eczema_ointment_jak_all <- rslt$eczema_ointment_jak %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept'),by=c('concept_id')) %>%
  mutate(category='Topical JAK Inhibitor') %>%
  compute_new(indexes=c('concept_id'))

rslt$eczema_ointment_cnis <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('1158739', '1164204')) %>%
  mutate(category='Topical Calcineurin Inhibitor') %>%
  left_join(vocabulary_tbl('concept'),by=c('concept_id'))

rslt$eczema_ointment_condensed <- rslt$topical_corticosteroids_all %>%
  dplyr::union(rslt$eczema_ointment_pde4_all) %>%
  dplyr::union(rslt$eczema_ointment_jak_all) %>%
  dplyr::union(rslt$eczema_ointment_cnis) %>%
  distinct(concept_id, category) %>%
  left_join(vocabulary_tbl('concept'),by=c('concept_id')) %>%
  compute_new(indexes=c('concept_id'))

rslt$eczema_ointment_condensed %>% collect() %>% as_data_frame() %>%
  write_csv('specs/outcome_codes/eczema_ointments_condensed_codes.csv')

##############################################################################
##############################################################################
##############################################################################
#' Eczema ointments: add the following:

#' Mild:
#' Alcomethasone
alcometasone <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('108088', '1154609', '378616', '2648469')) %>% 
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension'))

alcometasone_desc <- alcometasone %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  filter(!is.na(concept_id))

mild_codes <- alcometasone %>%
  dplyr::union(alcometasone_desc) %>%
  distinct(concept_id) %>%
  mutate(drug='alcometasone') %>%
  mutate(severity='mild')

#' Moderate:
#' Tacrolimus 0.1% and 0.03%: not sure how to differentiate
tacrolimus_topical <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('1164204','379082'))

tacrolimus_desc <- tacrolimus_topical %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  filter(!is.na(concept_id))

tacrolimus_topical_all <- tacrolimus_topical %>%
  dplyr::union(tacrolimus_desc) %>%
  mutate(drug='tacrolimus') %>%
  mutate(severity='moderate')

#' pimecrolimus-- includes Elidel
pimecrolimus_topical <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('562805','1158739'))

pimecrolimus_desc <- pimecrolimus_topical %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  filter(!is.na(concept_id))

pimecrolimus_topical_all <- pimecrolimus_topical %>%
  dplyr::union(pimecrolimus_desc) %>%
  mutate(drug='pimecrolimus') %>%
  mutate(severity='moderate')

#' Amcinonide
amcinonide_topical <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('1157140','377850','377851','374449'))

amcinonide_desc <- amcinonide_topical %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  filter(!is.na(concept_id))

amcinonide_topical_all <- amcinonide_topical %>%
  dplyr::union(amcinonide_desc) %>%
  mutate(drug='amcinonide') %>%
  mutate(severity='moderate')

#' Fluocinonide (can't limit to 0.05%)
fluocinonide_topical <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('378617','379077','372222','385001','1160817')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension'))

fluocinonide_desc <- fluocinonide_topical %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  filter(!is.na(concept_id))

fluocinonide_topical_all <- fluocinonide_topical %>%
  dplyr::union(fluocinonide_desc) %>%
  mutate(drug='fluocinonide') %>%
  mutate(severity='severe')

#' Halcinonide
halcinonide_topical <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('1160295','377707','372370','372371')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension'))

halcinonide_desc <- halcinonide_topical %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  filter(!is.na(concept_id))

halcinonide_topical_all <- halcinonide_topical %>%
  dplyr::union(halcinonide_desc) %>%
  mutate(drug='halcinonide') %>%
  mutate(severity='severe')

#' Betamethasone valerate
betval <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('227897','384499','2648726','2649020','1153804','2659633',
                             '19254','2471872', '1716093', '2646254',
                             '2646228', '2648710', '1153789', '2660029')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension'))

betval_desc <- betval %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  filter(!is.na(concept_id)) %>%
  filter(str_detect(lower(concept_name),'betamethasone')) %>%
  filter(str_detect(lower(concept_name),'valerate'))

betval_all <- betval %>%
  filter(str_detect(lower(concept_name),'betamethasone')) %>%
  filter(str_detect(lower(concept_name),'valerate')) %>%
  dplyr::union(betval_desc) %>%
  distinct(concept_id) %>%
  mutate(drug='betamethasone_valerate') %>%
  mutate(severity='moderate')

#' Desoximetasone
desoximetasone_topical <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('1154232','377670','385104','374373','1376337')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension'))

desoximetasone_desc <- desoximetasone_topical %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  filter(!is.na(concept_id))

desoximetasone_topical_all <- desoximetasone_topical %>%
  dplyr::union(desoximetasone_desc) %>%
  mutate(drug='desoximetasone') %>%
  mutate(severity='severe')

#' Fluticasone
fluticasone_topical <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('1165657','378620','577392','372254')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension'))

fluticasone_topical_desc <- fluticasone_topical %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  filter(!is.na(concept_id))

fluticasone_topical_all <- fluticasone_topical %>%
  dplyr::union(fluticasone_topical_desc) %>%
  distinct(concept_id) %>%
  mutate(drug='fluticasone') %>%
  mutate(severity='moderate')

#' Triamcinolone
triamcinolone_topical <- read_csv('specs/triamcinolone_topical.csv') %>%
  copy_to_new(dest=config('db_src'), name='triamcinolone_topical')
triamcinolone_paste <- read_csv('specs/triamcinolone_paste.csv') %>%
  copy_to_new(dest=config('db_src'), name='triamcinolone_paste')

triamcinolone_all <- triamcinolone_topical %>%
  dplyr::union(triamcinolone_paste) %>%
  distinct(concept_id) %>%
  mutate(drug='triamcinolone') %>%
  mutate(severity='moderate')

#' Mometasone
mometasone_topical <- read_csv('specs/mometasone_topical.csv') %>%
  copy_to_new(dest=config('db_src'), name='mometasone')

mometasone_all <- mometasone_topical %>%
  distinct(concept_id) %>%
  mutate(drug='mometasone') %>%
  mutate(severity='moderate')

#' Hydrocortisone butyrate
hydrocortisone_butyrate <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('103468','2648540','2645482','2646099','2647148','2662286')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension'))

hydrocortisone_butyrate_desc <- hydrocortisone_butyrate %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  filter(!is.na(concept_id))

hydrocortisone_butyrate_all <- hydrocortisone_butyrate %>%
  dplyr::union(hydrocortisone_butyrate_desc) %>%
  distinct(concept_id) %>%
  mutate(drug='hydrocortisone_butyrate') %>%
  mutate(severity='moderate')

moderate_codes <- tacrolimus_topical_all %>%
  dplyr::union(pimecrolimus_topical_all) %>%
  dplyr::union(amcinonide_topical_all) %>%
  dplyr::union(fluocinonide_topical_all) %>%
  dplyr::union(halcinonide_topical_all) %>%
  dplyr::union(betval_all) %>%
  dplyr::union(desoximetasone_topical_all) %>%
  dplyr::union(fluticasone_topical_all) %>%
  dplyr::union(triamcinolone_all) %>%
  dplyr::union(mometasone_all) %>%
  dplyr::union(hydrocortisone_butyrate_all)

#############################################################################
#' Severe:
#' betamethasone dipropionate

# betdip <- vocabulary_tbl('concept') %>%
#   filter(concept_code %in% c('19254','2471872', '1716093', '2646254',
#                              '2646228', '2648710', '1153789', '2660029')) %>%
#   filter(vocabulary_id=='RxNorm', 'RxNorm Extension')
# 
# betdip_desc <- betdip %>%
#   left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
#   distinct(descendant_concept_id) %>%
#   rename(concept_id=descendant_concept_id) %>%
#   left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
#   filter(!is.na(concept_id)) %>%
#   filter(str_detect(lower(concept_name),'betamethasone')) %>%
#   filter(str_detect(lower(concept_name),'dipropionate')) #' not sure if only bethametasone dipropionate should be used
# 
# betdip_all <- betdip %>%
#   dplyr::union(betdip_desc) %>%
#   distinct(concept_id) %>%
#   mutate(drug='betamethasone dipropionate')

betamethasone_topical <- read_csv('specs/betamethasone_topical.csv') %>%
  copy_to_new(dest=config('db_src'), name='betamethasone') %>%
  filter(!str_detect(lower(concept_name),'valerate')) %>%
  distinct(concept_id) %>%
  mutate(drug='betamethasone') %>%
  mutate(severity='severe')

#' Clobetasol
clobetasol <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('2590','1295100','1153083','438709','562031',
                             '379105','379076','379392','378982','544933',
                             '597751')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension'))

clobetasol_desc <- clobetasol %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  filter(!is.na(concept_id))

clobetasol_all <- clobetasol %>%
  dplyr::union(clobetasol_desc) %>%
  distinct(concept_id) %>%
  mutate(drug='clobetasol') %>%
  mutate(severity='severe')

#' Diflorasone
diflorasone <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('91311','378610','376493','1151507')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension'))

diflorasone_desc <- diflorasone %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  filter(!is.na(concept_id))

diflorasone_all <- diflorasone %>%
  dplyr::union(diflorasone_desc) %>%
  distinct(concept_id) %>%
  mutate(drug='diflorasone') %>%
  mutate(severity='severe')

#' Halobetasol
halobetasol <- vocabulary_tbl('concept') %>%
  filter(concept_code %in% c('41208','1158362','378637','2047645','1789961','372372')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension'))

halobetasol_desc <- halobetasol %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  filter(!is.na(concept_id))

halobetasol_all <- halobetasol %>%
  dplyr::union(halobetasol_desc) %>%
  distinct(concept_id) %>%
  mutate(drug='halobetasol') %>%
  mutate(severity='severe')

#' cyclosporine
cyclosporine_all <- read_csv('specs/cyclosporine.csv') %>%
  distinct(concept_id) %>%
  copy_to_new(dest=config('db_src'), name='cyclosporine') %>%
  mutate(drug='cyclosporine') %>%
  mutate(severity='severe')

#' methotrexate
methotrexate_injectable <- read_csv('specs/methotrexate_injectable.csv')
methotrexate_oral_liquid <- read_csv('specs/methotrexate_oral_liquid.csv')
methotrexate_oral <- read_csv('specs/methotrexate_oral.csv')
methotrexate_pill <- read_csv('specs/methotrexate_pill.csv')

methotrexate_all <- methotrexate_injectable %>%
  dplyr::union(methotrexate_oral_liquid) %>%
  dplyr::union(methotrexate_oral) %>%
  dplyr::union(methotrexate_pill) %>%
  distinct(concept_id) %>%
  copy_to_new(dest=config('db_src'), name='methotrexate') %>%
  mutate(drug='methotrexate') %>%
  mutate(severity='severe')


severe_codes <- betamethasone_topical %>%
  dplyr::union(clobetasol_all) %>%
  dplyr::union(diflorasone_all) %>%
  dplyr::union(halobetasol_all) %>%
  dplyr::union(cyclosporine_all) %>%
  dplyr::union(methotrexate_all) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id'))

##############################################################################

#' Atopic px:
mild_moderate_severe <- mild_codes %>%
  dplyr::union(moderate_codes) %>%
  dplyr::union(severe_codes) %>%
  select(concept_id, drug, severity)

atopic_px <- read_csv('specs/atopic_px.csv') %>%
  copy_to_new(dest=config('db_src'), name='atopic_px')

eczema_topical <- atopic_px %>%
  filter(drug_type=='eczema_topical')

eczema_topical_severe <- eczema_topical %>%
  filter(str_detect(lower(concept_name),'opzelura|ruxolitinib|desoximetasone|halcinonide|fluocinonide|clobetasol')) %>%
  mutate(drug=case_when(str_detect(lower(concept_name),'opzelura') ~ 'opzelura',
                        str_detect(lower(concept_name),'ruxolitinib') ~ 'ruxolitinib',
                        str_detect(lower(concept_name),'desoximetasone') ~ 'desoximetasone',
                        str_detect(lower(concept_name),'halcinonide') ~ 'halcinonide',
                        str_detect(lower(concept_name),'fluocinonide') ~ 'fluocinonide',
                        str_detect(lower(concept_name),'clobetasol') ~ 'clobetasol'
  )) %>%
  mutate(severity='severe') %>%
  select(concept_id, drug, severity)

eczema_topical_moderate <- eczema_topical %>%
  anti_join(eczema_topical_severe, by=c('concept_id')) %>%
  filter(str_detect(lower(concept_name),'flurandrenolide|crisaborole|budesonide|coal tar|tacrolimus|pimecrolimus|hydrocortisone valerate|amcinonide|westcort')) %>%
  mutate(drug=case_when(str_detect(lower(concept_name),'flurandrenolide') ~ 'flurandrenolide',
                        str_detect(lower(concept_name),'budesonide') ~ 'budesonide',
                        str_detect(lower(concept_name),'crisaborole') ~ 'crisaborole',
                        str_detect(lower(concept_name),'coal tar') ~ 'coal tar',
                        str_detect(lower(concept_name),'tacrolimus') ~ 'tacrolimus',
                        str_detect(lower(concept_name),'pimecrolimus') ~ 'pimecrolimus',
                        str_detect(lower(concept_name),'amcinonide') ~ 'amcinonide',
                        str_detect(lower(concept_name),'westcort') ~ 'westcort',
                        str_detect(lower(concept_name),'hydrocortisone valerate') ~ 'hydrocortisone_valerate'
  )) %>%
  mutate(severity='moderate') %>%
  select(concept_id, drug, severity)

eczema_topical_mild <- eczema_topical %>%
  anti_join(eczema_topical_moderate, by=c('concept_id')) %>%
  anti_join(eczema_topical_severe, by=c('concept_id')) %>%
  filter(str_detect(lower(concept_name),'desonide|hydrocortisone')) %>%
  mutate(drug=case_when(str_detect(lower(concept_name),'desonide') ~ 'desonide',
                        str_detect(lower(concept_name),'hydrocortisone') ~ 'hydrocortisone'
  )) %>%
  mutate(severity='mild') %>%
  select(concept_id, drug, severity)

eczema_other <- eczema_topical %>%
  anti_join(eczema_topical_mild, by=c('concept_id')) %>%
  anti_join(eczema_topical_moderate, by=c('concept_id')) %>%
  anti_join(eczema_topical_severe, by=c('concept_id')) %>%
  filter(!is.na(concept_id)) %>%
  mutate(drug=case_when(str_detect(lower(concept_name),'halobetasol') ~ 'halobetasol',
                        str_detect(lower(concept_name),'diflorasone') ~ 'diflorasone',
                        str_detect(lower(concept_name),'triamcinolone') ~ 'triamcinolone',
                        str_detect(lower(concept_name),'fluticasone') ~ 'fluticasone',
                        str_detect(lower(concept_name),'elidel') ~ 'pimecrolimus'
  )) %>%
  mutate(severity=case_when(drug %in% c('halobetasol','diflorasone') ~ 'severe',
                            drug %in% c('triamcinolone','fluticasone','pimecrolimus') ~ 'moderate'
  )) %>%
  filter(!is.na(severity)) %>%
  select(concept_id, drug, severity) %>%
  compute_new(indexes=c(list('concept_id')))





eczema_topical_combo <- mild_moderate_severe %>%
  dplyr::union(eczema_topical_severe) %>%
  dplyr::union(eczema_topical_moderate) %>%
  dplyr::union(eczema_other) %>%
  dplyr::union(eczema_topical_mild) %>%
  distinct(concept_id, drug, severity) %>%
  left_join(vocabulary_tbl('concept'),by=c('concept_id')) %>%
  filter(vocabulary_id %in% c('RxNorm', 'RxNorm Extension')) %>%
  compute_new(indexes=list(c('concept_id')))


write_csv(eczema_topical_combo %>% as.data.frame(),'specs/eczema_topical_codes_050324.csv')

#' Look to see counts of codes by drug type and severity level
drug_code_cts <- eczema_topical_combo %>%
  group_by(drug, severity) %>%
  summarise(n=n()) %>%
  ungroup()

#' Look to see if there are any codes classified as multiple severity levels
test <- mild_moderate_severe %>%
  dplyr::union(eczema_topical_severe) %>%
  dplyr::union(eczema_topical_moderate) %>%
  dplyr::union(eczema_other) %>%
  dplyr::union(eczema_topical_mild) %>%
  distinct(concept_id, drug, severity) %>%
  group_by(concept_id, drug) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  filter(n>1)

##############################################################################
##############################################################################
##############################################################################
#' Asthma severity: add the following:

asthma_severity <- read_xlsx('specs/asthma_meds_nodesc_severity.xlsx') %>%
  copy_to_new(dest=config('db_src'),name='asthma_severity_nodesc')

asthma_severity %>% group_by(category) %>% summarise(n=n())

albuterol <- load_codeset('rx_albuterol')

#' Try to make sure combinations are included in combination meds
asthma_combo_codes <- asthma_severity %>%
  filter(category=='ICS_and_LABA') %>%
  distinct(concept_id) %>%
  compute_new()

asthma_combo_codes_desc <- asthma_combo_codes %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  dplyr::union(asthma_combo_codes) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  mutate(has_combo=case_when(str_detect(lower(concept_name),'budesonide')==TRUE && str_detect(lower(concept_name),'formoterol')==TRUE ~ 1L,
                             str_detect(lower(concept_name),'fluticasone')==TRUE && str_detect(lower(concept_name),'salmeterol')==TRUE ~ 1L,
                             str_detect(lower(concept_name),'fluticasone')==TRUE && str_detect(lower(concept_name),'vilanterol')==TRUE ~ 1L,
                             str_detect(lower(concept_name),'formoterol')==TRUE && str_detect(lower(concept_name),'mometasone')==TRUE ~ 1L,
                             TRUE ~ 0L
  )) %>%
  filter(has_combo==1L) %>%
  select(-has_combo) %>%
  compute_new()

#' Do not include combination medications; ICS and LABA; do not include omalizumab-- will be included separately
asthma_categories <- list('zileuton','theophylline','LTRA','LABA','ICS','SABA','Cromolyn')

asthma_codes <- list()
for (a in 1:length(asthma_categories)) {
  
  asthma_category <- asthma_categories[[a]]
  
  asthma_severity_codes <- asthma_severity %>%
    filter(category==asthma_category) %>%
    distinct(concept_id) %>%
    compute_new()
  
  asthma_codes[[a]] <- asthma_severity_codes %>%
    left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
    distinct(descendant_concept_id) %>%
    rename(concept_id=descendant_concept_id) %>%
    dplyr::union(asthma_severity_codes) %>%
    left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
    mutate(category=asthma_category) %>%
    compute_new()
}

final_asthma_codes <- reduce(.x=asthma_codes,
                             .f=dplyr::union)

#' Moderate/severe asthma injection: omalizumab, mepolizumab, benralizumab, tezepelumab, and tralokinumab
asthma_injection <- vocabulary_tbl('concept') %>%
  filter(vocabulary_id %in% c('RxNorm')) %>%
  filter(concept_code %in% c('302379','1159566','2663919','1657208','2058942',
                             '1720597','2170989','1720600','2173819','1720599',
                             '1989100','2205631','1989103','1989102',
                             '2587789','2628919','2587793','2587802','2587792',
                             '2589225','2589377','2589376'))


asthma_injection_desc <- asthma_injection %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept_ancestor'), by=c('concept_id'='ancestor_concept_id')) %>%
  distinct(descendant_concept_id) %>%
  rename(concept_id=descendant_concept_id) %>%
  dplyr::union(asthma_injection %>% distinct(concept_id)) %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  compute_new()

#' Get all asthma codes with categories (including SABAs and albuterol) and descendants
all_asthma_codes_combined <- final_asthma_codes %>%
  filter(category=='SABA') %>%
  select(concept_id) %>%
  dplyr::union(albuterol %>%
                 select(concept_id)) %>%
  distinct(concept_id) %>%
  left_join(vocabulary_tbl('concept'), by=c('concept_id')) %>%
  mutate(category='SABA') %>%
  dplyr::union(final_asthma_codes %>%
                 filter(!category=='SABA')) %>%
  dplyr::union(asthma_combo_codes_desc %>%
                 mutate(category='combo')) %>%
  dplyr::union(asthma_injection_desc %>%
                 mutate(category='injection')) %>%
  compute_new()

all_asthma_codes_combined %>% as_data_frame() %>% write_csv('specs/all_asthma_codes_combined.csv')






