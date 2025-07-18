# ---- Section 1: Load and Clean Phase 2 Randomized Data ----
# Read raw Phase 2 randomized dataset and replace empty strings with NA for character columns
P2_Randomized <- haven::read_sas("Raw/p2randomized.sas7bdat") %>%
  mutate(across(
    .cols = where(is.character),
    .fns = ~ na_if(.x, "")
  ))

# ---- Section 2: Refine Randomized Data & Identify Treated Participants ----
# Read refined SAS file, select key columns, compute interval and treatment status
P2_Randomized_Refined <- haven::read_sas(
  "/Users/Macintosh/Documents/Work/Johns Hopkins /SKOAP/R_February/Analysis/p2randomized.sas7bdat"
) |>
  select(SubjectID, SiteIDchar, P2RANDOMDATE, P2InterventionDate, P2Direct) |>
  mutate(
    Interval_RandtoInter = time_length(
      P2RANDOMDATE - P2InterventionDate,
      "days"
    ),
    Treated = ifelse(is.na(P2InterventionDate), FALSE, TRUE)
  ) |>
  as_tibble()

# Subset to only treated participants
P2_Randomized_Refined_Treated <- P2_Randomized_Refined |>
  filter(Treated == TRUE)

# ---- Section 3: Distribution of Days from Randomisation to Treatment ----
# Plot histogram with red dashed line at 180 days and save to PDF
ggplot(P2_Randomized_Refined, aes(abs(Interval_RandtoInter))) +
  geom_histogram(binwidth = 10) +
  geom_vline(xintercept = 180, colour = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(
    title = "Distribution of Days from Randomisation to Treatment for SKOAP P2 Participants",
    subtitle = "Red line at 180 days",
    x = "Interval (Days)",
    y = "Count"
  )

# Count how many participants exceeded 180-day interval
P2_Randomized_Refined |>
  count(abs(Interval_RandtoInter) > 180)

# Export JHU site participants with non-missing intervention dates to CSV
P2_Randomized_Refined |>
  filter(str_detect(SubjectID, "JHU"), !is.na(P2InterventionDate)) |>
  select(
    SubjectID,
    P2RANDOMDATE,
    P2InterventionDate
  ) |>
  write_csv("Outputs/JHU.csv")

# ---- Section 4: Build "Interesting" Dataset for Correlation Analysis ----
# Read SAS file, clean, select demographic and outcome scores, compute age and recode sex
P2_Randomized_Interesting <-
  haven::read_sas("Analysis/p2randomized.sas7bdat") %>%
  mutate(across(
    .cols = where(is.character),
    .fns = ~ na_if(.x, "")
  )) |>
  select(
    SubjectID,
    P2RANDOMDATE,
    Sex,
    Ethnicity = Ethnic,
    Race,
    DOB,
    "PROMIS Fatigue Score" = P2PROMISFat4aTScore,
    "PROMIS Sleep Disturbance Score" = P2BLPROMISSlpDist6TScore,
    "PROMIS Physical Function Score" = P2BLPROMISPhysFx6TScore,
    PGIC = P2W12PGIC,
    "Baseline BPI Intensity" = P2BLBPIIntensScore,
    "Week 12 BPI Intensity" = P2W12BPIIntensScore_both,
    "Baseline BPI Interference" = P2BLBPIInterferScore,
    "Week 12 BPI Interference" = P2W12BPIInterferScore,
    "Baseline KOOS Summary Score" = P2BLKOOSSummScore,
    "Week 12KOOS Summary Score" = P2W12KOOSSummScore,
    "WPI Score" = P2WPIPast7DayScore,
    "Symptom Severity Score" = P2SympSevrPast6MoScore,
    "Diffused Pain Phenotype" = P2WPITotalScore,
    "Baseline Pain Catastrophising" = P2BLPCQTotalScore,
    "Week 12 Pain Catastrophising" = P2W12PCQTotalScore,
    "Chronic Pain Coping" = P2BLCPCIsummaryScore,
    "Self-Efficacy for Disease Management" = P2BLASESscore,
    Resiliency = P2CDRISCScore,
    "High Impact Chronic Pain" = P2BLHighImpChronPainYN,
    "Sleep Duration Score" = P2BLSleepDurationCat
  ) |>
  mutate(
    `Age at Randomisation` = round(
      time_length(interval(DOB, P2RANDOMDATE), "year"),
      1
    ),
    Sex = case_when(
      Sex == "Male" ~ 1,
      Sex == "Female" ~ 2,
      Sex == "Unknown" ~ 0
    )
  )

# ---- Section 5: Pearson Correlation Matrix & Heatmap ----
# 1. Compute absolute Pearson correlations
corr_mat <- P2_Randomized_Interesting %>%
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson") |>
  abs()

# 2. Melt matrix into long format
corr_df <- as_tibble(corr_mat) %>%
  tibble::rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "r")

# 3. Define manual ordering of variables
my_order <- c(
  "Baseline BPI Intensity",
  "Week 12 BPI Intensity",
  "Diffused Pain Phenotype",
  "High Impact Chronic Pain",
  "WPI Score",
  "Symptom Severity Score",
  "Baseline KOOS Summary Score",
  "Week 12KOOS Summary Score",
  "Week 12 BPI Interference",
  "Baseline BPI Interference",
  "PROMIS Physical Function Score",
  "PROMIS Fatigue Score",
  "PROMIS Sleep Disturbance Score",
  "Sleep Duration Score",
  "Baseline Pain Catastrophising",
  "Week 12 Pain Catastrophising",
  "Chronic Pain Coping",
  "Self-Efficacy for Disease Management",
  "Resiliency",
  "PGIC",
  "Race",
  "Ethnicity",
  "Age at Randomisation",
  "Sex"
)

# 4. Check for missing/extra variables and relevel
vars_present <- sort(unique(corr_df$Var1))
final_order <- c(
  intersect(my_order, vars_present),
  setdiff(vars_present, my_order)
)
corr_df <- corr_df %>%
  mutate(
    Var1 = fct_relevel(Var1, final_order),
    Var2 = fct_relevel(Var2, final_order)
  )

# 5. Plot and save heatmap
Corr_Plot <- ggplot(corr_df, aes(x = Var1, y = Var2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 1) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(0, 1),
    name = "Pearson\nr"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = "SKOAP Phase II Correlation Plot")

ggsave(
  Corr_Plot,
  filename = "~/Desktop/SKOAP Phase II Correlation Plot.pdf",
  width = 15,
  height = 15,
  units = "in"
)

# ---- Section 6: Summary Statistics of Numeric Variables ----
P2_Randomized_Interesting_Summary <- P2_Randomized_Interesting %>%
  summarise(
    across(
      .cols = where(is.numeric),
      .fns = list(
        mean = ~ mean(., na.rm = TRUE),
        sd = ~ sd(., na.rm = TRUE),
        min = ~ min(., na.rm = TRUE),
        q25 = ~ quantile(., 0.25, na.rm = TRUE),
        median = ~ median(., na.rm = TRUE),
        q75 = ~ quantile(., 0.75, na.rm = TRUE),
        max = ~ max(., na.rm = TRUE),
        IQR = ~ IQR(., na.rm = TRUE),
        q5 = ~ quantile(., 0.05, na.rm = TRUE),
        q95 = ~ quantile(., 0.95, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) |>
  pivot_longer(
    everything(),
    names_to = c("variable", "stat"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  arrange(variable, stat) |>
  mutate(value = round(value, 2)) |>
  pivot_wider(names_from = stat, values_from = value) |>
  relocate(variable, min, q5, q25, median, q75, q95, max, mean, sd, IQR)

# ---- Section 7: Load Additional Raw Datasets ----
# Read raw Phase 2 randomized from Downloads
P2_Randomized_New <- haven::read_sas(
  '/Users/Macintosh/Downloads/p2randomized.sas7bdat'
)

# ---- Section 8: Pain Impact (PI) Data Processing ----
PI <- haven::read_sas("Analysis/painimpact.sas7bdat") |>
  filter(
    SubjectID %in% P2_Randomized_Refined$SubjectID,
    StudyEvent %in%
      c(
        "P2 Online Direct Recruit Assessments",
        "P2 Online Sequential Recruits Assessments",
        "P2 Baseline"
      )
  ) |>
  select(
    SubjectID,
    StudyEvent,
    PainImpactDate,
    BPI_Intensity = BPIintensScore_Orig,
    mBPI = BPIintensScore,
    BPIinterferScore,
    KOOSPainScaleScore,
    KOOSFuncScaleScore,
    KOOSQOLScaleScore,
    KOOSSummScore,
    PCQTotalScore
  ) |>
  left_join(P2_Randomized_Refined, by = "SubjectID") |>
  mutate(
    Interval_PItoRand = time_length(P2RANDOMDATE - PainImpactDate, "days"),
    Interval_PItoInter = time_length(
      P2InterventionDate - PainImpactDate,
      "days"
    )
  ) |>
  select(
    SubjectID,
    P2RANDOMDATE,
    P2InterventionDate,
    Interval_RandtoInter,
    P2Direct,
    Treated,
    StudyEvent,
    PainImpactDate,
    Interval_PItoRand,
    Interval_PItoInter,
    everything()
  ) |>
  mutate(
    StudyEvent = if_else(
      StudyEvent %in%
        c(
          "P2 Online Direct Recruit Assessments",
          "P2 Online Sequential Recruits Assessments"
        ),
      "Randomisation",
      "Baseline"
    )
  ) |>
  slice(-23) |>
  pivot_wider(
    names_from = StudyEvent,
    values_from = c(
      PainImpactDate,
      Interval_PItoRand,
      Interval_PItoInter,
      BPI_Intensity,
      mBPI,
      BPIinterferScore,
      KOOSPainScaleScore,
      KOOSFuncScaleScore,
      KOOSQOLScaleScore,
      KOOSSummScore,
      PCQTotalScore
    ),
    values_fill = NA
  ) |>
  right_join(P2_Randomized_Refined |> select(SubjectID), by = "SubjectID") |>
  mutate(
    across(
      ends_with("_Randomisation"),
      ~ !is.na(.x),
      .names = "{.col}_Present"
    ),
    across(ends_with("_Baseline"), ~ !is.na(.x), .names = "{.col}_Present")
  ) |>
  filter(Treated == TRUE)

# ---- Section 9: Derive "Useful" Flags for PI Measures ----
PI <- PI |>
  mutate(
    across(
      matches("^(BPI_Intensity|mBPI|BPIinterferScore|KOOS.*Score)$"),
      list(
        Randomisation_Useful_For_Randomisation = ~ (!is.na(.x) &
          between(Interval_PItoRand_Randomisation, 0, 30)),
        Baseline_Useful_For_Baseline = ~ (!is.na(.x) &
          between(Interval_PItoInter_Baseline, 0, 30)),
        Randomisation_Useful_For_Baseline = ~ (!is.na(.x) &
          between(Interval_PItoRand_Baseline, 0, 30))
      ),
      .names = "{fn}"
    ),
    PCQTotalScore_Randomisation_Useful_For_Randomisation = if_else(
      !is.na(PCQTotalScore_Randomisation) &
        between(Interval_PItoRand_Randomisation, 0, 90),
      TRUE,
      FALSE
    ),
    PCQTotalScore_Baseline_Useful_For_Baseline = if_else(
      !is.na(PCQTotalScore_Baseline) &
        between(Interval_PItoInter_Baseline, 0, 90),
      TRUE,
      FALSE
    ),
    PCQTotalScore_Randomisation_Useful_For_Baseline = if_else(
      !is.na(PCQTotalScore_Randomisation) &
        between(Interval_PItoRand_Baseline, 0, 90),
      TRUE,
      FALSE
    )
  ) |>
  mutate(
    across(
      ends_with("_Useful_For_Baseline"),
      ~ .x == TRUE,
      .names = "{.col}_Any_Useful"
    )
  ) |>
  filter(Treated == TRUE)

# ---- Section 10: BMI Data: New and Old ----
AAH <- haven::read_sas("Analysis/painimpact.sas7bdat") |>
  filter(SubjectID %in% P2_Randomized_Refined$SubjectID) |>
  filter(
    StudyEvent %in%
      c(
        "P2 Online Direct Recruit Assessments",
        "P2 Online Sequential Recruits Assessments"
      )
  ) |>
  select(SubjectID, BMI_Date = PainImpactDate)

Opioid_P2_Screening <- haven::read_sas("Raw/_p2scr_p2screen.sas7bdat") |>
  filter(SubjectID %in% P2_Randomized_Refined$SubjectID) |>
  select(SubjectID, StudyEvent, `_TDATE`, `_TOTALOPIOIDCALC`) |>
  mutate(
    across(c(`_TDATE`), as_date),
    across(c(`_TOTALOPIOIDCALC`), as.double)
  ) |>
  filter(`_TOTALOPIOIDCALC` != 0) |>
  rename(Total_Daily_Opioid_Dose = `_TOTALOPIOIDCALC`) |>
  select(SubjectID, Total_Daily_Opioid_Dose)

Screening_Survey <- haven::read_sas("Raw/screeningsurvey.sas7bdat")

P2_Eligibility <- haven::read_sas("Raw/_p2eli_p2eligibility.sas7bdat") |>
  filter(SubjectID %in% P2_Randomized_Refined$SubjectID) |>
  select(SubjectID, record_id = `_P2ELIG_REDCAPID`)

# ---- Section 11: Combine Phase 2 Opioid Information ----
P2_Screeing <- P2_Eligibility |>
  left_join(Screening_Survey, by = "record_id") |>
  filter(!is.na(record_id), record_id != "") |>
  group_by(SubjectID) |>
  slice_max(redcap_repeat_instance) |>
  select(SubjectID, totalopioid_dailydose_calc) |>
  mutate(across(c(totalopioid_dailydose_calc), as.numeric)) |>
  rename(Total_Daily_Opioid_Dose = totalopioid_dailydose_calc) |>
  filter(Total_Daily_Opioid_Dose != 0)

Phase_2_Opioid_Information <- bind_rows(
  Opioid_P2_Screening,
  P2_Screeing
) |>
  mutate(Opioid = if_else(!is.na(Total_Daily_Opioid_Dose), TRUE, FALSE))

# ---- Section 12: Opioid Use Logs and Checklists ----
Opioid_Log <- haven::read_sas("Raw/_opioi_opioiduse.sas7bdat") |>
  filter(SubjectID %in% P2_Randomized_Refined$SubjectID) |>
  count(SubjectID) |>
  arrange(desc(n)) |>
  filter(n > 1)

Opi <- Opioid_Log |>
  left_join(Phase_2_Opioid_Information, by = "SubjectID") |>
  select(SubjectID, "Total_Daily_Opioid_Dose", Opioid)

Opioid_Use_Checklist_PreInt <- haven::read_sas(
  "Raw/_aecon_checklist.sas7bdat"
) |>
  filter(
    SubjectID %in% P2_Randomized_Refined$SubjectID,
    StudyEvent == "P2 Pre-Intervention"
  ) |>
  mutate(across(c("_OPIOIDUSEREPORT", "_OPIOIDUSERECORD"), as.integer)) |>
  filter(`_OPIOIDUSEREPORT` != 0)

Opioid_Use_Change_Week12 <- haven::read_sas("Raw/_p2kne_p2painq.sas7bdat") |>
  filter(
    SubjectID %in% P2_Randomized_Refined$SubjectID,
    StudyEvent == "P2 Week 12 Follow-up"
  ) |>
  select(SubjectID, `_P2QOPDUSECHNG`, `_P2QOPDUNRELATED`) |>
  mutate(across(c(`_P2QOPDUSECHNG`, `_P2QOPDUNRELATED`), as.integer)) |>
  filter(`_P2QOPDUSECHNG` %in% c(1, 2))

Opioid_Log <- haven::read_sas("Analysis/opioidlog.sas7bdat") |> count(SubjectID)

# ---- Section 13: Additional SAS Imports (Screening, Eligibility) ----
P2_Screening <- haven::read_sas('Raw/p2screened.sas7bdat')

# ---- Section 14: Insurance and BMI Data Joins ----
New_Data <- P2_Randomized_Refined |>
  left_join(
    haven::read_sas('/Users/Macintosh/Desktop/_insur_ins_bmi.sas7bdat') |>
      filter(SubjectID %in% P2_Randomized_Refined$SubjectID),
    by = "SubjectID"
  )

# New Insurance
New_Insurance_Data <- New_Data |>
  select(
    SubjectID,
    P2RANDOMDATE,
    P2InterventionDate,
    Interval_RandtoInter,
    `_INS_VERIFY_DATE`,
    `_MEDINSURANCEYN`,
    `_PARTICIPIANTPAY`,
    `_INSURANCETYPE`,
    `_MEDICARENAME`,
    `_MEDICAIDNAME`,
    `_PPONAME`,
    `_HMONAME`,
    `_OTHRINSNAME`
  ) |>
  mutate(across(`_INS_VERIFY_DATE`, as_date)) |>
  filter(!is.na(`_INS_VERIFY_DATE`)) |>
  distinct() |>
  group_by(SubjectID) |>
  slice_max(`_INS_VERIFY_DATE`) |>
  ungroup() |>
  mutate(
    Valid_for_Randomisation = time_length(
      interval(P2RANDOMDATE, `_INS_VERIFY_DATE`),
      "day"
    ) <=
      180,
    Valid_for_Treatment = time_length(
      interval(P2InterventionDate, `_INS_VERIFY_DATE`),
      "day"
    ) <=
      180
  ) |>
  select(-c(7:15))

# New BMI
New_BMI_Data <- New_Data |>
  select(
    SubjectID,
    P2RANDOMDATE,
    P2InterventionDate,
    Interval_RandtoInter,
    `_HEIGHT_WEIGHT_DATE`,
    `_HEIGHT_LENGTH`,
    `_WEIGHT`
  ) |>
  mutate(across(`_HEIGHT_WEIGHT_DATE`, as_date)) |>
  filter(!is.na(`_HEIGHT_WEIGHT_DATE`)) |>
  distinct() |>
  group_by(SubjectID) |>
  slice_max(`_HEIGHT_WEIGHT_DATE`) |>
  ungroup() |>
  mutate(
    Valid_for_Randomisation = time_length(
      interval(P2RANDOMDATE, `_HEIGHT_WEIGHT_DATE`),
      "day"
    ) <=
      180,
    Valid_for_Treatment = time_length(
      interval(P2InterventionDate, `_HEIGHT_WEIGHT_DATE`),
      "day"
    ) <=
      180
  )

# Old BMI
Old_BMI_Data <- P2_Randomized_Refined |>
  left_join(P2_Randomized |> select(SubjectID, BMI), by = "SubjectID") |>
  left_join(AAH, by = "SubjectID") |>
  mutate(
    Days_BMI_to_Random = time_length(interval(BMI_Date, P2RANDOMDATE), "day"),
    Viable_for_Randomisation = !is.na(BMI) &
      Days_BMI_to_Random >= 0 &
      Days_BMI_to_Random <= 180
  )

# ---- Section 15: Demographics & Education ----
Demographics <- P2_Randomized_Refined |>
  left_join(
    haven::read_sas('Raw/_demog_demographics.sas7bdat') |>
      filter(SubjectID %in% P2_Randomized_Refined$SubjectID),
    by = "SubjectID"
  ) |>
  slice_max(`_TDATE`, by = SubjectID, with_ties = FALSE)


Demographics_Filtered <- Demographics |>
  filter(
    (P2Direct == 1 & StudyEvent == "P2 Online Direct Recruit Assessments") |
      (P2Direct == 0 & StudyEvent == "P1 Pre-Intervention")
  )

Demographics_Anti <- anti_join(
  Demographics,
  Demographics_Filtered,
  by = "SubjectID"
)

Education_total <- bind_rows(Demographics_Filtered, Demographics_Anti)

# ---- Section 16: Chronic Opioid Use for Direct & Sequential Recruits ----
Chronic_Opioid_P2Seq <- P2_Randomized_Refined |>
  filter(P2Direct == 0) |>
  left_join(
    haven::read_sas('Raw/_bodyp_symptoms_problems.sas7bdat'),
    by = "SubjectID"
  ) |>
  mutate(across(`_OPIOIDUSEYN`, as.integer)) |>
  select(
    SubjectID,
    P2RANDOMDATE,
    P2InterventionDate,
    P2Direct,
    Treated,
    Date = `_TDATE`,
    `_OPIOIDUSEYN`
  )

Chronic_Opioid_P2Dir <- P2_Randomized_Refined |>
  filter(P2Direct == 1) |>
  left_join(
    haven::read_sas('Raw/_chron_crncpain.sas7bdat'),
    by = "SubjectID"
  ) |>
  filter(StudyEvent == "P2 Online Direct Recruit Assessments") |>
  mutate(across(`_OPIOIDUSEYN`, as.integer)) |>
  select(
    SubjectID,
    P2RANDOMDATE,
    P2InterventionDate,
    P2Direct,
    Treated,
    Date = `_TDATEP1BASELINE`,
    `_OPIOIDUSEYN`
  )

# ---- Section 17: Eligibility Joins & Recruit Types ----
P2_Eligibility <- haven::read_sas('Raw/_p2eli_p2eligibility.sas7bdat') |>
  filter(SubjectID %in% P2_Randomized_Refined$SubjectID) |>
  select(SubjectID, `_P2ELIG_REDCAPID`)
P1_Eligibility <- haven::read_sas('Raw/_p1eli_p1eligibility.sas7bdat') |>
  filter(SubjectID %in% P2_Randomized_Refined$SubjectID) |>
  select(SubjectID, `_REDCAPID`)

Direct_Recruits <- P2_Randomized_Refined |>
  filter(P2Direct == 1)
Sequential_Recruits <- P2_Randomized_Refined |>
  filter(P2Direct == 0)

# ---- Section 18: Screening Survey for Insurance (Direct & Sequential) ----
Screeing_Survey_Insurance_Direct <- Screening_Survey |>
  filter(record_id %in% P2_Eligibility$`_P2ELIG_REDCAPID`) |>
  select(record_id, screening_timestamp, 57:70) |>
  inner_join(
    P2_Eligibility |> select(SubjectID, `_P2ELIG_REDCAPID`),
    by = c("record_id" = "_P2ELIG_REDCAPID")
  ) |>
  relocate(SubjectID, .before = record_id) |>
  group_by(SubjectID) |>
  slice_max(screening_timestamp) |>
  filter(SubjectID %in% Direct_Recruits$SubjectID) |>
  right_join(Direct_Recruits, by = "SubjectID")

Screeing_Survey_Insurance_Sequential <- Screening_Survey |>
  filter(record_id %in% P1_Eligibility$`_REDCAPID`) |>
  select(record_id, screening_timestamp, 57:70) |>
  inner_join(
    P1_Eligibility |> select(SubjectID, `_REDCAPID`),
    by = c("record_id" = "_REDCAPID")
  ) |>
  relocate(SubjectID, .before = record_id) |>
  group_by(SubjectID) |>
  slice_max(screening_timestamp) |>
  filter(SubjectID %in% Sequential_Recruits$SubjectID) |>
  right_join(Sequential_Recruits, by = "SubjectID")

# ---- Section 19: Long Phase 2 Data for mBPI ----
Phase_2_Long <- haven::read_sas('Raw/phase2_long.sas7bdat') |>
  filter(
    P2RANDOMDATE == mBPI_Date,
    !is.na(mBPI),
    SubjectID %in% P2_Randomized_Refined$SubjectID
  ) |>
  distinct(SubjectID)

# ---- Section 20: Concomitant Medications ----
AECONMED <- haven::read_sas("Raw/_aecon_checklist.sas7bdat") |>
  filter(
    SubjectID %in% P2_Randomized_Refined$SubjectID,
    StudyEvent %in%
      c(
        "P2 Online Direct Recruit Assessments",
        "P2 Online Sequential Recruits Assessments",
        "P2 Pre-Intervention",
        "P2 Screening (Direct Recruits)",
        "P2 Screening (Sequential Recruits)",
        "P2 Week 12 Follow-up",
        "P2 Week 4 Follow-up",
        "P2 Week 8 Follow-up"
      )
  ) |>
  mutate(across(10:15, as.integer)) |>
  select(-c(2, 5:13))

Concomitant_Medications <- haven::read_sas('Raw/_conco_conmed.sas7bdat') |>
  filter(SubjectID %in% P2_Randomized_Refined$SubjectID) |>
  filter(
    str_detect(`_MEDNAME`, regex("fish", ignore_case = TRUE)) |
      str_detect(`_MEDNAME`, regex("glucosamine", ignore_case = TRUE)) |
      str_detect(`_MEDNAME`, regex("omega", ignore_case = TRUE))
  ) |>
  select(SubjectID, `_MEDNAME`) |>
  mutate(
    Concomitant_Medication = case_when(
      (str_detect(`_MEDNAME`, regex("fish|omega", ignore_case = TRUE)) &
        str_detect(`_MEDNAME`, regex("glucosamine", ignore_case = TRUE))) ~
        "Both",
      str_detect(`_MEDNAME`, regex("fish|omega", ignore_case = TRUE)) ~
        "Fish Oil",
      str_detect(`_MEDNAME`, regex("glucosamine", ignore_case = TRUE)) ~
        "Glucosamine"
    ),
    Concomitant_Medication = as_factor(Concomitant_Medication)
  ) |>
  group_by(SubjectID) |>
  summarise(
    Number_of_Medications = n_distinct(Concomitant_Medication),
    Medication_List = as_factor(paste(
      sort(unique(Concomitant_Medication)),
      collapse = ", "
    ))
  )

P2_Concomitant_Medication <- P2_Randomized_Refined |>
  left_join(Concomitant_Medications, by = "SubjectID") |>
  left_join(
    Randomisation_BPI |> select(SubjectID, "KOAPI" = `_AVERAGE`),
    by = "SubjectID"
  ) |>
  mutate(
    Medication_List = as_factor(Medication_List),
    KOAPI = as.double(KOAPI)
  ) |>
  group_by(Medication_List) |>
  summarise(
    Count = n(),
    KOAPI_Mean = mean(KOAPI, na.rm = TRUE),
    KOAPI_SD = sd(KOAPI, na.rm = TRUE),
    KOAPI_Min = min(KOAPI),
    KOAPI_Median = median(KOAPI, na.rm = TRUE),
    KOAPI_Max = max(KOAPI, na.rm = TRUE)
  )

# ---- Section 21: Screening Survey Joins & Exports ----
Screening_Survey_Direct <- Screening_Survey |>
  filter(record_id %in% P2_Eligibility$`_P2ELIG_REDCAPID`) |>
  inner_join(
    P2_Eligibility |> select(SubjectID, `_P2ELIG_REDCAPID`),
    by = c("record_id" = "_P2ELIG_REDCAPID")
  ) |>
  relocate(SubjectID, .before = record_id) |>
  group_by(SubjectID) |>
  slice_max(screening_timestamp) |>
  filter(SubjectID %in% Direct_Recruits$SubjectID) |>
  right_join(Direct_Recruits, by = "SubjectID")

Screening_Survey_Sequential <- Screening_Survey |>
  filter(record_id %in% P1_Eligibility$`_REDCAPID`) |>
  inner_join(
    P1_Eligibility |> select(SubjectID, `_REDCAPID`),
    by = c("record_id" = "_REDCAPID")
  ) |>
  relocate(SubjectID, .before = record_id) |>
  group_by(SubjectID) |>
  slice_max(screening_timestamp) |>
  filter(SubjectID %in% Sequential_Recruits$SubjectID) |>
  right_join(Sequential_Recruits, by = "SubjectID")

Screening_Survey_P2_Joined <- bind_rows(
  Screening_Survey_Direct,
  Screening_Survey_Sequential
)

# ---- Section 22: Final Status Export ----
P2_Final_Status <- P2_Randomized_Refined |>
  left_join(
    haven::read_sas('Raw/_final_status.sas7bdat') |>
      filter(SubjectID %in% P2_Randomized_Refined$SubjectID),
    by = "SubjectID"
  ) |>
  filter(`_FINALSTATUS` == 2, `_DCSTATUS` == 2) |>
  select(
    SubjectID,
    P2RANDOMDATE,
    P2InterventionDate,
    Interval_RandtoInter,
    P2Direct,
    Treated,
    `_FINALSTATUS`,
    `_FINALSTATUSDATE`,
    `_DCSTATUS`,
    `_DCINSURDENIALDATE`,
    `_DCINSURANCEREASON`,
    `_DCINSURANCEDETAILS`
  ) |>
  write_csv('Outputs/Final_Status.csv')

# ---- Section 23: Zip Codes Export ----
Zipcodes <- read_csv('Raw/SKOAP-Zip-Codes-Phase-1-and-Phase-2.csv') |>
  filter(is.na(zipcode), subjectid %in% P2_Randomized_Refined$SubjectID)

Zipcodes_Treated <- read_csv('Raw/SKOAP-Zip-Codes-Phase-1-and-Phase-2.csv') |>
  filter(is.na(zipcode), subjectid %in% P2_Randomized_Refined_Treated$SubjectID)

# ---- Section 24: Your Knee Pain Assessments ----
Your_Knee_Pain <- haven::read_sas("Raw/_yourk_bpi.sas7bdat") |>
  filter(
    SubjectID %in% P2_Randomized_Refined$SubjectID,
    StudyEvent %in% c("P2 Randomization BPI")
  ) |>
  left_join(P2_Randomized_Refined, by = "SubjectID") |>
  mutate(across(where(is.character), ~ na_if(.x, ""))) |>
  mutate(
    Yourk_Date = as.Date(`_TDATE`),
    Interval_YOURKtoRand = time_length(P2RANDOMDATE - Yourk_Date, "days"),
    Interval_YOURKtoInter = time_length(P2InterventionDate - Yourk_Date, "days")
  ) |>
  select(
    SubjectID,
    P2RANDOMDATE,
    P2InterventionDate,
    Interval_RandtoInter,
    P2Direct,
    Treated,
    StudyEvent,
    Yourk_Date,
    Interval_YOURKtoRand,
    Interval_YOURKtoInter,
    `_AVERAGE`
  )

Your_Knee_Pain_Online <- haven::read_sas("Raw/_yourk_bpi.sas7bdat") |>
  filter(
    SubjectID %in% P2_Randomized_Refined$SubjectID,
    StudyEvent %in% c("P2 BPI Online Assessment (ALL)")
  ) |>
  left_join(P2_Randomized_Refined, by = "SubjectID") |>
  mutate(across(where(is.character), ~ na_if(.x, ""))) |>
  mutate(
    Yourk_Date = as.Date(`_TDATE`),
    Interval_YOURKtoRand = time_length(P2RANDOMDATE - Yourk_Date, "days"),
    Interval_YOURKtoInter = time_length(P2InterventionDate - Yourk_Date, "days")
  ) |>
  select(
    SubjectID,
    P2RANDOMDATE,
    P2InterventionDate,
    Interval_RandtoInter,
    P2Direct,
    Treated,
    StudyEvent,
    Yourk_Date,
    Interval_YOURKtoRand,
    Interval_YOURKtoInter,
    `_AVERAGE`
  ) |>
  filter(Interval_YOURKtoInter > 0 | is.na(Interval_YOURKtoInter)) |>
  slice_min(Interval_YOURKtoInter, with_ties = FALSE, by = SubjectID) |>
  right_join(P2_Randomized_Refined |> select(SubjectID), by = "SubjectID") |>
  filter(
    (Interval_YOURKtoInter > 0 & Interval_YOURKtoInter < 30) |
      is.na(Interval_YOURKtoInter)
  ) |>
  filter(!is.na(`_AVERAGE`))

# ---- Section 25: Missing mBPI Data ----
No_mBPI_data <- haven::read_sas("Analysis/painimpact.sas7bdat") |>
  filter(
    SubjectID %in% P2_Randomized_Refined$SubjectID,
    StudyEvent %in%
      c(
        "P2 Online Direct Recruit Assessments",
        "P2 Online Sequential Recruits Assessments",
        "P2 Baseline"
      )
  ) |>
  select(SubjectID, StudyEvent, PainImpactDate, mBPI = BPIintensScore) |>
  left_join(P2_Randomized_Refined, by = "SubjectID") |>
  mutate(
    Interval_PItoRand = time_length(P2RANDOMDATE - PainImpactDate, "days"),
    Interval_PItoInter = time_length(
      P2InterventionDate - PainImpactDate,
      "days"
    )
  ) |>
  select(
    SubjectID,
    P2RANDOMDATE,
    P2InterventionDate,
    Interval_RandtoInter,
    P2Direct,
    Treated,
    StudyEvent,
    PainImpactDate,
    Interval_PItoRand,
    Interval_PItoInter,
    mBPI
  )
