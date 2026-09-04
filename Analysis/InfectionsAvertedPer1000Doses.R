rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 999)

cat("================================================================================\n")
cat("CORRECTED LA PrEP EFFICIENCY ANALYSIS (Infections Averted per 1000 Doses)\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================================\n\n")

# ==================== 1. LOAD PACKAGES ====================
required_packages <- c("dplyr", "tidyr", "purrr", "data.table", "ggplot2", "stringr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ==================== 2. CONFIGURATION ====================
BASE_PATH <- "/gpfs/scratch/mudime01/SouthAfricaPrEPTiming/250simulations_delayedPrEP_20Jan2026"
BASELINE_SCENARIO <- "MainAnalysisBaseline"
OUTPUT_DIR <- file.path(BASE_PATH, "LAPrEP_Efficiency_Corrected")

TEST_MODE <- FALSE
MAX_FILES_TEST <- 5

ANALYSIS_START_YEAR <- 2026
ANALYSIS_END_YEAR <- 2060

CENSUS_YEAR <- 2024.5
SA_CENSUS_POP <- 63015904

REPORT_SUBFOLDER <- "ReportHIVByAgeAndGender"

# Column names
DOSE_COLUMN <- "ReceivedLAPrEP"
INFECTION_COLUMN <- "Newly.Infected"  # Make sure this matches your data

# ==================== 3. CORE FUNCTIONS ====================

# Function to read and process a single file (identical to separate code)
read_emod_file_corrected <- function(file_path) {
  tryCatch({
    df <- data.table::fread(file_path, showProgress = FALSE, check.names = TRUE)
    data.table::setnames(df, names(df), make.names(names(df)))
    
    # Check required columns
    required_cols <- c("Year", "Gender", "Age", INFECTION_COLUMN, DOSE_COLUMN, "Population")
    missing_cols <- setdiff(required_cols, names(df))
    
    if (length(missing_cols) > 0) {
      return(NULL)
    }
    
    # Extract simulation ID
    filename <- basename(file_path)
    sim_id_match <- regmatches(filename, regexec("REP([0-9]+)", filename))
    if (length(sim_id_match[[1]]) > 1) {
      sim_id <- as.integer(sim_id_match[[1]][2])
    } else {
      sim_id_match2 <- regmatches(filename, regexec("TPI([0-9]+)", filename))
      if (length(sim_id_match2[[1]]) > 1) {
        sim_id <- as.integer(sim_id_match2[[1]][2])
      } else {
        sim_id <- as.integer(gsub("\\D", "", filename))
        if (is.na(sim_id)) sim_id <- 0
      }
    }
    
    # Population scaling (same as separate code)
    df_census <- df %>% dplyr::filter(Year == CENSUS_YEAR)
    
    if (nrow(df_census) == 0) {
      pop.scaling.factor <- 1
    } else {
      total.pop <- sum(df_census$Population, na.rm = TRUE)
      pop.scaling.factor <- ifelse(is.na(total.pop) || total.pop <= 0, 1, SA_CENSUS_POP / total.pop)
    }
    
    # Filter to analysis period
    df <- df %>%
      dplyr::filter(Year >= ANALYSIS_START_YEAR & Year <= ANALYSIS_END_YEAR)
    
    if (nrow(df) == 0) return(NULL)
    
    # Apply scaling
    df <- df %>%
      dplyr::mutate(
        Newly.Infected = Newly.Infected * pop.scaling.factor,
        ReceivedLAPrEP = ReceivedLAPrEP * pop.scaling.factor,
        Population = Population * pop.scaling.factor,
        Age_Group = ifelse(Age >= 15 & Age <= 49, "15-49", "0-99"),
        Gender_Label = ifelse(Gender == 0, "Male", "Female"),
        Simulation = sim_id,
        Pop_Scaling_Factor = pop.scaling.factor
      )
    
    # Summarize (matching separate code logic)
    result <- df %>%
      dplyr::group_by(Gender_Label, Age_Group, Simulation) %>%
      dplyr::summarise(
        Cumulative_Infections = sum(Newly.Infected, na.rm = TRUE),
        Cumulative_Doses = sum(ReceivedLAPrEP, na.rm = TRUE),
        Avg_Population = mean(Population, na.rm = TRUE),
        Pop_Scaling_Factor = dplyr::first(Pop_Scaling_Factor),
        .groups = "drop"
      )
    
    return(result)
    
  }, error = function(e) {
    return(NULL)
  })
}

# Process all files in a scenario
process_scenario_corrected <- function(scenario_name) {
  scenario_path <- file.path(BASE_PATH, scenario_name, REPORT_SUBFOLDER)
  
  if (!dir.exists(scenario_path)) return(NULL)
  
  csv_files <- list.files(scenario_path, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) return(NULL)
  
  csv_files <- sort(csv_files)
  if (TEST_MODE && length(csv_files) > MAX_FILES_TEST) {
    csv_files <- csv_files[1:MAX_FILES_TEST]
  }
  
  all_data <- list()
  
  for (i in seq_along(csv_files)) {
    result <- read_emod_file_corrected(csv_files[i])
    if (!is.null(result)) {
      all_data[[i]] <- result
    }
  }
  
  if (length(all_data) == 0) return(NULL)
  
  combined_data <- data.table::rbindlist(all_data, fill = TRUE)
  combined_data$Scenario <- scenario_name
  
  return(combined_data)
}

# Calculate efficiency (infections averted per 1000 doses)
calculate_efficiency_corrected <- function(baseline_data, intervention_data) {
  # Find common simulations
  common_sims <- intersect(unique(baseline_data$Simulation), 
                           unique(intervention_data$Simulation))
  
  if (length(common_sims) == 0) return(NULL)
  
  # Filter to common simulations
  baseline_common <- baseline_data %>% dplyr::filter(Simulation %in% common_sims)
  intervention_common <- intervention_data %>% dplyr::filter(Simulation %in% common_sims)
  
  # Calculate averted infections and doses delivered
  results <- baseline_common %>%
    dplyr::select(Simulation, Gender_Label, Age_Group, Cumulative_Infections) %>%
    dplyr::rename(Baseline_Infections = Cumulative_Infections) %>%
    dplyr::left_join(
      intervention_common %>%
        dplyr::select(Simulation, Gender_Label, Age_Group, 
                      Cumulative_Infections, Cumulative_Doses),
      by = c("Simulation", "Gender_Label", "Age_Group")
    ) %>%
    dplyr::rename(Intervention_Infections = Cumulative_Infections) %>%
    dplyr::mutate(
      Infections_Averted = Baseline_Infections - Intervention_Infections,
      # Efficiency: Infections averted per 1000 doses
      Efficiency = ifelse(Cumulative_Doses > 0, 
                          (Infections_Averted / Cumulative_Doses) * 1000, 
                          NA_real_)
    )
  
  # Summarize by group
  summary <- results %>%
    dplyr::group_by(Gender_Label, Age_Group) %>%
    dplyr::summarise(
      N_Simulations = dplyr::n(),
      Mean_Infections_Averted = mean(Infections_Averted, na.rm = TRUE),
      SD_Infections_Averted = sd(Infections_Averted, na.rm = TRUE),
      Mean_Doses_Delivered = mean(Cumulative_Doses, na.rm = TRUE),
      SD_Doses_Delivered = sd(Cumulative_Doses, na.rm = TRUE),
      Mean_Efficiency = mean(Efficiency, na.rm = TRUE),
      SD_Efficiency = sd(Efficiency, na.rm = TRUE),
      CI_Lower = Mean_Efficiency - qt(0.975, N_Simulations - 1) * (SD_Efficiency / sqrt(N_Simulations)),
      CI_Upper = Mean_Efficiency + qt(0.975, N_Simulations - 1) * (SD_Efficiency / sqrt(N_Simulations)),
      .groups = "drop"
    )
  
  # Add overall population-level summary
  overall <- results %>%
    dplyr::group_by(Simulation) %>%
    dplyr::summarise(
      Total_Infections_Averted = sum(Infections_Averted, na.rm = TRUE),
      Total_Doses = sum(Cumulative_Doses, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::summarise(
      Gender_Label = "Overall",
      Age_Group = "All",
      N_Simulations = dplyr::n(),
      Mean_Infections_Averted = mean(Total_Infections_Averted, na.rm = TRUE),
      SD_Infections_Averted = sd(Total_Infections_Averted, na.rm = TRUE),
      Mean_Doses_Delivered = mean(Total_Doses, na.rm = TRUE),
      SD_Doses_Delivered = sd(Total_Doses, na.rm = TRUE),
      Mean_Efficiency = mean((Total_Infections_Averted / Total_Doses) * 1000, na.rm = TRUE),
      SD_Efficiency = sd((Total_Infections_Averted / Total_Doses) * 1000, na.rm = TRUE),
      CI_Lower = Mean_Efficiency - qt(0.975, N_Simulations - 1) * (SD_Efficiency / sqrt(N_Simulations)),
      CI_Upper = Mean_Efficiency + qt(0.975, N_Simulations - 1) * (SD_Efficiency / sqrt(N_Simulations))
    )
  
  final_summary <- dplyr::bind_rows(summary, overall)
  
  return(list(summary = final_summary, detailed = results))
}

# ==================== 4. MAIN EXECUTION ====================
cat("\nProcessing baseline scenario...\n")
baseline_data <- process_scenario_corrected(BASELINE_SCENARIO)

if (is.null(baseline_data)) {
  stop("ERROR: Failed to process baseline scenario")
}

cat("\nBaseline summary (should match separate code):\n")
baseline_check <- baseline_data %>%
  dplyr::group_by(Gender_Label, Age_Group) %>%
  dplyr::summarise(
    Mean_Infections = mean(Cumulative_Infections, na.rm = TRUE),
    .groups = "drop"
  )
print(baseline_check)

# Get intervention scenarios
all_scenarios <- list.dirs(BASE_PATH, full.names = FALSE, recursive = FALSE)
intervention_scenarios <- sort(all_scenarios[grepl("^MainAnalysis.*VOICE[246]", all_scenarios)])

cat("\nFound", length(intervention_scenarios), "intervention scenarios\n")

all_efficiency_results <- list()

for (scenario in intervention_scenarios) {
  cat("\nProcessing:", scenario, "\n")
  
  intervention_data <- process_scenario_corrected(scenario)
  
  if (!is.null(intervention_data)) {
    # Check intervention doses (should match separate doses code)
    dose_check <- intervention_data %>%
      dplyr::filter(Cumulative_Doses > 0) %>%
      dplyr::group_by(Gender_Label, Age_Group) %>%
      dplyr::summarise(
        Mean_Doses = mean(Cumulative_Doses, na.rm = TRUE),
        .groups = "drop"
      )
    
    cat("  Mean doses delivered:\n")
    print(dose_check)
    
    # Calculate efficiency
    efficiency <- calculate_efficiency_corrected(baseline_data, intervention_data)
    
    if (!is.null(efficiency)) {
      efficiency$summary$Scenario <- scenario
      all_efficiency_results[[scenario]] <- efficiency
      cat("  ✓ Efficiency calculated\n")
      
      # Print results
      print(efficiency$summary %>% 
              dplyr::filter(Gender_Label == "Overall") %>%
              dplyr::select(Scenario, Mean_Efficiency, CI_Lower, CI_Upper))
    }
  }
}

# Save results
if (length(all_efficiency_results) > 0) {
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  
  # Combine all summaries
  all_summaries <- do.call(rbind, lapply(all_efficiency_results, function(x) x$summary))
  
  # Extract VOICE level
  all_summaries <- all_summaries %>%
    dplyr::mutate(
      VOICE_Level = case_when(
        grepl("VOICE2", Scenario) ~ "VOICE2",
        grepl("VOICE4", Scenario) ~ "VOICE4",
        grepl("VOICE6", Scenario) ~ "VOICE6",
        TRUE ~ NA_character_
      )
    )
  
  # Save main results
  write.csv(all_summaries, 
            file.path(OUTPUT_DIR, "LAPrEP_Efficiency_Corrected.csv"), 
            row.names = FALSE)
  
  # Save detailed results for each scenario
  for (scenario in names(all_efficiency_results)) {
    safe_name <- gsub("[^A-Za-z0-9]", "_", scenario)
    write.csv(all_efficiency_results[[scenario]]$detailed,
              file.path(OUTPUT_DIR, paste0("Detailed_", safe_name, ".csv")),
              row.names = FALSE)
  }
  
  cat("\n================================================================================\n")
  cat("RESULTS SUMMARY (Infections Averted per 1000 Doses)\n")
  cat("================================================================================\n")
  
  # Display overall results
  all_summaries %>%
    dplyr::filter(Gender_Label == "Overall") %>%
    dplyr::select(VOICE_Level, Mean_Efficiency, CI_Lower, CI_Upper, 
                  Mean_Infections_Averted, Mean_Doses_Delivered) %>%
    dplyr::arrange(VOICE_Level) %>%
    print()
  
  cat("\nResults saved to:", OUTPUT_DIR, "\n")
  cat("================================================================================\n")
}


####nothing is estimated correctly however the otput is good 
# rm(list = ls())
# options(stringsAsFactors = FALSE)
# options(scipen = 999)
# 
# cat("================================================================================\n")
# cat("STARTING EMOD LA PrEP DOSES DISTRIBUTED ANALYSIS (2026-2060)\n")
# cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
# cat("================================================================================\n\n")
# 
# # ==================== 1. LOAD PACKAGES ====================
# required_packages <- c("dplyr", "tidyr", "purrr", "data.table", "ggplot2")
# 
# for (pkg in required_packages) {
#   if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
#     install.packages(pkg, dependencies = TRUE, quiet = TRUE)
#     library(pkg, character.only = TRUE)
#   }
# }
# 
# # ==================== 2. CONFIGURATION ====================
# BASE_PATH <- "/gpfs/scratch/mudime01/SouthAfricaPrEPTiming/250simulations_delayedPrEP_20Jan2026"
# BASELINE_SCENARIO <- "MainAnalysisBaseline"
# OUTPUT_DIR <- "/gpfs/scratch/mudime01/SouthAfricaPrEPTiming/250simulations_delayedPrEP_20Jan2026/InfectionsAvertedPer1000Doses"
# 
# TEST_MODE <- FALSE
# MAX_FILES_TEST <- 5
# 
# ANALYSIS_START_YEAR <- 2026
# ANALYSIS_END_YEAR <- 2060
# ANALYSIS_PERIOD_YEARS <- ANALYSIS_END_YEAR - ANALYSIS_START_YEAR + 1
# 
# CENSUS_YEAR <- 2024.5
# SA_CENSUS_POP <- 63015904
# 
# REPORT_SUBFOLDER <- "ReportHIVByAgeAndGender"
# 
# DOSE_COLUMN_NAME <- "ReceivedLAPrEP"
# INFECTION_COLUMN_NAME <- "Newly.Infected"
# 
# # ==================== 3. FUNCTIONS ====================
# 
# read_emod_file_laprep_doses <- function(file_path) {
#   cat(paste("    Reading file:", basename(file_path), "\n"))
#   
#   tryCatch({
#     df <- data.table::fread(file_path, showProgress = FALSE, check.names = TRUE)
#     data.table::setnames(df, names(df), make.names(names(df)))
#     
#     dose_col <- make.names(DOSE_COLUMN_NAME)
#     infection_col <- make.names(INFECTION_COLUMN_NAME)
#     
#     required_cols <- c("Year", "Gender", "Age", dose_col, infection_col, "Population")
#     missing_cols <- setdiff(required_cols, names(df))
#     
#     if (length(missing_cols) > 0) {
#       warning(paste("Missing columns in", basename(file_path), ":",
#                     paste(missing_cols, collapse = ", ")))
#       cat("Available columns are:\n")
#       print(names(df))
#       return(NULL)
#     }
#     
#     filename <- basename(file_path)
#     sim_id_match <- regmatches(filename, regexec("REP([0-9]+)", filename))
#     
#     if (length(sim_id_match[[1]]) > 1) {
#       sim_id <- as.integer(sim_id_match[[1]][2])
#     } else {
#       sim_id_match2 <- regmatches(filename, regexec("TPI([0-9]+)", filename))
#       if (length(sim_id_match2[[1]]) > 1) {
#         sim_id <- as.integer(sim_id_match2[[1]][2])
#       } else {
#         sim_id <- as.integer(gsub("\\D", "", filename))
#         if (is.na(sim_id)) sim_id <- 0
#       }
#     }
#     
#     df_census <- df %>% dplyr::filter(Year == CENSUS_YEAR)
#     
#     if (nrow(df_census) == 0) {
#       pop.scaling.factor <- 1
#     } else {
#       total.pop <- sum(df_census$Population, na.rm = TRUE)
#       pop.scaling.factor <- ifelse(
#         is.na(total.pop) || total.pop <= 0,
#         1,
#         SA_CENSUS_POP / total.pop
#       )
#     }
#     
#     df <- df %>%
#       dplyr::filter(Year >= ANALYSIS_START_YEAR & Year <= ANALYSIS_END_YEAR)
#     
#     if (nrow(df) == 0) {
#       cat(paste("    WARNING: No data in analysis period\n"))
#       return(NULL)
#     }
#     
#     df[[dose_col]] <- df[[dose_col]] * pop.scaling.factor
#     df[[infection_col]] <- df[[infection_col]] * pop.scaling.factor
#     df$Population <- df$Population * pop.scaling.factor
#     
#     df <- df %>%
#       dplyr::mutate(
#         Gender_Label = ifelse(Gender == 0, "Male", "Female"),
#         Simulation = sim_id,
#         Pop_Scaling_Factor = pop.scaling.factor
#       )
#     
#     # Dataset 1: all ages
#     df_all <- df %>%
#       dplyr::mutate(Age_Group = "All")
#     
#     # Dataset 2: age 15-49 only
#     df_1549 <- df %>%
#       dplyr::filter(Age >= 15 & Age <= 49) %>%
#       dplyr::mutate(Age_Group = "15-49")
#     
#     df_age_groups <- dplyr::bind_rows(df_all, df_1549)
#     
#     result <- df_age_groups %>%
#       dplyr::group_by(Gender_Label, Age_Group, Simulation) %>%
#       dplyr::summarise(
#         Cumulative_Doses = sum(.data[[dose_col]], na.rm = TRUE),
#         Cumulative_Infections = sum(.data[[infection_col]], na.rm = TRUE),
#         Mean_Annual_Doses = mean(.data[[dose_col]], na.rm = TRUE),
#         Pop_Scaling_Factor = dplyr::first(Pop_Scaling_Factor),
#         .groups = "drop"
#       )
#     
#     cat(paste("    ✓ Processed", basename(file_path),
#               "- Sim:", sim_id,
#               "- Cum doses:", round(sum(result$Cumulative_Doses)), "\n"))
#     
#     return(result)
#     
#   }, error = function(e) {
#     cat(paste("    ✗ ERROR reading", basename(file_path), ":", e$message, "\n"))
#     return(NULL)
#   })
# }
# 
# 
# process_scenario_files_laprep_doses <- function(scenario_name, is_baseline = FALSE) {
#   scenario_path <- file.path(BASE_PATH, scenario_name, REPORT_SUBFOLDER)
#   
#   if (!dir.exists(scenario_path)) {
#     cat(paste("    ✗ ERROR: Directory not found:", scenario_path, "\n"))
#     return(NULL)
#   }
#   
#   csv_files <- list.files(scenario_path, pattern = "\\.csv$", full.names = TRUE)
#   if (length(csv_files) == 0) return(NULL)
#   
#   csv_files <- sort(csv_files)
#   if (TEST_MODE && length(csv_files) > MAX_FILES_TEST) {
#     csv_files <- csv_files[1:MAX_FILES_TEST]
#   }
#   
#   all_data <- list()
#   success_count <- 0
#   
#   for (i in seq_along(csv_files)) {
#     result <- read_emod_file_laprep_doses(csv_files[i])
#     if (!is.null(result)) {
#       all_data[[i]] <- result
#       success_count <- success_count + 1
#     }
#   }
#   
#   if (success_count == 0) return(NULL)
#   
#   combined_data <- data.table::rbindlist(all_data, fill = TRUE)
#   combined_data$Scenario <- scenario_name
#   combined_data$Scenario_Type <- ifelse(is_baseline, "Baseline", "Intervention")
#   return(combined_data)
# }
# 
# 
# calculate_laprep_efficiency <- function(baseline_data, intervention_data, intervention_scenario_name) {
#   
#   cat(paste("\n  Calculating infections averted per 1,000 doses for:", 
#             intervention_scenario_name, "\n"))
#   
#   if (is.null(baseline_data) || is.null(intervention_data)) return(NULL)
#   
#   baseline_sims <- unique(baseline_data$Simulation)
#   intervention_sims <- unique(intervention_data$Simulation)
#   common_sims <- intersect(baseline_sims, intervention_sims)
#   
#   cat(paste("    Baseline simulations:", length(baseline_sims), "\n"))
#   cat(paste("    Intervention simulations:", length(intervention_sims), "\n"))
#   cat(paste("    Common simulations:", length(common_sims), "\n"))
#   
#   if (length(common_sims) == 0) {
#     stop("No common simulation IDs between baseline and intervention.")
#   }
#   
#   baseline_matched <- baseline_data %>%
#     dplyr::filter(Simulation %in% common_sims) %>%
#     dplyr::mutate(Match_ID = Simulation)
#   
#   intervention_matched <- intervention_data %>%
#     dplyr::filter(Simulation %in% common_sims) %>%
#     dplyr::mutate(Match_ID = Simulation)
#   
#   # Baseline summary per simulation, gender, age group
#   baseline_summary <- baseline_matched %>%
#     dplyr::group_by(Match_ID, Gender_Label, Age_Group) %>%
#     dplyr::summarise(
#       Baseline_Infections = mean(Cumulative_Infections, na.rm = TRUE),
#       Baseline_Doses = mean(Cumulative_Doses, na.rm = TRUE),
#       .groups = "drop"
#     )
#   
#   # Intervention summary per simulation, gender, age group
#   intervention_summary <- intervention_matched %>%
#     dplyr::group_by(Match_ID, Gender_Label, Age_Group) %>%
#     dplyr::summarise(
#       Intervention_Infections = mean(Cumulative_Infections, na.rm = TRUE),
#       Intervention_Doses = mean(Cumulative_Doses, na.rm = TRUE),
#       .groups = "drop"
#     )
#   
#   combined_summary <- baseline_summary %>%
#     dplyr::inner_join(
#       intervention_summary,
#       by = c("Match_ID", "Gender_Label", "Age_Group")
#     ) %>%
#     dplyr::mutate(
#       Infections_Averted = Baseline_Infections - Intervention_Infections,
#       Doses_Delivered = Intervention_Doses - Baseline_Doses
#     )
#   
#   # =========================================================
#   # POPULATION-LEVEL EFFICIENCY:
#   # total infections averted across all ages/genders /
#   # total LEN doses delivered across all ages/genders
#   # =========================================================
#   
#   population_detailed <- combined_summary %>%
#     dplyr::group_by(Match_ID) %>%
#     dplyr::summarise(
#       Total_Baseline_Infections = sum(Baseline_Infections, na.rm = TRUE),
#       Total_Intervention_Infections = sum(Intervention_Infections, na.rm = TRUE),
#       Total_Infections_Averted = sum(Infections_Averted, na.rm = TRUE),
#       Total_Doses_Delivered = sum(Doses_Delivered, na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     dplyr::mutate(
#       Infections_Averted_per_1000_Doses =
#         dplyr::if_else(
#           Total_Doses_Delivered > 0,
#           (Total_Infections_Averted / Total_Doses_Delivered) * 1000,
#           NA_real_
#         )
#     )
#   
#   population_summary <- population_detailed %>%
#     dplyr::summarise(
#       N = dplyr::n(),
#       
#       Mean_Baseline_Infections =
#         mean(Total_Baseline_Infections, na.rm = TRUE),
#       
#       Mean_Intervention_Infections =
#         mean(Total_Intervention_Infections, na.rm = TRUE),
#       
#       Mean_Infections_Averted =
#         mean(Total_Infections_Averted, na.rm = TRUE),
#       
#       SD_Infections_Averted =
#         sd(Total_Infections_Averted, na.rm = TRUE),
#       
#       Infections_Averted_CI_Lower =
#         Mean_Infections_Averted -
#         qt(0.975, N - 1) *
#         (SD_Infections_Averted / sqrt(N)),
#       
#       Infections_Averted_CI_Upper =
#         Mean_Infections_Averted +
#         qt(0.975, N - 1) *
#         (SD_Infections_Averted / sqrt(N)),
#       
#       Mean_Doses_Delivered =
#         mean(Total_Doses_Delivered, na.rm = TRUE),
#       
#       Mean_Infections_Averted_per_1000_Doses =
#         mean(Infections_Averted_per_1000_Doses, na.rm = TRUE),
#       
#       SD_Efficiency =
#         sd(Infections_Averted_per_1000_Doses, na.rm = TRUE),
#       
#       Efficiency_CI_Lower =
#         Mean_Infections_Averted_per_1000_Doses -
#         qt(0.975, N - 1) *
#         (SD_Efficiency / sqrt(N)),
#       
#       Efficiency_CI_Upper =
#         Mean_Infections_Averted_per_1000_Doses +
#         qt(0.975, N - 1) *
#         (SD_Efficiency / sqrt(N)),
#       
#       Ratio_of_Means_per_1000 =
#         (Mean_Infections_Averted / Mean_Doses_Delivered) * 1000,
#       
#       .groups = "drop"
#     ) %>%
#     dplyr::mutate(
#       Scenario = intervention_scenario_name,
#       Analysis_Period = paste(ANALYSIS_START_YEAR, "-", ANALYSIS_END_YEAR),
#       Analysis_Years = ANALYSIS_PERIOD_YEARS,
#       Result_Type = "Population-level",
#       Gender_Label = "Overall",
#       Age_Group = "All"
#     ) %>%
#     dplyr::select(
#       Scenario,
#       Analysis_Period,
#       Analysis_Years,
#       Result_Type,
#       Gender_Label,
#       Age_Group,
#       dplyr::everything()
#     )
#   
#   return(list(
#     summary = population_summary,
#     detailed = combined_summary,
#     population_detailed = population_detailed
#   ))
# }
# 
# 
# save_laprep_dose_results <- function(results_list, output_dir) {
#   if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
#   
#   all_summaries <- do.call(rbind, lapply(results_list, function(x) x$summary))
#   
#   write.csv(
#     all_summaries,
#     file.path(output_dir, "Infections_Averted_per_1000_Doses_summary.csv"),
#     row.names = FALSE
#   )
#   
#   for (i in seq_along(results_list)) {
#     scenario_name <- names(results_list)[i]
#     safe_name <- gsub("[^A-Za-z0-9]", "_", scenario_name)
#     
#     write.csv(
#       results_list[[i]]$detailed,
#       file.path(output_dir, paste0("detailed_", safe_name, ".csv")),
#       row.names = FALSE
#     )
#     
#     write.csv(
#       results_list[[i]]$population_detailed,
#       file.path(output_dir, paste0("population_detailed_", safe_name, ".csv")),
#       row.names = FALSE
#     )
#   }
#   
#   return(all_summaries)
# }
# 
# # ==================== 4. MAIN EXECUTION ====================
# baseline_data_laprep_doses <- process_scenario_files_laprep_doses(
#   BASELINE_SCENARIO,
#   is_baseline = TRUE
# )
# 
# if (is.null(baseline_data_laprep_doses)) {
#   stop("ERROR: Failed to process baseline scenario.")
# }
# 
# all_scenarios <- list.dirs(BASE_PATH, full.names = FALSE, recursive = FALSE)
# 
# intervention_scenarios <- sort(
#   all_scenarios[grepl("^MainAnalysis.*VOICE[2468]", all_scenarios)]
# )
# 
# # For sensitivity analysis, use:
# # intervention_scenarios <- sort(all_scenarios[grepl("^Sens\\d{4}", all_scenarios)])
# 
# all_results_laprep_doses <- list()
# 
# for (intervention_scenario in intervention_scenarios) {
#   intervention_data <- process_scenario_files_laprep_doses(
#     intervention_scenario,
#     is_baseline = FALSE
#   )
#   
#   if (!is.null(intervention_data)) {
#     result <- calculate_laprep_efficiency(
#       baseline_data_laprep_doses,
#       intervention_data,
#       intervention_scenario
#     )
#     
#     if (!is.null(result)) {
#       all_results_laprep_doses[[intervention_scenario]] <- result
#     }
#   }
# }
# 
# if (length(all_results_laprep_doses) > 0) {
#   final_summary_laprep_doses <- save_laprep_dose_results(
#     all_results_laprep_doses,
#     OUTPUT_DIR
#   )
#   
#   cat("\n✓ LA PrEP infections averted per 1,000 doses analysis complete\n")
#   cat("Results saved to:", OUTPUT_DIR, "\n")
#   
# } else {
#   cat("\n✗ No results generated\n")
# }
# 




#####Estimates correct number of doses but not infections averted.

# rm(list = ls())
# options(stringsAsFactors = FALSE)
# options(scipen = 999)
# 
# cat("================================================================================\n")
# cat("STARTING EMOD LA PrEP DOSES DISTRIBUTED ANALYSIS (2026-2060)\n")
# cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
# cat("================================================================================\n\n")
# 
# # ==================== 1. LOAD PACKAGES ====================
# required_packages <- c("dplyr", "tidyr", "purrr", "data.table", "ggplot2")
# 
# for (pkg in required_packages) {
#   if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
#     install.packages(pkg, dependencies = TRUE, quiet = TRUE)
#     library(pkg, character.only = TRUE)
#   }
# }
# 
# # ==================== 2. CONFIGURATION ====================
# BASE_PATH <- "/gpfs/scratch/mudime01/SouthAfricaPrEPTiming/250simulations_delayedPrEP_20Jan2026"
# BASELINE_SCENARIO <- "MainAnalysisBaseline"
# OUTPUT_DIR <- "/gpfs/scratch/mudime01/SouthAfricaPrEPTiming/250simulations_delayedPrEP_20Jan2026/InfectionsAvertedPer1000Doses"
# 
# TEST_MODE <- FALSE
# MAX_FILES_TEST <- 5
# 
# ANALYSIS_START_YEAR <- 2026
# ANALYSIS_END_YEAR <- 2060
# ANALYSIS_PERIOD_YEARS <- ANALYSIS_END_YEAR - ANALYSIS_START_YEAR + 1
# 
# CENSUS_YEAR <- 2024.5
# SA_CENSUS_POP <- 63015904
# 
# # Change this if the event count lives in another report folder
# REPORT_SUBFOLDER <- "ReportHIVByAgeAndGender"
# 
# DOSE_COLUMN_NAME <- "ReceivedLAPrEP"
# 
# INFECTION_COLUMN_NAME <- "Newly.Infected"  # change this to the actual column name in your report
# 
# # ==================== 3. FUNCTIONS ====================
# 
# read_emod_file_laprep_doses <- function(file_path) {
#   cat(paste("    Reading file:", basename(file_path), "\n"))
# 
#   tryCatch({
#     df <- data.table::fread(file_path, showProgress = FALSE, check.names = TRUE)
#     data.table::setnames(df, names(df), make.names(names(df)))
# 
#     dose_col <- make.names(DOSE_COLUMN_NAME)
#     infection_col <- make.names(INFECTION_COLUMN_NAME)
# 
#     required_cols <- c("Year", "Gender", "Age", dose_col, infection_col, "Population")
#     missing_cols <- setdiff(required_cols, names(df))
# 
#     if (length(missing_cols) > 0) {
#       warning(paste("Missing columns in", basename(file_path), ":",
#                     paste(missing_cols, collapse = ", ")))
#       cat("Available columns are:\n")
#       print(names(df))
#       return(NULL)
#     }
# 
#     filename <- basename(file_path)
#     sim_id_match <- regmatches(filename, regexec("REP([0-9]+)", filename))
#     if (length(sim_id_match[[1]]) > 1) {
#       sim_id <- as.integer(sim_id_match[[1]][2])
#     } else {
#       sim_id_match2 <- regmatches(filename, regexec("TPI([0-9]+)", filename))
#       if (length(sim_id_match2[[1]]) > 1) {
#         sim_id <- as.integer(sim_id_match2[[1]][2])
#       } else {
#         sim_id <- as.integer(gsub("\\D", "", filename))
#         if (is.na(sim_id)) sim_id <- 0
#       }
#     }
# 
#     df_census <- df %>% dplyr::filter(Year == CENSUS_YEAR)
# 
#     if (nrow(df_census) == 0) {
#       pop.scaling.factor <- 1
#     } else {
#       total.pop <- sum(df_census$Population, na.rm = TRUE)
#       pop.scaling.factor <- ifelse(is.na(total.pop) || total.pop <= 0, 1, SA_CENSUS_POP / total.pop)
#     }
# 
#     df <- df %>%
#       dplyr::filter(Year >= ANALYSIS_START_YEAR & Year <= ANALYSIS_END_YEAR)
# 
#     if (nrow(df) == 0) {
#       cat(paste("    WARNING: No data in analysis period\n"))
#       return(NULL)
#     }
# 
#     df[[dose_col]] <- df[[dose_col]] * pop.scaling.factor
#     df[[infection_col]] <- df[[infection_col]] * pop.scaling.factor
#     df$Population <- df$Population * pop.scaling.factor
# 
#     df <- df %>%
#       dplyr::mutate(
#         Age_Group = ifelse(Age >= 15 & Age <= 49, "15-49", "0-99"),
#         Gender_Label = ifelse(Gender == 0, "Male", "Female"),
#         Simulation = sim_id,
#         Pop_Scaling_Factor = pop.scaling.factor
#       )
# 
#     result <- df %>%
#       dplyr::group_by(Gender_Label, Age_Group, Simulation) %>%
#       dplyr::summarise(
#         Cumulative_Doses = sum(.data[[dose_col]], na.rm = TRUE),
#         Cumulative_Infections = sum(.data[[infection_col]], na.rm = TRUE),
#         Mean_Annual_Doses = mean(.data[[dose_col]], na.rm = TRUE),
#         Pop_Scaling_Factor = dplyr::first(Pop_Scaling_Factor),
#         .groups = "drop"
#       )
# 
#     cat(paste("    ✓ Processed", basename(file_path),
#               "- Sim:", sim_id,
#               "- Cum doses:", round(sum(result$Cumulative_Doses)), "\n"))
# 
#     return(result)
# 
#   }, error = function(e) {
#     cat(paste("    ✗ ERROR reading", basename(file_path), ":", e$message, "\n"))
#     return(NULL)
#   })
# }
# 
# process_scenario_files_laprep_doses <- function(scenario_name, is_baseline = FALSE) {
#   scenario_path <- file.path(BASE_PATH, scenario_name, REPORT_SUBFOLDER)
# 
#   if (!dir.exists(scenario_path)) {
#     cat(paste("    ✗ ERROR: Directory not found:", scenario_path, "\n"))
#     return(NULL)
#   }
# 
#   csv_files <- list.files(scenario_path, pattern = "\\.csv$", full.names = TRUE)
#   if (length(csv_files) == 0) return(NULL)
# 
#   csv_files <- sort(csv_files)
#   if (TEST_MODE && length(csv_files) > MAX_FILES_TEST) csv_files <- csv_files[1:MAX_FILES_TEST]
# 
#   all_data <- list()
#   success_count <- 0
# 
#   for (i in seq_along(csv_files)) {
#     result <- read_emod_file_laprep_doses(csv_files[i])
#     if (!is.null(result)) {
#       all_data[[i]] <- result
#       success_count <- success_count + 1
#     }
#   }
# 
#   if (success_count == 0) return(NULL)
# 
#   combined_data <- data.table::rbindlist(all_data, fill = TRUE)
#   combined_data$Scenario <- scenario_name
#   combined_data$Scenario_Type <- ifelse(is_baseline, "Baseline", "Intervention")
#   return(combined_data)
# }
# 
# calculate_laprep_efficiency <- function(baseline_data, intervention_data, intervention_scenario_name) {
#   if (is.null(baseline_data) || is.null(intervention_data)) return(NULL)
# 
#   common_sims <- intersect(unique(baseline_data$Simulation),
#                            unique(intervention_data$Simulation))
# 
#   if (length(common_sims) == 0) {
#     stop("No common simulation IDs between baseline and intervention.")
#   }
# 
#   baseline_matched <- baseline_data %>%
#     dplyr::filter(Simulation %in% common_sims) %>%
#     dplyr::mutate(Match_ID = Simulation)
# 
#   intervention_matched <- intervention_data %>%
#     dplyr::filter(Simulation %in% common_sims) %>%
#     dplyr::mutate(Match_ID = Simulation)
# 
#   baseline_summary <- baseline_matched %>%
#     dplyr::group_by(Match_ID, Gender_Label, Age_Group) %>%
#     dplyr::summarise(
#       Baseline_Doses = sum(Cumulative_Doses, na.rm = TRUE),
#       Baseline_Infections = sum(Cumulative_Infections, na.rm = TRUE),
#       .groups = "drop"
#     )
# 
#   intervention_summary <- intervention_matched %>%
#     dplyr::group_by(Match_ID, Gender_Label, Age_Group) %>%
#     dplyr::summarise(
#       Intervention_Doses = sum(Cumulative_Doses, na.rm = TRUE),
#       Intervention_Infections = sum(Cumulative_Infections, na.rm = TRUE),
#       .groups = "drop"
#     )
# 
#   detailed <- baseline_summary %>%
#     dplyr::inner_join(
#       intervention_summary,
#       by = c("Match_ID", "Gender_Label", "Age_Group")
#     ) %>%
#     dplyr::mutate(
#       Doses_Delivered = Intervention_Doses - Baseline_Doses,
#       Infections_Averted = Baseline_Infections - Intervention_Infections,
#       Infections_Averted_per_1000_Doses = dplyr::if_else(
#         Doses_Delivered > 0,
#         (Infections_Averted / Doses_Delivered) * 1000,
#         NA_real_
#       )
#     )
# 
#   summary <- detailed %>%
#     dplyr::group_by(Gender_Label, Age_Group) %>%
#     dplyr::summarise(
#       N = dplyr::n(),
#       Mean_Infections_Averted = mean(Infections_Averted, na.rm = TRUE),
#       Mean_Doses_Delivered = mean(Doses_Delivered, na.rm = TRUE),
# 
#       Mean_Infections_Averted_per_1000_Doses =
#         mean(Infections_Averted_per_1000_Doses, na.rm = TRUE),
# 
#       SD_Efficiency =
#         sd(Infections_Averted_per_1000_Doses, na.rm = TRUE),
# 
#       CI_Lower =
#         Mean_Infections_Averted_per_1000_Doses -
#         qt(0.975, N - 1) * (SD_Efficiency / sqrt(N)),
# 
#       CI_Upper =
#         Mean_Infections_Averted_per_1000_Doses +
#         qt(0.975, N - 1) * (SD_Efficiency / sqrt(N)),
# 
#       Ratio_of_Means_per_1000 =
#         (Mean_Infections_Averted / Mean_Doses_Delivered) * 1000,
# 
#       .groups = "drop"
#     ) %>%
#     dplyr::mutate(
#       Scenario = intervention_scenario_name,
#       Analysis_Period = paste(ANALYSIS_START_YEAR, "-", ANALYSIS_END_YEAR),
#       Analysis_Years = ANALYSIS_PERIOD_YEARS
#     ) %>%
#     dplyr::select(
#       Scenario,
#       Analysis_Period,
#       Analysis_Years,
#       dplyr::everything()
#     )
# 
#   return(list(summary = summary, detailed = detailed))
# }
# 
# save_laprep_dose_results <- function(results_list, output_dir) {
#   if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# 
#   all_summaries <- do.call(rbind, lapply(results_list, function(x) x$summary))
#   write.csv(all_summaries, file.path(output_dir, "LAPrEP_doses_2026_2060_summary.csv"), row.names = FALSE)
# 
#   for (i in seq_along(results_list)) {
#     scenario_name <- names(results_list)[i]
#     safe_name <- gsub("[^A-Za-z0-9]", "_", scenario_name)
#     write.csv(results_list[[i]]$detailed,
#               file.path(output_dir, paste0("detailed_LAPrEP_doses_", safe_name, ".csv")),
#               row.names = FALSE)
#   }
# 
#   return(all_summaries)
# }
# 
# # ==================== 4. MAIN EXECUTION ====================
# baseline_data_laprep_doses <- process_scenario_files_laprep_doses(BASELINE_SCENARIO, is_baseline = TRUE)
# if (is.null(baseline_data_laprep_doses)) stop("ERROR: Failed to process baseline scenario.")
# 
# all_scenarios <- list.dirs(BASE_PATH, full.names = FALSE, recursive = FALSE)
# intervention_scenarios <- sort(all_scenarios[grepl("^MainAnalysis.*VOICE[2468]", all_scenarios)])
# ##intervention_scenarios <- sort(all_scenarios[grepl("^Sens\\d{4}", all_scenarios)])
# 
# all_results_laprep_doses <- list()
# 
# for (intervention_scenario in intervention_scenarios) {
#   intervention_data <- process_scenario_files_laprep_doses(intervention_scenario, is_baseline = FALSE)
#   if (!is.null(intervention_data)) {
#     result <- calculate_laprep_efficiency(
#       baseline_data_laprep_doses, intervention_data, intervention_scenario
#     )
#     if (!is.null(result)) all_results_laprep_doses[[intervention_scenario]] <- result
#   }
# }
# 
# if (length(all_results_laprep_doses) > 0) {
#   final_summary_laprep_doses <- save_laprep_dose_results(all_results_laprep_doses, OUTPUT_DIR)
#   cat("\n✓ LA PrEP doses distributed analysis complete\n")
# } else {
#   cat("\n✗ No results generated\n")
# }