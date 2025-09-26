library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/Users/kimberlymalloch/cplusplus/unchanged/MicroCOSMHPV/")
path <- ""
# Read and combine all 24 files
ProspCohort <- lapply(1:1, function(ii) {
	read.table(
		paste0(path, ii, "_0ProspectiveCohort.txt"),
		fill = TRUE,
		header = FALSE,
		sep = "",
		stringsAsFactors = FALSE
	)
}) %>%
	bind_rows() %>%
	as_tibble(.name_repair = "minimal") %>%
	setNames(paste0("V", seq_along(.)))

# HPV columns
hpv_cols <- paste0("V", 5:17)

# Rename some columns
ProspCohort <- ProspCohort %>%
	rename(
		Simulation = V1,
		Year = V2,
		ID = V3,
		ExactAge = V4
	)

# Type-specific prevalence
hpv_counts <- ProspCohort %>%
	pivot_longer(cols = all_of(hpv_cols),
							 names_to = "HPV_Type", values_to = "Stage") %>%
	filter(Stage %in% 1:5) %>%
	group_by(Year, HPV_Type) %>%
	summarise(Infected = n(), .groups = "drop")

total_per_year <- ProspCohort %>%
	group_by(Year) %>%
	summarise(Total = n(), .groups = "drop")

hpv_prevalence <- hpv_counts %>%
	left_join(total_per_year, by = "Year") %>%
	mutate(Prevalence = (Infected / Total) * 100)

# Map labels
hpv_labels <- c(
	"V5" = "HPV 16",
	"V6" = "HPV 18",
	"V7" = "HPV 31",
	"V8" = "HPV 33",
	"V9" = "HPV 35",
	"V10" = "HPV 39",
	"V11" = "HPV 45",
	"V12" = "HPV 51",
	"V13" = "HPV 52",
	"V14" = "HPV 56",
	"V15" = "HPV 58",
	"V16" = "HPV 59",
	"V17" = "HPV 68"
)

hpv_prevalence <- hpv_prevalence %>%
	mutate(HPV_Type = hpv_labels[HPV_Type])

# Plot type-specific prevalence
ggplot(hpv_prevalence, aes(x = Year, y = Prevalence, color = HPV_Type)) +
	geom_line(size = 1) +
	labs(
		title = "Type-Specific HPV Prevalence Over Years",
		y = "Prevalence (%)",
		x = "Year",
		color = "HPV Type"
	) +
	theme_minimal() +
	theme(legend.position = "right")

# Overall prevalence (infected with any HPV type)
overall_prevalence <- ProspCohort %>%
	mutate(AnyInfected = if_else(rowSums(across(all_of(hpv_cols), ~ .x %in% 1:5)) > 0, 1, 0)) %>%
	group_by(Year) %>%
	summarise(
		Infected = sum(AnyInfected),
		Total = n(),
		Prevalence = Infected / Total * 100,
		.groups = "drop"
	)

# Plot overall prevalence
ggplot(overall_prevalence, aes(x = Year, y = Prevalence)) +
	geom_line(size = 1.2, color = "darkred") +
	geom_point(color = "darkred") +
	labs(
		title = "Overall HPV Prevalence in 18-23 year old women",
		x = "Year",
		y = "Prevalence (%)"
	) +
	theme_minimal()

