#Multispectral Imaging and Machine Learning Can Accurately Predict Phenotypic Traits of Bluebunch Wheatgrass (Pseudoroegneria spicata)
#Brittany Johnson, Alexander Hernandez, Steve Larson, Kari E. Veblen, Efrain Duarte, Peter Porter, Holden Brecht, and Zayne Maughan


#### MAIN SCRIPT FOR THE ENTIRE ANALYSIS


## Libraries needed ##
library(Metrics)
library(valmetrics)
library(dplyr)
library(tidyr)
library(tidyplots)
library(ggplot2)
library(readr)
library(stringr)
library(doParallel)
library(randomForest)
library(e1071)
library(car)
library(creditmodel)
library(caret)
library(RColorBrewer)
library(patchwork)
library(stringr)
library(dplyr)
library(gridExtra)

# Data path for the entire project
# setwd("C:/data_1")

################################################################################
########################### FUNCTIONS ##########################################
metrics <- function(preds, actual, data) {
  rss <- sum((actual - preds)^2)
  tss <- sum((actual - mean(actual))^2)
  r_squared <- 1 - (rss / tss)
  n <- nrow(data)
  p <- ncol(data) - 1

  metric_list <- c(
    "R_squared" = r_squared,
    "adj_R_Squared" = 1 - ((1 - r_squared) * (n - 1)) / (n - p - 1),
    "mse" = mse(actual, preds),
    "rmse" = rmse(actual, preds),
    "mae" = mae(actual, preds),
    "lccc" = valmetrics::lccc(actual, preds)
  )
  metric_list
}

ph_plot <- function(data, avg_or_med) {
  if (avg_or_med == "average") {
    data |>
      rename(
        "UAV_measured" = mean_heights,
        "Infield_measured" = pth_tavg
      ) |>
      pivot_longer(
        cols = c(UAV_measured, Infield_measured),
        names_to = "variable",
        values_to = "value"
      ) |>
      group_by(variable, season) |>
      summarize(
        avg_height = mean(value),
        ymin = mean(value, na.rm = TRUE) - sd(value, na.rm = TRUE),
        ymax = mean(value, na.rm = TRUE) + sd(value, na.rm = TRUE)
      ) |>
      tidyplot(x = season, y = avg_height, color = variable) |>
      add_mean_bar(alpha = 0.4) |>
      add_data_points_beeswarm() |>
      add(ggplot2::geom_errorbar(aes(ymin = ymin, ymax = ymax),
                                 width = 0.2,
                                 position = position_dodge(width = 0.8)
      )) |>
      adjust_y_axis_title("Mean Height (cm)") |>
      adjust_x_axis_title("Season") |>
      adjust_colors(colors_discrete_metro[c(2, 4)]) |>
      adjust_size(NA, NA) |>
      add_title("Plant Height") |>
      adjust_font(fontsize = 12)
  } else {
    data |>
      rename(
        "UAV_measured" = heights,
        "Infield_measured" = pth_tmed
      ) |>
      pivot_longer(
        cols = c(UAV_measured, Infield_measured),
        names_to = "variable",
        values_to = "value"
      ) |>
      group_by(variable, season) |>
      summarize(
        avg_height = mean(value),
        ymin = mean(value, na.rm = TRUE) - sd(value, na.rm = TRUE),
        ymax = mean(value, na.rm = TRUE) + sd(value, na.rm = TRUE)
      ) |>
      tidyplot(x = season, y = avg_height, color = variable) |>
      add_mean_bar(alpha = 0.4) |>
      add_data_points_beeswarm() |>
      add(ggplot2::geom_errorbar(aes(ymin = ymin, ymax = ymax),
                                 width = 0.2,
                                 position = position_dodge(width = 0.8)
      )) |>
      adjust_y_axis_title("Median Height (cm)") |>
      adjust_x_axis_title("Season", fontsize = 10) |>
      adjust_colors(colors_discrete_metro[c(2, 4)]) |>
      adjust_size(NA, NA) |>
      add_title("Plant Height") |>
      adjust_font(fontsize = 12)
  }
}

create_scatterplots <- function(preds, actual, data, round, title) {
  metrics_vector <- metrics(preds, actual, data)

  test_data <- if (round == 3) {
    bind_cols(preds, actual) |>
      rename(
        preds = "...1",
        actual = "...2"
      )
  } else {
    bind_cols(preds, actual) |>
      rename(
        preds = "...1",
        actual = "...2"
      ) |>
      sample_frac(0.25)
  }

  scatterplot <- test_data |>
    tidyplot(x = actual, y = preds) |>
    add_data_points() |>
    add_curve_fit(se = TRUE) |>
    add_annotation_text(
      str_c(
        "R^2: ",
        formatC(
          metrics_vector["R_squared"],
          digits = 3,
          format = "f"
        )
      ),
      x = max(test_data$actual) - 0.5 * sd(test_data$actual),
      y = min(test_data$preds),
      fontsize = 8
    ) |>
    adjust_size(NA, NA) |>
    adjust_font(fontsize = 12) |>
    adjust_x_axis_title(title = "Observed", fontsize = 8) |>
    adjust_y_axis_title(title = "Prediction", fontsize = 8) |>
    add_title(as.character(title))


  scatterplot
}



##################### Initial Wrangling ##########################

# Read in the data
indiv_data <- read_csv("data/MergedSiteIndivPolys_Zstats_Join.csv")
agg_data <- read_csv("data/MergedSiteAggregatePolys_Zstats_Join.csv")

cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)


#### Begin wrangling
indiv_data |>
  #Filter out the Kernza plots
  filter(!(siteplot_id %in% c("MV778", "MV1556", "MV1654", "MV1834") |
             aggplot_id %in% c(195, 392, 414, 459))) |>
  # make more rows based on the predictors we want and it's date.
  pivot_longer(
    cols = `VOLUME06-05-2024`:`VOLUME10-30-2023`,
    names_to = c("variable", "date"),
    names_pattern = "(\\w+)([[:digit:]]{2}-[[:digit:]]{2}-[[:digit:]]{4}$)",
    values_to = "value"
  ) %>%
  # make the predictor value a column in itself
  pivot_wider(
    names_from = variable,
    values_from = value
  ) |>
  # separate the year and month column from the date
  mutate(
    year = str_extract(date, "[[:digit:]]+$"),
    month = str_extract(date, "^[[:digit:]]+")
  ) %>%
  # change the other columns that are also predictors to numeric.
  mutate(
    across(
      c(`BC_PTHTavg_24`:`MV_SDYLD_24`),
      ~ as.numeric(gsub("[^0-9.-]", "", .))
    ),
    na_count = rowSums(is.na(.) | sapply(., is.nan)) # Count NA/NaN per row
  ) %>%
  # make the other predictors have individual observations like before
  pivot_longer(
    cols = `BC_PTHTavg_24`:`MV_SDYLD_24`,
    names_to = c("variable", "year2"),
    names_pattern = "[[:alpha:]]_([[:alpha:]]+)_[[:alpha:]]{0,4}([[:digit:]]+$)",
    values_to = "value"
  ) |>
  # same as above make each predictor it's own column
  pivot_wider(
    names_from = variable,
    values_from = value,
    values_fn = ~ mean(.x, na.rm = TRUE)
  ) |>
  # focus on the frequency and fill in the missing areas with the other column
  mutate(
    year2 = str_c("20", as.character(year2)),
  ) |>
  # select which columns we want
  select(-c(
    "IndiviCts", "rpm",
    "FreqIndv", "SDYLD", "PTHTavg", "PTHTmed",
    "wetbiomass", "drybiomass", "Mort"
  )) |>
  # which groups we want the data to be seen
  group_by(
    fid, year2, year, plot_id, site, month, date
  ) |>
  # if there are multiple observations take the mean
  summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop"
  ) %>%
  # from the second group only take the rows where these are the same
  filter(year2 == year) %>%
  # filter where the observations are present (all the wanted
  #observations are present in the volume column so we filter on that)
  filter(
    !is.nan(VOLUME)
  ) %>%
  # drop these two additional columns
  select(-c(na_count, year2)) -> data

# filter the data into just millville and the selected dates for future analysis
data %>%
  filter(site == "Millville" &
           date %in% c("06-10-2022", "07-08-2022")) -> millville_freq

# filter the data into just bluecreek and the selected dates for future analysis
data %>%
  filter(site == "BlueCreek" &
           date %in% c("07-18-2024", "07-31-2024")) -> bluecreek_freq
# Free unused memory
gc()


# WRANGLE THE DATA INTO THE NEEDED FORMAT
agg_data |>
  pivot_longer(
    col = `VOLUME06-05-2024`:`VOLUME10-30-2023`,
    names_to = c("variable", "date"),
    names_pattern = "(\\w+)([[:digit:]]{2}-[[:digit:]]{2}-[[:digit:]]{4}$)",
    values_to = "value"
  ) |>
  pivot_wider(
    names_from = variable,
    values_from = value,
    values_fn = ~ mean(.x, na.rm = TRUE)
  ) |>
  mutate(
    across(c(`BC_aggLAI_24`:`MV_aggCCgrswd_24`), ~ as.numeric(gsub("[^0-9.-]", "", .)))
  ) |>
  pivot_longer(
    cols = `BC_aggLAI_24`:`MV_aggCCgrswd_24`,
    names_to = c("variable", "year2"),
    names_pattern = "[[:alpha:]]_[ag]{3}([[:alpha:]]+)_[[:alpha:]]{0,4}([[:digit:]]+$)",
    values_to = "value"
  ) |>
  pivot_wider(
    names_from = variable,
    values_from = value,
    values_fn = ~ mean(.x, na.rm = TRUE)
  ) |>
  mutate(
    year = str_extract(date, "[[:digit:]]+$"),
    year2 = str_c("20", as.character(year2))
  ) |>
  group_by(
    fid, year2, year, aggplot_id, site, date
  ) |>
  # if there are multiple observations take the mean
  summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  filter(site == "BlueCreek" | # Keep all rows where site is BlueCreek
           (site == "Millville" & !(year == 2022 & aggplot_id %in%
                                      c(190, 195, 288, 294, 328, 392, 414, 459)))) |>
  # from the second group only take the rows where these are the same
  filter(year2 == year) %>%
  # filter where the observations are present (all the wanted observations
  # are present in the volume column so we filter on that)
  filter(
    !is.nan(VOLUME) & !is.nan(LAI)
  ) -> agg_traits

indiv_data |>
  # make more rows based on the predictors we want and it's date.
  pivot_longer(
    cols = `VOLUME06-05-2024`:`VOLUME10-30-2023`,
    names_to = c("variable", "date"),
    names_pattern = "(\\w+)([[:digit:]]{2}-[[:digit:]]{2}-[[:digit:]]{4}$)",
    values_to = "value"
  ) %>%
  # make the predictor value a column in itself
  pivot_wider(
    names_from = variable,
    values_from = value
  ) |>
  # separate the year and month column from the date
  mutate(
    year = str_extract(date, "[[:digit:]]+$"),
    month = str_extract(date, "^[[:digit:]]+")
  ) %>%
  # make more rows based on the predictors we want and it's date.
  # change the other columns that are also predictors to numeric.
  mutate(
    across(c(`BC_PTHTavg_24`:`MV_SDYLD_24`), ~ as.numeric(gsub("[^0-9.-]", "",
                                                               .))),
  ) %>%
  # make the other predictors have individual observations like before
  pivot_longer(
    cols = `BC_PTHTavg_24`:`MV_SDYLD_24`,
    names_to = c("variable", "year2"),
    names_pattern = "[[:alpha:]]_([[:alpha:]]+)_[[:alpha:]]{0,4}([[:digit:]]+$)",
    values_to = "value"
  ) |>
  # same as above make each predictor it's own column
  pivot_wider(
    names_from = variable,
    values_from = value,
    values_fn = ~ mean(.x, na.rm = TRUE)
  ) |>
  mutate(
    year2 = str_c("20", as.character(year2)),
  ) |>
  filter(!(site == "Millville" &
             aggplot_id %in% c(190, 195, 288, 294, 328, 392, 414, 459) &
             siteplot_id %in% c(
               "MV759", "MV778", "MV1151", "MV1175", "MV1311",
               "MV1556", "MV1654", "MV1834"
             ) &
             year == "2022")) |>
  # from the second group only take the rows where these are the same
  filter(year2 == year) %>%
  filter(!is.nan(VOLUME) & !is.nan(drybiomass)) |>
  mutate(
    wetbiomass_adj = (wetbiomass / 1000) / 0.00005,
    drybiomass_adj = (drybiomass / 1000) / 0.00005
  ) |>
  # which groups we want the data to be seen
  group_by(
    fid, aggplot_id, plot_id, date, site, year, year2
  ) |>
  # if there are multiple observations take the mean
  summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  group_by(fid) %>%
  # summarize(
  #   across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
  #   across(where(is.character), ~ first(.x))
  # ) %>%
  ungroup() |>
  select(
    -c(FreqPercent, FreqPerc, IndiviCts, rpm, Mort, FreqIndv, PTHTavg, PTHTmed)
  ) -> indiv_traits


gc()
################################################################################
################### PLANT HEIGHTS ###################
ph_zstat <- read_csv("data/MergedSiteIndivPolys_PTHT_Zstats.csv")


ph_zstat |>
  pivot_longer(
    col = `mean-HEIGHTS06-05-2024`:`HEIGHTS10-30-2023`,
    names_to = c("variable", "date"),
    names_pattern = "(\\w+[-]{0,1}\\w+)([[:digit:]]{2}-[[:digit:]]{2}-[[:digit:]]{4}$)",
    values_to = "value"
  ) |>
  pivot_wider(
    names_from = variable,
    values_from = value,
    values_fn = ~ mean(.x, na.rm = TRUE)
  ) |>
  mutate(
    `mean-HEIGHTS` = `mean-HEIGHTS` * 100,
    HEIGHTS = HEIGHTS * 100,
    season = ifelse(site == "Millville" & str_extract(date, "\\d{4}$") == "2024",
                    "MVsum24",
                    ifelse(site == "Millville" & str_extract(
                      date,
                      "\\d{4}$"
                    ) == "2022",
                    "MVsum22",
                    ifelse(site == "Millville" & str_extract(date, "\\d{4}") == "2023",
                           ifelse(date %in% c(
                             "05-03-2023", "05-11-2023", "06-05-2023",
                             "06-21-2023", "07-06-2023", "07-27-2023"
                           ),
                           "MVsum23", "MVfall23"
                           ), "BCsum24"
                    )
                    )
    )
  ) |>
  group_by(
    fid, plot_id, aggplot_id, site, date, season
  ) |>
  # if there are multiple observations take the mean
  summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  filter(!is.nan(HEIGHTS)) -> observed_data

measured_data <- read_csv("data/MergedSiteIndivPolys_Zstats_join.csv") |>
  select(siteplot_id:MV_SDYLD_24) |>
  # change the other columns that are also predictors to numeric.
  mutate(
    across(c(`BC_PTHTavg_24`:`MV_SDYLD_24`), ~ as.numeric(gsub("[^0-9.-]", "", .))),
    # na_count = rowSums(is.na(.) | sapply(., is.nan))  # Count NA/NaN per row
  ) |>
  pivot_longer(
    col = `BC_PTHTavg_24`:MV_SDYLD_24,
    names_to = c("variable", "season"),
    names_pattern = "\\w{2}_(\\w+)_([[:alpha:]]{0,4}\\d+)$",
    values_to = "value"
  ) |>
  pivot_wider(
    names_from = variable,
    values_from = value,
    values_fn = ~ mean(.x, na.rm = TRUE)
  ) |>
  select(siteplot_id, season, PTHTavg, PTHTmed) |>
  filter(!is.na(PTHTavg)) |>
  mutate(
    year = as.character(str_extract(season, "\\d+$")),
    site = ifelse(str_detect(siteplot_id, "BC"), "BlueCreek", "Millville"),
    plot_id = as.numeric(str_extract(siteplot_id, "\\d+$"))
  ) |>
  mutate(
    season = ifelse(site == "Millville" & year %in% c("24", "22"), str_c("MVsum", year),
                    ifelse(site == "BlueCreek", "BCsum24",
                           ifelse(season %in% c("fall23", 'sum23'),
                                  str_c("MV", season), season))
    )
  )


plots_data <- inner_join(observed_data,
                         measured_data,
                         by = c("site", "plot_id", "season")
) |>
  janitor::clean_names()

plots_data$season <- factor(plots_data$season,
                            levels = c("MVsum22", "MVsum23",
                                       "MVfall23", "MVsum24", "BCsum24"))




test <- plots_data |> rename("UAV_measured" = mean_heights,
                            "Infield_measured" = pth_tavg) |>
  pivot_longer(
    cols = c(UAV_measured, Infield_measured),
    names_to = "variable",
    values_to = "value"
  )


plots_data |>
  rename(
    "UAV_measured" = heights,
    "Infield_measured" = pth_tmed
  ) |>
  pivot_longer(
    cols = c(UAV_measured, Infield_measured),
    names_to = "variable",
    values_to = "value"
  ) |>
  group_by(variable, season)|>
  summarize(avg_height = mean(value),
            ymin = mean(value, na.rm = TRUE) - sd(value, na.rm = TRUE),
            ymax = mean(value, na.rm = TRUE) + sd(value, na.rm = TRUE)) -> plot_test
plot_test |>
  tidyplot(x = season, y = avg_height, color = variable) |>
  add_mean_bar(alpha = 0.4)|>
  add_data_points_beeswarm() |>
  add(ggplot2::geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2,
                             position = position_dodge(width = 0.8))) |>
  adjust_y_axis_title("Average Height (cm)") |>
  adjust_x_axis_title("Season", fontsize = 10) |>
  adjust_colors(colors_discrete_metro[c(2,4)]) |>
  adjust_size(NA, NA) |>
  add_title("Plant Height") |>
  adjust_font(fontsize = 12)


sum(is.na(plot_test$ymin))

ph_plot(plots_data, 'average')
ph_plot(plots_data, 'median')

###############################################################################
############################### FREQUENCY MODELS ##############################
#############################  FILTER & WRANGLE ########################
#Combine the millville and blue creek selected datasets
freq_data <- bind_rows(millville_freq, bluecreek_freq)

#change the response variables so all are in one column and drop previous columns
freq_data %>%
  mutate(freq_perc = ifelse(site == 'Millville', FreqPerc, FreqPercent)) %>%
  select(-c("FreqPerc", "FreqPercent")) -> freq_data

#REMOVE THOSE VARIABLES
freq_data %>%
  select(-c("BGI", "GLA", "SI", "NGRDI")) -> freq_data

#scale the data (ensure that all predictors have same range and one
#cannot overpower or cause faulty analysis)
freq_data[,8:39] <- scale(freq_data[,8:39])


################################# VIF ###########################################
#FIT A SIMPLE MODEL TO GET THE VIF VALUES
freq_model <- glm(freq_perc ~., data = freq_data[,8:39])
#USING THIS FUNCTION SEE WHICH VALUES ARE TOO INTERDEPENDENT ON OTHER VARIABLES
alias(freq_model)


#FIT MODEL AGAIN USING PRE-SELECTION FOR VIF
freq_model <- glm(freq_perc ~., data= freq_data[,8:39])

#YOU DON"T REALLY NEED THE FOLLOWING ALIAS BUT CAN LOOK AGAIN
#alias(freq_model)

########START VIF SECTION
vif_mv <- vif(freq_model)

# Visualizing VIF
barplot(vif_mv, main = "Variance Inflation Factor (VIF)", las=2)

#DETERMINE WHICH VIF VALUES ARE GREATER THAN OR EQUAL TO FIVE
vif_mv[vif_mv >= 5]
############################ ROUND 1 #####################################
# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
iteration1_data <- train_test_split(freq_data[,8:39], prop = 0.7, seed = 43)

#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)

#FIT THE MODEL AND CROSS VALIDATE IT
rf_freq1 <- train(freq_perc ~ ., data = iteration1_data$train,
                  method='rf',
                  trControl = repeat_cv,
                  metric = 'Rsquared')

#PREDICT USING THE TEST SET AND THE RANDOM FOREST MODEL.
preds_freq1 <- predict(rf_freq1, newdata = iteration1_data$test)
metrics(preds = preds_freq1,
        actual = iteration1_data$test$freq_perc,
        iteration1_data$test)

############################ ROUND 2 #####################################
#FILTER the data
freq_data %>%
  select(1:7, "CIRE", "NDVI", "NDRE", "RVI", "TVI", "freq_perc") -> freq_data
# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
iteration1_2_data <- train_test_split(freq_data[,8:13], prop = 0.7, seed = 43)

#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)

#FIT THE MODEL AND CROSS VALIDATE IT
rf_freq1_2 <- train(freq_perc~., data = iteration1_2_data$train,
                    method='rf',
                    trControl = repeat_cv,
                    metric = 'Rsquared')

#PREDICT USING THE TEST SET AND THE RANDOM FOREST MODEL.
preds_freq1_2 <- predict(rf_freq1_2, newdata = iteration1_2_data$test)
metrics(preds = preds_freq1,
        actual = iteration1_2_data$test$freq_perc,
        iteration1_2_data$test)


############################ ROUND 3 #####################################
#FILTER the data
freq_data %>%
  select(1:7, "NDVI", "RVI", "TVI", "freq_perc") -> freq_data
# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
iteration1_3_data <- train_test_split(freq_data[,8:11], prop = 0.7, seed = 43)

#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)

#FIT THE MODEL AND CROSS VALIDATE IT
rf_freq1_3 <- train(freq_perc~., data = iteration1_3_data$train,
                    method='rf',
                    trControl = repeat_cv,
                    metric = 'Rsquared')

#PREDICT USING THE TEST SET AND THE RANDOM FOREST MODEL.
preds_freq1_3 <- predict(rf_freq1_3, newdata = iteration1_3_data$test)
metrics(preds = preds_freq1,
        actual = iteration1_3_data$test$freq_perc,
        iteration1_3_data$test)


########################## VARIABLE IMPORTANCE ################################

#LOOK AT THE VARIABLE IMPORTANCE
var_imp <- varImp(rf_freq1, scale=FALSE)$importance
var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)

## Create a plot of variable importance
var_imp %>%

  ## Sort the data by importance
  arrange(importance) %>%

  ## Create a ggplot object for aesthetic
  ggplot(aes(x=reorder(variables, importance), y=importance)) +

  ## Plot the bar graph
  geom_bar(stat='identity') +

  ## Flip the graph to make a horizontal bar plot
  coord_flip() +

  ## Add x-axis label
  xlab('Variables') +

  ## Add a title
  labs(title='Random forest variable importance') +

  ## Some layout for the plot
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20),
  )
############PREPARE DATA #############
#Combine the millville and blue creek selected datasets
freq_data <- bind_rows(millville_freq, bluecreek_freq)
#change the response variables so all are in one column and drop previous columns
freq_data %>%
  mutate(freq_perc = ifelse(site == 'Millville', FreqPerc, FreqPercent)) %>%
  select(-c("FreqPerc", "FreqPercent")) -> freq_data

#REMOVE THOSE VARIABLES
freq_data %>%
  select(-c("BGI", "GLA", "SI", "NGRDI")) -> freq_data
################################# FILTER ####################################
#Filter the data to try and build it using these dates
iteration2 <- freq_data |>
  filter(date %in% c("06-10-2022", "07-18-2024"))

#scale the data (ensure that all predictors have same range and one cannot
#overpower or cause faulty analysis), scaling the data (mean=0, sd=1)
iteration2[,8:38] <- scale(iteration2[,8:38])

############# FILTER VARIABLES #################
##Removed after initial VIF analysis
iteration2 %>%
  select(-c("VOLUME")) -> iteration2


########################## VIF ANALYSIS #####################################
#FIT A SIMPLE MODEL TO GET THE VIF VALUES
freq_model2 <- glm(freq_perc ~., data= iteration2[,8:38])
#USING THIS FUNCTION SEE WHICH VALUES ARE TOO INTERDEPENDENT ON OTHER VARIABLES
alias(freq_model2)

vif_mv2 <- vif(freq_model2)

# Visualizing VIF
barplot(vif_mv2, main = "Variance Inflation Factor (VIF)", las=2)

# Visualizing the logarithm of VIF
barplot(log(vif_mv2), main = "Log of Variance Inflation Factor (VIF)", las = 2)

#DETERMINE WHICH VIF VALUES ARE GREATER THAN OR EQUAL TO FIVE
# Check for high multicollinearity (threshold = 5) and highlight problematic variables
if (any(vif_mv2 > 5)) {
  # Extract variables with VIF > 5 and their corresponding values
  high_vif_vars <- vif_mv2[vif_mv2 > 5]  # Get VIF values for problematic variables

  # Sort the variables by VIF values in ascending order
  high_vif_vars <- sort(high_vif_vars)

  # Print the sorted high VIF values and the variable names
  print("High multicollinearity detected! Variables with VIF > 5 (sorted by VIF values):")
  print(high_vif_vars)

  # Print the number of variables with high VIF
  print(paste("Number of variables with VIF > 5:", length(high_vif_vars)))
} else {
  print("No significant multicollinearity detected.")
}

######################## TRAIN THE MODEL ######################################

######################## ROUND 1 ######################################
#SPLIT THE DATA INTO A TEST AND TRAIN SET
iteration2_data <- train_test_split(iteration2[,8:38], prop = 0.7)

#SET SEED FOR SIMILAR RESULTS
set.seed(24)
#CREATE A TEMPLATE FOR CROSS VALIDATION
repeat_cv <- trainControl(method = 'oob', number = 10)

#FIT THE MODEL
rf_2 <- train(freq_perc~., data = iteration2_data$train,
              method='rf',
              trControl = repeat_cv,
              metric = 'Rsquared')


###############PREDICT ON THE DATA###################
preds_freq2 <- predict(rf_2, newdata = iteration2_data$test)
metrics(preds = preds_freq2, actual = iteration2_data$test$freq_perc,
        data = iteration2_data$test)
################################## ROUND 2 #####################################
iteration2 %>%
  select(1:7, "HUE", "HI","TVI", "GLI",
         "RVI", "RGBVI", "NDVI", "freq_perc") -> iteration2

#SPLIT THE DATA INTO A TEST AND TRAIN SET
iteration2_2data <- train_test_split(iteration2[,8:15], prop = 0.7)

#SET SEED FOR SIMILAR RESULTS
set.seed(24)
#CREATE A TEMPLATE FOR CROSS VALIDATION
repeat_cv <- trainControl(method = 'oob', number = 10)

#FIT THE MODEL
rf_2_2 <- train(freq_perc~., data = iteration2_2data$train,
                method='rf',
                trControl = repeat_cv,
                metric = 'Rsquared')


###############PREDICT ON THE DATA###################
preds_freq2_2 <- predict(rf_2_2, newdata = iteration2_2data$test)
metrics(preds = preds_freq2, actual = iteration2_2data$test$freq_perc,
        data = iteration2_2data$test)

######################## ROUND 3 ######################################
iteration2 %>%
  select(1:7, "HUE", "HI", "RVI", "RGBVI","NDVI", "freq_perc") -> iteration2

#SPLIT THE DATA INTO A TEST AND TRAIN SET
iteration2_3data <- train_test_split(iteration2[,8:13], prop = 0.7)

#SET SEED FOR SIMILAR RESULTS
set.seed(24)
#CREATE A TEMPLATE FOR CROSS VALIDATION
repeat_cv <- trainControl(method = 'oob', number = 10)

#FIT THE MODEL
rf_2_3 <- train(freq_perc~., data = iteration2_3data$train,
                method='rf',
                trControl = repeat_cv,
                metric = 'Rsquared')


###############PREDICT ON THE DATA###################
preds_freq2_3 <- predict(rf_2_3, newdata = iteration2_3data$test)
metrics(preds = preds_freq2, actual = iteration2_3data$test$freq_perc,
        data = iteration2_3data$test)


##################### Post Model Selection ####################################


#CREATE THE DATA FOR VARIABLE IMPORTANCE
var_imp_it2 <- varImp(rf_2, scale=FALSE)$importance
var_imp_it2 <- data.frame(variables=row.names(var_imp_it2),
                          importance=var_imp_it2$Overall)

## Create a plot of variable importance
var_imp_it2 %>%

  ## Sort the data by importance
  arrange(importance) %>%

  ## Create a ggplot object for aesthetic
  ggplot(aes(x=reorder(variables, importance), y=importance)) +

  ## Plot the bar graph
  geom_bar(stat='identity') +

  ## Flip the graph to make a horizontal bar plot
  coord_flip() +

  ## Add x-axis label
  xlab('Variables') +

  ## Add a title
  labs(title='Random forest variable importance') +

  ## Some layout for the plot
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20),
  )

################################# FILTER THE DATA #########################
millville_freq |>
  rename(freq_percent = FreqPerc) |>
  select(-c("FreqPercent","BGI", "GLA", "SI", "NGRDI")) -> millville_freq

bluecreek_freq |>
  rename(freq_percent = FreqPercent) |>
  select(-c("FreqPerc", "BGI", "GLA", "SI", "NGRDI")) -> bluecreek_freq

millville_freq %>%
  select(-c("VOLUME")) -> millville_freq
bluecreek_freq %>%
  select(-c("VOLUME")) -> bluecreek_freq

######################## VIF ANALYSIS ###################################
###########VIF#############################
#FIT A SIMPLE MODEL TO GET THE VIF VALUES
freq_model3 <- glm(freq_percent ~., data= millville_freq[,8:38])
#USING THIS FUNCTION SEE WHICH VALUES ARE TOO INTERDEPENDENT ON OTHER VARIABLES
alias(freq_model3)

vif_mv3 <- vif(freq_model3)

# Visualizing VIF
barplot(vif_mv3, main = "Variance Inflation Factor (VIF)", las=2)

# Visualizing the logarithm of VIF
barplot(log(vif_mv3), main = "Log of Variance Inflation Factor (VIF)", las = 2)

#DETERMINE WHICH VIF VALUES ARE GREATER THAN OR EQUAL TO FIVE
# Check for high multicollinearity (threshold = 5) and highlight problematic variables
if (any(vif_mv3 > 5)) {
  # Extract variables with VIF > 5 and their corresponding values
  high_vif_vars <- vif_mv3[vif_mv3 > 5]  # Get VIF values for problematic variables

  # Sort the variables by VIF values in ascending order
  high_vif_vars <- sort(high_vif_vars)

  # Print the sorted high VIF values and the variable names
  print("High multicollinearity detected! Variables with VIF > 5 (sorted by VIF values):")
  print(high_vif_vars)

  # Print the number of variables with high VIF
  print(paste("Number of variables with VIF > 5:", length(high_vif_vars)))
} else {
  print("No significant multicollinearity detected.")
}





######################## TRAIN THE MODEL AND GET PREDICTIONS #############
#SET SEED
set.seed(24)
#TEMPLATE FOR CROSS_VALIDATION
repeat_cv <- trainControl(method = 'boot', number = 3)

# pretrain <- proc.time()
#FIT TEH MODEL
rf_3 <- train(freq_percent~., data = millville_freq[,8:38],
              method='rf',
              trControl = repeat_cv,
              metric = 'RMSE')
#proc.time() - pretrain



#PREDICT AND REPORT THE METRICS
predsfreq_3 <- predict(rf_3, newdata = bluecreek_freq[,8:38])
metrics(preds = predsfreq_3, actual = bluecreek_freq$freq_percent,
        bluecreek_freq)



############################VARIABLE IMPORTANCE ANALYSIS #####################
#CREATE THE DATA FOR THE VIF PLOT
var_imp_it3 <- varImp(rf_3, scale=FALSE)$importance
var_imp_it3 <- data.frame(variables=row.names(var_imp_it3),
                          importance=var_imp_it3$Overall)

## Create a plot of variable importance
var_imp_it3 %>%

  ## Sort the data by importance
  arrange(importance) %>%

  ## Create a ggplot object for aesthetic
  ggplot(aes(x=reorder(variables, importance), y=importance)) +

  ## Plot the bar graph
  geom_bar(stat='identity') +

  ## Flip the graph to make a horizontal bar plot
  coord_flip() +

  ## Add x-axis label
  xlab('Variables') +

  ## Add a title
  labs(title='Random forest variable importance') +

  ## Some layout for the plot
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20),
  )
################################## BIOMASS #####################################
indiv_traits |>
  filter((site == "Millville" & year == "2022" & date %in% c(
    "09-26-2022",
    "07-08-2022"
  )) |
    (site == "BlueCreek" & year == "2024" & date %in% c(
      "07-31-2024",
      "07-18-2024"
    ))) |>
  select(-c(fid, year, year2, aggplot_id, SDYLD)) -> bio_data

# SELECT COLUMNS WANT TO exCLUDE
bio_data %>%
  select(-c(GLA, BGI, SCI, NGRDI)) -> bio_data

### Do not scale response vars
bio_data[, 4:34] <- scale(bio_data[, 4:34])

bio_data %>%
  select(
    1:3, red, PSRI, GRVI,
    GLI, VARI, MGRVI, "HUE", "drybiomass_adj"
  ) -> bio_data
###################################Modeling###################################
# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
it1_dry <- train_test_split(bio_data[, 4:11], prop = 0.7, seed = 43)

# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = "oob", number = 10)


rf_dry <- train(drybiomass_adj ~ .,
                data = it1_dry$train,
                method = "rf",
                trControl = repeat_cv,
                metric = "Rsquared"
)


preds_dry_bio_global <- predict(rf_dry, newdata = it1_dry$test)
metrics(
  preds = preds_dry_bio_global, actual = it1_dry$test$drybiomass_adj,
  it1_dry$test
)

################################################################################
# VAR IMPORTANCE


# LOOK AT THE VARIABLE IMPORTANCE
var_imp_bio <- varImp(rf_dry, scale = FALSE)$importance
var_imp_bio <- data.frame(variables = row.names(var_imp_bio), importance = var_imp_bio$Overall)

## Create a plot of variable importance
var_imp %>%
  ## Sort the data by importance
  arrange(importance) %>%
  ## Create a ggplot object for aesthetic
  ggplot(aes(x = reorder(variables, importance), y = importance)) +

  ## Plot the bar graph
  geom_bar(stat = "identity") +

  ## Flip the graph to make a horizontal bar plot
  coord_flip() +

  ## Add x-axis label
  xlab("Variables") +

  ## Add a title
  labs(title = "Random forest variable importance") +

  ## Some layout for the plot
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 20),
  )

#############################PREPARE THE DATA###################################
indiv_traits |>
  filter(site == 'Millville' & year == '2023'&
           date %in% c('05-03-2023','05-11-2023',
                       '06-05-2023', '06-21-2023',
                       '07-06-2023', '07-26-2023', '10-30-2023')) |>
  select(-c(fid, year, year2, aggplot_id, SDYLD)) -> bio_data23

bio_data23 %>%
  select(-c(GLA, BGI, SCI, NGRDI)) -> bio_data23

###Do not scale response vars
bio_data23[,4:34] <- scale(bio_data23[,4:34])

bio_data23 %>%
  select(1:3, CVI, NGBDI, SI, RGBVI, "drybiomass_adj") -> bio_data23

#############################MODELING###########################################
# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
it2_dry23 <- train_test_split(bio_data23[,4:8], prop = 0.7, seed = 43)

#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)


rf_dry23 <- train(drybiomass_adj ~ ., data = it2_dry23$train,
                  method='rf',
                  trControl = repeat_cv,
                  metric = 'Rsquared')

preds_bio2_23 <- predict(rf_dry23, newdata = it2_dry23$test)
metrics(preds_bio2_23, it2_dry23$test$drybiomass_adj, it2_dry23$test)

################################2024############################################
indiv_traits |>
  filter(site == 'Millville' & year == '2024'&
           date %in% c('06-05-2024', '06-25-2024', '07-11-2024')) |>
  select(-c(fid, year, year2, aggplot_id, SDYLD)) -> bio_data24

bio_data24 %>%
  select(-c(GLA, BGI, SCI, NGRDI)) -> bio_data24

###Do not scale response vars
bio_data24[,4:34] <- scale(bio_data24[,4:34])

bio_data24 %>%
  select(1:3, CVI, green, blue, "drybiomass_adj") -> bio_data24

# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
it2_dry24 <- train_test_split(bio_data24[,4:7], prop = 0.7, seed = 43)

#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)


rf_dry24 <- train(drybiomass_adj~., data = it2_dry24$train,
                  method='rf',
                  trControl = repeat_cv,
                  metric = 'Rsquared')

preds_bio2_24 <- predict(rf_dry24, newdata = it2_dry24$test)
metrics(preds_bio2_24, it2_dry24$test$drybiomass_adj, it2_dry24$test)




###############################################################################
###################################WRANGLING ##################################

indiv_traits |> filter(site == 'Millville' & year != '2024')|>
  arrange(year, fid, plot_id) |>
  group_by(plot_id, fid, year) |>
  summarise(drybiomass_adj = mean(drybiomass_adj),
            wetbiomass_adj = mean(wetbiomass_adj)) |>
  rename(prev_dryadj = drybiomass_adj,
         prev_wetadj = wetbiomass_adj) |>
  ungroup()|>
  mutate(year = case_when(
    year == '2023' ~ '2024',
    year == '2022' ~ '2023',
    TRUE ~ '0'  # Handle other years
  )) -> prev_year


indiv_traits |> filter(site == 'Millville' & year != '2022') |>
  arrange(year, fid, plot_id) |>
  inner_join(prev_year, relationship = 'many-to-one') -> bio3_data23
################################################################################
#                                SPLIT DATA UP
bio3_data23 |>
  filter(year == '2023'&
           date %in% c('05-03-2023','05-11-2023',
                       '06-05-2023', '06-21-2023', '07-06-2023',
                       '07-26-2023','10-30-2023'))|>
  select(-c(fid, year, year2, SDYLD)) -> bio3_data23

bio3_data23 %>%
  select(-c(GLA, BGI, SCI, NGRDI)) -> bio3_data23

###Do not scale response vars
bio3_data23[,5:35] <- scale(bio3_data23[,5:35])

indiv_traits |> filter(site == 'Millville' & year != '2022') |>
  arrange(year, fid, plot_id) |>
  inner_join(prev_year, relationship = 'many-to-one') -> bio3_data24

bio3_data24 |>
  filter(site == 'Millville' & year == '2024' &
           date %in% c('06-05-2024', '06-25-2024','07-11-2024')) |>
  select(-c(fid, year, year2, SDYLD)) -> bio3_data24

bio3_data24 %>%
  select(-c(GLA, BGI, SCI, NGRDI)) -> bio3_data24

###Do not scale response vars
bio3_data24[,5:35] <- scale(bio3_data24[,5:35])

########################FILTER VARIABLES########################################
bio3_data24 %>%
  select(1:4, CIG, EVI, "prev_dryadj", "drybiomass_adj") -> bio3_data24

bio3_data23 |>
  select(1:4, NGBDI, prev_dryadj, drybiomass_adj) -> bio3_data23


###################################### MODELING ################################
it3_bio_23 <- train_test_split(bio3_data23[,5:7], prop = 0.7)
#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)


rf_bio_23 <- train(drybiomass_adj~., data = it3_bio_23$train,
                   method='rf',
                   trControl = repeat_cv,
                   metric = 'Rsquared')


preds_bio3_23 <- predict(rf_bio_23, newdata = it3_bio_23$test)
metrics(preds = preds_bio3_23, actual = it3_bio_23$test$drybiomass_adj,
        it3_bio_23$test)
##############################################################################
########################## 2024 ################################################
it3_bio_24 <- train_test_split(bio3_data24[,5:8], prop = 0.7)
#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)


rf_bio_24 <- train(drybiomass_adj~., data = it3_bio_24$train,
                   method='rf',
                   trControl = repeat_cv,
                   metric = 'Rsquared')


preds_bio3_24 <- predict(rf_bio_24, newdata = it3_bio_24$test)
metrics(preds = preds_bio3_24, actual = it3_bio_24$test$drybiomass_adj,
        it3_bio_24$test)

###############################################################################
###############################################################################
############################## CANOPY COVER ###################################
#### CANOPY COVER#####

### THE ONLY THING THAT CHANGES FROM OTHER METHODS OR ITERATIONS IS THE DATA

agg_traits |>
  # FILTER ON THE CONDITIONS WE WANT
  filter((site == "Millville" & year == "2022" &
            date %in% c("09-26-2022", "07-08-2022")) |
           (site == "BlueCreek" & year == "2024" &
              date %in% c("07-31-2024", "07-18-2024"))) |>
  # SELECT COLUMNS WE WANT TO INCLUDE
  select(-c(fid, LAI, GLA, BGI, SCI, NGRDI)) |> tidyr::drop_na() -> CC_data

### Do not scale response vars
CC_data[, 6:36] <- scale(CC_data[, 6:36])

######################## FILTER VARIABLES########################################
CC_data %>%
  select(1:5, PSRI, HUE, CCgrs) -> CC_data_grs


##############################################################################
# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
iteration1_data_grs <- train_test_split(CC_data_grs[, 6:8],
                                        prop = 0.7, seed = 43
)

# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = "oob", number = 10)

#train the model
rf_grass <- train(CCgrs ~ .,
                  data = iteration1_data_grs$train,
                  method = "rf",
                  trControl = repeat_cv,
                  metric = "Rsquared"
)


########### TEST #####################
preds_grs <- predict(rf_grass, newdata = iteration1_data_grs$test)
metrics(
  preds = preds_grs, actual = iteration1_data_grs$test$CCgrs,
  iteration1_data_grs$test
)



##################### GRSWD SECTION############################################
######################## FILTER VARIABLES########################################
CC_data %>%
  select(1:5, "PSRI", "HUE", "CCgrswd") -> CC_data
###############################################################################

# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
iteration1_data_grswd <- train_test_split(CC_data[, 6:8], prop = 0.7, seed = 43)

# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = "oob", number = 10)

#Train the model
rf_grasswd <- train(CCgrswd ~ .,
                    data = iteration1_data_grswd$train,
                    method = "rf",
                    trControl = repeat_cv,
                    metric = "Rsquared"
)

########### TEST #####################
preds_grswd <- predict(rf_grasswd, newdata = iteration1_data_grswd$test)
metrics(
  preds = preds_grswd, actual = iteration1_data_grswd$test$CCgrswd,
  iteration1_data_grswd$test
)

################################################################################
######################### Variable importance ##################################
# #LOOK AT THE VARIABLE IMPORTANCE
var_imp_cc <- varImp(rf_grass, scale=FALSE)$importance
var_imp_cc <- data.frame(variables=row.names(var_imp_cc), importance=var_imp_cc$Overall)
################################################################################


########################PREPARE / LOAD DATA ###################################
agg_traits |> filter(site == 'Millville' &
                       year == '2023'&
                       date %in% c('05-03-2023','05-11-2023',
                                   '06-05-2023', '06-21-2023', '07-06-2023',
                                   '07-26-2023','10-30-2023')) |>
  select(-c(fid, year, year2, LAI, GLA, BGI, SCI, NGRDI)) -> CC_data23

agg_traits |> filter(site == 'Millville' &
                       year == '2024'&
                       date %in% c('06-05-2024','06-25-2024','07-11-2024')) |>
  select(-c(fid, year, year2, LAI)) -> CC_data24

#SELECT COLUMNS WANT TO INCLUDE, changed depending on 23 or 24
CC_data24 %>%
  select(-c(GLA, BGI, SCI, NGRDI)) -> CC_data24

# #SCALE DATA
CC_data23[,4:34] <- scale(CC_data23[,4:34])
CC_data24[,4:34] <- scale(CC_data24[4:34])


#######################FILTER VARIABLES#######################################

CC_data23 %>%
  select(1:3, GRVI, HI, CVI, NGBDI, VOLUME, SI, blue, CCgrs) -> CC_data23

# CC_data24 %>%
#    select(-c(CCgrs)) -> CC_data24

CC_data24 %>%
  select(1:3,   "CIG", "GNDVI", "NGBDI", "CCgrs") -> CC_data24

#########################RANDOM FOREST MODEL##################################
# 24
it2_data_grs24 <- train_test_split(CC_data24[,4:7], prop = 0.7)


#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)


rf_grass24 <- train(CCgrs~., data = it2_data_grs24$train,
                    method='rf',
                    trControl = repeat_cv,
                    metric = 'Rsquared')

preds_grs_24 <- predict(rf_grass24, newdata = it2_data_grs24$test)
metrics(preds = preds_grs_24, actual = it2_data_grs24$test$CCgrs,
        it2_data_grs24$test)

#23
it2_data_grs23 <- train_test_split(CC_data23[,4:11], prop = 0.7)


#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)


rf_grass23 <- train(CCgrs~., data = it2_data_grs23$train,
                    method='rf',
                    trControl = repeat_cv,
                    metric = 'Rsquared')

preds_grs_23 <- predict(rf_grass23, newdata = it2_data_grs23$test)
metrics(preds = preds_grs_23, actual = it2_data_grs23$test$CCgrs,
        it2_data_grs23$test)

#############################SECTION FOR GRASS + WEEDS#########################

####################################### 2024 ##################################
# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
agg_traits |> filter(site == 'Millville' &
                       year == '2024'&
                       date %in% c('06-05-2024','06-25-2024','07-11-2024')) |>
  select(-c(fid, year, year2, LAI)) -> CC_data24

#SELECT COLUMNS WANT TO INCLUDE, changed depending on 23 or 24
CC_data24 %>%
  select(-c(GLA, BGI, SCI, NGRDI)) -> CC_data24

# #SCALE DATA
CC_data24[,4:34] <- scale(CC_data24[4:34])

CC_data24 %>%
  select(1:3, blue, CVI, BGI2, CIG, GNDVI, "CCgrswd") -> CC_data24

it2_data_grswd24 <- train_test_split(CC_data24[,4:9], prop = 0.7,
                                     seed = 43)

#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)


rf_grasswd24 <- train(CCgrswd ~., data = it2_data_grswd24$train,
                      method='rf',
                      trControl = repeat_cv,
                      metric = 'Rsquared')

preds_grswd24 <- predict(rf_grasswd24, newdata = it2_data_grswd24$test)
metrics(preds = preds_grswd24, actual = it2_data_grswd24$test$CCgrswd,
        it2_data_grswd24$test)

####################################### 2023 ##################################
agg_traits |> filter(site == 'Millville' &
                       year == '2023'&
                       date %in% c('05-03-2023','05-11-2023',
                                   '06-05-2023', '06-21-2023', '07-06-2023',
                                   '07-26-2023','10-30-2023')) |>
  select(-c(fid, year, year2, LAI, GLA, BGI, SCI, NGRDI)) -> CC_data23

CC_data23[,4:34] <- scale(CC_data23[,4:34])

CC_data23 %>%
  select(1:3,CVI, GR, BGI2, VOLUME, SI, CCgrswd) -> CC_data23

it2_data_grswd23 <- train_test_split(CC_data23[,4:9], prop = 0.7,
                                     seed = 43)
################################## Grass & Weed #################################
#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)


rf_grasswd23 <- train(CCgrswd ~., data = it2_data_grswd23$train,
                      method='rf',
                      trControl = repeat_cv,
                      metric = 'Rsquared')

preds_grswd23 <- predict(rf_grasswd23, newdata = it2_data_grswd23$test)
metrics(preds = preds_grswd23, actual = it2_data_grswd23$test$CCgrswd,
        it2_data_grswd23$test)

################################################################################
agg_traits |>
  filter(site == "Millville" & year != "2024") |>
  arrange(year, fid, aggplot_id) |>
  group_by(aggplot_id, fid, year) |>
  summarise(CCgrs = mean(CCgrs), CCgrswd = mean(CCgrswd)) |>
  rename(
    prev_CCgrs = CCgrs,
    prev_CCgrswd = CCgrswd
  ) |>
  ungroup() |>
  mutate(year = case_when(
    year == "2023" ~ "2024",
    year == "2022" ~ "2023",
    TRUE ~ "0" # Handle other years
  )) -> prev_year


agg_traits |>
  filter(site == "Millville" & year != "2022") |>
  arrange(year, fid, aggplot_id) |>
  inner_join(prev_year, relationship = "many-to-one") -> CC3_data24

agg_traits |>
  filter(site == "Millville" & year != "2022") |>
  arrange(year, fid, aggplot_id) |>
  inner_join(prev_year, relationship = "many-to-one") |>
  filter(year == "2023" &
           date %in% c(
             "05-03-2023", "05-11-2023",
             "06-05-2023", "06-21-2023", "07-06-2023",
             "07-26-2023", "10-30-2023"
           )) |>
  select(-c(fid, year, year2, LAI, GLA, BGI, SCI, NGRDI)) -> CC3_data23

CC3_data24 |>
  filter(site == "Millville" &
           year == "2024" &
           date %in% c("06-05-2024", "06-25-2024", "07-11-2024")) |>
  select(-c(fid, year, year2, LAI)) -> CC3_data24


CC3_data24 %>%
  select(-c(GLA, BGI, SCI, NGRDI)) -> CC3_data24

# #SCALE DATA
CC3_data24[, 4:34] <- scale(CC3_data24[, 4:34])
CC3_data23[, 4:34] <- scale(CC3_data23[, 4:34])

####################### FILTER VARIABLES#######################################

CC3_data23 %>%
  select(1:3, "GR", "prev_CCgrs", "CCgrs") -> CC3_data23_grs

CC3_data23 %>%
  select(1:3, "HI", "prev_CCgrswd", "CCgrswd") -> CC3_data23_grswd

CC3_data24 %>%
  select(
    1:3, "GNDVI", "CIG", green, EVI,
    blue, "NGBDI", BGI2, prev_CCgrs, "CCgrs"
  ) -> CC3_data24_grs

CC3_data24 %>%
  select(
    1:3, CVI, VOLUME, "GNDVI", blue, "CIG",
    EVI, NGBDI, prev_CCgrswd, "CCgrswd"
  ) -> CC3_data24_grswd


###################################### 2024 ####################################

############################# RANDOM FOREST######################################


it3_data_grs24 <- train_test_split(CC3_data24_grs[, 4:12],
                                   prop = 0.7,
                                   seed = 43
)
# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = "oob", number = 10)


rf3_grass24 <- train(CCgrs ~ .,
                     data = it3_data_grs24$train,
                     method = "rf",
                     trControl = repeat_cv,
                     metric = "Rsquared"
)

preds_3_grs <- predict(rf3_grass24, newdata = it3_data_grs24$test)
metrics(
  preds = preds_3_grs, actual = it3_data_grs24$test$CCgrs,
  it3_data_grs24$test
)

############################### Grass w Weeds ##################################



# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
it3_data_grswd24 <- train_test_split(CC3_data24_grswd[, 4:12],
                                     prop = 0.7,
                                     seed = 43
)

# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = "oob", number = 10)


rf3_grasswd24 <- train(CCgrswd ~ .,
                       data = it3_data_grswd24$train,
                       method = "rf",
                       trControl = repeat_cv,
                       metric = "Rsquared"
)

preds3_grswd <- predict(rf3_grasswd24, newdata = it3_data_grswd24$test)
metrics(
  preds = preds3_grswd, actual = it3_data_grswd24$test$CCgrswd,
  it3_data_grswd24$test
)

#################################### 2023 #####################################
############################# RANDOM FOREST######################################


it3_data_grs23 <- train_test_split(CC3_data23_grs[, 4:6],
                                   prop = 0.7,
                                   seed = 43
)
# #SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# #CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = "oob", number = 10)
#
#
rf3_grass23 <- train(CCgrs ~ .,
                     data = it3_data_grs23$train,
                     method = "rf",
                     trControl = repeat_cv,
                     metric = "Rsquared"
)

preds3_grs_23 <- predict(rf3_grass23, newdata = it3_data_grs23$test)
metrics(
  preds = preds3_grs_23, actual = it3_data_grs23$test$CCgrs,
  it3_data_grs23$test
)

################################################################################
################################## Grass w Weed ###############################

# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
it3_data_grswd23 <- train_test_split(CC3_data23_grswd[, 4:6],
                                     prop = 0.7,
                                     seed = 43
)

# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = "oob", number = 10)


rf3_grasswd23 <- train(CCgrswd ~ .,
                       data = it3_data_grswd23$train,
                       method = "rf",
                       trControl = repeat_cv,
                       metric = "Rsquared"
)

preds3_grswd_23 <- predict(rf3_grasswd23, newdata = it3_data_grswd23$test)
metrics(
  preds = preds3_grswd_23, actual = it3_data_grswd23$test$CCgrswd,
  it3_data_grswd23$test
)
################################################################################
####################################### LAI ####################################
agg_traits |>
  filter((site == 'Millville' & year == '2022' & date %in% c('09-26-2022',
                                                             '07-08-2022')) |
           (site == 'BlueCreek' & year == '2024' & date %in% c('07-31-2024',
                                                               '07-18-2024')))|>
  select(-c(fid, year, year2, CCgrs, CCgrswd, GLA, BGI, SCI, NGRDI)) -> LAI_data


#SCALE DATA
LAI_data[,4:34] <- scale(LAI_data[,4:34])

LAI_data <- LAI_data |> select(CVI, VOLUME, SI, PSRI, HUE, LAI)

################################# MODEL ########################################

# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
iteration1_data_lai <- train_test_split(LAI_data, prop = 0.7, seed = 43)

#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)


rf_lai_global <- train(LAI~., data = iteration1_data_lai$train,
                       method='rf',
                       trControl = repeat_cv,
                       metric = 'Rsquared')

preds_LAI1 <- predict(rf_lai_global, newdata = iteration1_data_lai$test)
metrics(preds = preds_LAI1, actual = iteration1_data_lai$test$LAI,
        iteration1_data_lai$test)

var_imp_lai <- varImp(rf_lai_global, scale=FALSE)$importance
var_imp_lai <- data.frame(variables=row.names(var_imp_lai), importance=var_imp_lai$Overall)
##############################################################################

agg_traits |>
  filter(site == 'Millville' & year == '2023'&
           date %in% c('05-03-2023','05-11-2023','06-05-2023', '06-21-2023',
                       '07-06-2023', '07-26-2023','10-30-2023')) |>
  select(-c(fid, year, year2, CCgrs, CCgrswd )) -> LAI_data23

LAI_data23 %>%
  select(-c(GLA, BGI, SCI, NGRDI)) -> LAI_data23

LAI_data23[,4:34] <- scale(LAI_data23[,4:34])
#######################FILTER VARIABLES#######################################

LAI_data23 %>%
  select(1:3, rededge, nir, green, CVI, BGI2, NGBDI, VOLUME, SI,
         LAI) -> LAI_data23

#########################RANDOM FOREST MODEL 2023##################################
# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
it2_23_lai_data <- train_test_split(LAI_data23[, 4:12],
                                    prop = 0.7, seed = 43)
#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)


rf_lai23 <- train(LAI~., data = it2_23_lai_data$train,
                  method='rf',
                  trControl = repeat_cv,
                  metric = 'Rsquared')
###################PREDICT TEST###########################################
preds_lai23 <- predict(rf_lai23, newdata = it2_23_lai_data$test)
metrics(preds_lai23, it2_23_lai_data$test$LAI,
        it2_23_lai_data$test)
##########################################################################
agg_traits |>
  filter(site == "Millville" & year == "2024" & date %in% c(
    "06-05-2024",
    "06-25-2024", "07-11-2024"
  )) |>
  select(-c(fid, year, year2, CCgrs, CCgrswd)) -> LAI_data24

# SELECT COLUMNS WANT TO INCLUDE, changed depending on 23 or 24
LAI_data24 %>%
  select(-c(GLA, BGI, SCI, NGRDI)) -> LAI_data24

# #SCALE DATA
LAI_data24[, 4:34] <- scale(LAI_data24[, 4:34])

####################### FILTER VARIABLES#######################################

LAI_data24 %>%
  select(1:3, GNDVI, CIG, CVI, BGI2, NGBDI, VOLUME,"LAI") -> LAI_data24
######################### RANDOM FOREST MODEL 2024 ##############################
##############################  2024 ###########################################
# SPLIT THE DATA INTO TRAIN AND TEST DATA SETS
iteration2_lai_data <- train_test_split(LAI_data24[, 4:10], prop = 0.7, seed = 43)

# SET SEED TO GET SIMILAR RESULTS
set.seed(24)
# CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = "oob", number = 10)


rf_lai2 <- train(LAI ~ .,
                 data = iteration2_lai_data$train,
                 method = "rf",
                 trControl = repeat_cv,
                 metric = "Rsquared"
)

################### PREDICT TEST###########################################
preds <- predict(rf_lai2, newdata = iteration2_lai_data$test)
metrics(
  preds = preds, actual = iteration2_lai_data$test$LAI,
  iteration2_lai_data$test
)
################################################################################
agg_traits |> filter(site == 'Millville' & year != '2024')|>
  arrange(year, fid, aggplot_id) |>
  group_by(aggplot_id, fid, year) |>
  summarise(LAI = mean(LAI)) |>
  rename(prev_LAI = LAI) |>
  ungroup()|>
  mutate(year = case_when(
    year == '2023' ~ '2024',
    year == '2022' ~ '2023',
    TRUE ~ '0'  # Handle other years
  )) -> prev_year


agg_traits |> filter(site == 'Millville' & year != '2022') |>
  arrange(year, fid, aggplot_id) |>
  inner_join(prev_year, relationship = 'many-to-one') |>
  filter(year == '2023'& date %in% c('05-03-2023','05-11-2023',
                                     '06-05-2023', '06-21-2023', '07-06-2023',
                                     '07-26-2023','10-30-2023')) -> LAI3_data23

LAI3_data23 %>%
  select(-c(GLA, BGI, SCI, NGRDI, CCgrs, CCgrswd)) -> LAI3_data23

LAI3_data23[,7:37] <- scale(LAI3_data23[,7:37])

#######################FILTER VARIABLES#######################################
LAI3_data23 %>%
  select(1:6, EVI, HUE, SI, prev_LAI, LAI) -> LAI3_data23

#############################RANDOM FOREST######################################

it3_lai23 <- train_test_split(LAI3_data23[,7:11], prop = 0.7)
#it3_data_lai24 <- train_test_split(LAI3_data24[ ,7:14], prop = 0.7)
#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)

rf3_LAI23 <- train(LAI~., data = it3_lai23$train,
                   method='rf',
                   trControl = repeat_cv,
                   metric = 'Rsquared')

#########################PREDICT TEST#########################################
preds_LAI_3_23 <- predict(rf3_LAI23, newdata = it3_lai23$test)
metrics(preds = preds_LAI_3_23, actual = it3_lai23$test$LAI,
        it3_lai23$test)

########################PREPARE DATA##########################################
agg_traits |> filter(site == 'Millville' & year != '2024')|>
  arrange(year, fid, aggplot_id) |>
  group_by(aggplot_id, fid, year) |>
  summarise(LAI = mean(LAI)) |>
  rename(prev_LAI = LAI) |>
  ungroup()|>
  mutate(year = case_when(
    year == '2023' ~ '2024',
    year == '2022' ~ '2023',
    TRUE ~ '0'  # Handle other years
  )) -> prev_year


agg_traits |> filter(site == 'Millville' & year != '2022') |> arrange(year, fid, aggplot_id) |>
  inner_join(prev_year, relationship = 'many-to-one') -> LAI3_data24

LAI3_data24 |>
  filter(site == 'Millville' & year == '2024'&
           date %in% c('06-05-2024','06-25-2024','07-11-2024')) -> LAI3_data24

LAI3_data24 %>%
  select(-c(GLA, BGI, SCI, NGRDI)) -> LAI3_data24

LAI3_data24 %>%
  select(-c(CCgrs, CCgrswd)) -> LAI3_data24
# #SCALE DATA
LAI3_data24[,7:37] <- scale(LAI3_data24[,7:37])

#######################FILTER VARIABLES#######################################
LAI3_data24 %>%
  select(1:6,  "VOLUME", "CVI", "GNDVI",
         "CIG", "NGBDI", "BGI2", "prev_LAI", "LAI") -> LAI3_data24

#############################RANDOM FOREST######################################

#it3_data_lai23 <- train_test_split(LAI3_data23[,7:11], prop = 0.7)
it3_data_lai24 <- train_test_split(LAI3_data24[ ,7:14], prop = 0.7)
#SET SEED TO GET SIMILAR RESULTS
set.seed(24)
#CREATE A FORMAT TO CROSS_VALIDATE THE MODEL
repeat_cv <- trainControl(method = 'oob', number = 10)

rf3_LAI24 <- train(LAI~., data = it3_data_lai24$train,
                   method='rf',
                   trControl = repeat_cv,
                   metric = 'Rsquared')

#########################PREDICT TEST#########################################
preds_LAI_3_24 <- predict(rf3_LAI24, newdata = it3_data_lai24$test)
metrics(preds = preds_LAI_3_24, actual = it3_data_lai24$test$LAI,
        it3_data_lai24$test)

#############################################################################
my_colors <- new_color_scheme(c("#006d2c","#238b45",
                                "#41ab5d", "#a1d99b", "#e5f5e0" ),
                              name = "my_color_scheme")

var_imp_cc <- var_imp_cc |>
  dplyr::arrange(importance) |>
  dplyr::mutate(variables = factor(variables, levels = unique(variables)))


var_imp_cc |>
  tidyplot(x = importance, y = variables, color = variables ) |>
  add_mean_bar() |>
  adjust_colors(RColorBrewer::brewer.pal(3, "Greens")) |>
  adjust_x_axis_title("Mean Decrease in Accuracy")|>
  adjust_y_axis_title("Spectral Index") |>
  adjust_font(10) |>
  adjust_size(NA, NA) |>
  remove_legend()|>
  add_title("Canopy Cover") -> var_cc_plot

var_imp_lai <- var_imp_lai |>
  dplyr::arrange(importance) |>
  dplyr::mutate(variables = factor(variables, levels = unique(variables)))

var_imp_lai |>
  tidyplot(x = importance, y = variables, color = variables ) |>
  add_mean_bar() |>
  adjust_colors(RColorBrewer::brewer.pal(5, "Greens")) |>
  adjust_x_axis_title("Mean Decrease in Accuracy")|>
  adjust_y_axis_title("Spectral Index") |>
  adjust_font(10)|>
  remove_legend()|>
  adjust_size(NA, NA) |>
  add_title("Leaf Area Index") -> var_lai_plot

var_imp_bio <- var_imp_bio |>
  dplyr::arrange(importance) |>
  dplyr::mutate(variables = factor(variables, levels = unique(variables)))

var_imp_bio |>
  tidyplot(x = importance, y = variables, color = variables ) |>
  add_mean_bar() |>
  adjust_x_axis(rotate_labels = TRUE) |>
  adjust_x_axis_title("Mean Decrease in Accuracy")|>
  adjust_y_axis_title("Spectral Index") |>
  adjust_colors(RColorBrewer::brewer.pal(7, "Greens")) |>
  adjust_font(10) |>
  adjust_size(NA, NA) |>
  remove_legend() |>
  add_title("Biomass") -> var_bio_plot


design <- "
  123
"

var_cc_plot + var_lai_plot + var_bio_plot +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

###############################################################################

round_1_global <- create_scatterplots(preds_LAI1, iteration1_data_lai$test$LAI,
                                      iteration1_data_lai$test, 1, "Global Model")

global_grs <- create_scatterplots(preds_grs, iteration1_data_grs$test$CCgrs,
                                  iteration1_data_grs$test, 1, "Global Model")
global_grswd <- create_scatterplots(preds_grswd, iteration1_data_grswd$test$CCgrswd,
                                    iteration1_data_grswd$test, 1, "Global Model")

round_1 <- create_scatterplots(
  preds_dry_bio_global, it1_dry$test$drybiomass_adj,
  it1_dry$test$drybiomass_adj, 1, "Global Model"
)


################################BIO MASS #######################################
bio_graph_data <- bind_cols(preds_dry_bio_global, it1_dry$test$drybiomass_adj) |>
  rename(
    preds = "...1",
    actual = "...2"
  ) |>
  mutate(
    residuals = actual - preds
  )


bio_metrics <- metrics(preds_dry_bio_global, it1_dry$test$drybiomass_adj,
                       it1_dry$test)

bio_graph_data |>
  tidyplot(x = actual, y = preds) |>
  add_data_points(color = "#238b45") |>
  add_curve_fit(se = TRUE, color = "forestgreen", fill= "#bae4b3") |>
  add_annotation_text(
    str_c(
      "R^2: ",
      formatC(
        bio_metrics["R_squared"],
        digits = 3,
        format = "f"
      )
    ),
    x = max(bio_graph_data$actual) - 0.5 * sd(bio_graph_data$actual),
    y = min(bio_graph_data$preds),
    fontsize = 10
  ) |>
  add_annotation_text(
    "CI: 95%",
    x = max(bio_graph_data$actual) - 0.5 * sd(bio_graph_data$actual),
    y = min(bio_graph_data$preds) - 0.5 * sd(bio_graph_data$actual),
    fontsize = 10
  )|>
  adjust_size(NA, NA) |>
  adjust_font(fontsize = 12) |>
  adjust_x_axis_title(title = "Observed", fontsize = 10) |>
  adjust_y_axis_title(title = "Predicted", fontsize = 10) |>
  add_title("Global Model: Biomass") -> global_bio
######################residuals plot
bio_graph_data |>
  tidyplot(x = actual, y = residuals) |>
  add_data_points(color = "#238b45") |>
  add_curve_fit(se = FALSE, color = "forestgreen") |>
  adjust_size(NA, NA) |>
  adjust_font(fontsize = 12) |>
  adjust_x_axis_title(title = "Observed", fontsize = 10) |>
  adjust_y_axis_title(title = "Residuals (Observed - Predicted)", fontsize = 10) |>
  add_title("Biomass")

################VARIABLE IMPORTANCE

var_imp_bio |>
  dplyr::arrange(importance) -> test1

test1$variables <- factor(test1$variables, levels = test1$variables)

test1 |>
  tidyplot(x = importance, y = variables, color = variables) |>
  add_mean_bar() |>
  adjust_colors(new_colors = brewer.pal(7, "Greens")) |>
  adjust_size(NA, NA) |>
  adjust_font(fontsize = 12) |>
  adjust_x_axis_title(title = "Mean Decrease in Accuracy", fontsize = 10) |>
  adjust_x_axis(limits = c(0, 20000000))|>
  adjust_y_axis_title(title = "Spectral index", fontsize = 10) |>
  add_title("Biomass") |> remove_legend() |>  add(
    ggplot2::theme(
      axis.title.x = element_text(size = 10),  # Adjust x-axis title font size
      axis.title.y = element_text(size = 10),  # Adjust y-axis title font size
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1)    # Adjust x-axis text font size
    )
  )



############################### CANOPY COVER ###################################
CC_graph_data <- bind_cols(preds_grs, iteration1_data_grs$test$CCgrs) |>
  rename(
    preds = "...1",
    actual = "...2"
  ) |>
  mutate(
    residuals = actual - preds
  )

CC_metrics <- metrics(preds_grs, iteration1_data_grs$test$CCgrs,
                      iteration1_data_grs$test)

CC_graph_data |>
  tidyplot(x = actual, y = preds) |>
  add_data_points(color = "#238b45") |>
  add_curve_fit(se = TRUE, color = "forestgreen", fill= "#bae4b3") |>
  add_annotation_text(
    str_c(
      "R^2: ",
      formatC(
        CC_metrics["R_squared"],
        digits = 3,
        format = "f"
      )
    ),
    x = max(CC_graph_data$actual) - 0.5 * sd(CC_graph_data$actual),
    y = min(CC_graph_data$preds),
    fontsize = 10
  ) |>
  add_annotation_text(
    "CI: 95%",
    x = max(CC_graph_data$actual) - 0.5 * sd(CC_graph_data$actual),
    y = 0,
    fontsize = 10
  )|>
  adjust_size(NA, NA) |>
  adjust_font(fontsize = 12) |>
  adjust_x_axis_title(title = "Observed", fontsize = 10) |>
  adjust_y_axis_title(title = "Predicted", fontsize = 10) |>
  add_title("Global Model: Canopy Cover") -> global_cc

######################################## RESIDUALS #############################
CC_graph_data |>
  tidyplot(x = actual, y = residuals) |>
  add_data_points(color = "#238b45") |>
  add_curve_fit(se = FALSE, color = "forestgreen") |>
  adjust_size(NA, NA) |>
  adjust_font(fontsize = 12) |>
  adjust_x_axis_title(title = "Observed", fontsize = 8) |>
  adjust_y_axis_title(title = "Residuals (Observed - Predicted)", fontsize = 8) |>
  add_title("Canopy Cover (grass)")

####################################Varaible importance ########################

var_imp_cc |>
  dplyr::arrange(importance) -> test_cc

test_cc$variables <- factor(test_cc$variables, levels = test_cc$variables)

test_cc|>
  tidyplot(x = importance, y = variables, color = variables) |>
  add_mean_bar() |>
  adjust_colors(new_colors = c( "#a1d99b",  "#74c476","#41ab5d","#238b45", "#005a32")) |>
  adjust_size(NA, NA) |>
  adjust_font(fontsize = 12) |>
  adjust_x_axis_title(title = "Mean Decrease in Accuracy", fontsize = 8) |>
  adjust_y_axis_title(title = "Spectral index", fontsize = 8) |>
  add_title("Canopy Cover")


################################LAI ############################################
LAI_graph_data <- bind_cols(preds_LAI1, iteration1_data_lai$test$LAI) |>
  rename(
    preds = "...1",
    actual = "...2"
  ) |>
  mutate(
    residuals = actual - preds
  )

LAI_metrics <- metrics(preds_LAI1, iteration1_data_lai$test$LAI, iteration1_data_lai$test)


LAI_graph_data |>
  tidyplot(x = actual, y = preds) |>
  add_data_points(color = "#238b45") |>
  add_curve_fit(se = TRUE, color = "forestgreen", fill= "#bae4b3") |>
  add_annotation_text(
    str_c(
      "R^2: ",
      formatC(
        LAI_metrics["R_squared"],
        digits = 3,
        format = "f"
      )
    ),
    x = max(LAI_graph_data$actual) - 0.5 * sd(LAI_graph_data$actual),
    y = min(LAI_graph_data$preds),
    fontsize = 8
  ) |>
  add_annotation_text(
    "CI: 95%",
    x = max(LAI_graph_data$actual) - 0.5 * sd(LAI_graph_data$actual),
    y = 0,
    fontsize = 8
  )|>
  adjust_size(NA, NA) |>
  adjust_font(fontsize = 12) |>
  adjust_x_axis_title(title = "Observed", fontsize = 8) |>
  adjust_y_axis_title(title = "Predicted", fontsize = 8) |>
  add_title("Global Model: LAI") -> global_lai
######################################## RESIDUALS #############################
LAI_graph_data |>
  tidyplot(x = actual, y = residuals) |>
  add_data_points(color = "#238b45") |>
  add_curve_fit(se = FALSE, color = "forestgreen") |>
  adjust_size(NA, NA) |>
  adjust_font(fontsize = 12) |>
  adjust_x_axis_title(title = "Observed", fontsize = 8) |>
  adjust_y_axis_title(title = "Residuals (Observed - Predicted)", fontsize = 8) |>
  add_title("LAI")

####################################Varaible importance ########################

var_imp_lai |>
  dplyr::arrange(importance) -> test_lai

test_lai$variables <- factor(test_lai$variables, levels = test_lai$variables)

test_lai |>
  tidyplot(x = importance, y = variables, color = variables) |>
  add_mean_bar() |>
  adjust_colors(new_colors = c( "#a1d99b",  "#74c476","#41ab5d","#238b45", "#005a32")) |>
  adjust_size(NA, NA) |>
  adjust_font(fontsize = 12) |>
  adjust_x_axis_title(title = "Mean Decrease in Accuracy", fontsize = 8) |>
  adjust_y_axis_title(title = "Spectral index", fontsize = 8) |>
  add_title("LAI")


##### CREATE THE TABLE

summary(LAI_graph_data$residuals)
summary(CC_graph_data$residuals)
names <- c("Canopy Cover", "Leaf Area Index", "Biomass(kg/ha)")
mean_residual <- c(mean(CC_graph_data$residuals), mean(LAI_graph_data$residuals),
                   mean(bio_graph_data$residuals))
sd_residual <- c(sd(CC_graph_data$residuals), sd(LAI_graph_data$residuals),
                 sd(bio_graph_data$residuals))
max_residual <- c(max(CC_graph_data$residuals), max(LAI_graph_data$residuals),
                  max(bio_graph_data$residuals))
min_residual <- c(c(min(CC_graph_data$residuals), min(LAI_graph_data$residuals),
                    min(bio_graph_data$residuals)))
residual_range <- c(max(CC_graph_data$residuals) - min(CC_graph_data$residuals),
                    max(LAI_graph_data$residuals) - min(LAI_graph_data$residuals),
                    max(bio_graph_data$residuals) - min(bio_graph_data$residuals))

residual_table <- bind_cols(names, mean_residual,
                            sd_residual, max_residual,
                            min_residual, residual_range) |>
  rename(
    Traits = "...1",
    `Mean Residual` = "...2",
    `SD Residual` = "...3",
    `Max Residual` = "...4",
    `Min Residual` = "...5",
    `Residual Range` = "...6"
  )

################################################################################


design <- "123"


global_cc + global_lai + global_bio +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")


################################################################################
# Create a scatterplot for the Frequency Percentages
freq_1 <- create_scatterplots(
  preds = preds_freq1,
  actual = iteration1_data$test$freq_perc,
  iteration1_data$test, 1, "Global Plot"
) |> adjust_y_axis(title = "Round 1 \n Predictions")

# Create a scatterplot for the secound round of frequency percentages
freq_2 <- create_scatterplots(
  preds = preds_freq1_2,
  actual = iteration1_2_data$test$freq_perc,
  iteration1_2_data$test, 2,
  ""
) |> adjust_y_axis(title = "Round 2 \n Predictions")
# Create a scatterplot of the third round of freq percentages
freq_3 <- create_scatterplots(
  preds = preds_freq1_3,
  actual = iteration1_3_data$test$freq_perc,
  iteration1_3_data$test, 2, ""
) |> adjust_y_axis(title = "Round 3 \n Predictions")
# Create a scatterplot for the Frequency Percentages
freq_4 <- create_scatterplots(
  preds = preds_freq2,
  actual = iteration2_data$test$freq_perc,
  iteration2_data$test, 1, "Global Plot 1st Flight Only"
)

# Create a scatterplot for the secound round of frequency percentages
freq_5 <- create_scatterplots(
  preds = preds_freq2_2,
  actual = iteration2_2data$test$freq_perc,
  iteration2_2data$test, 2,
  ""
)
# Create a scatterplot of the third round of freq percentages
freq_6 <- create_scatterplots(
  preds = preds_freq2_3,
  actual = iteration2_3data$test$freq_perc,
  iteration2_3data$test, 2, ""
)

design <- "
  14
  25
  36
"


#USING PATCHWORK PLACE THE GRAPHS NEXT TO EACH OTHER
freq_1 + freq_2 + freq_3 + freq_4 + freq_5 + freq_6 +
  plot_layout(design = design) +
  plot_annotation(title = "Seedling Frequency")
###############################################################################


stopCluster(cl)
foreach::registerDoSEQ()
