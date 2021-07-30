---
  title: "Experiment Data Exploration"
author: "Jesús Vélez Santiago"
date: "`r format(Sys.Date(), '%Y-%m')`"
output: 
  html_document:
  theme: readable
highlight: kate
toc: true
toc_float: true
toc_depth: 3
code_folding: show
self_contained: true
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  fig.align = "center",
  # dev = "svg",
  fig.retina = 2
)
```

## Libraries

```{r libraries, message=FALSE}
library(tidyverse)
library(ggpubr)
library(ggwaffle)
library(here)
library(glue)
```

## Load Data

```{r load_data}
lineages_files <- here("data", "processed", "lineages.tsv")
lineages_df <- read_tsv(lineages_files, show_col_types = FALSE) %>% 
  glimpse()
```

## Minimal preprocessing

```{r minimal_preprocessing}
processed_lineages_df <- lineages_df %>% 
  mutate(
    gfp = log10(gfp),
    ds_red = log10(ds_red),
    across(contains("filamentaded_"),~factor(.x, c(FALSE, TRUE), c("Not filamentaded", "Filamentaded")))
  ) %>% 
  add_count(experiment_id, trap_id, track_id, time) %>% 
  filter(n == 1) %>% 
  select(-n) %>% 
  glimpse()
```


## Exploratory Data Analysis

### Set default plot style

```{r default_plot_theme}
theme_set(theme_bw())
```


```{r donout_chart}
donout_df <- processed_lineages_df %>% 
  count(experiment_id, filamentaded_track) %>% 
  arrange(filamentaded_track) %>% 
  mutate(
    experiment_id = case_when(
      experiment_id == "Chromosome" ~ "C",
      TRUE ~ "P"
    ),
    percentage = n / sum(n) * 100,
    ymax = cumsum(percentage),
    ymin = c(0, head(ymax, -1)),
    label = glue("{experiment_id}: {format(percentage, digits=2)}%"),
    label_position = (ymax + ymin) / 2
  ) %>% 
  glimpse()

donout_total <- donout_df %>% 
  pull(n) %>% 
  sum()

donout_df %>% 
  ggplot(
    aes(
      ymax=ymax,
      ymin=ymin,
      xmax=4,
      xmin=3
    ),
  ) +
  geom_rect(
    size = 1.5,
    color = "white",
    aes(
      fill=filamentaded_track,
      group=experiment_id
    )
  ) +
  geom_label(x = 2, aes(y = label_position, label = label), size=3.5) +
  coord_polar(theta = "y") +
  xlim(c(-1, 4)) +
  labs(
    fill = "Cell status",
    caption = glue("Total: {format(donout_total, big.mark=',')} tracks")
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    plot.caption = element_text(face = "bold", hjust = 0.5)
  )
```


### Mean GFP
```{r gfp_distribution}
processed_lineages_df %>% 
  group_by(experiment_id, trap_id, track_id, filamentaded_track) %>% 
  summarize(mean_gfp = mean(gfp), .groups = "drop") %>%
  gghistogram(
    x = "mean_gfp",
    facet.by = "experiment_id",
    color = "filamentaded_track",
    fill = "filamentaded_track",
    alpha = 1/3,
    add = "mean",
    xlab = "Mean fluorescent intensity (log10)",
    ylab = "Count of cells"
  ) +
  labs(
    color = "Cell status",
    fill = "Cell status"
  )
```

```{r,area_chart}
processed_lineages_df %>% 
  count(experiment_id, filamentaded_at_frame, time) %>% 
  group_by(experiment_id, time) %>% 
  summarize(
    filamentaded_at_frame = filamentaded_at_frame,
    percentage = n / sum(n),
    .groups = "drop"
  ) %>% 
  ggplot(aes(x = time, y = percentage, fill = filamentaded_at_frame)) +
  geom_area(size = 0.5, alpha = 1/1) +
  geom_vline(xintercept = c(60, 140), linetype = "dashed") +
  geom_text(
    data = data.frame(
      x = 73,
      y = 0.8,
      label = "Start",
      experiment_id = "Plasmid"
    ),
    mapping = aes(x = x, y = y, label = label),
    size = 5.64,
    colour = "white",
    fontface = 2,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = data.frame(
      x = 151,
      y = 0.8,
      label = "End",
      experiment_id = "Plasmid"
    ),
    mapping = aes(x = x, y = y, label = label),
    size = 5.64,
    colour = "white",
    fontface = 2,
    inherit.aes = FALSE
  ) +
  facet_grid(experiment_id ~ .) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.spacing.y = unit(1, "lines")
  ) +
  labs(
    x = "Time (minutes)",
    y = "Percentage of cells",
    fill = "Cell status"
  )
```

```{r rain_plot}
library(ggdist)

p_1 <- processed_lineages_df %>% 
  filter(time == 0) %>% 
  ggplot(
    aes(
      x = filamentaded_track,
      y = length,
      group = filamentaded_track,
      fill = filamentaded_track,
      color = filamentaded_track
    )
  ) +
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2,
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15,
    outlier.shape = NA,
    alpha = 1/3
  ) +
  ggdist::geom_dots(side = "bottom", alpha = 1/10) +
  facet_wrap(experiment_id ~ .) +
  theme_bw() +
  theme(legend.position = "top") +
  labs(
    x = "Cell status",
    y = "Initial length",
    color = "Cell status",
    fill = "Cell status"
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  coord_flip() +
  stat_compare_means(label.y = 60, label.x = 1.5)

p_2 <- processed_lineages_df %>% 
  filter(time == 0) %>% 
  ggplot(
    aes(
      x = filamentaded_track,
      y = gfp,
      group = filamentaded_track,
      fill = filamentaded_track,
      color = filamentaded_track
    ),
    side = "bottom"
  ) +
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2,
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15,
    outlier.shape = NA,
    alpha = 1/3
  ) +
  ggdist::geom_dots(side = "bottom", alpha = 1/10) +
  facet_wrap(experiment_id ~ .) +
  theme_bw() +
  theme(legend.position = "top") +
  labs(
    x = "Cell status",
    y = "Initial GFP",
    color = "Cell status",
    fill = "Cell status"
  ) +
  coord_flip() +
  stat_compare_means(label.y = 2.5, label.x = 1.5)
```

```{r rain_cloud_2}
library(patchwork)
library(ggpubr)

(p_1 / p_2) +
  plot_layout(guides = 'collect')
```


```{r metric_charts}
processed_lineages_df %>% 
  select(experiment_id, filamentaded_track, time, length, gfp, ds_red) %>% 
  pivot_longer(
    cols = c(length, gfp, ds_red),
    names_to = "metric"
  ) %>%
  mutate(
    metric = case_when(
      metric == "ds_red" ~ "DsRed",
      metric == "gfp" ~ "GFP",
      metric == "length" ~ "Length"
    )
  ) %>% 
  group_by(experiment_id, filamentaded_track, time, metric) %>%
  summarise(
    ci = list(mean_cl_normal(value)),
    .groups = "drop"
  ) %>% 
  unnest(cols = c(ci)) %>% 
  ggplot(aes(x = time, y = y, ymin= ymin, ymax=ymax, color = filamentaded_track)) +
  annotate("rect", xmin=60, xmax=140, ymin=-Inf, ymax=Inf, alpha=1/2, color = "transparent", fill = "#FCB565") +
  geom_smooth(method = "loess") +
  facet_grid(metric ~ experiment_id, scales = "free_y") +
  labs(
    x = "Time (minutes)",
    y = "Value",
    color = "Cell status"
  ) +
  theme_bw() +
  theme(legend.position = "top")
```
```{r survival_probability}
gfp_control_hist <- processed_lineages_df %>% 
  ggplot(aes(x = gfp)) +
  geom_histogram(bins = 100)

gfp_hist_data <- gfp_control_hist %>% 
  ggplot_build() %>% 
  pluck("data", 1) %>% 
  select(count, x, xmin, xmax) %>% 
  as_tibble()

gfp_breaks <- gfp_hist_data %>% 
  {c(.$xmin, last(.$xmax))}

survival_probability_df <- processed_lineages_df %>%
  group_by(experiment_id, lineage_id, trap_id, filamentaded_track) %>% 
  summarise(
    initial_gfp = first(gfp),
    is_long_track = first(centered_frame) < unique(centered_antibiotic_start_frame) &&
      last(centered_frame) > unique(centered_antibiotic_end_frame),
    .groups = "drop"
  ) %>% 
  filter(is_long_track) %>%
  group_by(experiment_id, filamentaded_track) %>% 
  group_modify(~{
    tibble(
      plot = list(
        ggplot(data = .x, aes(x = initial_gfp)) +
          geom_histogram(breaks = gfp_breaks)
      )
    )
  }) %>% 
  mutate(
    counts = map(plot, ggplot_build),
    counts = map(counts, pluck, "data", 1),
    counts = map(counts, add_column, control_count = gfp_hist_data$count),
    counts = map(counts, select, gfp = x, control_count, count)
  ) %>% 
  unnest(counts) %>% 
  mutate(
    survival_probability = count / control_count,
    #survival_probability = survival_probability / max(survival_probability, na.rm = TRUE)
  ) %>% 
  filter(survival_probability != 0) %>% 
  glimpse()
```


```{r survival_probability_plot}
survival_probability_df %>% 
  ggplot(aes(x = gfp, y = survival_probability, color = filamentaded_track, linetype = experiment_id)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(
    legend.position = "top"
  ) +
  labs(
    x = "Initial GFP",
    y = "Survival probability",
    color = "Cell status",
    linetype = "Experiment"
  )
```

```{r survival_probability_plot_area}
survival_probability_df %>% 
  group_by(filamentaded_track, gfp) %>% 
  ggplot(aes(x = gfp, y = count, fill = filamentaded_track)) +
  geom_area(position = "fill", stat="identity") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  theme_bw() +
  theme(
    legend.position = "top"
  ) +
  labs(
    x = "Initial GFP (log10)",
    y = "Percentage of cells",
    fill = "Cell status"
  )
```

```{r}
status_points_df <- processed_lineages_df %>% 
  group_by(experiment_id, trap_id, track_id, filamentaded_track) %>% 
  summarize(
    first_false = which.min(filamentaded_at_frame),
    first_true = which.max(filamentaded_at_frame),
    initial_gfp = first(gfp),
    sos_gfp = gfp[first_true],
    end_gfp = last(gfp),
    #diff_sos_initial_gfp = sos_gfp -initial_gfp,
    #diff_end_sos_gfp = end_gfp - sos_gfp,
    diff_end_intial_gfp = end_gfp - initial_gfp,
    initial_length = first(length),
    sos_length = length[first_true],
    end_length = last(length),
    #diff_sos_initial_length = sos_length - initial_length,
    #diff_end_sos_length = end_length - sos_length,
    diff_end_intial_length = end_length - initial_length,
    initial_time = first(centered_frame) * 10,
    sos_time = centered_frame[first_true] * 10,
    end_time = last(centered_frame) * 10,
    life_time = end_time - initial_time,
    #diff_sos_intial_time = sos_time - initial_time,
    #diff_end_sos_time = end_time - sos_time,
    #diff_end_intial_time = end_time - initial_time,
    is_survivor = initial_time < unique(centered_antibiotic_end_frame) * 10 &&
      end_time > unique(centered_antibiotic_end_frame) * 10,
    is_survivor = ifelse(is_survivor, "Survived", "Dit not survive"),
    is_survivor = factor(is_survivor),
    .groups = "drop"
  ) %>% 
  #filter(initial_frame <= sos_frame, sos_frame <= end_frame) %>% 
  glimpse()  
```

```{r}
status_points_df %>% 
  count(experiment_id, filamentaded_track, life_time) %>% 
  ggplot(aes(x = as.factor(life_time), y = n, fill = filamentaded_track)) +
  geom_bar(position = "fill", stat="identity", width = 1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  facet_grid(experiment_id ~ .) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.spacing.y = unit(1, "lines")
  ) +
  labs(
    x = "Cell life time",
    y = "Percentage of cells",
    fill = "Cell status"
  )
```

```{r}
status_points_df %>% 
  group_by(experiment_id, filamentaded_track, life_time) %>% 
  summarize(
    n = n(),
    initial_length = median(initial_length),
    sos_length = median(sos_length),
    end_length = median(end_length),
    .groups = "drop"
  ) %>% 
  pivot_longer(
    cols = contains("length"),
    names_to = "length_type"
  ) %>% 
  mutate(
    length_type = factor(length_type, levels = c("initial_length", "sos_length", "end_length"), labels = c("Initial", "SOS", "End"))
  ) %>% 
  ggplot(aes(x = life_time, y = value)) +
  geom_line(aes(group = life_time)) +
  geom_point(aes(color = length_type), alpha = 1/1) +
  facet_grid(filamentaded_track ~ experiment_id) +
  theme_bw() +
  theme(
    legend.position = "top"
  ) +
  labs(
    x = "Cell life time",
    y = "Length value",
    color = "Length type"
  )
```

```{r}
status_points_df %>% 
  group_by(experiment_id, filamentaded_track, life_time) %>% 
  summarize(
    n = n(),
    initial_gfp = median(initial_gfp),
    sos_gfp = median(sos_gfp),
    end_gfp = median(end_gfp),
    .groups = "drop"
  ) %>% 
  pivot_longer(
    cols = contains("gfp"),
    names_to = "gfp_type"
  ) %>% 
  mutate(
    gfp_type = factor(gfp_type, levels = c("initial_gfp", "sos_gfp", "end_gfp"), labels = c("Initial", "SOS", "End"))
  ) %>% 
  ggplot(aes(x = life_time, y = value)) +
  geom_line(aes(group = life_time)) +
  geom_point(aes(color = gfp_type), alpha = 1/1) +
  facet_grid(filamentaded_track ~ experiment_id) +
  theme_bw() +
  theme(
    legend.position = "top"
  ) +
  labs(
    x = "Cell life time",
    y = "GFP value",
    color = "GFP type"
  )
```

```{r}
library(tidymodels)

model_data <- status_points_df %>% 
  select(is_survivor, filamentaded_track, contains("gfp"), contains("length"), -contains("diff"), -contains("sos")) %>% 
  mutate(
    out = interaction(is_survivor, filamentaded_track),
    out = as.character(out),
    out = as.factor(out)
  ) %>% 
  select(-is_survivor, -filamentaded_track) %>% 
  glimpse()

count(model_data, out)
```


```{r}
set.seed(123)
model_data_split <- initial_split(model_data, prop = 0.75, strata = out)

training_data <- training(model_data_split)
testing_data <- testing(model_data_split)

model_split
```

```{r}
data_folds <- vfold_cv(training_data, v = 20, strata = out)
data_folds
```


```{r}
library(themis)

data_recipe <- recipe(out ~ ., data = training_data) %>% 
  #step_corr(all_numeric(), threshold = 0.8) %>%
  step_normalize(all_numeric(), -contains("time")) %>% 
  step_zv(all_predictors()) %>% 
  step_dummy(all_nominal(), -all_outcomes()) %>% 
  step_downsample(out)

summary(data_recipe)
```

```{r}
dt_tune_model <- decision_tree(
  mode = "classification",
  engine = "rpart",
  cost_complexity = tune(),
  tree_depth = tune(),
  min_n = tune()
)

dt_tune_model <- rand_forest(
  mode = "classification",
  engine = "ranger",
  mtry = tune(),
  trees = tune(),
  min_n = tune()
)

dt_tune_model
```

```{r}
set.seed(123)
dt_grid <- grid_random(
  mtry() %>% range_set(c(2, 3)),
  trees(),
  min_n(),
  size = 10
)

dt_grid
```

```{r}
data_wkfl <- workflow() %>% 
  add_model(dt_tune_model) %>% 
  add_recipe(data_recipe)

data_wkfl
```

```{r}
dt_tuning <- data_wkfl %>% 
  tune_grid(
    resamples = data_folds,
    grid = dt_grid
  )

dt_tuning %>% 
  show_best(metric = "roc_auc", n = 5)
```
```{r}
best_dt_model <- dt_tuning %>% 
  select_best(metric = "roc_auc")

best_dt_model
```

```{r}
final_data_wkfl <- data_wkfl %>% 
  finalize_workflow(best_dt_model)

final_data_wkfl
```

```{r}
data_wf_fit <- final_data_wkfl %>% 
  fit(data = training_data)

tree_fit <- data_wf_fit %>% 
  extract_fit_parsnip()
```

```{r}
vip::vip(tree_fit)
```

```{r}
data_final_fit <- final_data_wkfl %>% 
  last_fit(split = model_data_split)

data_final_fit %>% 
  collect_metrics()
```


```{r}
data_final_fit %>% 
  collect_predictions() %>% 
  roc_curve(truth = is_survivor, .estimate = .pred_Survived) %>% 
  identity() %>% 
  autoplot()
```

```{r}
tree_predictions <- data_final_fit %>% collect_predictions()

conf_mat(tree_predictions, truth = is_survivor, estimate = .pred_class) %>% 
  autoplot()
```


```{r}
status_points_df %>% 
  count(is_survivor) %>% 
  identity()
```

