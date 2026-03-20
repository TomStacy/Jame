required_packages <- c("readxl", "dplyr", "plotly", "htmlwidgets")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(plotly)
  library(htmlwidgets)
})

workbook_path <- "fungal_mycorrhizal_tropical new_filtered_clean.xlsx"

if (!file.exists(workbook_path)) {
  stop(sprintf("Expected workbook at %s", normalizePath(workbook_path, winslash = "/", mustWork = FALSE)))
}

df <- read_excel(workbook_path)

required_columns <- c(
  "Sample",
  "Gen_sp",
  "MycorrhizalType",
  "RootWaterContent_pct",
  "RootVolume_cm3"
)

missing_columns <- setdiff(required_columns, names(df))
if (length(missing_columns) > 0) {
  stop(sprintf(
    "Missing expected columns: %s",
    paste(missing_columns, collapse = ", ")
  ))
}

plot_df <- df |>
  select(all_of(required_columns)) |>
  mutate(
    Sample = suppressWarnings(as.numeric(Sample)),
    RootWaterContent_pct = suppressWarnings(as.numeric(RootWaterContent_pct)),
    RootVolume_cm3 = suppressWarnings(as.numeric(RootVolume_cm3)),
    Gen_sp = if_else(is.na(Gen_sp), "NA", as.character(Gen_sp)),
    MycorrhizalType = if_else(is.na(MycorrhizalType), "NA", as.character(MycorrhizalType))
  ) |>
  filter(
    !is.na(Sample),
    !is.na(RootWaterContent_pct),
    !is.na(RootVolume_cm3)
  )

outlier_columns <- c("RootWaterContent_pct", "RootVolume_cm3")
outlier_mask <- rep(FALSE, nrow(plot_df))
outlier_bounds <- vector("list", length(outlier_columns))
names(outlier_bounds) <- outlier_columns

for (column in outlier_columns) {
  q1 <- quantile(plot_df[[column]], probs = 0.25, na.rm = TRUE)
  q3 <- quantile(plot_df[[column]], probs = 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr

  outlier_bounds[[column]] <- c(lower = lower_bound, upper = upper_bound)
  outlier_mask <- outlier_mask |
    plot_df[[column]] < lower_bound |
    plot_df[[column]] > upper_bound
}

outlier_df <- plot_df |>
  filter(outlier_mask) |>
  arrange(MycorrhizalType, Sample)

plot_df <- plot_df |>
  filter(!outlier_mask) |>
  arrange(MycorrhizalType, RootVolume_cm3) |>
  mutate(
    RootWaterContent_zscore = (
      RootWaterContent_pct - mean(RootWaterContent_pct, na.rm = TRUE)
    ) / sd(RootWaterContent_pct, na.rm = TRUE)
  )

group_levels <- sort(unique(plot_df$MycorrhizalType))

pvalue_records <- lapply(group_levels, function(mycorrhizal_type) {
  group_values <- plot_df |>
    filter(MycorrhizalType == mycorrhizal_type) |>
    pull(RootWaterContent_pct)

  other_values <- plot_df |>
    filter(MycorrhizalType != mycorrhizal_type) |>
    pull(RootWaterContent_pct)

  p_value <- NA_real_
  if (length(group_values) >= 2 && length(other_values) >= 2) {
    p_value <- t.test(group_values, other_values, var.equal = FALSE)$p.value
  }

  data.frame(
    MycorrhizalType = mycorrhizal_type,
    n = length(group_values),
    p_value = p_value,
    stringsAsFactors = FALSE
  )
})

pvalue_df <- bind_rows(pvalue_records) |>
  arrange(MycorrhizalType)

format_pvalue <- function(value) {
  if (is.na(value)) {
    return("NA")
  }

  format(value, digits = 4)
}

pvalue_annotation <- paste(
  vapply(seq_len(nrow(pvalue_df)), function(index) {
    sprintf(
      "%s: p=%s (n=%d)",
      pvalue_df$MycorrhizalType[[index]],
      format_pvalue(pvalue_df$p_value[[index]]),
      pvalue_df$n[[index]]
    )
  }, character(1)),
  collapse = "<br>"
)

cat(sprintf(
  "Removed %d outlier rows using the 1.5xIQR rule across RootWaterContent_pct and RootVolume_cm3.\n",
  nrow(outlier_df)
))

print(
  pvalue_df |>
    mutate(p_value = vapply(p_value, format_pvalue, character(1)))
)

if (nrow(outlier_df) > 0) {
  print(
    outlier_df |>
      select(
        Sample,
        Gen_sp,
        MycorrhizalType,
        RootWaterContent_pct,
        RootVolume_cm3
      )
  )
}

sample_hover <- with(
  plot_df,
  paste0(
    "Sample: ", Sample,
    "<br>Gen_sp: ", Gen_sp,
    "<br>Mycorrhizal Type: ", MycorrhizalType,
    "<br>Root Volume (cm^3): ", sprintf("%.4f", RootVolume_cm3),
    "<br>Root Water Content (%): ", RootWaterContent_pct,
    "<br>Root Water Content z-score: ", sprintf("%.2f", RootWaterContent_zscore)
  )
)

sample_figure <- plot_ly(
  data = plot_df,
  x = ~RootVolume_cm3,
  y = ~RootWaterContent_pct,
  color = ~MycorrhizalType,
  type = "scatter",
  mode = "lines+markers",
  text = sample_hover,
  hoverinfo = "text"
) |>
  layout(
    title = list(text = "Root Water Content by Root Volume"),
    template = "plotly_dark",
    hovermode = "x unified",
    xaxis = list(title = "Root Volume (cm^3)"),
    yaxis = list(title = "Root Water Content (%)"),
    legend = list(
      title = list(text = "Mycorrhizal Type"),
      x = 1.02,
      y = 1,
      xanchor = "left",
      yanchor = "top"
    ),
    margin = list(r = 220),
    annotations = list(list(
      xref = "paper",
      yref = "paper",
      x = 1.02,
      y = 0.02,
      xanchor = "left",
      yanchor = "bottom",
      showarrow = FALSE,
      align = "left",
      text = paste0(
        "Welch t-test p-values<br>vs remaining types<br><br>",
        pvalue_annotation
      ),
      bgcolor = "rgba(20, 20, 20, 0.9)",
      bordercolor = "rgba(255, 255, 255, 0.25)",
      borderwidth = 1,
      font = list(color = "#f5f5f5")
    ))
  )

box_hover <- with(
  plot_df,
  paste0(
    "Sample: ", Sample,
    "<br>Gen_sp: ", Gen_sp,
    "<br>Mycorrhizal Type: ", MycorrhizalType,
    "<br>Root Volume (cm^3): ", sprintf("%.4f", RootVolume_cm3),
    "<br>Root Water Content (%): ", RootWaterContent_pct
  )
)

box_figure <- plot_ly(
  data = plot_df,
  x = ~MycorrhizalType,
  y = ~RootWaterContent_pct,
  color = ~MycorrhizalType,
  type = "box",
  boxpoints = "all",
  jitter = 0.3,
  pointpos = 0,
  text = box_hover,
  hoverinfo = "text"
) |>
  layout(
    title = list(text = "Distribution of Root Water Content by Mycorrhizal Type"),
    template = "plotly_dark",
    showlegend = FALSE,
    xaxis = list(title = "Mycorrhizal Type"),
    yaxis = list(title = "Root Water Content (%)"),
    margin = list(r = 220),
    annotations = list(list(
      xref = "paper",
      yref = "paper",
      x = 1.02,
      y = 0.02,
      xanchor = "left",
      yanchor = "bottom",
      showarrow = FALSE,
      align = "left",
      text = paste0(
        "Welch t-test p-values<br>vs remaining types<br><br>",
        pvalue_annotation
      ),
      bgcolor = "rgba(20, 20, 20, 0.9)",
      bordercolor = "rgba(255, 255, 255, 0.25)",
      borderwidth = 1,
      font = list(color = "#f5f5f5")
    ))
  )

if (interactive()) {
  print(sample_figure)
  print(box_figure)
}

output_html <- "root_water_content_plot.html"
box_output_html <- "root_water_content_boxplot.html"
output_png <- "root_water_content_plot.png"

saveWidget(sample_figure, output_html, selfcontained = TRUE)
cat(sprintf("Saved interactive chart to %s\n", normalizePath(output_html, winslash = "/", mustWork = FALSE)))

saveWidget(box_figure, box_output_html, selfcontained = TRUE)
cat(sprintf("Saved box plot to %s\n", normalizePath(box_output_html, winslash = "/", mustWork = FALSE)))

tryCatch(
  {
    plotly::save_image(sample_figure, output_png, width = 1200, height = 700, scale = 2)
    cat(sprintf("Saved static chart to %s\n", normalizePath(output_png, winslash = "/", mustWork = FALSE)))
  },
  error = function(error) {
    cat("Static export skipped. Install or update kaleido if needed.\n")
    message(error$message)
  }
)