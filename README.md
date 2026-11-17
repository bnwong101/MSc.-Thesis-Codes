# MSc.-Thesis-Codes

This repository contains the R scripts, statistical models, and spatial data workflows developed for my Master of Science (M.Sc.) thesis in Epidemiology and Biostatistics at Western University. 

The project investigated the association between short-term exposure to high ambient temperatures/heatwaves and emergency department (ED) visits for stroke across Ontario, Canada (2011-2021) from the National Ambulatory Care Reporting System (NACRS), linked to Daymet Version 4 temperature data using distributed lag non-linear models (DLNMs) to assess 6-day cumulative effects.

---

## 🔬 Study Overview & Methodology
* **Study Design:** Population-based time-stratified case-crossover design matching each case day to control days on the same day of the week within the same calendar month.
* **Study Population:** 110,880 stroke-related ED visits across Ontario during the warm season (May to September) from 2011 to 2021.
* **Exposure Assessment:** Daily meteorological variables (temperature and relative humidity) extracted at a 1 km $\times$ 1 km spatial resolution from the **Daymet Version 4** dataset, linked to patient Forward Sortation Areas (FSAs).
* **Statistical Modeling:** Two-stage analytical approach utilizing **Distributed Lag Non-Linear Models (DLNMs)** via the `dlnm` package to model non-linear exposure-response functions and cumulative lag effects (lags 0–6 days), combined via multivariate random-effects meta-regression (`mixmeta`) across 35 Ontario health regions.

---

## 📂 Repository Contents
* **`Sample.R`**: Core script and two-stage DLNM implementation for total stroke risk (Objective 1) and subtype-stratified analyses for ischemic, hemorrhagic, and unspecified strokes (Objective 2).
* **`Objective3_Brandonthesis.R`**: Stratified models and meta-regressions examining effect modification by age categories and sex, complete with Wald interaction tests.
* **`Objective4_BrandonThesis.R`**: Heatwave definition algorithms (evaluating $\ge 2$ consecutive days above the regional 95th temperature percentile) and region-specific conditional logistic regression scripts with pooled meta-analysis and forest plot visualization.
* **`Sensitivityanalysis_BrandonThesis.R`**: Robustness checks testing alternative modeling specifications, including substituting daily maximum temperature for mean temperature, omitting relative humidity, and redefining extreme heat thresholds at the 95th percentile.
* **`DataCleaning1.R`**: Data cleaning scripts, descriptive statistics generation, case-weighted environmental distribution summaries, and spatial choropleth map generation using `sf` and `ggplot2`.

---

## 🛠️ Key Technical Skills & Packages
* **Languages & Software:** R (v4.2.2)
* **Epidemiological & Statistical Methods:** Case-crossover designs, Distributed Lag Non-Linear Models (DLNM), multivariate random-effects meta-analysis, conditional logistic regression (`clogit`), natural cubic splines, log-knots positioning.
* **Key R Packages:** `dlnm`, `mixmeta`, `survival`, `splines`, `sf`, `dplyr`, `ggplot2`, `lubridate`.
* **Spatial & Environmental Data Integration:** High-resolution gridded meteorological processing (Daymet v4), FSA-to-Health Region correspondence mapping, and geospatial visualization.

---

## 📊 Summary of Principal Findings
* **Overall Stroke Risk:** Short-term exposure to high ambient temperatures (at the 99th percentile versus the minimum-risk temperature) showed no statistically significant association with overall stroke-related ED visits (Cumulative OR: 0.967, 95% CI: 0.847–1.104).
* **Subtype & Subgroup Consistency:** Null findings remained consistent across stroke subtypes (ischemic, hemorrhagic, unspecified) and age stratifications. While effect modification by sex was significant ($p = 0.048$), individual stratifications for males and females both yielded non-significant associations.
* **Heatwaves:** Regional heatwaves were associated with a small decrease in stroke ED visits (OR: 0.936, 95% CI: 0.906–0.967), potentially reflecting behavioral adaptations, heat warnings, or competing pre-hospital mortality risks.
* **Robustness:** Results remained stable across all sensitivity checks regarding alternative temperature metrics and humidity adjustments.

---

## 🔒 Data Privacy & Reproducibility Note

* **Data Governance & Privacy:** Due to data privacy regulations and institutional data-sharing agreements, the underlying administrative health datasets (National Ambulatory Care Reporting System [NACRS] and proprietary spatial linkage files (`stroke_pm25_temp_imputed.csv`, `Canada_FSA_to_HR18_Correspondence.xlsx`) cannot be shared publicly in this repository. All analyses were conducted on de-identified, population-level records in full compliance with provincial health privacy standards.
* **Computational Reproducibility:** To run these scripts locally on authorized datasets, update the working directory paths (`setwd(...)`) and ensure all required R packages (`dlnm`, `mixmeta`, `survival`, `splines`, `sf`, `dplyr`, `ggplot2`, `lubridate`, `readxl`, `timeDate`, `data.table`, `broom`, `viridis`, `lme4`, `boot`) are installed under R $\ge 4.2.2$[cite: 50]. Publicly accessible meteorological and boundary inputs can be obtained directly from NASA ORNL DAAC (**Daymet Version 4**) and **Statistics Canada** (2016 FSA Boundary Files & Census Population Counts).
* * **Thesis Publication:** The complete M.Sc. thesis detailing these findings and methodologies is permanently archived and publicly accessible on the Western University Institutional Repository: [Western Scholaris Repository Link](https://uwo.scholaris.ca/items/363c18c8-8f85-4f34-a0aa-335b63e3fec0).
