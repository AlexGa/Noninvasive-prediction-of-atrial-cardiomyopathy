Noninvasive prediction of atrial cardiomyopathy characterized by
multipolar high-density contact mapping
================

- [Import Data](#import-data)
- [General Stats](#general-stats)
  - [Gender](#gender)
  - [Median (IQR)](#median-iqr)
  - [Mean ± SD](#mean--sd)
  - [Paroxysmal AF bzw. Persistent AF](#paroxysmal-af-bzw-persistent-af)
- [Gaussian Mixture Model (GMM)](#gaussian-mixture-model-gmm)
  - [GMM - LVA 0.5 (cm²)](#gmm---lva-05-cm²)
  - [GMM - LVA 1.0 (cm²)](#gmm---lva-10-cm²)
  - [Summarizing GMM experiments (Two
    classes)](#summarizing-gmm-experiments-two-classes)
  - [Proof of concept](#proof-of-concept)
- [Echocardiography-dependent
  Classification](#echocardiography-dependent-classification)
  - [Feature-Selection (Boruta
    Algorithm)](#feature-selection-boruta-algorithm)
  - [Feature importance:](#feature-importance)
  - [SVM - Klassifikator](#svm---klassifikator)
  - [SVM - Klassifikator lineare
    Kernel](#svm---klassifikator-lineare-kernel)
  - [Paarweise Vergleich zwischen Ac und noAC basierend auf selektierten
    Parametern](#paarweise-vergleich-zwischen-ac-und-noac-basierend-auf-selektierten-parametern)
- [Analyse Voltageverteilungen](#analyse-voltageverteilungen)
  - [PCA: Area Ratio](#pca-area-ratio)
  - [PCA: Area\[mm^2\]](#pca-areamm2)
- [Validation Outcome after PVI](#validation-outcome-after-pvi)
  - [Recurrence AF ECG ≤ 12 months (0=no,
    1=yes)](#recurrence-af-ecg--12-months-0no-1yes)
  - [Recurrence AF ECG FU (0=no, 1=yes)](#recurrence-af-ecg-fu-0no-1yes)
  - [Recurrence AF symptoms +/- ECG FU (0=no,
    1=yes)](#recurrence-af-symptoms---ecg-fu-0no-1yes)
  - [Recurrence AF symptoms +/- ECG ≤ 12 months (0=no,
    1=yes)](#recurrence-af-symptoms---ecg--12-months-0no-1yes)
- [Tables](#tables)
  - [Table 1](#table-1)
  - [Table 1’ Median: All patients (missing
    entries)](#table-1-median-all-patients-missing-entries)
  - [Table 1’ Mean ± SD: All patients (missing
    entries)](#table-1-mean--sd-all-patients-missing-entries)
  - [Table 2: Comparison of selected echocardiographic parameters
    between groups mild AC and severe
    AC.](#table-2-comparison-of-selected-echocardiographic-parameters-between-groups-mild-ac-and-severe-ac)
  - [Table 3: Features evaluated with the Boruta algorithm. Display of
    all features evaluated with the Boruta algorithm as a function of
    feature
    importance.](#table-3-features-evaluated-with-the-boruta-algorithm-display-of-all-features-evaluated-with-the-boruta-algorithm-as-a-function-of-feature-importance)
  - [Table 4: Prediction of SVM for each CV
    run](#table-4-prediction-of-svm-for-each-cv-run)

Moritz T. Huttelmaier MD<sup>a</sup>; Alexander Gabel, PhD
<sup>b,c</sup>; Stefan Störk, MD PhD <sup>a,d</sup>; Stefan Frantz, MD
<sup>a,d</sup>; Caroline Morbach, MD <sup>a,d</sup>; Thomas H. Fischer,
MD <sup>a/sup\>

<sup>a</sup> Dept. of Internal Medicine I, University Hospital Würzburg
(UKW), Germany

<sup>b</sup> Institute of Medical Virology, Goethe-University Frankfurt,
60596 Frankfurt am Main, Germany

<sup>c</sup> Infection Control and Antimicrobial Stewardship Unit, UKW,
Würzburg, Germany

<sup>d</sup> Dept. Clinical Research & Epidemiology, Comprehensive Heart
Failure Centre Würzburg, University Hospital Würzburg, Germany

<br />

# Import Data

``` r
library(dplyr)
library(ggplot2)
library(dendextend)
library(Boruta)
library(Amelia)
library(caret)
library(yardstick)
library(plotly)
library(mclust)

source("R/helper_functions.R")

input_dir <- "~/Projekte/Huttelmaier_Atriale_Kardiomyopathie/"  
plot.dir <- "plots"

if(!file.exists(plot.dir)){
  dir.create(plot.dir)
}

data.dir <- "data"

if(!file.exists(data.dir)){
  dir.create(data.dir)
}

if(! dir.exists("data")){
  dir.create("data")
}

compl_df <- readr::read_delim("data/input_dataset.csv")

set.seed(1)
```

# General Stats

## Gender

``` r
compl_df %>% mutate(gender = if_else(`gender (0=female, 1=male)` == 0, true = "female", false = "male")) %>% 
             group_by(gender) %>% summarise(n = n(), .groups = "drop") %>% 
             mutate(comb =  paste0(n, " (", round(n/sum(n) * 100, 2),"%)")) %>% dplyr::select(-n) %>%
             tidyr::pivot_wider(names_from = gender, values_from = "comb") %>% 
  mutate(across(where(is.double), 
                ~ format(.x, digits = 2, scientific = T))) %>% 
  knitr::kable()
```

| female   | male     |
|:---------|:---------|
| 19 (38%) | 31 (62%) |

``` r
# 
# %>% 
#   DT::datatable(extensions = 'Buttons', options = list(
#     dom = 'Blfrtip',
#     buttons = c('copy', 'csv', 'excel', 'pdf'),
#     lengthMenu = list(c(10,30, 50, -1), 
#                       c('10', '30', '50', 'All')),
#     paging = F))
```

## Median (IQR)

``` r
compl_df %>% dplyr::summarise(Age = med_iqr(age_examination),
                       Height = med_iqr(`height (cm)`),
                       Weight = med_iqr(`weight (kg)`),
                       BMI = med_iqr(`BMI (kg/m²)`),
                       `NT_pro_BNP (pg/ml)` = med_iqr(`NT_pro_BNP (pg/ml)`),
                       `CHA2DS2-VASc` = med_iqr(`CHA2DS2-VASc`),
                        `AA_duration (months)` = med_iqr(`AA_duration (months)`),
                       `BSADuBois (m²)` = med_iqr(`BSADuBois (m²)`),
                       `LVEF biplan (%)` = med_iqr(`LVEF biplan (%)`),
                       `GLPS AFI average (%)` = med_iqr(`GLPS AFI average (%)`),
                       `Map Points` = med_iqr(map_points),
                       `LVA_0.5 (cm²)` = med_iqr(`LVA_0.5 (cm²)`),
                       `LVA_1.0 (cm²)` = med_iqr(`LVA_1.0 (cm²)`),
                       ) %>% 
              tidyr::pivot_longer(cols = matches("."), names_to = "Property", values_to = "Median (IQR)") %>% 
  mutate(across(where(is.double), 
                ~ format(.x, digits = 2, scientific = T))) %>% 
  knitr::kable() 
```

| Property             | Median (IQR)           |
|:---------------------|:-----------------------|
| Age                  | 64 (56 - 70)           |
| Height               | 174.5 (168.2 - 180.8)  |
| Weight               | 85 (75 - 99.8)         |
| BMI                  | 28 (24.6 - 30.8)       |
| NT_pro_BNP (pg/ml)   | 189.5 (70.8 - 391.2)   |
| CHA2DS2-VASc         | 2 (1 - 3)              |
| AA_duration (months) | 0 (0 - 5)              |
| BSADuBois (m²)       | 2 (1.9 - 2.2)          |
| LVEF biplan (%)      | 61 (56.2 - 63)         |
| GLPS AFI average (%) | -17.6 (-18.9 - -16.2)  |
| Map Points           | 5711 (5217.2 - 6987.5) |
| LVA_0.5 (cm²)        | 1.8 (0.5 - 5.9)        |
| LVA_1.0 (cm²)        | 7.8 (3.9 - 22.9)       |

``` r
  # DT::datatable(extensions = 'Buttons', options = list(
  #   dom = 'Blfrtip',
  #   buttons = c('copy', 'csv', 'excel', 'pdf'),
  #   lengthMenu = list(c(10,30, 50, -1), 
  #                     c('10', '30', '50', 'All')),
  #   paging = F))
```

## Mean ± SD

``` r
compl_df %>% summarise(Age = mean_sd(age_examination),
                       Height = mean_sd(`height (cm)`),
                       Weight = mean_sd(`weight (kg)`),
                       BMI = mean_sd(`BMI (kg/m²)`),
                       `NT_pro_BNP (pg/ml)` = mean_sd(`NT_pro_BNP (pg/ml)`),
                       `CHA2DS2-VASc` = mean_sd(`CHA2DS2-VASc`),
                        `AA_duration (months)` = mean_sd(`AA_duration (months)`),
                       `BSADuBois (m²)` = mean_sd(`BSADuBois (m²)`),
                       `LVEF biplan (%)` = mean_sd(`LVEF biplan (%)`),
                       `GLPS AFI average (%)` = mean_sd(`GLPS AFI average (%)`),
                       `Map Points` = mean_sd(map_points),
                       `LVA_0.5 (cm²)` = mean_sd(`LVA_0.5 (cm²)`),
                       `LVA_1.0 (cm²)` = mean_sd(`LVA_1.0 (cm²)`),
                       ) %>% 
              tidyr::pivot_longer(cols = matches("."), names_to = "Property", values_to = "Mean ± SD") %>% 
  DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1), 
                      c('10', '30', '50', 'All')),
    paging = F))
```

![](analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Paroxysmal AF bzw. Persistent AF

``` r
compl_df %>% mutate(AF = if_else(`AF (0=paroxysmal, 1=persistent)` == 0, true = "paroxysmal", false = "persistent")) %>% 
             group_by(AF) %>% summarise(n = n(), .groups = "drop") %>% 
             mutate(comb =  paste0(n, " (", round(n/sum(n) * 100, 2),"%)")) %>% dplyr::select(-n) %>%
             tidyr::pivot_wider(names_from = AF, values_from = "comb") %>% 
  DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1), 
                      c('10', '30', '50', 'All')),
    paging = F))
```

![](analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
compl_df %>% dplyr::select(`LVA_0.5 (cm²)`, `LVA_1.0 (cm²)`) %>% 
             dplyr::rename(`LVA 0.5` = `LVA_0.5 (cm²)`,
                           `LVA 1.0` = `LVA_1.0 (cm²)`) %>%
             tidyr::pivot_longer(cols = c(`LVA 0.5`, `LVA 1.0`),names_to = "area (cm²)", values_to = "LVA") %>%
  ggplot2::ggplot(ggplot2::aes(y = `area (cm²)`, x = LVA, fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(bandwidth = 3, jittered_points = TRUE,
                                         position = ggridges::position_points_jitter(width = 0.05, height = 0),
                                          point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
  ggplot2::theme_void() +
  ggplot2::xlab("") +
  ggplot2::ylab("") +
  scale_fill_viridis_c(option = "C", name ="[cm²]") + 
  ggplot2::scale_y_discrete(expand = c(0,0,0,1.25)) +
  ggplot2::scale_x_continuous(limits = c(0, 70)) +
  ggplot2::theme(legend.position = "right",
                 # legend.title = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(colour = "black", size = 12),
                 strip.background = ggplot2::element_rect(fill = "white"),
                 axis.title = ggplot2::element_text(face = "bold", size = 12),
                 strip.text = ggplot2::element_text(face = "bold", size = 12), panel.border = ggplot2::element_blank())
```

![](analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggplot2::ggsave(filename = file.path(plot.dir, "Density_plots_LVA.pdf"), width = 8, height = 4, device = "pdf", dpi = 300)
ggplot2::ggsave(filename = file.path(plot.dir, "Density_plots_LVA.svg"), width = 8, height = 4, device = "svg", dpi = 300)
```

# Gaussian Mixture Model (GMM)

## GMM - LVA 0.5 (cm²)

### GMM cluster parameter LVA 0.5 (cm²)

``` r
x_obs <- as.matrix(compl_df %>% pull(`LVA_0.5 (cm²)`))
set.seed(10)
gmm_obj05 <- max.em.gmm(x = x_obs, K = 2, epsilon = 1e-3, do.plot = F)

log_ll <- gmm_obj05$log.ll[length(gmm_obj05$log.ll)]

clusters <- data.frame(ID = compl_df %>% pull(ID), cluster_number = apply(gmm_obj05$gamma, 1, FUN = which.max))

orderAC <- order(gmm_obj05$mu[,1], decreasing = F)

ordered_clusters <- data.frame(cluster_order = orderAC, 
                               gmm05_1_2 = factor( c("mild", "severe"), levels =  c("mild", "severe")))


gmm_cluster_vec <- clusters %>% 
                   dplyr::inner_join(ordered_clusters, by = c("cluster_number" = "cluster_order"))

compl_df <- dplyr::inner_join(compl_df, gmm_cluster_vec, by = "ID")


gmm_obj05[[1]] <- as.matrix(gmm_obj05[[1]][orderAC,])
gmm_obj05[[2]] <- as.matrix(gmm_obj05[[2]][orderAC,])
gmm_obj05[[3]] <- gmm_obj05[[3]][orderAC]

rownames(gmm_obj05[[1]]) <- rownames(gmm_obj05[[2]]) <- names(gmm_obj05[[3]]) <- c("mild", "severe")

para_mat <- do.call("cbind", gmm_obj05[1:2])
if(ncol(para_mat) == 2){
      colnames(para_mat) <- c("mu", "var")
}else{
      colnames(para_mat) <- paste0(rep(c("mu.", "var."), each = 2), 1:2)
}

DT::datatable(apply(para_mat, 2, format, digits = 2), 
                extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1), 
                      c('10', '30', '50', 'All')),
    paging = F)) %>% 
      htmltools::tagList() 
```

<div class="datatables html-widget html-fill-item" id="htmlwidget-484ed42468d86de4a421" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-484ed42468d86de4a421">{"x":{"filter":"none","vertical":false,"extensions":["Buttons"],"data":[["mild","severe"],[" 0.96","14.80"],["  0.76","178.66"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>mu<\/th>\n      <th>var<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Blfrtip","buttons":["copy","csv","excel","pdf"],"lengthMenu":[[10,30,50,-1],["10","30","50","All"]],"paging":false,"columnDefs":[{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"mu","targets":1},{"name":"var","targets":2}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

### Shapiro-Wilk test for cluster LVA 0.5 (cm²)

``` r
compl_df %>% group_by(gmm05_1_2) %>% dplyr::rename(LVA_0.5 = `LVA_0.5 (cm²)`) %>% rstatix::shapiro_test(LVA_0.5)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["gmm05_1_2"],"name":[1],"type":["fct"],"align":["left"]},{"label":["variable"],"name":[2],"type":["chr"],"align":["left"]},{"label":["statistic"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["p"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"mild","2":"LVA_0.5","3":"0.8879104","4":"0.003091940"},{"1":"severe","2":"LVA_0.5","3":"0.8129042","4":"0.002324572"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

## GMM - LVA 1.0 (cm²)

``` r
x_obs <- as.matrix(compl_df %>% pull(`LVA_1.0 (cm²)`))
set.seed(10)
gmm_obj <- max.em.gmm(x = x_obs, K = 2, epsilon = 1e-3, do.plot = F)

log_ll <- gmm_obj$log.ll[length(gmm_obj$log.ll)]

clusters <- data.frame(ID = compl_df %>% pull(ID), cluster_number = apply(gmm_obj$gamma, 1, FUN = which.max))

orderAC <- order(gmm_obj$mu[,1], decreasing = F)

ordered_clusters <- data.frame(cluster_order = orderAC, 
                               gmm_1_2 = factor( c("mild", "severe"), levels =  c("mild", "severe")))


gmm_cluster_vec <- clusters %>% 
                   dplyr::inner_join(ordered_clusters, by = c("cluster_number" = "cluster_order"))

compl_df <- dplyr::inner_join(compl_df, gmm_cluster_vec, by = "ID")


gmm_obj[[1]] <- as.matrix(gmm_obj[[1]][orderAC,])
gmm_obj[[2]] <- as.matrix(gmm_obj[[2]][orderAC,])
gmm_obj[[3]] <- gmm_obj[[3]][orderAC]

rownames(gmm_obj[[1]]) <- rownames(gmm_obj[[2]]) <- names(gmm_obj[[3]]) <- c("mild", "severe")

x <- seq(from = 0, to = 80, length.out = 10000)
ll_mat <- matrix(nc = 2, nr = length(x))
for(k in 1:2){
  ll_mat[,k] <- dnorm(x, gmm_obj$mu[k], sd = sqrt(gmm_obj$sigma[k]))
}
lls <- rowSums(ll_mat)

org_points <- data.frame(lva = x_obs, logll = -0.008, AC = compl_df$gmm_1_2) 

gmm_plot <- data.frame(lva = x, logll = lls) %>%
  ggplot2::ggplot(ggplot2::aes(y = logll, x = lva, fill = lva)) +
  ggplot2::geom_bar(stat = "identity", width = 0.1) +
  ggplot2::geom_line(col ="black", size = 1) +
  ggplot2::geom_point(data = org_points, aes(shape = AC),col = "black", size = 4, alpha = 0.5) + 
  scale_shape_manual(values=c(21, 24)) +
  scale_color_viridis_d() +
  scale_fill_viridis_c(option = "C", name ="[cm²]") +
  ggplot2::xlab("") +
  ggplot2::ylab("") +
  ggplot2::scale_y_continuous(limits= c(-0.008, 0.17)) +
  ggplot2::scale_x_continuous(expand = c(0,0,0,0)) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "right",
                 axis.text = ggplot2::element_text(colour = "black", size = 12),
                 strip.background = ggplot2::element_rect(fill = "white"),
                 axis.title = ggplot2::element_text(face = "bold", size = 12),
                 strip.text = ggplot2::element_text(face = "bold", size = 12))

vis_arrow_df <- data.frame(lva = gmm_obj$mu) 
log_ll_tmp <- matrix(nc = 2, nr = 2)
for(k in 1:2){
  log_ll_tmp[,k] <- dnorm(vis_arrow_df$lva, gmm_obj$mu[k], sd = sqrt(gmm_obj$sigma[k]))
}
vis_arrow_df$logll <- rowSums(log_ll_tmp)
vis_arrow_df$sigma <- gmm_obj$sigma
vis_arrow_df$mu_label <- c("mu[mild]", "mu[severe]")

gmm_plot <- gmm_plot + ggplot2::geom_point(data = vis_arrow_df, pch = 3, size = 4) + 
           ggplot2::geom_errorbarh(data = vis_arrow_df, aes(xmax = lva + sqrt(sigma), xmin = lva - sqrt(sigma), height = 0)) +
           ggplot2::geom_segment(data = vis_arrow_df, aes(x = lva - sqrt(sigma), xend = lva + sqrt(sigma), y = logll, yend = logll),
                                 arrow = arrow(angle = 45, length = unit(0.25, "cm"), ends = "both")) +
           ggplot2::geom_text(data = vis_arrow_df, aes(label = mu_label), parse = T, nudge_y = 0.01)



ggplot2::ggsave(gmm_plot, filename = file.path(plot.dir, "Density_plots_LVA_GMM.pdf"), width = 8, height = 4, device = "pdf", dpi = 300)
ggplot2::ggsave(gmm_plot, filename = file.path(plot.dir, "Density_plots_LVA_GMM.svg"), width = 8, height = 4, device = "svg", dpi = 300)
```

### GMM cluster parameter 1.0 LVA

``` r
para_mat <- do.call("cbind", gmm_obj[1:2])
if(ncol(para_mat) == 2){
      colnames(para_mat) <- c("mu", "var")
}else{
      colnames(para_mat) <- paste0(rep(c("mu.", "var."), each = 2), 1:2)
}

DT::datatable(apply(para_mat, 2, format, digits = 2), 
                extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1), 
                      c('10', '30', '50', 'All')),
    paging = F)) %>% 
      htmltools::tagList() 
```

<div class="datatables html-widget html-fill-item" id="htmlwidget-697d58abce2053a24549" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-697d58abce2053a24549">{"x":{"filter":"none","vertical":false,"extensions":["Buttons"],"data":[["mild","severe"],[" 4.6","29.6"],["  7.1","333.3"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>mu<\/th>\n      <th>var<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Blfrtip","buttons":["copy","csv","excel","pdf"],"lengthMenu":[[10,30,50,-1],["10","30","50","All"]],"paging":false,"columnDefs":[{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"mu","targets":1},{"name":"var","targets":2}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

### Shapiro-Wilk test for cluster LVA 1.0 (cm²)

``` r
compl_df %>% group_by(gmm_1_2) %>% 
             dplyr::rename(LVA_1.0 = `LVA_1.0 (cm²)`) %>% 
             rstatix::shapiro_test(LVA_1.0) %>% 
              mutate(across(where(is.double), ~ round(.x, digits = 2)))
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["gmm_1_2"],"name":[1],"type":["fct"],"align":["left"]},{"label":["variable"],"name":[2],"type":["chr"],"align":["left"]},{"label":["statistic"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["p"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"mild","2":"LVA_1.0","3":"0.95","4":"0.19"},{"1":"severe","2":"LVA_1.0","3":"0.90","4":"0.03"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
gmm_plot
```

![](analysis_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

## Summarizing GMM experiments (Two classes)

### Numbers of patients in mild and severe AC

``` r
rbind(compl_df %>% 
        group_by(gmm05_1_2) %>% 
        count() %>% 
        tidyr::pivot_wider(names_from = gmm05_1_2, values_from = n) %>% 
        mutate(name = "LVA 0.5"),
      compl_df %>% 
        group_by(gmm_1_2) %>% 
        count() %>% 
        tidyr::pivot_wider(names_from = gmm_1_2, values_from = n) %>% 
        mutate(name = "LVA 1.0")) %>% 
  dplyr::relocate(name) %>% 
DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1), 
                      c('10', '30', '50', 'All')),
    paging = F)) %>% 
      htmltools::tagList() 
```

<div class="datatables html-widget html-fill-item" id="htmlwidget-29a56e939b13f92b7856" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-29a56e939b13f92b7856">{"x":{"filter":"none","vertical":false,"extensions":["Buttons"],"data":[["1","2"],["LVA 0.5","LVA 1.0"],[32,28],[18,22]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>name<\/th>\n      <th>mild<\/th>\n      <th>severe<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Blfrtip","buttons":["copy","csv","excel","pdf"],"lengthMenu":[[10,30,50,-1],["10","30","50","All"]],"paging":false,"columnDefs":[{"className":"dt-right","targets":[2,3]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"name","targets":1},{"name":"mild","targets":2},{"name":"severe","targets":3}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

## Proof of concept

**Wilcoxon-Test for “CHA2DS2-VASc” and “NT_pro_BNP (pg/ml)”**

``` r
comp_features <- c("CHA2DS2-VASc", "NT_pro_BNP (pg/ml)")

res_list <- sapply(comp_features, function(sel_feature){
  if(!all(is.na(compl_df %>% dplyr::pull(!!sym(sel_feature))))){
    
     return(wilcox_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>% 
              mutate(Feature = sel_feature))
  }else{
    return(NA)
  }
}, simplify = FALSE)

all_tests_cluster <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>% 
  relocate(Feature)

cbind(all_tests_cluster, "p.adjusted" = p.adjust(all_tests_cluster$`P-value`, method = "BH")) %>% as_tibble() %>%
  mutate(`P-adj. (short)` = dplyr::case_when(
                                                 `p.adjusted` >= 0.01 ~ as.character(round(`p.adjusted`, 2)),
                                                   `p.adjusted` < 0.0001 ~ "< 0.0001",
                                                   `p.adjusted` < 0.001 ~ "< 0.001",
                                                   `p.adjusted` < 0.01 ~ "< 0.01"),
         p.adjusted = format(p.adjusted, digits = 2),
         `P-value` = format(`P-value`, digits = 2)) %>% 
  dplyr::arrange(p.adjusted) %>% 
   DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1), 
                      c('10', '30', '50', 'All')),
    paging = F))
```

![](analysis_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

**Fisher Test for AF (0=paroxysmal, 1=persistent)**

``` r
fisher_test_table(df = compl_df, gr_row = "gmm_1_2", gr_col = "AF (0=paroxysmal, 1=persistent)", as_kable = FALSE) %>% 
   DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1), 
                      c('10', '30', '50', 'All')),
    paging = F))
```

    ## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
    ## ℹ Please use `all_of()` or `any_of()` instead.
    ##   # Was:
    ##   data %>% select(gr_col)
    ## 
    ##   # Now:
    ##   data %>% select(all_of(gr_col))
    ## 
    ## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning in data.frame(..., check.names = FALSE): row names were found from a
    ## short variable and have been discarded

![](analysis_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

# Echocardiography-dependent Classification

## Feature-Selection (Boruta Algorithm)

``` r
two_cl_df_set <- compl_df %>% dplyr::select(ID, age_examination, `gender (0=female, 1=male)`, `BMI (kg/m²)`, 
                                            `LAVI/a'2D ((ml/m²)/(m/s))`, `LAVI/a' 4D ((ml/m²)/(m/s))`, 
                                            `LAVImax 4D (ml/m²)`, `total LA-EF (%)`, `LA-EF total 4D (%)`,
                                            `LASr R av. (%)`,
                                            `LAScd R av. (%)`, `LASct R av. (%)`,
                                            `LASr 4D (%)`, `LAScd 4D (%)`, `LASct 4D (%)`,
                                            `LASr P av. (%)`,
                                            `LASr_c 4D (%)`, `LAScd_c 4D (%)`,  `LASct_c 4D (%)`,
                                            `LASr_c 4D (%)`,   `LAScd_c 4D (%)`,  `LASct_c 4D (%)`,
                                            `MV e' sept. (cm/s)`,
                                            `MV e' lat. (cm/s)`, `MV E/e' average`,
                                            gmm_1_2, `LVA_1.0 (cm²)`, `LVA_0.5 (cm²)`, `TR Vmax (m/s)`)

imputed_dataset <- c()

if(!file.exists(file.path(data.dir, "imputed_finale_dataset.rds"))){
  
  boruta_test <- two_cl_df_set %>% as_tibble() %>% 
               mutate(gender = factor(`gender (0=female, 1=male)`),
                      y = factor(gmm_1_2)) %>% 
               dplyr::select(-gmm_1_2, -ID, -`gender (0=female, 1=male)`, -`LVA_1.0 (cm²)`, -`LVA_0.5 (cm²)`)

boruta_test_compl <- Amelia::amelia(boruta_test, m=5, parallel = "multicore",
                            noms=c('gender','y'))

imputed_dataset <- data.frame(ID = two_cl_df_set$ID, boruta_test_compl$imputations$imp3, check.names = F)
saveRDS(object = imputed_dataset, file = file.path(data.dir, "imputed_finale_dataset.rds"))

}else{
  
  imputed_dataset <- readRDS(file.path(data.dir, "imputed_finale_dataset.rds"))
}
```

## Feature importance:

``` r
boruta.af_df <- c()
if(file.exists(file.path(data.dir, "boruta_af_dataset.rds"))){
  boruta.af_df <- readRDS(file.path(data.dir, "boruta_af_dataset.rds"))

}else{
  boruta.df_train <- Boruta(y~., data = imputed_dataset %>% dplyr::select(-ID), doTrace = 2)
  #take a call on tentative features
  boruta.af_df <- TentativeRoughFix(boruta.df_train)

  saveRDS(object = boruta.af_df, file = file.path(data.dir, "boruta_af_dataset.rds"))
}

df_importance <- reshape2::melt(as.data.frame(boruta.af_df$ImpHistory)) %>% dplyr::filter(is.finite(value))
  

df_imp_order <- df_importance %>% dplyr::group_by(variable) %>% 
                  dplyr::summarize(med_imp = median(value)) %>% 
                  ungroup() %>% mutate(med_rank = rank(med_imp))

df_importance <- df_importance %>%
                 left_join(data.frame(variable = names(boruta.af_df$finalDecision), 
                                       decision = boruta.af_df$finalDecision)) %>% 
                 mutate(decision = factor(if_else(is.na(decision), true = "Tentative", false = decision), levels = c("Tentative", "Rejected", "Confirmed")),
                        variable = factor(variable, levels = df_imp_order$variable[order(df_imp_order$med_rank)]))
  
df_importance %>% dplyr::filter(decision != "Tentative") %>% 
  ggplot2::ggplot(ggplot2::aes(x = variable, y = value, fill = decision)) + 
  ggplot2::stat_boxplot(geom = "errorbar", width = 0.5) + 
  ggplot2::geom_boxplot() +
  ggplot2::scale_fill_manual(values = c( "#F37651FF","#40B7ADFF")) +
  ggplot2::xlab("") + 
  ggplot2::ylab("Importance score") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom",
                 axis.text = ggplot2::element_text(colour = "black", size = 12),
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
                 strip.background = ggplot2::element_rect(fill = "white"),
                 axis.title = ggplot2::element_text(face = "bold", size = 12),
                 strip.text = ggplot2::element_text(face = "bold", size = 12))
```

![](analysis_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
ggplot2::ggsave(file.path(plot.dir, "Boruta_feature_importance.pdf"), width = 12, height = 7, device = "pdf", dpi = 300)
ggplot2::ggsave(file.path(plot.dir, "Boruta_feature_importance.svg"), width = 12, height = 7, device = "svg", dpi = 300)

boruta_attributes <- getSelectedAttributes(boruta.af_df, withTentative = F)

imputed_dataset <- imputed_dataset[c(boruta_attributes, "y")]
final_df <- imputed_dataset[c(boruta_attributes, "y")]

final_df <- final_df %>% dplyr::rename(age = "age_examination",
                                       LAVI_a_2D = "LAVI/a'2D ((ml/m²)/(m/s))",  
                                       LAVI_a_4D = "LAVI/a' 4D ((ml/m²)/(m/s))", 
                                       LASr_R_av = "LASr R av. (%)",             
                                       LASct_R_av = "LASct R av. (%)",            
                                       LASr_4D = "LASr 4D (%)",
                                       LASr_P_av = "LASr P av. (%)",
                                       LASct_c_4D = "LASct_c 4D (%)",             
                                       MV_e_sep = "MV e' sept. (cm/s)") %>% 
    mutate(y = factor(if_else(y == 1, true = "mild", false = "severe"), 
                      levels = c("mild", "severe")))
```

## SVM - Klassifikator

- zur Bestimmung der Vorhersagefähigkeit der ausgewählten Features
  (Echoparameter) wird ein SVM-Klassifikator trainiert.
- Training und Evaluation erfolgt mittels 5-facher Kreuzvalidierung
- Güte mittels ROC-Kurve (Receiver operating characteristic curve) und
  der Fläche unter der ROC (AUC) bestimmt

## SVM - Klassifikator lineare Kernel

``` r
# calculation of f1 score

f1 <- function (data, lev = NULL, model = NULL) {
  
  precision <- posPredValue(data$pred, data$obs, positive = "pass")
  recall  <- sensitivity(data$pred, data$obs, postive = "pass")
  f1_val <- (2 * precision * recall) / (precision + recall)
  names(f1_val) <- c("F1")
  
  return(f1_val)
}

plot_ROC_5cv <- function(roc_merged_df){
  
  ggplot2::ggplot(data = roc_merged_df) +
  ggplot2::geom_line(ggplot2::aes(x = fpr, 
                      y = tpr, 
                      colour = CV_fold,
                      linetype = CV_fold,
                      size = CV_fold,
                      alpha = CV_fold)) + 
  ggplot2::scale_alpha_manual(values = c(rep(0.8, 5), 1)) + 
  ggplot2::scale_size_manual(values = c(rep(1, 5), 1.5)) + 
  ggplot2::scale_linetype_manual(values = c(8, 6:3, 1)) +
  ggplot2::scale_color_manual(values = c("#DF8F44FF", "#00A1D5FF","#B24745FF", "#79AF97FF", "#6A6599FF", "#000000")) +
  ggplot2::scale_x_continuous(expand = expansion(add = c(0.01, 0.01))) + 
  ggplot2::scale_y_continuous(expand = expansion(add = c(0.01, 0.01))) + 
  ggplot2::geom_abline(slope = 1, linetype = "dashed") +
  ggplot2::xlab("False positive rate") +
  ggplot2::ylab("True positive rate") + 
  ggplot2::theme_bw() + 
  ggplot2::theme(legend.title = ggplot2::element_blank(),
                 legend.position = "bottom",
                 axis.text = ggplot2::element_text(colour = "black", size = 12),
                 strip.background = ggplot2::element_rect(fill = "white"),
                 axis.title = ggplot2::element_text(face = "bold", size = 12),
                 strip.text = ggplot2::element_text(face = "bold", size = 12)) +
  ggplot2::guides(col = ggplot2::guide_legend(nrow = 2))
}

# Training and validation of SVM model 
#-------------------------------------
# - 5-fold cross validation for model performance
# - 10-fold cross validation for parameter tuning within each cross validation

svm_cross_validation <- function(data_set, svm_method = "svmLinear", number_cv = 5){
  
  cv_index_mat <- matrix(nrow = number_cv, data = sample(seq_along(data_set[, 1])))
  
  roc_plt_list <- list()
  auc_list <- list()
  
  svmModels <- list()
  svmFitProbs_list <- list()
  svmFitClasses_list <- list()
  svmFitPredictions_list <- list()
  overall_perf <- c()
  
  approx_tpr <- list()
  mean_fpr <- seq(0, 1, length.out = 100)
  
  for (i in 1:nrow(cv_index_mat)) {
    
    train_df <- data_set[as.vector(cv_index_mat[-i, ]), ]
    test_df <- data_set[as.vector(cv_index_mat[i, ]), ]
    
    ctrl <- trainControl(method = "cv",
                         savePred = T,
                         classProb = T)
    
    svmFit <- c()
    
    if(svm_method == "svmLinear"){
      
      svmFit <- train(
      y ~ .,
      data = train_df,
      method = "svmLinear",
      trControl = ctrl,
      preProcess = c("center", "scale"),
      tuneGrid = expand.grid(C = seq(0, 2, length = 20)))
      
    }else if(svm_method == "svmRadial"){
      
      svmFit <- train(y ~ ., 
                      data = train_df,
                      method = "svmRadial", 
                      trControl = ctrl, 
                      preProcess = c("center","scale"), 
                      tuneLength = 10)
      
    }else if(svm_method == "svmPoly"){
      
      svmFit <- train(y ~ ., 
                      data = train_df, 
                      method = "svmPoly", 
                      trControl = ctrl, 
                      preProcess = c("center","scale"), 
                      tuneLength = 4)
  
    }else{
      
      stop("Please choose for svm_method between \"svmLinear\", \"svmRadial\", and \"svmPoly\"")
      
    }
    
    svmModels[[i]] <- svmFit
    
    svmFitProbs_list[[i]] <- predict(svmFit, test_df, type = "prob")
    svmFitClasses_list[[i]] <- predict(svmFit, test_df, type = "raw")
    svmFitPredictions_list[[i]] <- cbind(test_df$y, svmFitClasses_list[[i]], svmFitProbs_list[[i]])
    
    names(svmFitPredictions_list[[i]]) <- c("Observed", "Predicted", "probAC", "probnoAC")
    
    auc_list[[i]] <- svmFitPredictions_list[[i]] %>% roc_auc(truth = Observed, probAC)
    
    roc_plt_list[[i]] <- svmFitPredictions_list[[i]] %>%
      roc_curve(truth = Observed, probAC) %>%
      dplyr::mutate(CV_fold = paste0(i, ". run (auROC = ", format(auc_list[[i]]$.estimate, digits = 2), ")"),
                    fpr = 1 - specificity) %>%
      dplyr::rename(tpr = sensitivity) %>%
      dplyr::select(-specificity) %>% 
      dplyr::arrange(-.threshold)  %>% # rm threshold column for later merge with mean ROC coordinates
      dplyr::relocate(fpr, tpr, .threshold)
    
    overall_perf[i] <- caret::confusionMatrix(svmFitClasses_list[[i]], test_df$y)$overall["Accuracy"]
    
  }
  
  overall_auc <- do.call("rbind", auc_list) %>% 
                  summarize(.estimate = mean(.estimate)) %>% 
                  pull(.estimate)
  
  average_ROC <- genscore::avgrocs(lapply(roc_plt_list, function(tab) 
                                          tab %>% 
                                            dplyr::select(tpr, fpr) %>% 
                                            as.matrix()), 12, 8
                                   )
  
  # Add [0, 0] starting point to avg ROC
  average_ROC <- rbind(c(0,0), average_ROC)
  
  colnames(average_ROC) <- c("tpr", "fpr")
  roc_mean_df <- average_ROC %>% as_tibble() %>% 
                                 mutate(CV_fold = paste0("approx. mean ROC (mean auROC = ",
                                                          round(overall_auc, digits = 2),")"))
  
  rm(average_ROC)
  
  roc_single_df <- do.call("rbind", roc_plt_list) %>% dplyr::select(fpr, tpr, CV_fold)
  roc_merged_single_mean <- rbind(roc_single_df, roc_mean_df) %>% 
                            dplyr::mutate(CV_fold = as.factor(CV_fold))
  
  svmFitPredictions_df <- data.frame(do.call("rbind", svmFitPredictions_list),
                                     run = rep(
                                     seq_along(svmFitPredictions_list),
                                     each = nrow(svmFitPredictions_list[[1]])
                                     ))
  
  perf_plot <- plot_ROC_5cv(roc_merged_df = roc_merged_single_mean)
  
  return(
    list(roc_plt_list= roc_plt_list,
         auc_list = auc_list,
         svmFitPredictions_df = svmFitPredictions_df,
         svmFitProbs_list = svmFitProbs_list,
         svmFitClasses_list = svmFitClasses_list,
         svmFitPredictions_list = svmFitPredictions_list,
         overall_perf = overall_perf,
         approx_tpr = approx_tpr,
         mean_fpr = mean_fpr,
         roc_merged_single_mean = roc_merged_single_mean,
         perf_plot = perf_plot
    )
  )
}

svm_linear_result_list <- svm_cross_validation(final_df, svm_method = "svmLinear", number_cv = 5)
svm_linear_result_list$perf_plot 
```

![](analysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
ggplot2::ggsave(svm_linear_result_list$perf_plot, filename = file.path(plot.dir, "SVM_ROC_plot.pdf"), 
                dpi = 300, device = "pdf", height = 6, width = 8)

ggplot2::ggsave(svm_linear_result_list$perf_plot, filename = file.path(plot.dir, "SVM_ROC_plot.svg"), 
                dpi = 300, device = "svg", height = 6, width = 8)
```

### SVM Klassifikator radiale Kernel

``` r
svm_radial_result_list <- svm_cross_validation(final_df, svm_method = "svmRadial", number_cv = 5)
svm_radial_result_list$perf_plot
```

![](analysis_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
ggplot2::ggsave(svm_radial_result_list$perf_plot, filename = file.path(plot.dir, "SVM_radial_ROC_plot.pdf"), 
                dpi = 300, device = "pdf", height = 4, width = 5)
ggplot2::ggsave(svm_radial_result_list$perf_plot, filename = file.path(plot.dir, "SVM_radial_ROC_plot.svg"), 
                dpi = 300, device = "svg", height = 4, width = 5)
```

### SVM Klassifikator polynomiale Kernel

``` r
svm_poly_result_list <- svm_cross_validation(final_df, svm_method = "svmRadial", number_cv = 5)
svm_poly_result_list$perf_plot
```

![](analysis_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
ggplot2::ggsave(svm_poly_result_list$perf_plot, filename = file.path(plot.dir, "SVM_poly_ROC_plot.pdf"), 
                dpi = 300, device = "pdf", height = 4, width = 5)
ggplot2::ggsave(svm_poly_result_list$perf_plot, filename = file.path(plot.dir, "SVM_poly_ROC_plot.svg"), 
                dpi = 300, device = "svg", height = 4, width = 5)
```

## Paarweise Vergleich zwischen Ac und noAC basierend auf selektierten Parametern

``` r
comp_features <- boruta_attributes

res_list <- sapply(comp_features, function(sel_feature){
  if(!all(is.na(compl_df %>% dplyr::pull(!!sym(sel_feature))))){
    
     return(wilcox_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>% 
              mutate(Feature = sel_feature))
  }else{
    return(NA)
  }
}, simplify = FALSE)

all_tests_cluster <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>% 
  relocate(Feature)

all_tests_cluster_comp <- cbind(all_tests_cluster, "p.adjusted" = p.adjust(all_tests_cluster$`P-value`, method = "BH")) %>% as_tibble() %>%
  mutate(`P-adj. (short)` = dplyr::case_when(
                                                 `p.adjusted` >= 0.01 ~ as.character(round(`p.adjusted`, 2)),
                                                   `p.adjusted` < 0.0001 ~ "< 0.0001",
                                                   `p.adjusted` < 0.001 ~ "< 0.001",
                                                   `p.adjusted` < 0.01 ~ "< 0.01"),
         p.adjusted = format(p.adjusted, digits = 2),
         `P-value` = format(`P-value`, digits = 2))

all_tests_cluster_comp %>% 
  dplyr::arrange(p.adjusted) %>% 
   DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1), 
                      c('10', '30', '50', 'All')),
    paging = F))
```

![](analysis_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
plt_list <- list()

plot_df <- compl_df[c(boruta_attributes, "gmm_1_2")] %>%
            rename(group = gmm_1_2) %>%
          dplyr::rename(age = `age_examination`,
                        LAVI_a_2D = "LAVI/a'2D ((ml/m²)/(m/s))",  
                        LAVI_a_4D = "LAVI/a' 4D ((ml/m²)/(m/s))", 
                        LASr_R_av = "LASr R av. (%)",             
                        LASct_R_av = "LASct R av. (%)",            
                        LASr_4D = "LASr 4D (%)",
                        LASr_P_av = "LASr P av. (%)",
                        LASct_c_4D = "LASct_c 4D (%)",             
                        MV_e_sep = "MV e' sept. (cm/s)")

all_tests_cluster_comp <- all_tests_cluster_comp %>% dplyr::mutate(Feature = dplyr::if_else(Feature == "age_examination", true = "Age", false = Feature))

feature_print_df <- data.frame(print_name = c("Age", "LAVI/a'2D ((ml/m²)/(m/s))", "LAVI/a' 4D ((ml/m²)/(m/s))", "LASr R av. (%)", 
                                              "LASct R av. (%)", "LASr 4D (%)", "LASr P av. (%)", "LASct_c 4D (%)",  "MV e' sept. (cm/s)"),
                              feature = c("age", "LAVI_a_2D", "LAVI_a_4D", "LASr_R_av", "LASct_R_av", "LASr_4D", "LASr_P_av", "LASct_c_4D", "MV_e_sep"))

for(i in seq_along(feature_print_df$feature)){

  sel_feature <- feature_print_df$feature[i]
  print_name_feature <- feature_print_df$print_name[i]
  
  pval_anno <- all_tests_cluster_comp %>% filter(Feature == print_name_feature) %>% pull(`P-adj. (short)`)

  plt_list[[i]] <- plot_df %>%
  ggplot2::ggplot(ggplot2::aes_string(x = "group", y = sel_feature, fill = "group")) +
  ggplot2::stat_boxplot(geom = "errorbar", lwd = 1, width = 0.5, show.legend = FALSE) +
  ggplot2::geom_boxplot(color = "black", lwd = 1, show.legend = FALSE, outlier.size = 0) +
  ggbeeswarm::geom_quasirandom(alpha = 0.8, width = 0.4, pch = 20, size = 1) +
  ggplot2::theme_bw() +
  ggsci::scale_fill_jama() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 14, colour = "black"),
                 axis.title = ggplot2::element_text(face = "bold", size = 12),
                 legend.title = ggplot2::element_blank(),
                 legend.text = ggplot2::element_text(size = 12),
                 legend.position = "none") +
  ggplot2::xlab("") +
  ggplot2::ylab(print_name_feature) +
  ggsignif::geom_signif(comparisons = list(c("mild", "severe")),
                        annotation = pval_anno,
                        tip_length = c(0.2, 0.04)) + 
  ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, .1)))
  
  ggplot2::ggsave(plt_list[[i]], filename = file.path(plot.dir, paste0("Boxplot_comparison_MWU_", sel_feature, ".pdf")), 
                  dpi = 300, height = 5, width = 6, device = "pdf")
  ggplot2::ggsave(plt_list[[i]], filename = file.path(plot.dir, paste0("Boxplot_comparison_MWU_", sel_feature, ".svg")), 
                  dpi = 300, height = 5, width = 6, device = "svg")
}

cowplot::plot_grid(plotlist = plt_list, nrow = 3)
```

![](analysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
ggplot2::ggsave(filename = file.path(plot.dir, paste0("Grid_Boxplots_comparison_MWU_", sel_feature, ".pdf")), dpi = 300, height = 10, width = 15, device = "pdf")
ggplot2::ggsave(filename = file.path(plot.dir, paste0("Grid_Boxplots_comparison_MWU_", sel_feature, ".svg")), dpi = 300, height = 10, width = 15, device = "svg")
```

# Analyse Voltageverteilungen

<!-- Analyse Voltage-Verteilung (Datensätze 1-3) -->
<!-- In den Datensätzen 1-2 findest du die Voltage-Verteilung, die ich aus dem Software-Tool (CARTO-Net) exportiert habe. Hier sind die Pulmonalvenen und alle anderen Strukturen entfernt worden. In diesen Datensätzen fehlen die IDs 26 u. 46. Wie besprochen ist hier das Ziel zu prüfen, ob die Voltageverteilung (Mittelwert, Standardabweichung,…) in den beiden Gruppen (mild AC und severe AC) unterschiedlich ist bzw. die Gruppen so vorhergesagt werden können. -->
<!-- Der Datensetz 3 umfasst den Export aller Mapping-Daten aus CARTO-Net. Hier fehlen leider die IDs 3, 29 u. 45. Diese wurden leider bis eben nicht erfolgreich aus CARTO-Net exportiert. Dieser Datensatz hat den Vorteil, dass wir nicht an die maximal 20 Unterteilungen der Voltage-Verteilung gebunden sind, hat jedoch den Nachteil, dass alle Mapping-Punkte (auch die der Pulmonalvenen und des Ventrikels) als potentielle Fehl-Informationen mit eingehen. Bei der Analyse beschränken wir uns auf alle Punkte des 1. Maps (MapId=1, Spalte B). Relevant bzgl. der Voltage-Verteilung ist VoltageBipolar (Spalte F). Die Spalten T bis Y geben die Koordinaten im dreidimensionalen Raum an. -->

``` r
val_ds_volt_dir <- "~/Projekte/Huttelmaier_Voltage_Verteilungen/data/Histogramm_modified 21-12-23/"

hist_df <- data.frame()
hist_df_full_area <- hist_df_full_ratio <- list()

for(xl_file in dir(val_ds_volt_dir, full.names = T, pattern = "xlsx")){
  id_val <- gsub(x = basename(xl_file), pattern = "([0-9]*)?_(.*)", repl = "\\1")
  
  hist_data <- readxl::read_excel(xl_file, sheet = 1, col_types = "numeric") 
  ratio <- as.numeric(na.omit(hist_data$`Area Ratio`))
  area <- as.numeric(na.omit(hist_data$`Area[mm^2]`))
  
  hist_df_full_ratio[[id_val]] <- ratio
  hist_df_full_area[[id_val]] <- area
  
  hist_data <- hist_data  %>% summarize(sd_area = sd(`Area[mm^2]`, na.rm = T),
                                                                      mean_area = mean(`Area[mm^2]`, na.rm = T),
                                                                      median_area = median(`Area[mm^2]`, na.rm = T),
                                                                      median_area = mad(`Area[mm^2]`, na.rm = T),
                                                                      sd_ratio = sd(`Area Ratio`, na.rm = T),
                                                                      mean_ratio = mean(`Area Ratio`, na.rm = T),
                                                                      median_ratio = median(`Area Ratio`, na.rm = T),
                                                                      median_ratio = mad(`Area Ratio`, na.rm = T),
                                                                      )
  
  hist_df <- rbind(hist_df, hist_data %>% mutate(ID = id_val))
}
```

## PCA: Area Ratio

### PCA - 1. und 2. Komponente

``` r
ratio_df <- do.call("cbind", hist_df_full_ratio)
id_names_ratio <- colnames(ratio_df)
ratio_df <- t(ratio_df)

prin_comp <- prcomp(ratio_df, rank. = 3)
pc_vars <- round(prin_comp$sdev / sum(prin_comp$sdev) * 100, 2)

components <- prin_comp[["x"]]
components <- data.frame(components)
components <- cbind(components, ID = id_names_ratio)

components %>% dplyr::inner_join(compl_df %>% dplyr::select(ID, gmm_1_2) %>% mutate(ID = as.character(ID)), by = "ID") %>%
  ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, fill = gmm_1_2)) + 
  ggplot2::geom_point(aes(shape = gmm_1_2), size = 4) + 
  scale_shape_manual(values=c(21, 24)) +
  scale_fill_manual(values = c("#0D0887FF", "#FA9E3BFF")) + 
  ggplot2::xlab(paste0("PC1 (", pc_vars[1],"%)")) +
  ggplot2::ylab(paste0("PC2 (", pc_vars[2],"%)")) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "right",
                 legend.title = ggplot2::element_blank(),
                 legend.text = ggplot2::element_text(colour = "black", size = 12),
                 axis.text = ggplot2::element_text(colour = "black", size = 12),
                 strip.background = ggplot2::element_rect(fill = "white"),
                 axis.title = ggplot2::element_text(face = "bold", size = 12),
                 strip.text = ggplot2::element_text(face = "bold", size = 12))
```

![](analysis_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
ggplot2::ggsave(filename = file.path(plot.dir, "PCA_PC12_voltageverteilung_area_ratio.pdf"), width = 8, height = 6, device = "pdf", dpi = 300)
ggplot2::ggsave(filename = file.path(plot.dir, "PCA_PC12_voltageverteilung_area_ratio.svg"), width = 8, height = 6, device = "svg", dpi = 300)
```

### PCA - ersten drei Komponenten

``` r
fig <- components %>% dplyr::inner_join(compl_df %>% dplyr::select(ID, gmm_1_2) %>% mutate(ID = as.character(ID)), by = "ID") %>% 
       plot_ly(x = ~PC1, y = ~PC2, z = ~PC3, color = ~gmm_1_2, colors = c("#0D0887FF", "#FA9E3BFF"))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = paste0("PC1 (", pc_vars[1],"%)")),
                                   yaxis = list(title = paste0("PC2 (", pc_vars[2],"%)")),
                                   zaxis = list(title = paste0("PC3 (", pc_vars[3],"%)"))))

fig
```

![](analysis_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

## PCA: Area\[mm^2\]

``` r
area_df <- do.call("cbind", hist_df_full_area)
id_names_area <- colnames(area_df)
area_df <- t(area_df)

prin_comp <- prcomp(area_df, rank. = 3)
pc_vars <- round(prin_comp$sdev / sum(prin_comp$sdev) * 100, 2)
components <- prin_comp[["x"]]
components <- data.frame(components)
components <- cbind(components, ID = id_names_area)


components %>% dplyr::inner_join(compl_df %>% dplyr::select(ID, gmm_1_2) %>% mutate(ID = as.character(ID)), by = "ID") %>%
  ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, fill = gmm_1_2)) + 
  ggplot2::geom_point(aes(shape = gmm_1_2), size = 4) + 
  scale_shape_manual(values=c(21, 24)) +
  scale_fill_manual(values = c("#0D0887FF", "#FA9E3BFF")) + 
  ggplot2::xlab(paste0("PC1 (", pc_vars[1],"%)")) +
  ggplot2::ylab(paste0("PC2 (", pc_vars[2],"%)")) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "right",
                 legend.title = ggplot2::element_blank(),
                 legend.text = ggplot2::element_text(colour = "black", size = 12),
                 axis.text = ggplot2::element_text(colour = "black", size = 12),
                 strip.background = ggplot2::element_rect(fill = "white"),
                 axis.title = ggplot2::element_text(face = "bold", size = 12),
                 strip.text = ggplot2::element_text(face = "bold", size = 12))
```

![](analysis_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
ggplot2::ggsave(filename = file.path(plot.dir, "PCA_PC12_voltageverteilung_area.pdf"), width = 8, height = 6, device = "pdf", dpi = 300)
ggplot2::ggsave(filename = file.path(plot.dir, "PCA_PC12_voltageverteilung_area.svg"), width = 8, height = 6, device = "svg", dpi = 300)
```

### PCA - ersten drei Komponenten

``` r
fig <- components %>% dplyr::inner_join(compl_df %>% dplyr::select(ID, gmm_1_2) %>% mutate(ID = as.character(ID)), by = "ID") %>% 
       plot_ly(x = ~PC1, y = ~PC2, z = ~PC3, color = ~gmm_1_2, colors = c("#0D0887FF", "#FA9E3BFF"))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = paste0("PC1 (", pc_vars[1],"%)")),
                                   yaxis = list(title = paste0("PC2 (", pc_vars[2],"%)")),
                                   zaxis = list(title = paste0("PC3 (", pc_vars[3],"%)"))))

fig
```

![](analysis_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

# Validation Outcome after PVI

<!-- Die Outcome-Daten nach PVI sind auf Blatt 4 der Excel-Tabelle erfasst.  -->
<!-- Relevant sind die letzten 4 Spalten:  -->
<!-- a) Recurrence AF ECG ≤12 months,  -->
<!-- b) Recurrence AF ECG FU,  -->
<!-- c) Recurrence AF symptoms +/- ECG ≤12 months und  -->
<!-- d) Recurrence AF symptoms +/- ECG FU.  -->
<!-- Diese vier verschiedenen Variablen haben jeweils Stärken und Schwächen.  -->
<!-- Die Variable a) beachtet nur AF-Rezidive in einem Zeitraum von 12 Monaten nach der PVI die im EKG dokumentiert wurden.  -->
<!-- Die Variable b) beachtet alle mittels eines EKGs dokumentiertem AF-Rezidive des gesamten Follow-Ups (FU). Der Zeitraum des Fus ist jedoch bei den Patienten unterschiedlich lang.  -->
<!-- Die Variable c) inkludiert noch typische Symptome. Da ist infolge der fehlenden EKG-Dokumentation ein gewisser Fehler möglich. Die Variable bietet aber die Chance, dass besonders kranke Patienten, die sich mit ihrer Krankheit abfinden und nicht mehr zum Arzt gehen, nicht verpasst werden. Der Unterschied zwischen c) und d) ist jeweils wieder der Zeitraum des Follow-Ups. -->
<!-- Unsere Hypothesen sind, dass  -->
<!-- i) Patienten mit hohen LVAs auch vermehrt ein AF-Rezidiv zeigen und dass  -->
<!-- ii) unser Klassifikator bzw. ein neuer Klassifikator (vgl. unser letztes Telefonat) das AF-Rezidiv vorhersagen kann. -->

``` r
val_data <- readxl::read_excel("~/Projekte/Huttelmaier_Voltage_Verteilungen/data/Atriale_Kardiomyopathie_Outcome_pseudonymisiert 08-05-24.xlsx", sheet = 1)
val_data <- val_data %>% dplyr::select(ID,
                                       `Recurrence AF ECG ≤ 12 months (0=no, 1=yes)`, 
                                       `Recurrence AF ECG FU (0=no, 1=yes)`,
                                       `Recurrence AF symptoms +/- ECG  FU (0=no, 1=yes)`, 
                                       `Recurrence AF symptoms +/- ECG  ≤ 12 months (0=no, 1=yes)`)

val_compl <- compl_df %>% inner_join(val_data, by = "ID")

vars2check <- c("Recurrence AF ECG ≤ 12 months (0=no, 1=yes)",
  "Recurrence AF ECG FU (0=no, 1=yes)", 
  "Recurrence AF symptoms +/- ECG  FU (0=no, 1=yes)", 
  "Recurrence AF symptoms +/- ECG  ≤ 12 months (0=no, 1=yes)")
```

## Recurrence AF ECG ≤ 12 months (0=no, 1=yes)

``` r
fisher_test_table(df = val_compl, gr_row = "gmm_1_2", gr_col = "Recurrence AF ECG ≤ 12 months (0=no, 1=yes)", as_kable = FALSE) %>% 
   DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1), 
                      c('10', '30', '50', 'All')),
    paging = F))
```

![](analysis_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

## Recurrence AF ECG FU (0=no, 1=yes)

``` r
fisher_test_table(df = val_compl, gr_row = "gmm_1_2", gr_col = vars2check[2], as_kable = FALSE) %>% 
   DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1), 
                      c('10', '30', '50', 'All')),
    paging = F))
```

![](analysis_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

## Recurrence AF symptoms +/- ECG FU (0=no, 1=yes)

``` r
fisher_test_table(df = val_compl, gr_row = "gmm_1_2", gr_col = vars2check[3], as_kable = FALSE) %>% 
   DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1), 
                      c('10', '30', '50', 'All')),
    paging = F))
```

![](analysis_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

## Recurrence AF symptoms +/- ECG ≤ 12 months (0=no, 1=yes)

``` r
fisher_test_table(df = val_compl, gr_row = "gmm_1_2", gr_col = vars2check[4], as_kable = FALSE) %>% 
   DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1), 
                      c('10', '30', '50', 'All')),
    paging = F))
```

![](analysis_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

# Tables

## Table 1

Vergleiche:

- Fisher’s Exakt Test für Vergleiche von binären Variablen (bspw. AF
  (0=paroxysmal, 1=persistent))
  - pro Test zwei Zeilen
  - Spalte Feature gibt die Eigenschaft (Tabllenename aus original
    Excel)
  - Spalte value gibt den Wert der Eigenschaft bspw. **paroxysmal** oder
    **persistent** für Feature **AF**
- Wilcoxon-Rangsummen-Test für stetige Werte durchgeführt
  - pro Test nur eine Zeile
  - Spalte Odds Ratio nur eingefügt, damit die Zeilen der verschiedenen
    Tests übereinandergelegt werden können
  - es gibt im Wilcoxon-Test keine Odds-Ratio

``` r
comp_features <- c("age_examination", "gender (0=female, 1=male)", "AF (0=paroxysmal, 1=persistent)",
                   "BMI (kg/m²)", "time_since_diagnosis (years)", "CHA2DS2-VASc", "AA_therapy (0=no AA, 1=amiodarone, 2=flecainide)",
                   "AA_duration (months)", "aHT    (0=no aHT, 1=aHT)", "stroke (0=no stroke, 1=stroke)",
                   "Diabetes_mellitus (0=no diabetes, 1=type 1, 2=type 2)", "CRP (mg/dl)", 
                   "Creatinin (mg/dl)", "GFR (MDRD)", "NT_pro_BNP (pg/ml)", "HbA1c (%)")


table_df <- compl_df[,c("ID", "gmm_1_2", comp_features)] %>% 
            mutate(gender = factor(if_else(`gender (0=female, 1=male)` == 0, 
                                                               true = "female",
                                                               false = "male"), 
                                 levels = c("female", "male")),
                   `AF` = factor(if_else(`AF (0=paroxysmal, 1=persistent)` == 0, 
                                                               true = "paroxysmal AF",
                                                               false = "persistent AF"), 
                                 levels = c("paroxysmal AF", "persistent AF")),
                   `AA_therapy` = factor(if_else(`AA_therapy (0=no AA, 1=amiodarone, 2=flecainide)` == 0, 
                                                               true = "no AA therapy",
                                                               false = "amiodarone or flecainide"), 
                                 levels = c("no AA therapy", "amiodarone or flecainide")),
                   hypertension = factor(if_else(`aHT    (0=no aHT, 1=aHT)` == 0, 
                                                               true = "no hypertension",
                                                               false = "hypertension"), 
                                 levels = c("no hypertension", "hypertension")),
                   stroke = factor(if_else(`stroke (0=no stroke, 1=stroke)` == 0, 
                                                               true = "no stroke",
                                                               false = "stroke"), 
                                 levels = c("no stroke", "stroke")),
                   diabetes_mellitus = factor(if_else(`Diabetes_mellitus (0=no diabetes, 1=type 1, 2=type 2)` == 0, 
                                                               true = "no diabetes",
                                                               false = "diabetes (type 1 or 2)"), 
                                 levels = c("no diabetes", "diabetes (type 1 or 2)")),
                   ) %>% 
                  dplyr::rename(age = age_examination) %>%
                  dplyr::select(-`gender (0=female, 1=male)`, -`AF (0=paroxysmal, 1=persistent)`,
                                -`AA_therapy (0=no AA, 1=amiodarone, 2=flecainide)`, 
                                -`aHT    (0=no aHT, 1=aHT)`, -`stroke (0=no stroke, 1=stroke)`,
                                -`Diabetes_mellitus (0=no diabetes, 1=type 1, 2=type 2)`)

res_list <- sapply(colnames(table_df)[-(1:2)], function(sel_feature){
  
  if(!all(is.na(table_df %>% dplyr::pull(!!sym(sel_feature))))){

    if(!is.factor(table_df %>% dplyr::pull(!!sym(sel_feature)))){

      return(wilcox_test_table(df = table_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>%
             mutate(`Odds-Ratio` = "-",
                    Feature = sel_feature,
                    value = sel_feature) %>%
             relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`)) 
             
    }else{
      return(fisher_test_table(df = table_df, gr_col = "gmm_1_2", gr_row = sel_feature, as_kable = FALSE) %>% 
               rename(value = !!sym(sel_feature)) %>%
               mutate(Feature = sel_feature) %>%
               relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`))
    }
    
  }else{
    
    return(NA)
    
  }
}, simplify = FALSE)

na_comps_idx <- unlist(lapply(res_list, function(elem) all(is.na(elem))))

all_tests_cluster_list <- list()
all_tests_cluster_list[[1]] <- do.call("rbind", res_list[!na_comps_idx]) %>% relocate(Feature)
```

<!-- ## Table 1 - Sheet 2 -->

``` r
comp_features <- c("IVSd M (cm)", "LVIDd M (cm)", "LVPWd M (cm)", 
                  "LVEDV biplan (ml)", "LVEF biplan (%)", "GLPS AFI average (%)",
                  "GLPS 4D (%)")


res_list <- sapply(seq_along(comp_features), function(i){

  sel_feature <- comp_features[i]
      return(wilcox_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>%
            mutate(`Odds-Ratio` = "-",
                    Feature = sel_feature,
                    value = sel_feature) %>%
             relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`)) 
}, simplify = FALSE)

all_tests_cluster_list[[2]] <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>%
  relocate(Feature)
```

<!-- ## Table 1 - Sheet 2: LA 2D Biplan -->

``` r
comp_features <- c("LAEDV MOD Biplan 2D (ml)", "LAESV MOD Biplan 2D (ml)", 
                   "LAESV Mod Biplan Index BSA 2D (ml/m²)")

res_list <- sapply(seq_along(comp_features), function(i){

  sel_feature <- comp_features[i]
      return(wilcox_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>%
            mutate(`Odds-Ratio` = "-",
                    Feature = sel_feature,
                    value = sel_feature) %>%
             relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`)) 
}, simplify = FALSE)

all_tests_cluster_list[[3]] <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>%
  relocate(Feature)
```

<!-- ## Table 1 - Sheet 2: RA 2D Monoplan -->

``` r
comp_features <- c("RAA (s) (cm²)", "RAESVMod (ml)", "RAESV AL (ml)")

res_list <- sapply(seq_along(comp_features), function(i){

  sel_feature <- comp_features[i]
      return(wilcox_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>%
             mutate(`Odds-Ratio` = "-",
                    Feature = sel_feature,
                    value = sel_feature) %>%
             relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`)) 
}, simplify = FALSE)

all_tests_cluster_list[[4]] <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>%
  relocate(Feature)
```

<!-- ## Table 1 - Sheet 2: 4D -->

``` r
comp_features <- c("LA Vol.  min 4D (ml)", "LA Vol. max 4D (ml)", "LA Vol pre A (ml)", "LAVImax 4D (ml/m²)")

res_list <- sapply(seq_along(comp_features), function(i){

  sel_feature <- comp_features[i]
      return(wilcox_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>%
             mutate(`Odds-Ratio` = "-",
                    Feature = sel_feature,
                    value = sel_feature) %>%
             relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`)) 
}, simplify = FALSE)

all_tests_cluster_list[[5]] <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>%
  relocate(Feature)
```

<!-- ## Table 1 - Sheet 2: LAEF 2D Biplan -->

``` r
comp_features <- c("total LA-EF (%)")

res_list <- sapply(seq_along(comp_features), function(i){

  sel_feature <- comp_features[i]
      return(wilcox_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>%
             mutate(`Odds-Ratio` = "-",
                    Feature = sel_feature,
                    value = sel_feature) %>%
             relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`)) 
}, simplify = FALSE)

all_tests_cluster_list[[6]] <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>%
  relocate(Feature)
```

<!-- ## Table 1 - Sheet 2: 4D -->

``` r
comp_features <- c("LA-EF total 4D (%)", "LA-EF passive 4D (%)", "LA-EF active 4D (%)", 
                  "MV a' av. (m/s)", "LAVI/a'2D ((ml/m²)/(m/s))", "LAVI/a' 4D ((ml/m²)/(m/s))")

res_list <- sapply(seq_along(comp_features), function(i){

  sel_feature <- comp_features[i]
      return(wilcox_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>%
             mutate(`Odds-Ratio` = "-",
                    Feature = sel_feature,
                    value = sel_feature) %>%
             relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`)) 
}, simplify = FALSE)

all_tests_cluster_list[[7]] <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>%
  relocate(Feature)
```

<!-- ## Table 1 - Sheet 2: Distolic function -->

``` r
comp_features <- c("MV e' sept. (cm/s)", "MV e' lat. (cm/s)","MV E/e' average", "TR Vmax (m/s)")

res_list <- sapply(seq_along(comp_features), function(i){

  sel_feature <- comp_features[i]
      return(wilcox_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>%
             mutate(`Odds-Ratio` = "-",
                    Feature = sel_feature,
                    value = sel_feature) %>%
             relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`)) 
}, simplify = FALSE)

all_tests_cluster_list[[8]] <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>%
  relocate(Feature)
```

<!-- ## Table 1 - Sheet 2: 2D strain -->

``` r
comp_features <- c("LASr R av. (%)", "LAScd R av. (%)", "LASct R av. (%)", 
                   "LASr P av. (%)", "LAScd P av. (%)", "LASct P av. (%)")

res_list <- sapply(seq_along(comp_features), function(i){

  sel_feature <- comp_features[i]
      return(wilcox_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>%
             mutate(`Odds-Ratio` = "-",
                    Feature = sel_feature,
                    value = sel_feature) %>%
             relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`)) 
}, simplify = FALSE)

all_tests_cluster_list[[9]] <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>%
  relocate(Feature)
```

<!-- ## Table 1 - Sheet 2: 4D strain -->

``` r
comp_features <- c("LASr 4D (%)", "LAScd 4D (%)", "LASct 4D (%)", 
                  "LASr_c 4D (%)", "LAScd_c 4D (%)", "LASct_c 4D (%)")

res_list <- sapply(seq_along(comp_features), function(i){

  sel_feature <- comp_features[i]
      return(wilcox_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>%
             mutate(`Odds-Ratio` = "-",
                    Feature = sel_feature,
                    value = sel_feature) %>%
             relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`)) 
}, simplify = FALSE)

all_tests_cluster_list[[10]] <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>%
  relocate(Feature)
```

<!-- ## Table 1: -->

``` r
comp_features <- c("map_points", "P_wave_duration (ms)", "PTFV1 (µV*ms)", 
                   "LVA_1.0 (cm²)", "LVA_0.5 (cm²)")

res_list <- sapply(seq_along(comp_features), function(i){

  sel_feature <- comp_features[i]
      return(wilcox_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>%
             mutate(`Odds-Ratio` = "-",
                    Feature = sel_feature,
                    value = sel_feature) %>%
             relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`)) 
}, simplify = FALSE)

all_tests_cluster_list[[11]] <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>%
  relocate(Feature)

 all_tests_cluster_table1 <- do.call("rbind",all_tests_cluster_list)

cbind(all_tests_cluster_table1, "p.adjusted" = p.adjust(all_tests_cluster_table1$`P-value`, method = "BH")) %>% as_tibble() %>%
  mutate(`P-adj. (short)` = dplyr::case_when(`p.adjusted` >= 0.01 ~ as.character(round(`p.adjusted`, 2)),
                                             `p.adjusted` < 0.0001 ~ "< 0.0001",
                                             `p.adjusted` < 0.001 ~ "< 0.001",
                                             `p.adjusted` < 0.01 ~ "< 0.01"),
         p.adjusted = format(as.numeric(p.adjusted), digits = 3, scientific = TRUE),
         `P-value` = format(as.numeric(`P-value`), digits = 3, scientific = TRUE)) %>%
   DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1),
                      c('10', '30', '50', 'All')),
    paging = F))
```

![](analysis_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

### Hinzufügen LVA Welch-Test (Mittelwerte ± SD): 1.0 (cm²) und 0.5 (cm²)

``` r
comp_features <- c("LVA_1.0 (cm²)", "LVA_0.5 (cm²)")

res_list <- sapply(seq_along(comp_features), function(i){

  sel_feature <- comp_features[i]
      return(welch_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>%
             mutate(`Odds-Ratio` = "Welch-Test",
                    Feature = sel_feature,
                    value = sel_feature) %>%
             relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`)) 
}, simplify = FALSE)

all_tests_cluster_list[[12]] <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>%
  relocate(Feature)

 all_tests_cluster_table1 <- do.call("rbind",all_tests_cluster_list)

cbind(all_tests_cluster_table1, "p.adjusted" = p.adjust(all_tests_cluster_table1$`P-value`, method = "BH")) %>% as_tibble() %>%
  mutate(`P-adj. (short)` = dplyr::case_when(`p.adjusted` >= 0.01 ~ as.character(round(`p.adjusted`, 2)),
                                             `p.adjusted` < 0.0001 ~ "< 0.0001",
                                             `p.adjusted` < 0.001 ~ "< 0.001",
                                             `p.adjusted` < 0.01 ~ "< 0.01"),
         p.adjusted = format(as.numeric(p.adjusted), digits = 3, scientific = TRUE),
         `P-value` = format(as.numeric(`P-value`), digits = 3, scientific = TRUE)) %>%
   DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1),
                      c('10', '30', '50', 'All')),
    paging = F))
```

![](analysis_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

## Table 1’ Median: All patients (missing entries)

``` r
med_iqr_feats <- c("CRP (mg/dl)", "GFR (MDRD)", "NT_pro_BNP (pg/ml)", "HbA1c (%)", "LVEF biplan (%)", 
                   "LAVI/a'2D ((ml/m²)/(m/s))", "MV E/e' average", "TR Vmax (m/s)","LVA_1.0 (cm²)", "LVA_0.5 (cm²)")

apply(compl_df[,med_iqr_feats], 2, med_iqr, get_n = TRUE) %>% as_tibble(rownames = "Feature") %>% dplyr::rename("all patients" = value)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["Feature"],"name":[1],"type":["chr"],"align":["left"]},{"label":["all patients"],"name":[2],"type":["chr"],"align":["left"]}],"data":[{"1":"CRP (mg/dl)","2":"0.1 (0.1 - 0.2) n=50"},{"1":"GFR (MDRD)","2":"75 (63 - 85) n=50"},{"1":"NT_pro_BNP (pg/ml)","2":"189.5 (70.8 - 391.2) n=50"},{"1":"HbA1c (%)","2":"5.7 (5.4 - 6) n=50"},{"1":"LVEF biplan (%)","2":"61 (56.2 - 63) n=50"},{"1":"LAVI/a'2D ((ml/m²)/(m/s))","2":"733.3 (523.8 - 890) n=45"},{"1":"MV E/e' average","2":"7.6 (6 - 9) n=49"},{"1":"TR Vmax (m/s)","2":"2.5 (2.4 - 2.7) n=31"},{"1":"LVA_1.0 (cm²)","2":"7.8 (3.9 - 22.9) n=50"},{"1":"LVA_0.5 (cm²)","2":"1.8 (0.5 - 5.9) n=50"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

## Table 1’ Mean ± SD: All patients (missing entries)

``` r
med_iqr_feats <- c("CRP (mg/dl)", "GFR (MDRD)", "NT_pro_BNP (pg/ml)", "HbA1c (%)", "LVEF biplan (%)", 
                   "LAVI/a'2D ((ml/m²)/(m/s))", "MV E/e' average", "TR Vmax (m/s)","LVA_1.0 (cm²)", "LVA_0.5 (cm²)")

apply(compl_df[,med_iqr_feats], 2, mean_sd, get_n = TRUE) %>% as_tibble(rownames = "Feature") %>% dplyr::rename("all patients" = value)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["Feature"],"name":[1],"type":["chr"],"align":["left"]},{"label":["all patients"],"name":[2],"type":["chr"],"align":["left"]}],"data":[{"1":"CRP (mg/dl)","2":"0.27 ± 0.32 n=50"},{"1":"GFR (MDRD)","2":"75.9 ± 18.2 n=50"},{"1":"NT_pro_BNP (pg/ml)","2":"292.28 ± 308.52 n=50"},{"1":"HbA1c (%)","2":"5.84 ± 0.75 n=50"},{"1":"LVEF biplan (%)","2":"60.4 ± 5.35 n=50"},{"1":"LAVI/a'2D ((ml/m²)/(m/s))","2":"965.04 ± 650.9 n=45"},{"1":"MV E/e' average","2":"7.87 ± 2.64 n=49"},{"1":"TR Vmax (m/s)","2":"2.63 ± 0.33 n=31"},{"1":"LVA_1.0 (cm²)","2":"16.34 ± 17.98 n=50"},{"1":"LVA_0.5 (cm²)","2":"6.32 ± 10.84 n=50"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

## Table 2: Comparison of selected echocardiographic parameters between groups mild AC and severe AC.

``` r
comp_features <- c("LASr R av. (%)", "LASr P av. (%)", "LAVI/a' 4D ((ml/m²)/(m/s))", 
                   "LASct R av. (%)", "LASr 4D (%)", "MV e' sept. (cm/s)", "LAVI/a'2D ((ml/m²)/(m/s))", 
                   "LASct_c 4D (%)")

res_list <- sapply(seq_along(comp_features), function(i){

  sel_feature <- comp_features[i]
      return(wilcox_test_table(df = compl_df, group = "gmm_1_2", feature = sel_feature, as_kable = FALSE) %>%
             mutate(`Odds-Ratio` = "-",
                    Feature = sel_feature,
                    value = sel_feature) %>%
             relocate(Feature, value, mild, severe, `Odds-Ratio`, `P-value`)) 
}, simplify = FALSE)

all_tests_cluster <- do.call("rbind", res_list[!unlist(lapply(res_list, function(elem) all(is.na(elem))))]) %>%
  relocate(Feature)

cbind(all_tests_cluster, "p.adjusted" = p.adjust(all_tests_cluster$`P-value`, method = "BH")) %>% as_tibble() %>%
  mutate(`P-adj. (short)` = dplyr::case_when(
                                                 `p.adjusted` >= 0.01 ~ as.character(round(`p.adjusted`, 2)),
                                                   `p.adjusted` < 0.0001 ~ "< 0.0001",
                                                   `p.adjusted` < 0.001 ~ "< 0.001",
                                                   `p.adjusted` < 0.01 ~ "< 0.01"),
         p.adjusted = format(as.numeric(p.adjusted), digits = 3, scientific = TRUE),
         `P-value` = format(as.numeric(`P-value`), digits = 3, scientific = TRUE)) %>%
   DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1),
                      c('10', '30', '50', 'All')),
    paging = F))
```

![](analysis_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

## Table 3: Features evaluated with the Boruta algorithm. Display of all features evaluated with the Boruta algorithm as a function of feature importance.

``` r
df_importance %>% dplyr::filter(decision != "Tentative") %>% 
  group_by(variable) %>% summarise(`Importance score (Median)` = format(median(value), digits = 3, 
                                                   scientific = F),
                                   `Importance score (IQR)` = paste0("[", 
                                                format(quantile(value, 0.25), 
                                                       digits = 3, scientific = F), 
                                                ", ", 
                                                format(quantile(value, 0.75), 
                                                       digits = 3, scientific = F), "]"),
                                   Decision = unique(decision)) %>%
   DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1),
                      c('10', '30', '50', 'All')),
    paging = F))
```

![](analysis_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

## Table 4: Prediction of SVM for each CV run

``` r
svm_linear_result_list$svmFitPredictions_df %>%
   DT::datatable(extensions = 'Buttons', options = list(
    dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf'),
    lengthMenu = list(c(10,30, 50, -1),
                      c('10', '30', '50', 'All')),
    paging = F))
```

![](analysis_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->
