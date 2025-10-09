Warm-up mini-Report: Mosquito Blood Hosts in Salt Lake City, Utah
================
Dawson Smith
2025-10-09

- [ABSTRACT](#abstract)
- [BACKGROUND](#background)
- [STUDY QUESTION and HYPOTHESIS](#study-question-and-hypothesis)
  - [Questions](#questions)
  - [Hypothesis](#hypothesis)
  - [Prediction](#prediction)
- [MOTHODS](#mothods)
  - [Fill in first analysis](#fill-in-first-analysis)
- [put code for analysis here](#put-code-for-analysis-here)
  - [Fill in second analysis/plot](#fill-in-second-analysisplot)
- [Load libraries](#load-libraries)
- [Prepare data](#prepare-data)
- [Calculate map extent with
  padding](#calculate-map-extent-with-padding)
- [Handle edge case: if all values are
  identical](#handle-edge-case-if-all-values-are-identical)
- [Create the hotspot map](#create-the-hotspot-map)
- [GLM: Generalized Linear Modeling](#glm-generalized-linear-modeling)
- [DISCUSSION](#discussion)
  - [Interpretation - fill in
    analysis](#interpretation---fill-in-analysis)
  - [Interpretation - fill in
    analysis/plot](#interpretation---fill-in-analysisplot)
- [CONCLUSION](#conclusion)
- [REFERENCES](#references)

# ABSTRACT

\#West Nile virus (WNV) is maintained in a bird–mosquito transmission
cycle, with certain bird species acting as amplifying hosts. To
investigate which species play this role in Salt Lake City, Utah,
mosquitoes were collected from gravid and other trap types, and
bloodmeal DNA was analyzed using PCR and BLAST to identify avian hosts.
Statistical analyses and graphs were used to compare feeding frequencies
among bird species and to examine how long WNV can persist in avian
hosts. Our results showed that mosquitoes frequently fed on house
finches (Haemorhous mexicanus) and that these locations corresponded
with higher WNV infection rates in mosquito pools. These findings
support the hypothesis that house finches act as important amplifying
hosts of WNV in urban areas of Salt Lake City.

# BACKGROUND

Fill in some text here that provides background info on the WNV system,
the blood meal DNA extractions, PCR, sequencing, etc. and the foundation
for our question/hypothesis.

For example, we can use the viremia duration (Kumar et al., 2003) bar
plot (make sure to reference sources!!!) to illustrate the potential
importance of house finches in WNV transmission and as the logical
foundation for our hypothesis that house finches serve as amplifying
hosts for WNV… and the prediction that locations with more house finches
in our blood host analysis are also the same locations with higher
positive tests for WNV in mosquito pools…

NOTE: Examples of data you can plot for the background info at
<https://github.com/saarman/BIOL3070/>

``` r
# Manually transcribe duration (mean, lo, hi) from the last table column
duration <- data.frame(
  Bird = c("Canada Goose","Mallard", 
           "American Kestrel","Northern Bobwhite",
           "Japanese Quail","Ring-necked Pheasant",
           "American Coot","Killdeer",
           "Ring-billed Gull","Mourning Dove",
           "Rock Dove","Monk Parakeet",
           "Budgerigar","Great Horned Owl",
           "Northern Flicker","Blue Jay",
           "Black-billed Magpie","American Crow",
           "Fish Crow","American Robin",
           "European Starling","Red-winged Blackbird",
           "Common Grackle","House Finch","House Sparrow"),
  mean = c(4.0,4.0,4.5,4.0,1.3,3.7,4.0,4.5,5.5,3.7,3.2,2.7,1.7,6.0,4.0,
           4.0,5.0,3.8,5.0,4.5,3.2,3.0,3.3,6.0,4.5),
  lo   = c(3,4,4,3,0,3,4,4,4,3,3,1,0,6,3,
           3,5,3,4,4,3,3,3,5,2),
  hi   = c(5,4,5,5,4,4,4,5,7,4,4,4,4,6,5,
           5,5,5,7,5,4,3,4,7,6)
)

# Choose some colors
cols <- c(rainbow(30)[c(10:29,1:5)])  # rainbow colors

# horizontal barplot
par(mar=c(5,12,2,2))  # wider left margin for names
bp <- barplot(duration$mean, horiz=TRUE, names.arg=duration$Bird,
              las=1, col=cols, xlab="Days of detectable viremia", xlim=c(0,7))

# add error bars
arrows(duration$lo, bp, duration$hi, bp,
       angle=90, code=3, length=0.05, col="black", xpd=TRUE)
```

<img src="Warm-up-mosquitoes-TEMPLATE-COMPLETE-FINAL_files/figure-gfm/viremia-1.png" style="display: block; margin: auto auto auto 0;" />

# STUDY QUESTION and HYPOTHESIS

## Questions

Fill in here, the question we want to answer… e.g. What bird species is
acting as WNV amplifying host in Salt Lake City?

## Hypothesis

Fill in hypothesis… e.g. House finches are acting as important
amplifying hosts of WNV in Salt Lake City.

## Prediction

Fill in prediction… e.g. If house finches are acting as important
amplifying hosts, we predict that trapping locations where mosquitoes
feed on house finches will also have higher rates of confirmed WNV in
tested mosquito pools.

# MOTHODS

Fill in here, including overview of procedure and methods used for this
project.

## Fill in first analysis

# put code for analysis here

Horizontal plots:

``` r
## import counts_matrix: data.frame with column 'loc_positives' (0/1) and host columns 'host_*'
counts_matrix <- read.csv("./bloodmeal_plusWNV_for_BIOL3070.csv")

## 1) Identify host columns
host_cols <- grep("^host_", names(counts_matrix), value = TRUE)

if (length(host_cols) == 0) {
  stop("No columns matching '^host_' were found in counts_matrix.")
}

## 2) Ensure loc_positives is present and has both levels 0 and 1 where possible
counts_matrix$loc_positives <- factor(counts_matrix$loc_positives, levels = c(0, 1))

## 3) Aggregate host counts by loc_positives
agg <- stats::aggregate(
  counts_matrix[, host_cols, drop = FALSE],
  by = list(loc_positives = counts_matrix$loc_positives),
  FUN = function(x) sum(as.numeric(x), na.rm = TRUE)
)

## make sure both rows exist; if one is missing, add a zero row
need_levels <- setdiff(levels(counts_matrix$loc_positives), as.character(agg$loc_positives))
if (length(need_levels)) {
  zero_row <- as.list(rep(0, length(host_cols)))
  names(zero_row) <- host_cols
  for (lv in need_levels) {
    agg <- rbind(agg, c(lv, zero_row))
  }
  ## restore proper type
  agg$loc_positives <- factor(agg$loc_positives, levels = c("0","1"))
  ## coerce numeric host cols (they may have become character after rbind)
  for (hc in host_cols) agg[[hc]] <- as.numeric(agg[[hc]])
  agg <- agg[order(agg$loc_positives), , drop = FALSE]
}

## 4) Decide species order (overall abundance, descending)
overall <- colSums(agg[, host_cols, drop = FALSE], na.rm = TRUE)
host_order <- names(sort(overall, decreasing = TRUE))
species_labels <- rev(sub("^host_", "", host_order))  # nicer labels

## 5) Build count vectors for each panel in the SAME order
counts0 <- rev(as.numeric(agg[agg$loc_positives == 0, host_order, drop = TRUE]))
counts1 <- rev(as.numeric(agg[agg$loc_positives == 1, host_order, drop = TRUE]))

## 6) Colors: reuse your existing 'cols' if it exists and is long enough; otherwise generate
if (exists("cols") && length(cols) >= length(host_order)) {
  species_colors <- setNames(cols[seq_along(host_order)], species_labels)
} else {
  species_colors <- setNames(rainbow(length(host_order) + 10)[seq_along(host_order)], species_labels)
}

## 7) Shared x-limit for comparability
xmax <- max(c(counts0, counts1), na.rm = TRUE)
xmax <- if (is.finite(xmax)) xmax else 1
xlim_use <- c(0, xmax * 1.08)

## 8) Plot: two horizontal barplots with identical order and colors
op <- par(mfrow = c(1, 2),
          mar = c(4, 12, 3, 2),  # big left margin for species names
          xaxs = "i")           # a bit tighter axis padding

## Panel A: No WNV detected (loc_positives = 0)
barplot(height = counts0,
        names.arg = species_labels, 
        cex.names = .5,
        cex.axis = .5,
        col = rev(unname(species_colors[species_labels])),
        horiz = TRUE,
        las = 1,
        xlab = "Bloodmeal counts",
        main = "Locations WNV (-)",
        xlim = xlim_use)

## Panel B: WNV detected (loc_positives = 1)
barplot(height = counts1,
        names.arg = species_labels, 
        cex.names = .5,
        cex.axis = .5,
        col = rev(unname(species_colors[species_labels])),
        horiz = TRUE,
        las = 1,
        xlab = "Bloodmeal counts",
        main = "Locations WNV (+)",
        xlim = xlim_use)
```

<img src="Warm-up-mosquitoes-TEMPLATE-COMPLETE-FINAL_files/figure-gfm/horiz-plot-1.png" style="display: block; margin: auto auto auto 0;" />

``` r
par(op)

## Keep the colors mapping for reuse elsewhere
host_species_colors <- species_colors
```

## Fill in second analysis/plot

# Load libraries

library(dplyr) library(tidyr) library(ggplot2) library(maps) \# Needed
for the borders() layer library(viridis) \# For the viridis color scale

# Prepare data

hm \<- counts_matrix %\>% filter(!is.na(long), !is.na(lat)) %\>%
transmute( long, lat, `WNV-positive pools (rate)` = loc_rate,
`House finch blood meals` = host_House_finch ) %\>% pivot_longer(
-c(long, lat), names_to = “metric”, values_to = “value” ) %\>%
mutate(value = as.numeric(ifelse(is.na(value), 0, value)))

# Calculate map extent with padding

xpad \<- diff(range(hm$long)) * 0.05
ypad <- diff(range(hm$lat)) \* 0.05

# Handle edge case: if all values are identical

if (xpad == 0) xpad \<- 0.05 if (ypad == 0) ypad \<- 0.05

xlim \<- range(hm$long) + c(-xpad, xpad)
ylim <- range(hm$lat) + c(-ypad, ypad)

# Create the hotspot map

ggplot(hm, aes(long, lat)) + borders(“state”, regions = “utah”, colour =
“grey75”, fill = NA) + \# KDE heatmap layer stat_density_2d_filled(
aes(weight = value, fill = after_stat(level)), contour_var = “ndensity”,
geom = “polygon”, bins = 12, alpha = 0.9, adjust = 0.8 ) + \# Optional
contour outlines stat_density_2d( aes(weight = value, color =
after_stat(level)), contour_var = “ndensity”, bins = 12, linewidth =
0.15, show.legend = FALSE ) + \# Overlay sampling sites geom_point(
aes(size = pmax(value, 0.0001)), alpha = 0.5, shape = 21, stroke = 0.2,
fill = NA, color = “grey30”, show.legend = FALSE ) + coord_quickmap(xlim
= xlim, ylim = ylim) + facet_wrap(~ metric, nrow = 1) +
scale_fill_viridis_c(option = “C”, direction = -1, name = “Hotspot
Intensity”) + guides(color = “none”) + theme_minimal(base_size = 11) +
labs( x = NULL, y = NULL, title = “Hotspot Maps for WNV-positive Pools
and House Finch Blood Meals”, subtitle = “Density estimates across
sampling sites in Utah” )

# GLM: Generalized Linear Modeling

House finch GLM:

Does house finch have an effect by location positive +/- (binary)? by
location positivity rate?

``` r
#glm with house finch alone against binary +/_
glm1 <- glm(loc_positives ~ host_House_finch,
            data = counts_matrix,
            family = binomial)
summary(glm1)
```

    ## 
    ## Call:
    ## glm(formula = loc_positives ~ host_House_finch, family = binomial, 
    ##     data = counts_matrix)
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)       -0.1709     0.1053  -1.622   0.1047  
    ## host_House_finch   0.3468     0.1586   2.187   0.0287 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 546.67  on 394  degrees of freedom
    ## Residual deviance: 539.69  on 393  degrees of freedom
    ## AIC: 543.69
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
#glm with house-finch alone against positivity rate
glm2 <- glm(loc_rate ~ host_House_finch,
            data = counts_matrix)
summary(glm2)
```

    ## 
    ## Call:
    ## glm(formula = loc_rate ~ host_House_finch, data = counts_matrix)
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      0.054861   0.006755   8.122 6.07e-15 ***
    ## host_House_finch 0.027479   0.006662   4.125 4.54e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.01689032)
    ## 
    ##     Null deviance: 6.8915  on 392  degrees of freedom
    ## Residual deviance: 6.6041  on 391  degrees of freedom
    ##   (2 observations deleted due to missingness)
    ## AIC: -484.56
    ## 
    ## Number of Fisher Scoring iterations: 2

\#Answer: YES. Statistically significant association.

\#Possible caveats? Limitations?

# DISCUSSION

## Interpretation - fill in analysis

## Interpretation - fill in analysis/plot

# CONCLUSION

# REFERENCES

1.  Komar N, Langevin S, Hinten S, Nemeth N, Edwards E, Hettler D, Davis
    B, Bowen R, Bunning M. Experimental infection of North American
    birds with the New York 1999 strain of West Nile virus. Emerg Infect
    Dis. 2003 Mar;9(3):311-22. <https://doi.org/10.3201/eid0903.020628>

2.  ChatGPT. OpenAI, version Jan 2025. Used as a reference for functions
    such as plot() and to correct syntax errors. Accessed 2025-10-09.
