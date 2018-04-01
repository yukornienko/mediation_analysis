# mediation_analysis
**Student:**

Julia Kornienko

**Supervisors:**

Oleg Sergeev, Yulia Medvedeva

**What is already known?**

• Peripubertal exposure to TCDD is associated with poorer semen quality (RCS, Minguez-Alarcon et al., 2017)

• 52 differentially methylated regions (DMRs) were identified that distinguished lowest and highest peripubertal serum TCDD concentrations (RCS, Pilsner et al., 2018)

**What is needed to be known?**

• How do other factors influence the methylation level of the human sperm? • Which factors can mediate the effect of the peripubertal exposure to TCDD?
 
**Aim of the project:**
Mediation analysis using regression models and longitudinal design
Predictor: peripubertal TCDD concentration
Mediators: smoking and tempo of puberty
Outcomes: sperm methylation profiles (all CpGs with coverage >=10 using RRBS

**What are the data to analyse:**

• 34 samples with different concentration of TCDD in the boy’s serum at enrollment in the study (8-9 years old) 
• Methylation level of 2 611 773 CpGs with coverage >=10 presented in at least one of 34 samples
• Data of lifestyle and pubertal tempo of each of 34 selected participants

**What is going to be done with these data:**

1. Building the linear regression model of TCDD concentration and methylation level
2. Selection only those CpGs, which methylation level is significantly dependent on the TCDD concentration.
3. Regress the smoking and pubertal tempo (mediators) on the TCDD (predictor) to confirm that the TCDD is a significant predictor of the mediator. If the mediator is not associated with the TCDD level, then it couldn’t possibly mediate anything
4. Regress the CpGs methylation on both the mediators (smoking and pubertal tempo) and TCDD level to confirm that the mediators are significant predictor of the CpG methylation, and the previously significant predictor in Step #1 is now greatly reduced, if not - nonsignificant
5. To collect a) CpGs with significant mediation by pubertal tempo and/or smoking and b) non-mediated CpGs
6. Biological sense analysis (functional enrichment analyses)




###############
*Samples selection:*

51 out of 516 samples were chosen by the following criterias:

•	Prepubertal TCDD exposure at 8-9 years old
•	Semen quality at 18-19 years old
•	Frozen semen samples at 18-19 years old
•	Buffy coat samples at 18-19 years old closest to date of semen sampling


To start with, 10 IDs with either the highest or the lowest TCDD concentration were selected. Then were selected one ID with the highest semen quality and one with the lowest. Later, 39 more IDs were chosen by random selection 13 IDs belonging to each tercile in terms of TCDD concentration. Therefore the data set of 51 samples was created. Finally, 34 samples with satisfactory sequencing quality were selected for this project .



*Files description:*

all_cpgs_10x.txt - all CpGs that are present in at least 1 sample out of 40 - 7 132 072 CpGs CpGs (not filtered, contains empty rows).
all_cpgs_10x_filtered.csv - all CpGs at least 10X coverage that are present in at least 1 sample out of 40 - 2 647 787 CpGs.

34_CpGs_10x.csv - all CpGs with at least 10X coverage that are present in at least 1 sample out of 34 - 2 611 773 CpGs.
34_CpGs_10x_present_everywhere.csv - all CpGs with at least 10X coverage that are present in all samples out of 34 - 29 161 CpGs.

34_all_CpGs_10x_with_regression_and_CpGs_info.csv - updated 34_CpGs_10x.csv file with regression parameters and descriptive statistics calculated for each CpG.

34_CpGs_10x_with_regression_and_CpGs_info_sign.csv - a subset of the 34_all_CpGs_10x_with_regression_and_CpGs_info.csv with p_value <0.05 and |R2| > 0.80 for CpGs presented in at least 3 samples. Moreover to this file was added a column with nstq2378d variable and its descriptive statistics.

TCDD_info.csv - descriptive statistics calculated for the predictor variable nstq2378d.

################

From the file all_cpgs_10x.txt containing information regarding all CpGs presented in at least 1 of all 40 samples (with CpGs coverage > 10x), I extracted methylation values for the selected 34 samples (the list of selected samples was given in a file RRBS_TCDD_34 selected subjects.xls), therefore the file named 34_CpGs_10x.csv (tab-separated csv file without empty rows) was created. Moreover, for the selected 34 samples, I extracted methylation values for only those CpGs which are presented in all samples into the another file named 34_CpGs_10x_present_everywhere.csv.
The script is saved under the name Creating_two_files_with_CpGs.py.

##########

Next the linear regression models (y = Ax +B, where y was the vector of methylation level for the one of 2 611 773 CpGs from file 34_CpGs_10x.csv and x was the vector of nstq2378d values) were built. Moreover, this file was updated with the mean, median, st. deviation and interquartile range for each of CpGs.

The output file with the regression coefficients (TCDD_A - intercept, TCDD_B - slope), and regression parameters (TCDD_R2 - R2, TCDD_P - p-value, st_err - standard error) is named 34_CpGs_10x_with_regression_and_CpGs_info.csv and downloaded on the google drive. For the CpGs which were analyzed in less than 3 samples, in all the new positions (TCDD_A, TCDD_B, TCDD_R2, TCDD_P, st_err) "-" were written. Comparative statistics for each of  2 611 773 CpGs were saved in the columns named Mean, St_d, Min, 1st_q, Median, 3rd_q, Max.

This was performed with pre-processed data from files 34_CpGs_10x.csv and RRBS_TCDD_34_selected_subjects.csv. The code for the data pre-procession, further regression building and descriptive statistics calculating was written on the Python3. Linear regression was built with the use of scipy.stats library (this script is downloaded on google-drive under the name Dioxines_CpGs_linear_regression.py and can be run).

Moreover, the same statistics were calculated for the predictor nstq2378d and downoaded under the name TCDD_info.csv.

Plots showing the distribution of independent variable (nstq2378d, predictor) and all regression parameters were downloaded to the Figures folder. 
Plots showing the distribution of all regression parameters were also downloaded to the Figures folder.

Significant correlations (p < 0.05 and |R2| > 0.80) between CpG methylation level and TCDD concentration were then found and a new data frame containing only information regarding these TCDD-dependent CpGs was created and named 34_CpGs_10x_with_regression_and_CpGs_info_sing.csv. 26 395 CpGs were found to be strongly associated with the TCDD concentration.




  
