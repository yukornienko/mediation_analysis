# Mediation of effects of persistent chemicals on the human sperm epigenome

**Student:**

Julia Kornienko

**Supervisors:**

Oleg Sergeev, Yulia Medvedeva

## Project description 

The Russian Children’s Study is a prospective cohort of 516 boys who were enrolled at 8–9 years of age and provided semen samples at 18–19 years of age. RRBS of sperm was conducted to identify the methylation level of CpG denucleotides. 
At the moment of enrollment into the study, the TCDD dioxine (which is one of the most harmful endocrine disrupting chemicals) concentration in the blood of each boy was measured for further longitudinal study of its influence on the reproductive health. Moreover, each boy visited the clinic biennially - for blood sampling; annually - for urine sampling, follow up of growth and puberty  and interviewing => 20 000+ sample aliquots and 1000+ analyzing parameters were collected in total for further analysis. 

**What is already known?**

• Peripubertal exposure to TCDD is associated with poorer semen quality (RCS, Minguez-Alarcon et al., 2017)

• 52 differentially methylated regions (DMRs) were identified that distinguished lowest and highest peripubertal serum TCDD concentrations (RCS, Pilsner et al., 2018)

**What is needed to be known?**

• How do other factors influence the methylation level of the human sperm?

• Which factors can mediate the effect of the peripubertal exposure to TCDD?
 
**Aim of the project:**

Mediation analysis using regression models and longitudinal design

Predictors: peripubertal TCDD concentration and smoking

Outcomes: sperm methylation profiles (all CpGs with coverage >=10 using RRBS)

**What was the data to analyse:**

• 34 samples with different concentrations of TCDD in the boy’s serum at enrollment in the study (8-9 years old) (you can see the histogram of TCDD concentration destribution on the figure 1).

• Methylation level of 2 611 773 CpGs with coverage >=10 presented in at least one of 34 samples

• Data regarding lifestyle habits of each of 34 chosen participants (you can see the list of questions based on which we evaluated smoking within 6 months before sperm collection as a range variable with 6 categories on the figure 2 and bar chart with the count of partisipants belonging to each category on the figure 3).

**What was going to be done with these data:**

1. Building the linear regression model of TCDD concentration and methylation level
2. Regress the smoking on the TCDD (predictor) to recognize whether the TCDD is a significant predictor of the mediator. If the mediator is not associated with the TCDD level, then it could not possibly mediate anything and then these two variables can be used only as independent predictors in multi-factorial regression model.
3. Regress the CpGs methylation on smoking to confirm that it is a significant predictor of the CpG methylation.
4. Building multi-factorial regression model of TCDD concentration in the boy’s serum at enrollment in the study and smoking within 6 months before sperm collection either as two independent predictors or as a predictor and mediator (based on the result of n.2) and methylation level of sperm as a dependent variable.
5. Collect CpGs whose methylation level is significantly dependent on TCDD concentration and/or smoking (three models in total).
6. Mapping of all significant CpGs found in each model to the human genome
7. Biological sense analysis of significant CpGs found (functional enrichment analyses)

**Breif summary of the results:**

1. Linear regression models for the influences of TCDD concentration in prepubertal age and smoking within 6 months before collection of sperm on the young adults semen methylation level were built (both for separate predictors and for their combination).
2. Predictors were considered to be independent on each other.
3. Only CpGs whose methylation levels were significantly dependent on the TCDD concentration and/or smoking, were chosen for further analysis. 
4. All significant CpGs were mapped to either gene, promoter or enhancer regions.
5. Enrichment analysis for all genes associated with CpGs, whose methylation levels were significantly dependent on one or both of predictors, was performed.
6. Some interesting findings are needed to be analyzed in more details.


## Detailed description of the results:

*Samples selection:*

51 out of 516 samples were chosen by the following criterias:

•	Prepubertal TCDD exposure at 8-9 years old

•	Semen quality at 18-19 years old

•	Frozen semen samples at 18-19 years old

•	Buffy coat samples at 18-19 years old closest to date of semen sampling

To start with, 10 IDs with either the highest or the lowest TCDD concentration were selected. Then were selected one ID with the highest semen quality and one with the lowest. Later, 39 more IDs were chosen by random selection 13 IDs belonging to each tercile in terms of TCDD concentration. Therefore the data set of 51 samples was created. Finally, 34 samples with satisfactory sequencing quality were selected for this project .


*Data preprocessing:*

From the file containing information regarding all CpGs presented in at least 1 of all 40 samples (with CpGs coverage > 10x), I extracted methylation values for the selected 34 samples (the list of selected samples was given in the another file)  
The script is saved under the name Creating_two_files_with_CpGs.py.

Then for further analysis the dataset with 2 611 773 was restricted according to the following criterias:
- Methilation range of CpG is more than 20 (histogram of methylation range distribution - figure 4)
- This CpG is presented in at least 10 of 34 samples (histogram of sample number distribution - figure 5)
So the subset of initial dataframe was created (it contains 307 538 CpGs which satisfy described criterias)
The script is saved under the name histograms.py.

*Regression model of TCDD concentration (predictor) and methylation level (dependent variable):*

Next the linear regression models were built. This was performed with the Dioxines_CpGs_linear_regression.py script, which takes file with CpG methilation level and TCDD concentration as inputs and gives table with all regression coefficients and parameters for each CpG position, as well as the descriptive statistics for this CpG methylation level destribution. Here the scipy.stats python package was used for the linear regression model with one continuous predictor. 


*Regression model of TCDD concentration (predictor) and smoking (possible mediator):*

To answer the question whether the smoking habit can act as a mediator of TCDD exposure effects on human sperm epigenome, we performed regression (TCDD concentration - predictor and smoking - dependent variable) to recognize if smoking habit is dependent on the TCDD exposure. It was done in R with simple functions (the code is somewhere in Stuff.R file).
It appeared that:
When smoking is a binary variable (did the person smoke within 6 months before semen collection), there is significant (pvalue = 0.03) but slight (R2 = 0.14) dependence of smoking habbit on TCDD concentration
When smoking is a range variable (did and how much did the person smoke within 6 months before semen collection) – no significant dependence was observed. To confirm that, the histogram of TCDD concentration, filled with the category of smoking variable, was plotted (figure 6).
For further analysis the range variable was chosen as more informative one, so TCDD concentration at the moment of enrollment and smoking within 6 months before semen collection were used as two independent predictors of methylation level in sperm.


*Regression model of smoking (predictor) and methylation level and multi-factor regression model with two predictors (smoking and TCDD concentration):*

To perform the regression analysis of smoking influence on CpG methylation level, the statsmodels.formula.api python package was used (as smoking is a categorial variable, the scipy.stats python package can not be used). 
For the multi-factorial regression analysis with two predictors, the same package was exploited.
It was done with the All_regressions_statsmodels.py script, which takes the file with CpG methilation level, TCDD concentration and smoking as inputs and gives table with all regression coefficients and parameters for each CpG position, as well as descriptive statistics of this CpG methylation level destribution as an output. 
Moreover, we also built three additional regression models with another representation of the same predictors:
- TCDD as a categorial variable (belonging to one of three TCDD concentration tertiles)
- Smoking as a categorial variable with reduced number of categories (three instead of six)
- Multi-factor regression model with these two predictors (they were also concidered to be independent on each other).
However, for further analysis we chose the first set of predictors representation (TCDD as continuous variable and smoking with 6 categories) as they are more informative ones. 


*Selection of only CpGs whose methylation level is significantly dependent on one of or both of predictors (three models):*

For each model, the p-value fdr correction was performed (Stuff.R). Then, based on the results of regressions, selection of CpGs, whose methylation level was significantly dependent on the predictors, was performed according to the following criterias:
- P value of the overall model with fdr correction is less than 0.05
- Absolute value of R2 is more than 0.7
 
 The table with the total amount of significat CpGs in each model is on the figure 7.
 
*Mapping of all significant CpGs found in each model to the human genome:*

For mapping of all significant CpGs found in each model to the human genome regions, the following data were used:
 - Promoter regions: 2000 upstream:1000 downstream from TSS (Ensembl annotation)
 - Intraenhancer regions (Fantom5 annotation)
 - Intragenic regions (hg19 annotation)

For mapping, positions of all significant CpGs in each model were saved as bed-files. The mapping itself was performed with the use of BEDOPS v2.4.35 tool (closest-features command, example: closest-features --closest --delim "\t" --dist only_smoke6.bed  enhancers_hg19.bed  > smoke_enhansers.bed).
Then, all found features were filtered according to the distance to the closest feature (not more than 500bp upstream/downstream far from the promoter regions, and only intra-enhancer and intra-genic regions) and merged in the one file containing information about all genes, associated with significant CpGs found in each model.
These sets of genes were then used for the enrichment analysis in each model.
The script of data processing - Mapping_and_enrichment.R.
Amount of all promoters, enhancers and intragenic found to be associated with significant CpGs in each model is presented on the figure 8. 

*Enrichment analysis for all genes, associated with significant CpGs found in each model:*

Finally, enrichment analysis for all genes, associated with significant CpGs found in each model was performed. It was done with the use of R packages (see Mapping_and_enrichment.R) and three databases: KEGG, GO, Reactome. 
For each model, the p-value<0.05 cutoff of enrichment result was chosen. The list of most enriched pathways in each model is shown on the figure 9. The example of cnetplot for the top 5 enriched GO categories for the genes significantly associated with  TCDD and Smoking  is shown on the figure 10. 
However, the more detailed analysis of obtained results in each model is needed to be done in the nearest future.




