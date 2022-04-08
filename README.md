# Imputation-LaboratoryValues-EHR_v2.0
#### This repository has all the codes related to the manuscript - <i> "Imputation of Missing Values for EHR Laboratory Data" </i>
#### <ins>Objective</ins>: 
Laboratory measurements from Electronic Health Records (EHR) are increasingly used in machine learning, however, missingness of lab values are rarely considered. Mishandling of missingness could lead to biased estimation. We investigated patterns of missingness in laboratory variables, and evaluated performance of commonly imputation algorithms based on lab data from two distinct healthcare systems.

#### <ins>Methods</ins>: 
We assessed the missingness pattern for lab measures and applied two commonly-used imputation methods (Multi-level (2l.pan), single-level imputation methods) in combination with 2 imputation technology (monotone imputation, fully conditional specification (FCS)) to impute missing values. We evaluated the performance of imputation methods using normalized RMSE (nRMSE). We further conducted a case study to illustrate imputed lab value (i.e. HbA1c) has improved model performance in prediction.

#### <ins>Results</ins>: 
The pattern of missingness was not at random and was highly associated with patients’ comorbidity data. Multi-level imputation (2l.pan) showed smaller nRMSE for most variables compared to other methods. In the case study of HbA1c lab result, we further evaluated how the imputed values impacts on predicting microvascular outcome. Univariate imputation using multi-level model with FCS, which took comorbidity as latent variables in the imputation, has superior performance compared to other methods.

#### <ins>Conclusion</ins>: 
Overall performance of multi-level method (2l.pan) is superior to cross-sectional pmm method. Multi-level univariate imputation using latent variables derived from comorbidity showed better performance for variables with high missingness. The better imputed value could potentially further improve model performance in prediction.

#### <ins>File name description</ins>:
1. "<i>03012021_imputation_codes.R</i>"- has the codes related to statistical anlayses as well as data visualization
2. "<i>03012021_imputation_function.R</i>"- has the codes related to all imputation algorithms used in this manuscript.

# Publication
[Imputation of missing values for electronic health record laboratory data](https://www.nature.com/articles/s41746-021-00518-0)
