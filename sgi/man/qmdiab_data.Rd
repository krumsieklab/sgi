% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qmdiab_data.R
\docType{data}
\name{qmdiab_clin}
\alias{qmdiab_clin}
\alias{qmdiab_plasma}
\alias{qmdiab_urine}
\alias{qmdiab_saliva}
\title{QMDiab: Qatar Metabolomics Study on Diabetes}
\format{







}
\source{
\url{https://figshare.com/articles/dataset/Qatar_Metabolomics_Study_on_Diabetes/5904022}
}
\usage{
qmdiab_clin

qmdiab_plasma

qmdiab_urine

qmdiab_saliva
}
\description{
The Qatar Metabolomics Study on Diabetes (QMDiab) is a type 2 diabetes case-control study, which was conducted in
2012 at the Dermatology Department of Hamad Medical Corporation and the
Weill Cornell Medical College in Doha, Qatar.

Untargeted metabolomics measurements (LC/MS+, LC/MS-, and GC/MS) from plasma, urine, and saliva samples of 374 participants were performed by Metabolon Inc.

The preprocessed dataset is missing values treated, normalized, transformed, winsorized,
and scaled. Rows correspond to participants and columns correspond to metabolites.
}
\details{
\describe{
\item{\code{qmdiab_clin}}{clinical variables:
AGE, SEX, BMI, HbA1c (Hemoglobin A1c), Diabetes (0=No, 1=Yes), Sodium, Hemoglobin,
HDL-C (High-density lipoprotein cholesterol), Cholesterol, AF score (skin Auto-Fluorescence)
}

\item{\code{qmdiab_plasma}}{plasma metabolomics}

\item{\code{qmdiab_urine}}{urine metabolomics}

\item{\code{qmdiab_urine}}{urine metabolomics}
}
}
\references{
\describe{
\item{Mook-Kanamori, Dennis O., et al. "1, 5-Anhydroglucitol in saliva is a noninvasive marker of
short-term glycemic control." The Journal of Clinical Endocrinology & Metabolism 99.3 (2014): E479-E483.}{
\url{https://doi.org/10.1210/jc.2013-3596}}

\item{Yousri, N.A., Mook-Kanamori, D.O., El-Din Selim, M.M. et al.
Erratum to: A systems view of type 2 diabetes-associated metabolic perturbations in saliva,
blood and urine at different timescales of glycaemic control. Diabetologia 58, 2199 (2015).}{
\url{https://doi.org/10.1007/s00125-015-3682-9}}

\item{Do, K.T., Pietzner, M., Rasp, D.J. et al.
Phenotype-driven identification of modules in a hierarchical map of multifluid metabolic correlations.
npj Syst Biol Appl 3, 28 (2017).
}{\url{https://doi.org/10.1038/s41540-017-0029-9}}
}
}
\keyword{datasets}
