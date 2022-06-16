\name{cmef}
\alias{cmef}
\docType{data}
\title{Performance and subsidies of EU farms.}
\description{Aggregated data on farm performance and public subsidies for EU countries in the period 2004-2019.
Source: Common Monitoring and Evaluation Framework (CMEF) for the CAP 2014-2020 (European Commission, 2020).}
\usage{data(cmef)}
\format{
  A data.frame with a total of 433 observations on the following 15 variables:
  \describe{
    \item{\code{Country}}{Country name.}
    \item{\code{Country_code}}{Country code.}
    \item{\code{Year}}{Time of measurement (year).}
    \item{\code{Labour}}{Total labour input (annual work unit), item SE010.}
    \item{\code{Land}}{Total utilised agricultural area (ha), SE025.}
    \item{\code{Capital}}{Total intermediate consumption (euro), item SE275.}
    \item{\code{Int_cons}}{Depreciation of fixed capital, excluding land (euro), SE360.}
    \item{\code{Inputs_total}}{Total inputs (euro), item SE270.}
    \item{\code{Output_total}}{Total output (euro), item SE131.}
    \item{\code{TFP}}{Total factor productivity, computed as the ratio of total output to total input, item SE132.}
    \item{\code{Net_income}}{Farm net income (euro), item SE420.}
    \item{\code{Subs_prod}}{Subsidies on production (euro), item SE605.}
    \item{\code{Subs_inv}}{Subsidies on investments (euro), item SE406.}
    \item{\code{Subs_rur}}{Total support for rural development (euro), item SE624.}
    \item{\code{Subs_dec}}{Decoupled payments (euro), item SE630.}
  }
}
%\note{}
\references{
European Commission (2020). Common Monitoring and Evaluation Framework (CMEF) for the CAP 2014-2020. https://agridata.ec.europa.eu/extensions/DataPortal/cmef_indicators.html
}