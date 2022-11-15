\name{fadn}
\alias{fadn}
\docType{data}
\title{Farm performance and subsidies in EU countries}
\description{Data on farm performance and public subsidies in EU countries in the period 1989-2019 from the Farm Accountancy Data Network (FADN, European Commission, 2020). Values refer to nationally representative farms.}
\usage{data(fadn)}
\format{
  A data.frame with a total of 668 observations on the following 15 variables:
  \describe{
    \item{\code{Country}}{Country name.}
    \item{\code{Country_code}}{Country code.}
    \item{\code{Year}}{Time of measurement (year).}
    \item{\code{Labour}}{Total labour input (annual work unit), item \code{SE010}.}
    \item{\code{Land}}{Total utilised agricultural area (ha), \code{SE025}.}
    \item{\code{Capital}}{Total intermediate consumption (euro), item \code{SE275}.}
    \item{\code{Int_cons}}{Depreciation of fixed capital, excluding land (euro), item \code{SE360}.}
    \item{\code{Inputs_total}}{Total inputs (euro), item \code{SE270}.}
    \item{\code{Output_total}}{Total output (euro), item \code{SE131}.}
    \item{\code{TFP}}{Total factor productivity, computed as the ratio of total output to total input, item \code{SE132}.}
    \item{\code{Net_income}}{Farm net income (euro), item \code{SE420}.}
    \item{\code{Subs_prod}}{Subsidies on production (euro), item \code{SE610+SE615}.}
    \item{\code{Subs_inv}}{Subsidies on investments (euro), item \code{SE406}.}
    \item{\code{Subs_rur}}{Total support for rural development (euro), item \code{SE624}.}
    \item{\code{Subs_dec}}{Decoupled payments (euro), item \code{SE630}.}
  }
}
%\note{}
\references{
European Commission (2020). Farm Accountancy Data Network (FADN) public database. Accessed: 2022/06/21.
https://agridata.ec.europa.eu/extensions/FarmEconomyFocus/FADNDatabase.html
}