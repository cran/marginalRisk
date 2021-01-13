\name{bias.factor}
\alias{bias.factor}
\alias{E.value}
\alias{controlled.risk.bias.factor}
\title{
E-value and Controlled Risk Curve
}
\description{
Compute E-values (Equation 4 in Gilbert et al.) and controlled risk curve bias factor (Equation 6 in Gilbert et al.).
}
\usage{

bias.factor(RRud, RReu)

E.value(rr)

controlled.risk.bias.factor(ss, s.cent, s1, s2, RRud)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{RRud}{
%%     ~~Describe \code{RRud} here~~
}
  \item{RReu}{
%%     ~~Describe \code{RReu} here~~
}

  \item{rr}{
%%     ~~Describe \code{rr} here~~
}

  \item{ss}{
A vector of marker values
}
  \item{s.cent}{
Central marker value
}
  \item{s1}{
s1 and s2 are a pair of marker values for which we set a RRud.
}
  \item{s2}{
%%     ~~Describe \code{s2} here~~
}
}
\details{
These three functions constitute an implementation of the core functionality in Gilbert et al. (2020). For examples on how to use these functions, see the code for Gilbert et al. at https://github.com/youyifong/CoPveryhighVE \cr
\cr

}
\value{
controlled.risk.bias.factor returns a vector of bias factors corresponding to the vector of marker values in ss.
}
\references{
Gilbert, Fong, Carone (2020) Assessment of Immune Correlates of Protection via Controlled Riskof Vaccine Recipients
}
