{smcl}
{* updated 12jun2023}{* Ian White}{...}
{* updated 16may2022}{* Ian White}{...}
{* 7sep2007}{* Ian White}{...}
{vieweralsosee "metan (if installed)" "metan"}{...}
{vieweralsosee "metamiss2 (if installed)" "metamiss2"}{...}
{viewerjumpto "Description" "metamiss##Description"}{...}
{viewerjumpto "Data" "metamiss##Data"}{...}
{viewerjumpto "Syntax" "metamiss##Syntax"}{...}
{viewerjumpto "Meta-analysis options" "metamiss##Metaanalysisoptions"}{...}
{viewerjumpto "Imputation options" "metamiss##Imputationoptions"}{...}
{viewerjumpto "Imputation using reasons" "metamiss##Imputationusingreasons"}{...}
{viewerjumpto "Imputation weight options" "metamiss##Imputationweightoptions"}{...}
{viewerjumpto "Imputation debugging options" "metamiss##Imputationdebuggingoptions"}{...}
{viewerjumpto "Bayes options" "metamiss##Bayesoptions"}{...}
{viewerjumpto "Saved variables" "metamiss##Savedvariables"}{...}
{viewerjumpto "Examples" "metamiss##Examples"}{...}
{viewerjumpto "Note" "metamiss##Note"}{...}
{viewerjumpto "References" "metamiss##References"}{...}
{viewerjumpto "Please cite this program" "metamiss##Pleasecitethisprogram"}{...}
{viewerjumpto "Authors and Updates" "metamiss##Authors"}{...}
{viewerjumpto "See also" "metamiss##Seealso"}{...}
{title:Title}

{phang}
{bf:metamiss} {hline 2} Meta-analysis with missing data


{title:Description}{marker Description}

{pstd}
{cmd:metamiss} performs a meta-analysis of studies comparing a binary outcome between two groups (such as in a randomised controlled trial) where the outcome may be missing. Three methods are available:

{pstd}
(i) imputation methods as described in {help metamiss##Higgins++08:Higgins et al (2008)}. 
Missing values may be imputed as failures, successes, same as control group, same as 
experimental group, same as own group, or using Informative Missingness Odds Ratios (IMORs).
When reasons for missingness are known, a mixture of the methods may be used.

{pstd}The IMOR is defined as the odds ratio between the outcome (success/failure) and the missingness of the outcome (missing/observed), within groups.
A value of 1 indicates that the data are missing at random (MAR), while IMOR = 0 means that missing values are all failures, and IMOR = infinity means that missing values are all successes.

{pstd}or 

{pstd}
(ii) methods allowing for uncertainty about how informative the missing data are, as in {help metamiss##White++08:White et al (2008)} based on {help metamiss##ForsterSmith:Forster and Smith (1998)}. 
These use the following model for outcome Y, missingness M and group X:
{p_end}
          M|X=i       ~ Bernoulli(alpha_i)
          Y|M=0,X=i   ~ Bernoulli(pi_i)
          OR(Y,M|X=i) = IMOR_i (informative missingness odds ratio)

{pstd}
The values of logIMOR_i in the two groups can be specified exactly or through a Normal prior distribution. 

{pstd}or 

{pstd}(iii) the "uncertainty method" of {help metamiss##GambleHollis:Gamble and Hollis (2005)}.

{pstd}
With all three methods, {cmd:metamiss} only prepares the data for each study: it then calls {help metan} to 
perform the meta-analysis.
A note "(Calling metan ...)" is printed at this point 
to help the user identify the source of any error messages.
We have tried to make {cmd:metamiss} compatible with all versions of {help metan}.


{title:Data}{marker Data}

{pstd}
{cmd:metamiss} requires 6 variables specifying the numbers of successes, 
failures and missing values in each group. The first group is considered to be the experimental group and the second group the control group, although applications outside randomised controlled trials are also possible.
The 6 variables are described as {it:rE fE mE rC fC mC} in the syntax descriptions below.


{title:Syntax}{marker Syntax}

Simple imputation:

{phang}
   {cmd:metamiss} {it:rE fE mE rC fC mC}, {it:imputation_type} {cmd:w1|w2|w3|w4}
                      [{it:meta_options} {it:imor_spec}]
   
{phang}
   where {it:imputation_type} is {cmd:aca|ica0|ica1|icapc|icape|icap|icab|icaw|icaimor}

Imputation using reasons:

{phang}
   {cmd:metamiss} {it:rE fE mE rC fC mC}, {it:imputation_type_1}{cmd:(#|var [#|var])}
                      [{it:imputation_type_2}{cmd:(#|var [#|var])} [...]] 
                      {cmd:w1|w2|w3|w4} [{it:meta_options} {it:imor_spec}] 

{phang}
   where {it:imputation_type_1}, {it:imputation_type_2} etc. are {cmd:ica0|ica1|icapc|icape|icap|icaimor}.

Bayesian analysis using priors:

{phang}
   {cmd:metamiss} {it:rE fE mE rC fC mC}, {cmdab:sd:logimor(#|var [#|var])}
                      [{it:meta_options} {it:imor_spec} {cmdab:corr:logimor(#|var)} 
                      {cmd:method(GH|Taylor|MC)} {cmd:reps(#)} {cmd:nip(#)} {cmdab:miss:prior(# # [# #])} 
                      {cmdab:resp:prior(# # [# #])} {cmd:details} {cmd:nodots}]

Gamble-Hollis analysis:

{phang}
   {cmd:metamiss} {it:rE fE mE rC fC mC}, {cmdab:gamble:hollis} [{it:meta_options}]

{phang}
{it:meta_options} are meta-analysis options including {cmd:or rr rd log id(varname)}
   and any valid option for {help metan} including {cmd:random, by(), xlab()}. 

{phang}
{it:imor_spec} is {cmd:imor(#|var [#|var])} or {cmd:logimor(#|var [#|var])}.


{title:Meta-analysis options}{marker Metaanalysisoptions}

{phang}
{cmd:rr, or, rd} specify the measures to be analysed. 
Usually, only one measure may be specified; the default is {cmd:rr}.
The exception is using {cmd:method(mc)}, when all 3 measures may be obtained for no extra effort, 
so any combination is allowed, and the default is all 3. However, the formal meta-analysis
is not performed if more than 1 measure is specified.

{phang}
{cmd:log} has the results reported on the log RR or log OR scale.

{phang}
{cmd:id(var)} specifies a trial identifier for the forest plot.

{phang}
All other options allowed with {help metan} are allowed: for example, {cmd:by()}, {cmd:random}.


{title:Imputation options}{marker Imputationoptions}

{phang}{cmd:aca} performs an available cases analysis.

{phang}{cmd:ica0} imputes missing values as zeroes.

{phang}{cmd:ica1} imputes missing values as ones.

{phang}{cmd:icab} performs best-case analysis, which imputes missing values as ones 
in the experimental group and zeroes in the control group.
Note that if {it:rE} and {it:rC} count adverse events, 
not beneficial events, then {cmd:icab} will yield a worst-case analysis.

{phang}{cmd:icaw} performs worst-case analysis, which imputes missing values as
zeroes in the experimental group and ones in the control group.
Note that if {it:rE} and {it:rC} count adverse events, 
not beneficial events, then {cmd:icaw} will yield a best-case analysis.

{phang}{cmd:icape} imputes missing values using the observed probability in the experimental group.

{phang}
{cmd:icapc} imputes missing values using the observed probability in the control group.

{phang}
{cmd:icap} imputes missing values using the observed probability within groups.

{phang}
{cmd:icaimor} imputes missing values using the IMORs specified by {cmd:imor()} or {cmd:logimor()} within groups.

{phang}
The default is {cmd:icaimor} if {cmd:imor()} or {cmd:logimor()} is specified, and otherwise {cmd:aca}.


{title:Imputation using reasons}{marker Imputationusingreasons}

{pstd}
Alternatively, different subgroups of missing values may be assigned using
different reasons. Thus for example ica0(mfE mfC) icap(mpE mpC) indicates 
that mfE individuals in group E and mfC individuals in group C are imputed using 
ICA0, while mpE individuals in group E and mpC individuals in group C are imputed 
using ICAp. If, for some group,  the total over all reasons does not equal the number of missing observations (e.g. if mfE+mpE does not equal mE), then 
the missing observations are shared between imputation types in the given ratio. 
If the total is zero for some group then the missing observations are shared 
between imputation types in the ratio formed by summing over all studies. 
If the total is zero for all groups then an error is returned.
Numerical values may also be given: e.g. ica0(50 50) icap(50 50) indicates that 50% 
missing values in each group are imputed using ICA0 and the rest are imputed 
using ICAp.


{title:Imputation weight options}{marker Imputationweightoptions}

{phang}
{cmd:w1} specifies that standard errors be computed treating the imputed 
values as if they were observed. This is included for didactic purposes 
and should not be used in real analyses.

{phang}
{cmd:w2} specifies that standard errors from the available cases analysis should be used.

{phang}
{cmd:w3} specifies that standard errors be computed by scaling the imputed data down to the number 
of available cases in each group, and treating these data as if they were observed.

{phang}
{cmd:w4} specifies that standard errors be computed algebraically, conditional
on the IMORs. Note that conditioning on the IMORs is not strictly correct for
schemes including ICA-pE or ICA-pC, but the conditional standard errors appear
to be more realistic than the unconditional standard errors in this setting.


{title:Imputation debugging options}{marker Imputationdebuggingoptions}

{phang}
{cmd:listnum} lists the reason counts for each study implied by the {cmd:ica0} etc. options.

{phang}
{cmd:listall} lists the reason counts for each study after scaling to match the number 
of missing values and imputing missing for studies with no reasons.

{phang}
{cmd:listp} lists the imputed probabilities for each study.


{title:Bayes options}{marker Bayesoptions}

{phang}
{cmd:logimor()} sets the prior means for the log IMOR in the experimental and control groups. 
Both values default to 0. 
The arguments must be expressions, not variables.

{phang}
{cmd:imor()} is an alternative to {cmd:logimor()}.
It sets the exponentiated values of the prior means for the log IMOR in the experimental and control groups. 
Both values default to 1. 
The arguments must be expressions, not variables.

{phang}
{cmd:sdlogimor()} sets the prior standard deviation for the log IMOR in the experimental and
control groups. Both values default to 0.
The arguments must be expressions, not variables.

{phang}
{cmd:corrlogimor()} sets the prior correlation between the logIMORs in the experimental and
control groups. Default is 0.
The argument must be an expression, not a variable.

{phang}
{cmd:method(MC|Taylor|GH)}. {cmd: method(GH)} uses 2-dimensional Gauss-Hermite quadrature to integrate over the distribution of the IMORs and is the recommended 
method (and the default). {cmd: method(MC)} performs a full Bayesian analysis by sampling 
directly from the posterior. This is time-consuming, so dots display progress, and you can request 
more than one of the measures RR, OR, RD. {cmd: method(Taylor)} uses a Taylor series approximation 
as in section 4 of {help metamiss##ForsterSmith:Forster and Smith (1998)} and is faster than the default but inaccurate 
for wide log IMOR distributions. 

{phang}
{cmd:nip(#)} specifies the number of integration points under {cmd:method(gh)}. Default is 10.

{phang}
{cmd:reps(#)} specifies the number of MC draws under {cmd:method(mc)}. Default is 100.

{phang}
{cmd:missprior(# # [# #])} and {cmd:respprior(# # [# #])} specify beta(#,#) priors for alpha and pi under {cmd:method(mc)}. 
The 3rd and 4th arguments, if present, apply to the control group; otherwise the 1st and 2nd arguments are applied to both groups. Defaults are beta(1,1).

{phang}
{cmd:nodots} suppresses the dots that are displayed to mark the number of MC draws completed.


{title:Saved variables}{marker Savedvariables}

{pstd}
{cmd:metamiss} saves variables in the same way as {help metan}: _ES _selogES etc. _SS, 
the sample size, excludes the missing values, but an additional variable _SSmiss 
gives the total number of missing values.

{pstd}
When {cmd:method(mc)} is run, the following variables are saved for each {it:measure}:
the available cases analysis estimate ESTRAW_{it:measure}, 
the available cases analysis variance VARRAW_{it:measure}, the corrected 
estimate ESTSTAR_{it:measure}, and the corrected variance VARSTAR_{it:measure}.


{title:Examples}{marker Examples}

We analyse the haloperidol data as in {help metamiss##citethis:White and Higgins (2009)}.

{phang}{stata use haloperidol, clear}

The data set contains the key variables: author r1 f1 m1 r2 f2 m2

Available case analysis (three equivalent commands):

{phang}{stata metan r1 f1 r2 f2, rr fixedi label(namevar=author)}

{phang}{stata metamiss r1 f1 m1 r2 f2 m2, rr id(author) aca}

{phang}{stata metamiss r1 f1 m1 r2 f2 m2, rr id(author) imor(1)}

ICA-0, impute missing as zeroes (three equivalent commands):

{phang}{stata metamiss r1 f1 m1 r2 f2 m2, rr id(author) ica0 w4}

{phang}{stata metamiss r1 f1 m1 r2 f2 m2, rr id(author) ica0(m1 m2) w4}

{phang}{stata metamiss r1 f1 m1 r2 f2 m2, rr id(author) imor(0)}

Impute using reasons for missingness:

{phang}{stata metamiss r1 f1 m1 r2 f2 m2, fixed id(author) ica0(df1 df2) ica1(ds1 ds2) icapc(dc1 dc2) icap(dg1 dg2) w4}

Fixed equal IMORs (two equivalent commands):

{phang}{stata metamiss r1 f1 m1 r2 f2 m2, rr id(author) imor(2) nograph}

{phang}{stata metamiss r1 f1 m1 r2 f2 m2, rr id(author) logimor(log(2)) nograph}

Fixed opposite IMORs:

{phang}{stata metamiss r1 f1 m1 r2 f2 m2, rr id(author) imor(2 1/2) nograph}

Random equal IMORs:

{phang}{stata metamiss r1 f1 m1 r2 f2 m2, rr id(author) logimor(0) sdlogimor(1) corrlogimor(1)}

Random uncorrelated IMORs:

{phang}{stata metamiss r1 f1 m1 r2 f2 m2, rr id(author) logimor(0) sdlogimor(1) corrlogimor(0)}

{p}Possible ways to improve - unlikely to make much difference in practice:

{phang}{stata metamiss r1 f1 m1 r2 f2 m2, rr id(author) logimor(0) sdlogimor(1) corrlogimor(1) method(mc) reps(10000)}

{phang}{stata metamiss r1 f1 m1 r2 f2 m2, rr id(author) logimor(0) sdlogimor(1) corrlogimor(0) method(gh) nip(50)}


{title:Note}{marker Note}

{pstd}Different methods may give slightly different answers due to 
different handling of zero cells. For example, 

{phang}{cmd:. metamiss r1 f1 m1 r2 f2 m2, ica1 w4}

{pstd}and 

{phang}{cmd:. metamiss r1 f1 m1 r2 f2 m2, logimor(99) w4}

{pstd} differ slightly in the haloperidol data: the logimor analysis adds 1/2 to r1, f1, r2 and f2 for 6 studies with r2==0, whereas the ica0 analysis 
only does this for 3 studies with r2+m2==0. 


{title:References}{marker References}

{pstd}{marker ForsterSmith}Forster JJ, Smith PWF (1998).  Model-based inference for categorical survey data 
subject to non-ignorable non-response.  Journal of the Royal Statistical 
Society (B) 60: 57-70. 
{browse "https://doi.org/10.1111/1467-9868.00108"}

{pstd}{marker GambleHollis}Gamble C, Hollis S (2005).  Uncertainty method improved on best-worst case analysis 
in a binary meta-analysis.  Journal of Clinical Epidemiology 58: 579-588. 
{browse "https://doi.org/10.1016/j.jclinepi.2004.09.013"}

{pstd}{marker Higgins++08}Higgins JPT, White IR, Wood A (2008).  
Imputation methods for missing outcome data in meta-analysis of clinical trials.  
Clinical Trials 2008; 5: 225-239. 
{browse "https://doi.org/10.1177/1740774508091600"}

{pstd}{marker White++08}White IR, Higgins JPT, Wood AM (2008).
Allowing for uncertainty due to missing data in meta-analysis. 
1. Two-stage methods.  Statistics in Medicine 27: 711-727.
{browse "http://doi.wiley.com/10.1002/sim.3008"}

{pstd}{marker WhiteHiggins}White IR, Higgins JPT (2016). Meta-Analysis with missing data. 
In: Palmer TM, Sterne JAC, editors. Meta-Analysis in Stata: An Updated Collection from the Stata Journal. 2nd ed. 
{browse "https://www.stata.com/bookstore/meta-analysis-in-stata/"}

{pstd}{marker Chaimani++}Chaimani A, Mavridis D, Higgins J, Salanti G, White IR (2018). 
Allowing for informative missingness in aggregate data meta-analysis with continuous or binary outcomes: extensions to metamiss. 
Stata J. 3:716–40. 
{browse "https://www.stata-journal.com/article.html?article=st0540"}


{title:Please cite this program}{marker Pleasecitethisprogram}

{phang}If you find this command useful, please cite it as below (and please check for updates): 

{pstd}{marker citethis}White IR, Higgins JPT. Meta-analysis with missing data. Stata J. 2009;9:57–69. 
{browse "http://www.stata-journal.com/article.html?article=st0157"}

{phang}The lead author (below) would love to hear from you about how you are using the package.


{title:Authors and Updates}{marker Authors}

{pstd}Ian White, MRC Clinical Trials Unit at UCL, London, UK. 

{pstd}Julian Higgins, School of Social and Community Medicine, University of Bristol, UK.

{pstd}We thank David Fisher for help in making this version compatible with version 4 of {help metan},
and Ella Marley-Zagar for help with testing.

{pstd}You can get the latest version of this package by visiting {browse "https://github.com/UCL/metamiss"}
or within Stata by running 
{stata "net from https://raw.githubusercontent.com/UCL/metamiss/master"}.

{pstd}You can browse my other Stata software using 
{stata "net from http://www.homepages.ucl.ac.uk/~rmjwiww/stata/"}.

{pstd}For assistance please email {browse "mailto:ian.white@ucl.ac.uk":ian.white@ucl.ac.uk}.


{title:See also}{marker Seealso}

{help metan} if installed.
{help metamiss2} if installed.
