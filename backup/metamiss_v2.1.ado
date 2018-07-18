*! version 2.1   7feb2007 - tidier output
* version 2.0   10jul2006
* version 1.1   28oct2002
* version 1.2   02dec2002 - uses new _sensmiss.ado
* version 1.3   03jan2003 - weights disabled in metamiss, debug enabled
* version 1.3.1 09jan03 - replace option
* version 1.4   5feb2003 - can specify measures, prefix and suffix
* version 1.5   9may2005 - upgrade to stata 8; rename phi (logIMOR) as delta
* version 2.0   10jul2006 - options specified as or|rr|rd not measures(logRR|logOR|RD); variables not saved by default; with 1 measure, finally invokes metan; >1 measure allowed only with method(mc) (makes computational sense); replace option deleted - potentially confusing; help file updated

prog def metamiss
version 8
* metamiss rE fE mE rC fC mC, deltacorr(0.5) deltasd(1.22)
* currently assumes delta's are uncorrelated between studies
* mode 1: only 1 measure
* mode 2: method(mc), more than 1 measure
/* TO DO LIST:
1. SORT OUT NEW VARIABLES FROM SAVE
SORT OUT NEW VARIABLES FROM METAN
    NEW VARIABLES SHOULD BE THE SAME AS FOR METAN: 
         _SS             Sample size
         _ES             RD
         _seES           se(RD)
         _LCI            Lower CI (RD)
         _UCI            Upper CI (RD)
         _WT             M-H weight
2. Allow for zeroes - I think this should be as -metan- does
3. Edit _sensmiss to run over all studies at once
4. Don't print method when sd=0
NOTES
log - gives metan results on logRR/logOR scale
agrees with metan, fixedi
can use metan options including random[i], by()
*/

syntax varlist(min=6 max=6), [id(varname) debug or rr rd method(string) PREfix(string) SUFfix(string) PRInt save(string) log                 /* _sensmiss options start here
*/ DELTAMean(passthru) DELTACorr(passthru) DELTASd(passthru) /* nonresponse options
*/ MISSprior(passthru) RESPprior(passthru)          /* priors for nuisance parameters
*/ reps(int -1) nip(int 10)    /* method options
*/ details PROGress TIMEr table *]
local sensmissoptions `deltamean' `deltacorr' `deltasd' `missprior' `respprior' method(`method') reps(`reps') nip(`nip') `details' `progress' `timer' `table' 

tokenize "`varlist'"
local rE `1'
local fE `2'
local mE `3'
local rC `4'
local fC `5'
local mC `6'
* SORT OUT METHOD
local method = lower("`method'")
if "`method'"=="" {
   if `reps'~=-1 local method mc
   else local method gh
}
if `reps'==-1 local reps 100
if "`method'" != "mc" & "`method'" != "taylor" & "`method'" != "gh" {
   di as error "Please specify method(MC|Taylor|GH)"
}
* SORT OUT MEASURES
if "`method'"=="mc" & "`rr'`or'`rd'"=="" {
    local rr rr
    local or or
    local rd rd
}
if "`rr'"=="rr" local measures logRR
if "`or'"=="or" local measures `measures' logOR
if "`rd'"=="rd" local measures `measures' RD
local measlist = upper("`rr' `or' `rd'")
local nmeasures = ("`or'"=="or") + ("`rr'"=="rr") + ("`rd'"=="rd")
if `nmeasures'==0 {
    di in red "Please specify or, rr or rd"
    exit 498
}
if "`method'"=="gh" di as text "Method: Gauss-Hermite quadrature (`nip' integration points)"
if "`method'"=="mc" di as text "Method: Monte Carlo integration (`reps' replications)"
if "`method'"=="taylor" di as text "Method: Taylor series approximation"

if `nmeasures'==1 {
    local mode 1
    di as text "Measure: `measlist'"
}
if `nmeasures'>1 {
    if "`method'"=="mc" {
        if "`save'"=="" local save star varstar
        local mode 2
        di as text "Measures: `measlist'"
        di as text "Estimates and standard errors will be saved for future meta-analysis"
    }
    else {
       di in red "Please specify only one of or, rr or rd"
       exit 498
   }
}
local typelist raw varraw star varstar varimor varobs
* SORT OUT OTHER OPTIONS
* note save can contain any or all of raw varraw star varstar varobs varimor 
if "`save'"=="" & "`prefix'`suffix'"~="" local save raw varraw star varstar varobs varimor 
if "`debug'"=="debug" local dicmd dicmd
if "`log'"=="log" & (`mode'==2 | "`measlist'"=="RD") {
    di as text "(log option ignored)"
}
* END OF PARSING

* SET UP
foreach meas in `measures' {
   foreach type in `save' {
      confirm new var `prefix'`meas'`type'`suffix'
   }
   foreach type in `typelist' {
      tempvar `meas'`type'
      qui gen ``meas'`type'' = .
      label var ``meas'`type'' "`meas'`type'"
   }
}

* RUN SENSMISS FOR EACH STUDY
local n = _N
forvalues i = 1/`n' {
   if "`method'"=="mc" di "Study `i'... " _c
   local printprior = cond(`i'==1,"printprior","")
   local values
   foreach var of var `rE' `fE' `mE' `rC' `fC' `mC' {
      local val = `var'[`i']
      local values `values' `val'
   }
   if "`print'"=="print" | "`debug'"=="debug" {
      if "`id'"~="" di _new in green "`id' = " in yellow `id'[`i']
      else          di _new in green "i = " in yellow `i'
   }
   if "`print'"!="print" local print noprint
   `dicmd' _sensmiss `values', measures(`measures') `sensmissoptions' `print' `printprior'
   foreach meas in `measures' {
      foreach type in `typelist' {
         qui replace ``meas'`type'' = r(`meas'`type') in `i'
      }
   }
}
if "`method'"=="mc" di

* RUN METAN
if `mode'==1 {
    tempvar stderr
    qui gen `stderr'=.
    cap drop _SS
    cap drop _SSmiss
    qui gen _SS = `rE'+`fE'+`rC'+`fC'
    qui gen _SSmiss = `mE'+`mC'
    label var _SS "Sample size (observed)"
    label var _SSmiss "Sample size (missing)"
    * NB only 1 measure
    qui replace `stderr'=sqrt(``measures'varstar')
    if "`measures'"!="RD" & "`log'"~="log" local eform eform
    local idopt = cond("`id'"~="","label(namevar=`id')","")
    if "`debug'"=="debug" di
    `dicmd' metan ``measures'star' `stderr', `idopt' `options' `eform'
    qui gen _ES = ``measures'star'
    if "`eform'"=="eform" qui replace _ES=exp(_ES)
    label var _ES "`log'`measlist'"
    qui gen _se`log'ES = sqrt(``measures'varstar')
    label var _se`log'ES "se(`measures')"
    move _ES _LCI
    move _se`log'ES _LCI
    `pause'
}

* SAVE VARIABLES
foreach meas in `measures' {
   foreach type in `save' {
      gen `prefix'`meas'`type'`suffix' = ``meas'`type''
   }
}

end



prog def gauher
syntax newvarlist(min=2 max=2), n(int) [scaled maxit(int 1000) tol(int 10) trace]
tokenize "`varlist'"
local genx `1'
local genw `2'
* output variables w, x
* n=number of points
* w=weights
* x=gamma

qui count
if r(N)<`n' {
   set obs `n'
}
local pim4 .751125544464925
tempvar x w
qui gen `x'=.
qui gen `w'=.
local m=(`n'+1)/2

forvalues i=1/`m' {

   if "`trace'"=="trace" {
      di "Starting i=`i'"
   }
   if `i'==1 {
      local z=sqrt(float(2*`n'+1))-1.85575*(2*`n'+1)^(-0.16667)
   }
   else if `i'==2 {
      local z=`z'-1.14*`n'^0.426/`z'
   }
   else if `i'==3 {
      local z=1.86*`z'-0.86*`x'[1]
   }
   else if `i'==4 {
      local z=1.91*`z'-0.91*`x'[2]
   }
   else {
      local z=2*`z'-`x'[`i'-2]
   }

   local stop 0
   local iter 0
   while `stop'<1 {
      local ++iter
      local p1=`pim4'
      local p2=0
      forvalues j=1/`n' {
         local p3=`p2'
         local p2=`p1'
         local p1=`z'*sqrt(2/`j')*`p2'-sqrt((`j'-1)/`j')*`p3'
      }
      local pp=sqrt(2*`n')*`p2'
      local z1=`z'
      local z=`z1'-`p1'/`pp'
      local stop = abs(`z'-`z1')<10^(-`tol')
      if "`trace'"=="trace" {
         di "Iteration `iter', change = " `z'-`z1'
      }
      if `iter'>`maxit' {
         di as error "Too many iterations for i=`i'"
         exit 498
      }
   }
   if "`trace'"=="trace" {
      di "i=`i' required `iter' iterations"   
   }
   qui replace `x'=`z' if _n==`i'
   qui replace `x'=-`z' if _n==`n'+1-`i'
   qui replace `w'=2/(`pp'*`pp') if _n==`i' | _n==`n'+1-`i'
      
}

qui replace `x'=0 if _n==`m' & abs(`x')<1E-10

if "`scaled'"=="scaled" {
   gen `genx'=`x'*sqrt(2)
   gen `genw'=`w'/sqrt(_pi)
   label var `genx' "GH quadrature point, scaled to SD 1"
   label var `genw' "GH quadrature weight, scaled to sum 1"
}
else {
   gen `genx'=`x'
   gen `genw'=`w'
   label var `genx' "GH quadrature point, original scale"
   label var `genw' "GH quadrature weight, original scale"
}
sort `genx'

end


prog def random
* version 2.6  06jul2005 - lp option for dist(bernoulli)
version 7

/* examples of use:
random x, mean(6) dist(poisson)
random x, mean(1) var(0.25) dist(gamma)
random x, mean(1) var(0.25) dist(normal)
random x, mean(1) var(0.25) dist(binomial)

To do:
allow 2nd syntax like dist(gamma a l)
hence also allow dist(beta a b)
*/
local opts Mean(string) Variance(string) SD(string) Dist(string) LP(string) CHeck 
cap syntax varlist(max=1) [if] [in], [REPlace `opts']
if _rc==0 { /* existing variable */
  if "`replace'"=="replace" {
    local origvar `varlist'
    tempvar varlist
  }
  else {
    di in red "`varlist' already defined"
    exit 110
  }
}
else {
  syntax newvarlist(max=1) [if] [in], [`opts']
}

tempvar touse meanval varval
mark `touse' `if' `in'

foreach distribution in normal poisson bernoulli gamma binomial {
   if upper(substr("`dist'",1,3))==upper(substr("`distribution'",1,3)) {
      local dist `distribution'
   }
}

if "`dist'"=="" { 
   if "`lp'"~="" {
      local dist bernoulli
   }
   else {
      local dist normal
   }
   di as text "Assuming dist(`dist')"
}

if "`mean'"=="" {
   if "`dist'"=="normal" { 
      gen `meanval' = 0
      di as text "Taking mean = 0"   
   }
   else if "`dist'"=="bernoulli" { 
      cap gen `meanval' = 1/(1+exp(-(`lp')))
      if _rc {
         di as error "Can't use lp(`lp')"
         exit 498
      }
   }
   else {
      di as error "mean() must be specified"
      exit 498
   }
}
else {gen `meanval' = `mean'}

if "`variance'"=="" {
    if "`sd'"~="" { 
        cap assert `sd'>=0
        if _rc>0 {
            di in red "SD must be non-negative"
            exit 498
        }
        gen `varval'=(`sd')^2
    }
    else if "`dist'"=="normal" { 
        gen `varval' = 1
        di as text "Taking variance 1"
    }
    else if "`dist'"=="bernoulli" { 
       gen `varval'=`meanval'*(1-`meanval')
    }
    else {
        di in red "var() must be specified"
        exit 498
    }
}
else {
   cap assert `variance'>=0
   if _rc>0 {
      di in red "Variance must be non-negative"
      exit 498
   }
   gen `varval' = `variance'
}

local mean `meanval'
local variance `varval'

if "`dist'"=="normal" {
   qui gen `varlist' = `mean' + sqrt(`variance')*invnorm(uniform()) if `touse'
}
   
else if "`dist'"=="poisson" {
   if "`variance'"~="" {di "var() ignored"}
   cap assert `mean'>=0
   if _rc>0 {
     di in red "Mean must be non-negative"
     exit 498
   }
   tempvar uniform 
   gen `uniform' = uniform() if `touse'
   local i = 0
   qui gen `varlist'=. if `touse'
   while "`stop'"~="stop" {
      qui replace `varlist' = `i' if (`varlist'==.) & (1-gammap(`i'+1,`mean') > `uniform') & `touse'
      local i = `i'+1
      qui count if `varlist'==.
      if r(N)==0 {local stop stop}
      }
}
   
else if "`dist'"=="gamma" {
   cap assert `mean'>=0
   if _rc>0 {
     di in red "Mean must be non-negative"
     exit 498
   }
   * mean = alpha/beta, var = alpha/beta^2
   tempvar alpha beta
   gen `beta' = `mean'/`variance' if `touse'
   gen `alpha' = `mean'*`beta' if `touse'
   qui gen `varlist' = invgammap(`alpha',uniform())/`beta' if `touse'
}

else if "`dist'"=="binomial" {
   cap assert `mean'>=0
   if _rc>0 {
     di in red "Mean must be non-negative"
     exit 498
   }
   local p = 1-`variance'/`mean'
   local n = `mean'/`p'
   tempvar uniform 
   gen `uniform' = uniform() if `touse'
   local i = 0
   qui gen `varlist'=. if `touse'
   while "`stop'"~="stop" {
      qui replace `varlist' = `i' if (`varlist'==.) & (1-Binomial(`n',`i'+1,`mean'/`n') > `uniform') & `touse' 
      local i = `i'+1
      qui count if `varlist'==.
      if r(n)==0 {local stop stop}
      }
}
   
else if "`dist'"=="bernoulli" {
   cap assert `mean'>=0
   if _rc>0 {
     di in red "Mean must be non-negative"
     exit 498
   }
   qui gen `varlist'=uniform()<`mean' if `touse'
}
   
else {
   di in red "Must specify dist(normal|poisson|gamma|binomial|bernoulli)"
   exit 498
}

if "`replace'"=="replace" {
  cap drop `origvar'
  gen `origvar' = `varlist' if `touse'
}



if "`check'"=="check" {
    summ `varlist' if `touse'
    gra `varlist' if `touse'
    }

end


prog def dicmd
* works for version (6?) 7 or 8
noi di as text "command is: " as result `"`0'"'
`0'
end
