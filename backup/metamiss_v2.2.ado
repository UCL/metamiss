*! version 2.2   9feb2007 - uses _sensmiss2
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
2. Allow for zeroes - I think this should be as -metan- does
3. Edit _sensmiss to run over all studies at once
4. Don't print method when sd=0
NOTES
log - gives metan results on logRR/logOR scale
agrees with metan, fixedi
can use metan options including random[i], by()
*/

syntax varlist(min=6 max=6), [                          ///
    or rr rd                                            /// meta-analysis options
    DELTAMean(string) DELTACorr(string) DELTASd(string) /// nonresponse options
    id(varname) log                                     /// reporting options
    debug PROGress TIMEr                                /// debugging options
    method(string) reps(int -1) nip(int 10)             /// method options
    MISSprior(string) RESPprior(string)                 /// priors for nuisance parameters 
    PREfix(string) SUFfix(string) PRInt save(string)    /// output options
    *]

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
if "`method'"=="mc" {
    if "`rd'"=="rd" local measure RD
    if "`rr'"=="rr" local measure `measure' logRR
    if "`or'"=="or" local measure `measure' logOR
    if "`measure'"=="" local measure RD logRR logOR
}
else {
   local measure = upper("`or'`rr'`rd'")
   if length("`measure'")>2 {
      di as error "Please specify only one of or, rr, rd"
      exit 498
   }
   if length("`measure'")==0 local measure RR
}
local nmeasures = wordcount("`measure'")

* PARSE PRIORS FOR INFORMATIVENESS PARAMETERS delta
tokenize "`deltamean'"
if "`1'"=="" local deltamE 0
else local deltamE `1'
if "`2'"=="" local deltamC `deltamE'
else local deltamC `2'
tokenize "`deltasd'"
if "`1'"=="" local deltasE 0
else local deltasE `1'
if "`2'"=="" local deltasC `deltasE'
else local deltasC `2'
if "`deltacorr'"=="" local deltacorr 1
assert `deltacorr'>=-1 & `deltacorr'<=1

* END OF PARSING

* SORT OUT DATA INCLUDING ZEROES
tokenize "`varlist'"
tempvar anyzero rE fE mE rC fC mC
gen `anyzero' = min(`1', `2', `4', `5')==0
gen `rE' = `1'+`anyzero'/2
gen `fE' = `2'+`anyzero'/2
gen `mE' = `3'+`anyzero'/2
gen `rC' = `4'+`anyzero'/2
gen `fC' = `5'+`anyzero'/2
gen `mC' = `6'+`anyzero'/2

* PRINT DETAILS
noi di in green "Meta-analysis allowing for missing data."
di as text "Priors used:  Arm 1: N(`deltamE',`deltasE'^2). Arm 2: N(`deltamC',`deltasC'^2). Corr: `deltacorr'."

***************** TAYLOR AND GH METHODS **********************
qui {
* BASICS AND VAROBS
if "`method'"=="taylor" | "`method'"=="gh" {
   tempvar eststar varobs varimor sestar 
   foreach R in E C {
      * "star" indicates values deltamE, deltamC
      tempvar  imor`R' alpha`R' piobs`R' pimiss`R' pistar`R' dpistar`R'_dpi`R' dpistar`R'_dalpha`R' varpiobs`R' varalpha`R' varpistar`R' K`R' 
      gen `imor`R'' = exp(`deltam`R'')
      gen `alpha`R'' = `m`R''/(`r`R''+`f`R''+`m`R'')
      gen `piobs`R'' = `r`R''/(`r`R''+`f`R'')
      gen `pimiss`R'' = `imor`R''*`piobs`R''/(`imor`R''*`piobs`R''+1-`piobs`R'')
      gen `pistar`R'' = (1-`alpha`R'')*`piobs`R'' + `alpha`R''*`pimiss`R''
      gen `dpistar`R'_dpi`R'' = 1-`alpha`R''+`alpha`R''*`imor`R''/(`imor`R''*`piobs`R''+1-`piobs`R'')^2
      gen `dpistar`R'_dalpha`R'' = (`imor`R''-1)*`piobs`R''*(1-`piobs`R'')/(`imor`R''*`piobs`R''+1-`piobs`R'')
      gen `varpiobs`R'' = `piobs`R''*(1-`piobs`R'')/(`r`R''+`f`R'')
      gen `varalpha`R'' = `alpha`R''*(1-`alpha`R'')/(`r`R''+`f`R''+`m`R'')
      gen `varpistar`R'' = (`dpistar`R'_dpi`R'')^2 * `varpiobs`R'' + (`dpistar`R'_dalpha`R'')^2 * `varalpha`R'' 
      gen `K`R'' = `alpha`R''*`imor`R''/(`imor`R''*`piobs`R''+1-`piobs`R'')^2
   }

   if "`measure'"=="RD" gen `varobs' = `varpistarE' + `varpistarC'
   if "`measure'"=="RR" gen `varobs' = `varpistarE'/`pistarE'^2 + `varpistarC'/`pistarC'^2 
   if "`measure'"=="OR" gen `varobs' = `varpistarE'/(`pistarE'*(1-`pistarE'))^2 + `varpistarC'/(`pistarC'*(1-`pistarC'))^2
}

* TAYLOR METHOD
if "`method'"=="taylor" {
      if "`print'"~="noprint" {
         noi di in green "Method: Taylor series approximation."
      }
      if "`measure'"=="RD" {
        gen `eststar' = `pistarE' - `pistarC'
        gen `varimor' = (`KE'*`piobsE'*(1-`piobsE')*`deltasE')^2 + (`KC'*`piobsC'*(1-`piobsC')*`deltasC')^2 - 2*`KE'*`KC'*`piobsE'*(1-`piobsE')*`piobsC'*(1-`piobsC')*`deltasE'*`deltasC'*`deltacorr'
      }
      if "`measure'"=="RR" {
        gen `eststar' = log(`pistarE' / `pistarC')
        gen `varimor' = (`KE'*(1-`piobsE')*`deltasE')^2 + (`KC'*(1-`piobsC')*`deltasC')^2 - 2*`KE'*`KC'*(1-`piobsE')*(1-`piobsC')*`deltasE'*`deltasC'*`deltacorr'
      }
      if "`measure'"=="OR" {
        gen `eststar' = log(`pistarE'/(1-`pistarE')) - log(`pistarC'/(1-`pistarC'))
        gen `varimor' = (`KE'*`deltasE')^2 + (`KC'*`deltasC')^2 - 2*`KE'*`KC'*`deltasE'*`deltasC'*`deltacorr'
      }
}

* GAUSS-HERMITE METHOD 
if "`method'"=="gh" {
   if "`print'"~="noprint" {
      noi di in green "Method: Gauss-Hermite quadrature (`nip' integration points)."
   }
   if "`progress'"=="progress" noi di as text "Performing GH quadrature"
   tempvar a1 a2
   gen `a1'=sqrt((1+`deltacorr')/2)
   gen `a2'=sqrt((1-`deltacorr')/2)
   local imorE exp(`deltamE'+`deltasE'*(`a1'*Z1+`a2'*Z2))
   local imorC exp(`deltamC'+`deltasC'*(`a1'*Z1-`a2'*Z2))
   foreach R in E C {
      local pimiss`R' (`piobs`R''*`imor`R'' / (1-`piobs`R''+`piobs`R''*`imor`R''))
      local pistar`R' ((1-`alpha`R'')*`piobs`R'' + `alpha`R''*`pimiss`R'')
   }
   if "`measure'"=="RD" local formula `pistarE'-`pistarC'
   if "`measure'"=="RR" local formula log(`pistarE') - log(`pistarC')
   if "`measure'"=="OR" local formula log(`pistarE'/(1-`pistarE')) - log(`pistarC'/(1-`pistarC'))
   expectn `formula', gen(`eststar') genvar(`varimor') token(Z1 Z2) `debug'
}

gen `sestar' = sqrt(`varobs' + `varimor')
}

***************** MONTE CARLO METHOD **********************
qui {
if "`method'"=="mc" {
   noi di in green "Method: Monte Carlo (`reps' draws)."
   foreach meas in `measure' {
      foreach adj in RAW STAR {
         foreach out in EST VAR {
            if "`replace'"=="" confirm new var `out'`adj'_`meas' 
            else cap drop `out'`adj'_`meas' 
         }
      }
   }
   noi di "Measures used: `measure'"

   if `reps'==-1 local reps 1000
 
   * PARSE PRIORS FOR P(MISSING)
   * note beta(c,d) with mean c/(c+d) is like c 1s and d 0s
   tokenize "`missprior'"
   if "`1'"=="" {
      local aE1 1
      local aE0 1
   }
   else {
      local aE1 `1'
      local aE0 `2'
   }   
   if "`3'"=="" {
      local aC1 `aE1'
      local aC0 `aE0'
   }
   else {
      local aC1 `3'
      local aC0 `4'
   }
 
   * PARSE PRIORS FOR P(OUTCOME|NON-MISSING)
   tokenize "`respprior'"
   if "`1'"=="" {
      local bE1 1
      local bE0 1
   }
   else {
      local bE1 `1'
      local bE0 `2'
   }   
   if "`3'"=="" {
      local bC1 `bE1'
      local bC0 `bE0'
   }
   else {
      local bC1 `3'
      local bC0 `4'
   }   
 
   * INITIALISE
   if "`timer'"=="timer" {
      now
      local start=r(datetime)
   }
 
   * COMPUTE DRAWS FROM POSTERIOR
   qui {
   tempvar deltaE deltaC x0 x1 alphaE alphaC piobsE piobsC pistarE pistarC RAW_RD STAR_RD RAW_logRR STAR_logRR RAW_logOR STAR_logOR
      
   foreach meas in `measure' {
      foreach adj in RAW STAR {
         tempvar sum`adj'_`meas' sum`adj'_`meas'2 
         gen `sum`adj'_`meas''  = 0
         gen `sum`adj'_`meas'2' = 0
      }
   }
   forvalues i=1/`reps' {
      noi dotter `i' `reps'
      if "`progress'"=="progress" {
         noi di "Drawing deltaE " _c
      }
      random `deltaE', dist(normal) 
      if "`progress'"=="progress" {
         noi di " `deltaC' " _c
      }
      random `deltaC', dist(normal) mean(`deltacorr'*`deltaE') var(1-`deltacorr'^2) 
      foreach i in E C {
         replace `delta`i'' = `deltam`i'' + `delta`i''*`deltas`i''

         if "`progress'"=="progress" {
            noi di " alpha`i' " _c
         }
         random `x1', dist(gamma) mean(`a`i'1'+`m`i'')        var(`a`i'1'+`m`i'')        
         random `x0', dist(gamma) mean(`a`i'0'+`f`i''+`r`i'') var(`a`i'0'+`f`i''+`r`i'') 
         gen `alpha`i''=`x1'/(`x0'+`x1')
         drop `x0' `x1'
     
         if "`progress'"=="progress" {
            noi di " piobs`i' " _c
         }
         random `x1', dist(gamma) mean(`b`i'1'+`r`i'') var(`b`i'1'+`r`i'') 
         random `x0', dist(gamma) mean(`b`i'0'+`f`i'') var(`b`i'0'+`f`i'') 
         gen `piobs`i''=`x1'/(`x0'+`x1')
         drop `x0' `x1'
   
         gen `pistar`i''=(1-`alpha`i'')*`piobs`i'' + `alpha`i''*exp(`delta`i'')*`piobs`i''/(1-`piobs`i''+exp(`delta`i'')*`piobs`i'')
         drop `delta`i'' `alpha`i'' 
      }
  
      gen `RAW_RD'=`piobsE'-`piobsC'
      gen `STAR_RD'=`pistarE'-`pistarC'
  
      gen `RAW_logRR'=log(`piobsE'/`piobsC')
      gen `STAR_logRR'=log(`pistarE'/`pistarC')
  
      gen `RAW_logOR' = log(`piobsE'/(1-`piobsE')) - log(`piobsC'/(1-`piobsC'))
      gen `STAR_logOR'= log(`pistarE'/(1-`pistarE')) - log(`pistarC'/(1-`pistarC'))
      
      drop `piobsE' `piobsC' `pistarE' `pistarC'
      
      foreach meas in `measure' {
         foreach adj in RAW STAR {
            replace `sum`adj'_`meas''  = `sum`adj'_`meas'' + ``adj'_`meas''
            replace `sum`adj'_`meas'2' = `sum`adj'_`meas'2' + ``adj'_`meas''^2
         }
      }
      
      drop `RAW_RD' `STAR_RD' `RAW_logRR' `STAR_logRR' `RAW_logOR' `STAR_logOR'

      if "`progress'"=="progress" {
         noi di " done." _newline
      }
   }
   }
   * SUMMARISE POSTERIOR
   foreach meas in `measure' {
      foreach adj in RAW STAR {
         gen EST`adj'_`meas' = `sum`adj'_`meas''/`reps'
         gen VAR`adj'_`meas' = `sum`adj'_`meas'2'/`reps' - EST`adj'_`meas'^2
      }
   }
   if "`timer'"=="timer" {
      now
      local elapse=r(datetime)-`start'
      noi di in green _newline "Elapsed time: " in yellow `elapse' in green " seconds"
   }
} 

} /* end of quietly */

********************************** END OF MONTE CARLO PART **********************************

* RUN METAN
if `nmeasures'==1 {
    cap drop _SS
    cap drop _SSmiss
    qui gen _SS = `rE'+`fE'+`rC'+`fC'
    qui gen _SSmiss = `mE'+`mC'
    label var _SS "Sample size (observed)"
    label var _SSmiss "Sample size (missing)"
    if "`measure'"=="RD" local semeasurename `measure'
    else local semeasurename log`measure'
    if "`measure'"!="RD" & "`log'"~="log" local eform eform
    if "`measure'"!="RD" & "`log'"=="log" local measurename log`measure'        
    else local measurename `measure'        
    local idopt = cond("`id'"~="","label(namevar=`id')","")
    metan `eststar' `sestar', `idopt' `options' `eform' nograph
    qui gen _ES = `eststar'
    if "`eform'"=="eform" qui replace _ES=exp(_ES)
    qui gen _se`log'ES = `sestar'
    label var _ES        "`measurename'"
    label var _LCI       "Lower CI (`measurename')"
    label var _UCI       "Upper CI (`measurename')"
    label var _se`log'ES "se(`semeasurename')"
    move _ES _LCI
    move _se`log'ES _LCI
}

end

******************************************************************************

* NB also requries random.ado, expectn.ado
