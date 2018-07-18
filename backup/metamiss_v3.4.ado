*! version 3.4  4sep2007 - new syntax allowing ica0 or ica0(var var) etc.; icab/icaw the right way round; big rewrite of Julian's part
* version 3.3  24apr2007 - (Julian) changed adjustments for empty cells
* version 3.2   24apr2007 - selects by `measure' not by `rr', `or', `rd'
* version 3.1   24apr2007 - graphs fixed
* version 3.0   9feb2007 - combined with Julian's program
* version 2.2   9feb2007 - uses _sensmiss2
* version 2.1   7feb2007 - tidier output
* version 2.0   10jul2006
* version 1.1   28oct2002
* version 1.2   02dec2002 - uses new _sensmiss.ado
* version 1.3   03jan2003 - weights disabled in metamiss, debug enabled
* version 1.3.1 09jan03 - replace option
* version 1.4   5feb2003 - can specify measures, prefix and suffix
* version 1.5   9may2005 - upgrade to stata 8; rename phi (logIMOR) as delta
* version 2.0   10jul2006 - options specified as or|rr|rd not measures(logRR|logOR|RD); variables not saved by default; with 1 measure, finally invokes metan; >1 measure allowed only with method(mc) (makes computational sense); replace option deleted - potentially confusing; help file updated

prog def metamiss
version 9

/* TO DO LIST:
Check it works when all reason variables are zero.
NOTES
log - gives metan results on logRR/logOR scale
agrees with metan, fixedi
can use metan options including random[i], by()
*/

noi di as text "*******************************************************************"
noi di as text "******** METAMISS: meta-analysis allowing for missing data ********"
noi di as text "*******************************************************************"

************************************ COMMON PARSING ************************************

* DETERMINE IF OPTIONS ICA0 ETC. HAVE ARGUMENTS (SYNTAX 2) OR NOT (SYNTAX 1)
syntax varlist(min=6 max=6) [if] [in], [aca ica0(string) ica1(string) icapc(string) icape(string) icap(string) icaimor(string) *]
if "`ica0'`ica1'`icapc'`icape'`icap'`icaimor'"!="" {
    if "`aca'"=="aca" exit498 Can't have ACA and ICA
    local syntax 2
    local JHoptions ica0(string) ica1(string) icapc(string) icape(string) icap(string) icaimor(string)
}
syntax varlist(min=6 max=6) [if] [in], [ica0 ica1 icapc icape icap icaimor aca icab icaw *]
local words = wordcount("`ica0' `ica1' `icap' `icape' `icapc' `icaimor' `aca' `icab' `icaw'")
if `words' > 1 exit498 State only one ACA or ICA analysis, or several ICA analyses with arguments
if `words' == 1 {
    if "`syntax'"=="2" exit498 Can't mix syntaxes 1 and 2
    local syntax 1
    local JHoptions aca ica0 ica1 icapc icape icap icaimor 
}

* MAIN PARSING
syntax varlist(min=6 max=6) [if] [in], [                  ///
///Common options:
    or rr rd                                              /// meta-analysis options
    id(varname) log nograph                               /// reporting options
    logimor(string) imor(string)                          /// nonresponse options
///Ian's options:
    SDlogimor(string) CORRlogimor(string)                 /// nonresponse options
    method(string) reps(int -1) nip(int 10)               /// method options
    MISSprior(string) RESPprior(string) replace           /// full-bayes options
    debug TIMEr                                           /// debugging options
///Julian's options:
    `JHoptions' icab icaw                                 /// main nonresponse options
    w1 w2 w3 w4                                           /// std error options
    listp listx listnum listall                           /// debugging options
*]

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
if "`debug'"=="debug" di "Measure: `measure'"

if "`sdlogimor'"!="" local syntax 3

* default to ACA or ICA-IMOR:
if "`syntax'"=="" {
    di as text "(No method specified: " _c
    if "`logimor'`imor'"=="" {
       di as text "assuming ACA.)"
       local aca aca
    }
    else {
       di as text "assuming ICAimor.)"
       local icaimor icaimor
    }
    local syntax 1
}

* parse imor() to imorE, imorC, etc.
local parselist logimor imor 
if `syntax'==3 local parselist `parselist' sdlogimor
foreach thing in `parselist' {
   if "``thing''"!="" {
       tokenize "``thing''"
       local `thing'E `1'
       if "`2'"=="" local `thing'C `1'
       else local `thing'C `2'
       if "`3'"!="" {
           di in red "Too many arguments in `thing'"
           exit 498
       }
   }
}
if `syntax'==2 {
   foreach thing in ica0 ica1 icap icapc icape icaimor {
      if "``thing''"!="" {
          tokenize "``thing''"
          local `thing'E `1'
          if "`2'"=="" local `thing'C 0
          else local `thing'C `2'
          if "`3'"!="" {
              di in red "Too many arguments in `thing'"
              exit 498
          }
      }
   }
}

if "`imor'"!="" & "`logimor'"!="" {
    di as error "Can't specify both logimor() and imor()"
    exit 498
}
foreach arm in E C {
   if "`imor'"!="" {
      local logimor`arm'=log(`imor`arm'')
      if `imor`arm''==0 local logimor`arm' -50
   }
   if "`logimor'"!="" {
      if `logimor`arm''>50 & `logimor`arm''!=. local logimor`arm' 50
      if `logimor`arm''<-50 local logimor`arm' -50
      local imor`arm'=exp(`logimor`arm'')
      if `logimor`arm''<-50 local imor`arm' 0
   }
}

********************************* END OF COMMON PARSING *********************************

tempvar eststar sestar

********************************* IAN'S PARSING *********************************

if `syntax'==3 {

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

* PARSE PRIORS FOR INFORMATIVENESS PARAMETERS logimor
di as text "Approach: allowing for uncertainty in IMORs."
if "`logimorcorr'"=="" local logimorcorr 0
assert `logimorcorr'>=-1 & `logimorcorr'<=1

********************************* END OF IAN'S PARSING *********************************

********************************* IAN'S CALCULATION *********************************

if "`debug'"=="debug" di as text "Starting Ian's calculation"

* SORT OUT DATA INCLUDING ZEROES
tokenize "`varlist'"
tempvar anyzero rE fE mE rC fC mC
qui {
   gen `anyzero' = min(`1', `2', `4', `5')==0
   gen `rE' = `1'+`anyzero'/2
   gen `fE' = `2'+`anyzero'/2
   gen `mE' = `3'
   gen `rC' = `4'+`anyzero'/2
   gen `fC' = `5'+`anyzero'/2
   gen `mC' = `6'
   count if `anyzero'==1
}
if r(N)>0 di as text "Zero cells detected: adding 1/2 to " as result r(N) as text " studies."

* PRINT DETAILS
di as text "Priors used:  Arm 1: N(" as result `logimorE' as text "," as result `sdlogimorE' as text "^2). Arm 2: N(" as result `logimorC' as text "," as result `sdlogimorC' as text "^2). Corr: " as result `logimorcorr' as text "."

***************** TAYLOR AND GH METHODS **********************
qui {
   * BASICS AND VAROBS
   if "`method'"=="taylor" | "`method'"=="gh" {
      tempvar varobs varimor
      foreach R in E C {
         * "star" indicates values logimorE, logimorC
         tempvar imor`R' alpha`R' piobs`R' pimiss`R' pistar`R' dpistar`R'_dpi`R' dpistar`R'_dalpha`R' varpiobs`R' varalpha`R' varpistar`R' K`R'
         gen `imor`R'' = exp(`logimor`R'')
         replace `imor`R'' = 0 if `logimor`R'' < -50
         replace `imor`R'' = 9999 if `logimor`R'' > 50
         label var `imor`R'' "imor`R'"
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

      * TAYLOR METHOD
      if "`method'"=="taylor" {
            noi di in green "Method: Taylor series approximation."
            if "`measure'"=="RD" {
              gen `eststar' = `pistarE' - `pistarC'
              gen `varimor' = (`KE'*`piobsE'*(1-`piobsE')*`sdlogimorE')^2 + (`KC'*`piobsC'*(1-`piobsC')*`sdlogimorC')^2 - 2*`KE'*`KC'*`piobsE'*(1-`piobsE')*`piobsC'*(1-`piobsC')*`sdlogimorE'*`sdlogimorC'*`logimorcorr'
            }
            if "`measure'"=="RR" {
              gen `eststar' = log(`pistarE' / `pistarC')
              gen `varimor' = (`KE'*(1-`piobsE')*`sdlogimorE')^2 + (`KC'*(1-`piobsC')*`sdlogimorC')^2 - 2*`KE'*`KC'*(1-`piobsE')*(1-`piobsC')*`sdlogimorE'*`sdlogimorC'*`logimorcorr'
            }
            if "`measure'"=="OR" {
              gen `eststar' = log(`pistarE'/(1-`pistarE')) - log(`pistarC'/(1-`pistarC'))
              gen `varimor' = (`KE'*`sdlogimorE')^2 + (`KC'*`sdlogimorC')^2 - 2*`KE'*`KC'*`sdlogimorE'*`sdlogimorC'*`logimorcorr'
            }
      }

      * GAUSS-HERMITE METHOD
      if "`method'"=="gh" {
         noi di as text "Method: Gauss-Hermite quadrature (" as result `nip' as text " integration points)."
         if "`debug'"=="debug" noi di as text "Performing GH quadrature"
         tempvar a1 a2
         gen `a1'=sqrt((1+`logimorcorr')/2)
         gen `a2'=sqrt((1-`logimorcorr')/2)
         local imorE exp(`logimorE'+`sdlogimorE'*(`a1'*Z1+`a2'*Z2))
         local imorC exp(`logimorC'+`sdlogimorC'*(`a1'*Z1-`a2'*Z2))
         foreach R in E C {
            local pimiss`R' (`piobs`R''*`imor`R'' / (1-`piobs`R''+`piobs`R''*`imor`R''))
            local pistar`R' ((1-`alpha`R'')*`piobs`R'' + `alpha`R''*`pimiss`R'')
         }
         if "`measure'"=="RD" local formula `pistarE'-`pistarC'
         if "`measure'"=="RR" local formula log(`pistarE') - log(`pistarC')
         if "`measure'"=="OR" local formula log(`pistarE'/(1-`pistarE')) - log(`pistarC'/(1-`pistarC'))
         qui expectn `formula', gen(`eststar') genvar(`varimor') token(Z1 Z2) `debug'
      }
    gen `sestar' = sqrt(`varobs' + `varimor')
    }
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
   tempvar logimorE logimorC x0 x1 alphaE alphaC piobsE piobsC pistarE pistarC RAW_RD STAR_RD RAW_logRR STAR_logRR RAW_logOR STAR_logOR

   foreach meas in `measure' {
      foreach adj in RAW STAR {
         tempvar sum`adj'_`meas' sum`adj'_`meas'2
         gen `sum`adj'_`meas''  = 0
         gen `sum`adj'_`meas'2' = 0
      }
   }
   forvalues i=1/`reps' {
      noi dotter `i' `reps'
      if "`debug'"=="debug" {
         noi di "Drawing logimorE " _c
      }
      random `logimorE' `if' `in', dist(normal)
      if "`debug'"=="debug" {
         noi di " `logimorC' " _c
      }
      random `logimorC' `if' `in', dist(normal) mean(`logimorcorr'*`logimorE') var(1-`logimorcorr'^2)
      foreach i in E C {
         replace `logimor`i'' = `logimor`i'' + `logimor`i''*`sdlogimor`i''

         if "`debug'"=="debug" {
            noi di " alpha`i' " _c
         }
         random `x1' `if' `in', dist(gamma) mean(`a`i'1'+`m`i'')        var(`a`i'1'+`m`i'')
         random `x0' `if' `in', dist(gamma) mean(`a`i'0'+`f`i''+`r`i'') var(`a`i'0'+`f`i''+`r`i'')
         gen `alpha`i''=`x1'/(`x0'+`x1')
         drop `x0' `x1'

         if "`debug'"=="debug" {
            noi di " piobs`i' " _c
         }
         random `x1' `if' `in', dist(gamma) mean(`b`i'1'+`r`i'') var(`b`i'1'+`r`i'')
         random `x0' `if' `in', dist(gamma) mean(`b`i'0'+`f`i'') var(`b`i'0'+`f`i'')
         gen `piobs`i''=`x1'/(`x0'+`x1')
         drop `x0' `x1'

         gen `pistar`i''=(1-`alpha`i'')*`piobs`i'' + `alpha`i''*exp(`logimor`i'')*`piobs`i''/(1-`piobs`i''+exp(`logimor`i'')*`piobs`i'')
         drop `logimor`i'' `alpha`i''
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

      if "`debug'"=="debug" {
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
   if `nmeasures'==1 {
      gen `eststar' = ESTSTAR_`measure'
      gen `sestar' = sqrt(VARSTAR_`measure')
   }

   if "`timer'"=="timer" {
      now
      local elapse=r(datetime)-`start'
      noi di in green _newline "Elapsed time: " in yellow `elapse' in green " seconds"
   }
}

} /* end of quietly */

********************************** END OF MONTE CARLO PART **********************************
}
********************************** END OF IAN'S CALCULATION **********************************

else {
di as text "Approach: imputation."
********************************* JULIAN'S PARSING *********************************

if `syntax'==1 {
if "`aca'"=="aca" local method "Available case analysis"
else if "`ica0'"! = ""  local method "imputed case analysis ICA-0 (impute zeros)"
else if "`ica1'"! = ""  local method "imputed case analysis ICA-1 (impute ones)"
else if "`icape'"!= ""  local method "imputed case analysis ICA-pE (impute experimental group risk)"
else if "`icapc'"!= ""  local method "imputed case analysis ICA-pC (impute control group risk)"
else if "`icap'"! = ""  local method "imputed case analysis ICA-p (impute group-specific risk)"
else if "`icab'"! = ""  local method "imputed case analysis ICA-b ('best'-case [0 in exp, 1 in control])"
else if "`icaw'"! = ""  local method "imputed case analysis ICA-w ('worst'-case [1 in exp, 0 in control])"
else if "`icaimor'"!="" local method "imputed case analysis ICA-IMOR (impute using IMORs `imorE' `imorC')"
}
else if `syntax'==2 local method "ICA-r (use specified reasons for missingness)"
di as text "Method: `method'."

* quit if imor reasons indicated but IMOR not given
if ("`icaimor'" != "" & "`imorE'" == "") {
    di as error "icaimor specified but no IMOR provided"
    exit 498
}
* now if user has not entered imors, replace with ones to avoid errors in later code (coeff will always be zero)
foreach imorR in imorE imorC {
    if "``imorR''"=="" local `imorR' = 1
}

** warning if user gives imors as well as reasons without providing reason=imor
*if ("`imorE'" ~= "") & (`numimorE' == 0) {
*    di in ye "Ignoring imorE and using reasons (no icaimorE)"
*}
*if ("`imorC'" ~= "") & (`numimorC' == 0) {
*    di in ye "Ignoring imorC and using reasons (no icaimorC)"
*}

** if user has specified IMOR and nothing else, will set number of imors to number of missing
** DELETE?
*tempvar useimorEforall useimorCforall
*quiet gen `useimorEforall' = 0
*quiet gen `useimorCforall' = 0
*if ("`imorC'" ~= "") & ("`aca'" == "" & "`ica0'" == "" & "`ica1'" == "" & "`icap'" == "" & "`icape'" == "" & "`icapc'" == "" & "`icab'" == "" & "`icaw'" == "" & "`ica0C'" == "" & "`ica1C'" == "" & "`icapeC'" == "" & "`icapcC'" == "" & "`icapC'" == "" & "`icaimorC'" == "") {
*    di in ye "Applying imorC to all missing in control group"
*    quiet replace `useimorCforall' = 1
*}
*if ("`imorE'" ~= "") & ("`aca'" == "" & "`ica0'" == "" & "`ica1'" == "" & "`icap'" == "" & "`icape'" == "" & "`icapc'" == "" & "`icab'" == "" & "`icaw'" == "" & "`ica0E'" == "" & "`ica1E'" == "" & "`icapeE'" == "" & "`icapcE'" == "" & "`icapE'" == "" & "`icaimorE'" == "") {
*    di in ye "Applying imorE to all missing in experimental group"
*    quiet replace `useimorEforall' = 1
*}

* check only one weight specified
local weight `w1'`w2'`w3'`w4'
if length("`weight'")>2 {
    di as error "Please specify only one of w1/w2/w3/w4"
    exit 498
}
if length("`weight'")==0 {
    if "`aca'"=="aca" {
       local w2 w2
       local weight w2
    }
    else {
       di as error "Please specify w1/w2/w3/w4"
       exit 498
    }
}

*local aca = cond("`aca'" =="", "*", "quiet")
*local ica0 = cond("`ica0'" =="", "*", "quiet")
*local ica1 = cond("`ica1'" =="", "*", "quiet")
*local icap = cond("`icap'" =="", "*", "quiet")
*local icape = cond("`icape'" =="", "*", "quiet")
*local icapc = cond("`icapc'" =="", "*", "quiet")
*local icab = cond("`icab'" =="", "*", "quiet")
*local icaw = cond("`icaw'" =="", "*", "quiet")
*local icaimor = cond("`icaimor'" =="", "*", "quiet")


*if user has not entered reasons for a whole column, replace with zeros (i.e. assumes no people missing for that particular reason)
*foreach reasoncount in ica0E ica1E icapeE icapcE icapE icaimorE ica0C ica1C icapeC icapcC icapC icaimorC {
*    if "``reasoncount''"=="" {
*        tempvar `reasoncount'
*        quiet gen ``reasoncount'' = 0
*    }
*
*}


************************ END OF JULIAN'S PARSING ****************************

************************ JULIAN'S CALCULATION ****************************
if "`debug'"=="debug" di as text "Starting Julian's calculation"

*tempvar uica0E uica1E uicapeE uicapcE uicapE uicaimorE uica0C uica1C uicapeC uicapcC uicapC uicaimorC uimorE uimorC umE umC
** make copy of variables to avoid re-writing existing variables
*quiet{
*    g `umE' = `mE'
*    g `umC' = `mC'
*    g `uica0E' = `ica0E'
*    g `uica1E' = `ica1E'
*    g `uicapeE' = `icapeE'
*    g `uicapcE' = `icapcE'
*    g `uicapE' = `icapE'
*    g `uicaimorE' = `icaimorE'
*    g `uica0C' = `ica0C'
*    g `uica1C' = `ica1C'
*    g `uicapeC' = `icapeC'
*    g `uicapcC' = `icapcC'
*    g `uicapC' = `icapC'
*    g `uicaimorC' = `icaimorC'
*    g `uimorE' = `imorE'
*    g `uimorC' = `imorC'
*}

* SORT OUT DATA NOT INCLUDING ZEROES
tokenize "`varlist'"
tempvar rE fE mE rC fC mC
qui {
   gen `rE' = `1'
   gen `fE' = `2'
   gen `mE' = `3'
   gen `rC' = `4'
   gen `fC' = `5'
   gen `mC' = `6'
}
foreach arm in E C {
   foreach reason in 0 1 pe pc p imor {
      tempvar num`reason'`arm'
      if `syntax'==1 {
         if "`ica`reason''"=="ica`reason'" qui gen `num`reason'`arm'' = `m`arm''
         else gen `num`reason'`arm'' = 0
      }
      if `syntax'==2 {
         if "`ica`reason'`arm''"!="" qui gen `num`reason'`arm'' = `ica`reason'`arm''
         else gen `num`reason'`arm'' = 0
      }
      char `num`reason'`arm''[varname] num`reason'`arm'
      local numlist `numlist' `num`reason'`arm''
   }
}
if "`listnum'"=="listnum" l `numlist', subvarname


* for studies with missing data, but no reasons, assign number for each reason to overall proportion (across studies) with that reason

foreach arm in E C {
   tempvar denom`arm'
   egen `denom`arm'' = sum(`num0`arm''+`num1`arm''+`numpe`arm''+`numpc`arm''+`nump`arm''+`numimor`arm'')
   foreach reason in 0 1 pe pc p imor {
      tempvar numer`reason'`arm' 
      egen `numer`reason'`arm'' = sum(`num`reason'`arm'')
      qui replace `num`reason'`arm'' = `numer`reason'`arm'' / `denom`arm'' if `num`reason'`arm''==.
* if user specifies aca, ica0, ica1, icapc, icape, icap, icab, icaw, then do that and override other imputation arguments
      if `syntax'==1 {
         qui replace `num`reason'`arm'' = 0
         if "`ica`reason''"!="" qui replace `num`reason'`arm'' = 1
      }
   }
*   quiet replace `numimor`arm'' = 1 if `useimor`arm'forall' == 1
}
if "`icab'"!="" qui replace `num1E' = 1
if "`icab'"!="" qui replace `num0C' = 1
if "`icaw'"!="" qui replace `num0E' = 1
if "`icaw'"!="" qui replace `num1C' = 1


* calculate basic sample sizes and probabilities

tempvar nE nC aE aC pE pC numallE numallC


quiet{
*    gen `nE' = `rE' + `fE' + `mE'
*    gen `nC' = `rC' + `fC' + `mC'
*    gen `aE' = `mE' / `nE'
*    gen `aC' = `mC' / `nC'
*    gen `pE' = `rE' / (`rE' + `fE')
*    gen `pC' = `rC' / (`rC' + `fC')


* pretend no missing data when ACA required

*    `aca' replace `nE' = `rE' + `fE'
*    `aca' replace `nC' = `rC' + `fC'
    if "`aca'"!="" replace `mE' = 0
    if "`aca'"!="" replace `mC' = 0
*    `aca' replace `aE' = 0
*    `aca' replace `aC' = 0
}


* calculate numbers that indicate whether there are zero cells after imputations of 0s and 1s

tempvar temprE tempfE temprC tempfC addpointfive

quiet{
    gen `addpointfive' = 0
    gen `temprE' = `rE' + `num1E'*`mE'
    gen `tempfE' = `fE' + `num0E'*`mE'
    gen `temprC' = `rC' + `num1C'*`mC'
    gen `tempfC' = `fC' + `num0C'*`mC'
    replace `addpointfive' = 1 if (`temprE'==0 & `tempfE'~=0 & `temprC'~=0 & `tempfC'~=0 | `temprE'~=0 & `tempfE'==0 & `temprC'~=0 & `tempfC'~=0 | `temprC'==0 & `tempfC'~=0 & `temprE'~=0 & `tempfE'~=0  | `temprC'~=0 & `tempfC'==0 & `temprE'~=0 & `tempfE'~=0 )
    count if `addpointfive'==1
}
if r(N)>0 di as text "Zero cells detected: adding 1/2 to " as result r(N) as text " studies."

* create new data set with 0.5s added where necessary

quiet{
foreach arm in E C {
    replace `r`arm'' = `r`arm''+`addpointfive'/2
    replace `f`arm'' = `f`arm''+`addpointfive'/2
    gen `n`arm'' = `r`arm'' + `f`arm'' + `m`arm''
    gen `p`arm'' = `r`arm'' / (`r`arm'' + `f`arm'')
    gen `a`arm'' = `m`arm'' / `n`arm''
    foreach start in r f n p a n r f m num0 num1 numpe numpc nump numimor {
       char ``start'`arm''[varname] `start'`arm'
    }
}
}
if "`listx'"=="listx" {
    l `id' `rE' `fE' `nE' `pE' `aE' `rC' `fC' `nC' `pC' `aC', subvarname
}

* calculate probability of imputing each specific reason
* asup*E's sum to aE, asup*C's sum to aC
quiet {
   foreach arm in E C {
      gen `numall`arm'' = `num0`arm''+`num1`arm''+`numpe`arm''+`numpc`arm''+`nump`arm''+`numimor`arm''
      foreach reason in 0 1 pe pc p imor {
         tempvar asup`reason'`arm'
         gen `asup`reason'`arm'' = `a`arm''* `num`reason'`arm'' / `numall`arm'' if `numall`arm''~=0
         replace `asup`reason'`arm'' = 0 if `numall`arm''==0
         replace `num`reason'`arm'' = `n`arm''*`asup`reason'`arm''
      }
   }
}

if "`listall'"=="listall" {
    foreach arm in E C {
        l `id' `n`arm'' `r`arm'' `f`arm'' `m`arm'' `num0`arm'' `num1`arm'' `numpe`arm'' `numpc`arm'' `nump`arm'' `numimor`arm'', subvarname
    }
}

* calculate risks among missing participants using:
*       p`arm' - AC success fraction
*       a`arm' - fraction of missing data
*       asup`reason'`arm' - fraction of missing data by reason 

foreach arm in E C {
    tempvar pfromimor`arm' p`arm'star
    if "`icaimor'"!="" qui gen `pfromimor`arm'' = (`p`arm''*`imor`arm'') / (`p`arm''*`imor`arm'' + 1 - `p`arm'')
    else gen `pfromimor`arm'' = 0
    qui gen `p`arm'star' = `p`arm''*(1-`a`arm'') + `asup1`arm'' + `asuppc`arm''*`pC' + `asuppe`arm''*`pE' + `asupp`arm''*`p`arm'' + `asupimor`arm''*`pfromimor`arm''
    char `pfromimor`arm''[varname] pfromimor`arm'
    char `p`arm'star'[varname] p`arm'star   
}

* calculate new 'effective 2x2 table based on new overall risks, incorporating missing data imputations, with original sample size

tempvar erE efE enE erC efC enC

quiet {
    gen `erE' = `nE' * `pEstar'
    gen `efE' = `nE' * (1-`pEstar')
    gen `enE' = `erE' + `efE'
    gen `erC' = `nC' * `pCstar'
    gen `efC' = `nC' * (1-`pCstar')
    gen `enC' = `erC' + `efC'
}

if "`listp'"=="listp" {
    list `id' `pE' `pEstar' `pC' `pCstar', subvarname
}

tempvar RR lnRR W1selnRR W2selnRR OR lnOR W1selnOR W2selnOR VE VC W4varlnRR W4selnRR


if "`debug'"=="debug" di "Calculating estimates ..."
qui {
   if "`measure'"=="RD" gen `eststar' = `pEstar' - `pCstar'
   else if "`measure'"=="RR" gen `eststar' = log(`pEstar'/`pCstar')
   else if "`measure'"=="OR" gen `eststar' = log((`pEstar'/(1-`pEstar'))/(`pCstar'/(1-`pCstar')))
}

* calculate std errors
if "`debug'"=="debug" di "Calculating standard errors ..."
if "`weight'"=="w1" {
   if "`measure'"=="RD" gen `sestar' = sqrt(`erE'*(`enE'-`erE')/(`enE')^3 + `erC'*(`enC'-`erC')/(`enC')^3 )
   if "`measure'"=="RR"  gen `sestar' = sqrt(1/`erE' - 1/`enE' + 1/`erC' - 1/`enC')
   if "`measure'"=="OR" gen `sestar' = sqrt(1/`erE' + 1/(`enE'-`erE') + 1/`erC' + 1/(`enC'-`erC'))
}
else if "`weight'"=="w2" {
   * ACA fixed-effect analysis for standard errors
   quiet metan `rE' `fE' `rC' `fC', fixedi `or' `rr' `rd' nointeger nograph
   if "`measure'"=="RD" gen `sestar' = _seES
   else quiet gen `sestar' = _selogES
}
else if "`weight'"=="w3" {
   * apply pstar probabilities to original ACA sample sizes (after zero cell adjustment)
   tempvar W3rE W3rC W3fE W3fC W3selnRR W3selnOR
   quiet {
      gen `W3rE' = `pEstar' * (`rE'+`fE')
      gen `W3rC' = `pCstar' * (`rC'+`fC')
      gen `W3fE' = (1-`pEstar') * (`rE'+`fE')
      gen `W3fC' = (1-`pCstar') * (`rC'+`fC')
      if "`measure'"=="RD" gen `sestar' = sqrt(`W3rE'*`W3fE'/(`rE'+`fE')^3 + `W3rC'*`W3fC'/(`rC'+`fC')^3)
      if "`measure'"=="RR" gen `sestar' = sqrt(1/`W3rE' - 1/(`rE'+`fE') + 1/`W3rC' - 1/(`rC'+`fC'))
      if "`measure'"=="OR" gen `sestar' = sqrt(1/`W3rE' + 1/`W3fE' + 1/`W3rC' + 1/`W3fC')
   }
}
else if "`weight'"=="w4" {
   * from formula derived from Taylor series
   * var of pstar adds over multiple reasons for missingness - here 0, 1, E, C, P and imor
   * formulae use 0-cell adjusted PSTARs, Ns, Ps
   tempvar OR_CE OR_EC varpEstar varpCstar firstsumE lastsumE firstsumC lastsumC

   quiet {
      gen `OR_CE' = (`pC' / (1-`pC')) / (`pE' / (1-`pE'))
      gen `OR_EC' = (`pE' / (1-`pE')) / (`pC' / (1-`pC'))

      replace `OR_CE' = 9999999 if (`pC'==1 | `pE'==0)
      replace `OR_EC' = 9999999 if (`pE'==1 | `pC'==0)

      gen `firstsumE' = - `asup0E' - `asup1E' + `asuppcE' * (`OR_CE' / ((`pE'*`OR_CE' + 1 - `pE' )^2) - 1) + `asupimorE' * (`imorE' / ((`pE'*`imorE' + 1 - `pE' )^2) - 1)
      gen `firstsumC' = - `asup0C' - `asup1C' + `asuppeC' * (`OR_EC' / ((`pC'*`OR_EC' + 1 - `pC' )^2) - 1) + `asupimorC' * (`imorC' / ((`pC'*`imorC' + 1 - `pC' )^2) - 1)

      gen `lastsumE' = `asup0E' * (`pEstar'^2 - (`pE'-`pEstar')^2) + `asup1E' * ((1-`pEstar')^2 - (`pE'-`pEstar')^2) + `asuppcE' * ((`pC'-`pEstar')^2 - (`pE'-`pEstar')^2) + `asupimorE' * ((`pfromimorE'-`pEstar')^2 - (`pE'-`pEstar')^2)
      gen `lastsumC' = `asup0C' * (`pCstar'^2 - (`pC'-`pCstar')^2) + `asup1C' * ((1-`pCstar')^2 - (`pC'-`pCstar')^2) + `asuppeC' * ((`pE'-`pCstar')^2 - (`pC'-`pCstar')^2) + `asupimorC' * ((`pfromimorC'-`pCstar')^2 - (`pC'-`pCstar')^2)

      gen `varpEstar' = (`pE' * (1-`pE') / (`nE' - `mE') ) * (1 + `firstsumE')^2 + (1/`nE') * ((`pE'-`pEstar')^2 + `lastsumE')
      gen `varpCstar' = (`pC' * (1-`pC') / (`nC' - `mC') ) * (1 + `firstsumC')^2 + (1/`nC') * ((`pC'-`pCstar')^2 + `lastsumC')

      if "`measure'"=="RD" gen `sestar' = sqrt( cond( `pEstar'==0, 0, `varpEstar') +  cond(`pCstar'==0, 0, `varpCstar'))

      if "`measure'"=="RR"  gen `sestar' = sqrt( cond( `pEstar'==0, 0, `varpEstar'/(`pEstar'^2)) +  cond(`pCstar'==0, 0, `varpCstar'/(`pCstar'^2)))

      if "`measure'"=="OR" gen `sestar' = sqrt( cond( `pEstar'==0, 0, `varpEstar'/((`pEstar'*(1-`pEstar'))^2) ) + cond(`pCstar'==0, 0, `varpCstar'/((`pCstar'*(1-`pCstar'))^2) ) )
   }
}

************************ END OF JULIAN'S CALCULATION ****************************
}

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
    if "`debug'"=="debug" di "metan `eststar' `sestar' `if' `in', `idopt' `options' `eform' `graph'"
*mac dir
*pause
    metan `eststar' `sestar' `if' `in', `idopt' `options' `eform' `graph'
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

*! version 2, 13 Oct 99
*! version 2.1, 15may2002-11oct2002 - update to version 7, cope with missing mean
*! version 2.2, 24oct2002 - allow if and in
*! version 2.3, 26jun2003 - allow abbreviated option names
*! version 2.4, 08mar2004 - dist defaults to normal
*! version 2.5, 26aug2004 - allow sd() option
*! version 2.6  06jul2005 - lp option for dist(bernoulli)
* note: random c fails if a variable starting with c exists!
prog def random
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

****************************************************************************************

prog def expectn
syntax anything, gen(name) [genvar(name) nip(int 10) maxit(int 1000) tol(int 10) trace token(string) list debug]
confirm new var `gen'
if "`genvar'"!="" confirm new var `genvar'
if "`token'"=="" local token @
di as text "Finding expectation of `anything', assuming `token'~N(0,1) ..."
di as text "(using Gauss-Hermite quadrature with `nip' integration points)"
* GENERATE INTEGRATION POINTS AND WEIGHTS
local pim4 .751125544464925
local m=(`nip'+1)/2
forvalues i=1/`m' {
   if "`trace'"=="trace" di "Starting i=`i'"
   if `i'==1 local z=sqrt(float(2*`nip'+1))-1.85575*(2*`nip'+1)^(-0.16667)
   else if `i'==2 local z=`z'-1.14*`nip'^0.426/`z'
   else if `i'==3 local z=1.86*`z'-0.86*`zlag2'
   else if `i'==4 local z=1.91*`z'-0.91*`zlag2'
   else local z=2*`z'-`zlag2'

   local stop 0
   local iter 0
   while `stop'<1 {
      local ++iter
      local p1=`pim4'
      local p2=0
      forvalues j=1/`nip' {
         local p3=`p2'
         local p2=`p1'
         local p1=`z'*sqrt(2/`j')*`p2'-sqrt((`j'-1)/`j')*`p3'
      }
      local pp=sqrt(2*`nip')*`p2'
      local z1=`z'
      local z=`z1'-`p1'/`pp'
      local stop = abs(`z'-`z1')<10^(-`tol')
      if "`trace'"=="trace" di "Iteration `iter', change = " `z'-`z1'
      if `iter'>`maxit' {
         di as error "Too many iterations for i=`i'"
         exit 498
      }
   }
   if "`trace'"=="trace" di "i=`i' required `iter' iterations"

   if `i'==`m' & abs(`z')<1E-10 local z 0

   local j=`nip'+1-`i'
   local w`i' = 2/(`pp'*`pp'*sqrt(_pi))
   local x`i' = -`z'*sqrt(2)
   local w`j' = `w`i''
   local x`j' = -`x`i''
   local zlag2 `zlag1'
   local zlag1 `z'

}

if "`list'"=="list" {
   di as text _col(10) "Integration point" _col(30) "Weight"
   forvalues i=1/`nip' {
      di as result `i' _col(10) `x`i'' _col(30) `w`i''
   }
}

* DO THE INTEGRATION
tempvar argument sumarg2 sumarg sumw
gen `sumarg2' = 0
gen `sumarg' = 0
gen `sumw' = 0
qui gen `argument' = .
local wordcount=wordcount("`token'")
forvalues j=1/`wordcount' {
   local token`j'=word("`token'",`j')
}
if `wordcount'==1 {
   forvalues i1=1/`nip' {
      local arg : subinstr local anything "`token1'" "(`x`i1'')", all
      qui replace `argument' = `arg'
      qui count if missing(`argument')
      if r(N)>0 {
          di /*as error*/ "expectn: missing expression (`argument') in " r(N) " cases at `token1'=`x`i1''"
          l if missing(`argument')
*          exit 498
      }
      qui replace `sumarg2' = `sumarg2' + `w`i1''*(`argument')^2
      qui replace `sumarg' = `sumarg' + `w`i1''*(`argument')
      qui replace `sumw' = `sumw' + `w`i1''
if "`debug'"=="debug" noi di `arg', `sumw', `sumarg', `sumarg2'
   }
}
else if `wordcount'==2 {
   forvalues i1=1/`nip' {
      forvalues i2=1/`nip' {
         local arg : subinstr local anything "`token1'" "(`x`i1'')", all
         local arg : subinstr local arg "`token2'" "(`x`i2'')", all
*        note: the following syntax fails as it truncates strings
*        local arg = subinstr("`anything'","`token1'","(`x`i1'')",.)
         qui replace `argument' = `arg'
         qui count if missing(`argument')
         if r(N)>0 {
             di /*as error*/ "expectn: missing expression `argument' in " r(N) " cases at `token1'=`x`i1'', `token2'=`x`i2''"
             l if missing(`argument')
*             exit 498
         }
         local w=`w`i1''*`w`i2''
         qui replace `sumarg2' = `sumarg2' + `w'*(`argument')^2
         qui replace `sumarg' = `sumarg' + `w'*`argument'
         qui replace `sumw' = `sumw' + `w'
if "`debug'"=="debug" noi di `arg', `sumw', `sumarg', `sumarg2'
      }
   }
}
else if `wordcount'==3 {
   forvalues i1=1/`nip' {
      forvalues i2=1/`nip' {
         forvalues i3=1/`nip' {
            local arg : subinstr local anything "`token1'" "(`x`i1'')", all
            local arg : subinstr local arg "`token2'" "(`x`i2'')", all
            local arg : subinstr local arg "`token3'" "(`x`i3'')", all
            qui replace `argument' = `arg'
            qui count if missing(`argument')
            if r(N)>0 {
                di /*as error*/ "expectn: missing expression `argument' in " r(N) " cases at `token1'=`x`i1'', `token2'=`x`i2'', `token3'=`x`i3''"
                l if missing(`argument')
*                exit 498
            }
            local w=`w`i1''*`w`i2''*`w`i3''
            qui replace `sumarg2' = `sumarg2' + `w'*(`argument')^2
            qui replace `sumarg' = `sumarg' + `w'*(`argument')
            qui replace `sumw' = `sumw' + `w'
if "`debug'"=="debug" noi di `arg', `sumw', `sumarg', `sumarg2'
         }
      }
   }
}
else {
    di as error "This program hasn't yet been adapted to handle >3 Normal deviates"
    di as error "Please fiddle with the code"
    exit 498
}
cap assert abs(`sumw'-1)<1E-6
if _rc {
    di as error "Weights failed to sum to 1 in following cases:"
    l `sumw' if abs(`sumw'-1)>=1E-6
}
gen `gen' = `sumarg'/`sumw'
label var `gen' "Expectation of `anything'"
if "`genvar'"!="" {
    gen `genvar' = `sumarg2'/`sumw'-`gen'^2
    label var `genvar' "Variance of `anything'"
}
di as text "Variable `gen' created"
end

prog def exit498
di as error `"`0'"'
exit 498
end
