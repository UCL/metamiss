*! version 1.1   Ian White   28oct2002
*! version 1.2   Ian White   02dec2002 - uses new _sensmiss.ado
*! version 1.2.1 09jan03 - replace option

prog def metamiss

* metamiss rE fE mE rC fC mC, phicorr(0.5) phisd(1.22)
* currently assumes phi's are uncorrelated between studies

syntax varlist(min=6 max=6) [fweight], [id(varname) *]
tokenize "`varlist'"
local rE `1'
local fE `2'
local mE `3'
local rC `4'
local fC `5'
local mC `6'

foreach meas in RD logRR logOR {
   foreach adj in raw star {
      if "`replace'"=="replace" {
        cap drop `meas'`adj'
        cap drop `meas'var`adj'
      }
     qui gen `meas'`adj'=.
     qui gen `meas'var`adj' = .
   }
}

local n = _N
forvalues i = 1/`n' {
local values
foreach var of var `rE' `fE' `mE' `rC' `fC' `mC' {
  local val = `var'[`i']
  local values `values' `val'
}
if "`id'"~="" {di _new in green "`id' = " in yellow `id'[`i']}
else {di _new in green "i = " in yellow `i'}
_sensmiss `values' [weight][ex], `options'
foreach meas in RD logRR logOR {
  foreach adj in raw star {
    qui replace `meas'`adj' = r(`meas'`adj') in `i'
    qui replace `meas'var`adj' = r(`meas'var`adj') in `i'
  }
}

end

