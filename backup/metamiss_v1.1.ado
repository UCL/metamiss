*! version 1.1   Ian White   28oct2002
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
  qui gen `meas'`adj' = .
  qui gen `meas'var`adj' = .
}
}

local n = _N
forvalues i = 1/`n' {
local values
foreach var of var rE fE mE rC fC mC {
  local val = `var'[`i']
  local values `values' `val'
}
if "`id'"~="" {di _new "`id' = " `id'[`i']}
else {di _new "i = `i'"}
_forstersmith `values' [weight][ex], `options' rd rr or
foreach meas in RD logRR logOR {
  foreach adj in raw star {
    qui replace `meas'`adj' = r(`meas'`adj') in `i'
    qui replace `meas'var`adj' = r(`meas'var`adj') in `i'
  }
}

end

