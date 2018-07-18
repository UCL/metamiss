*! version 1.1   Ian White   28oct2002
*! version 1.2   Ian White   02dec2002 - uses new _sensmiss.ado
*! version 1.3   Ian White   03jan2003 - weights disabled in metamiss, debug enabled
*! version 1.3.1 Ian         09jan03 - replace option
*! version 1.4   5feb2003 - can specify measures, prefix and suffix

prog def metamiss
version 7
* metamiss rE fE mE rC fC mC, phicorr(0.5) phisd(1.22)
* currently assumes phi's are uncorrelated between studies

syntax varlist(min=6 max=6), [id(varname) debug replace MEASures(string) PREfix(string) SUFfix(string) *]
tokenize "`varlist'"
local rE `1'
local fE `2'
local mE `3'
local rC `4'
local fC `5'
local mC `6'
if "`measures'"=="" {local measures RD logRR logOR}

* SET UP
foreach meas in `measures' {
   if "`meas'"~="RD" & "`meas'"~="logRR" & "`meas'"~="logOR" {
      di in red "measures() must contain any or all of RD logRR logOR"
      exit 498
   }
   foreach adj in raw star {
      if "`replace'"=="replace" {
        cap drop `prefix'`meas'`adj'`suffix'
        cap drop `prefix'`meas'var`adj'`suffix'
      }
     qui gen `prefix'`meas'`adj'`suffix'=.
     qui gen `prefix'`meas'var`adj'`suffix' = .
   }
}

* RUN SENSMISS FOR EACH STUDY
local n = _N
forvalues i = 1/`n' {
   local values
   foreach var of var `rE' `fE' `mE' `rC' `fC' `mC' {
      local val = `var'[`i']
      local values `values' `val'
   }
   if "`id'"~="" {di _new in green "`id' = " in yellow `id'[`i']}
   else {di _new in green "i = " in yellow `i'}
   if "`debug'"=="debug" {di "Running command: _sensmiss `values', measures(`measures') `options'"}
   _sensmiss `values', measures(`measures') `options'
   foreach meas in `measures' {
      foreach adj in raw star {
         qui replace `prefix'`meas'`adj'`suffix' = r(`meas'`adj') in `i'
         qui replace `prefix'`meas'var`adj'`suffix' = r(`meas'var`adj') in `i'
      }
   }
}
end

