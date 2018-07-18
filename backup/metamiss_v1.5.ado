*! version 1.1   Ian White   28oct2002
*! version 1.2   Ian White   02dec2002 - uses new _sensmiss.ado
*! version 1.3   Ian White   03jan2003 - weights disabled in metamiss, debug enabled
*! version 1.3.1 Ian         09jan03 - replace option
*! version 1.4   5feb2003 - can specify measures, prefix and suffix
*! version 1.5   9may2005 - upgrade to stata 8; rename phi (logIMOR) as delta

prog def metamiss
version 8
* metamiss rE fE mE rC fC mC, deltacorr(0.5) deltasd(1.22)
* currently assumes delta's are uncorrelated between studies

syntax varlist(min=6 max=6), [id(varname) debug replace MEASures(string) PREfix(string) SUFfix(string) PRInt *]
tokenize "`varlist'"
local rE `1'
local fE `2'
local mE `3'
local rC `4'
local fC `5'
local mC `6'
if "`measures'"=="" {
   local measures RD logRR logOR
}
if "`print'"=="" {
   local print noprint
}
else {
   local print
}

* SET UP
foreach meas in `measures' {
   if "`meas'"~="RD" & "`meas'"~="logRR" & "`meas'"~="logOR" {
      di in red "measures() must contain any or all of RD logRR logOR"
      exit 498
   }
   foreach type in raw star varraw varobs varimor varstar {
      if "`replace'"=="replace" {
        cap drop `prefix'`meas'`type'`suffix'
      }
     qui gen `prefix'`meas'`type'`suffix'=.
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
   if "`id'"~="" & "`print'"!="noprint" {
      di _new in green "`id' = " in yellow `id'[`i']
   }
   else if "`print'"!="noprint" {
      di _new in green "i = " in yellow `i'
   }
   if "`debug'"=="debug" {
      di "Running command: _sensmiss `values', measures(`measures') `options'"
   }
   _sensmiss `values', measures(`measures') `options' `print'
   foreach meas in `measures' {
      foreach type in raw star varraw varobs varimor varstar {
         qui replace `prefix'`meas'`type'`suffix' = r(`meas'`type') in `i'
      }
   }
}
end

