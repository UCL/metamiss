/*
Test script for metamiss
Added dicmd 12/5/2022
Set seed for reproducibility 12/5/2022
Extended 11apr2022 to loop over all three measures
Extended 27sep2018 to include simple tests under MAR and M=F
*/

// PRELIMINARIES
cd c:\ado\ian\metamiss
set linesize 79
cap log close
log using metamisstest.log, replace

cscript metamiss
which metamiss
which metan
di c(stata_version)
set seed 3467145

// UTILITIES
prog def dicmd
noi di as input `"`0'"'
`0'
end

prog def mycompare
if mi("`3'") local 3 1E-7
if reldif(`1',`2')>`3' {
	if !mi("`4'") di as error "Error in `4'"
	di as error "Results differ relatively by " reldif(`1',`2') " > `3'"
	di as error "Result 1: `1'"
	di as error "Result 2: `2'"
	exit 9
}
else di as text "Results differ relatively by " reldif(`1',`2') " <= `3'"
end

// CREATE DATA
clear
input study r1 f1 m1 r2 f2 m2
1 10 80 10 10 80 0
2 20 70 10 10 50 40
end
gen logimor=_n
gen zero = 0
gen fm1=f1+m1
gen fm2=f2+m2

// RUN TESTS
foreach measure in rd rr or {

// compare with metan under MAR
dicmd metamiss r1 f1 zero r2 f2 zero, logimor(0) nograph `measure'
local result_metamiss = r(ES)
dicmd metan r1 f1 r2 f2, nograph fixedi `measure'
local result_metan = r(ES)
mycompare `result_metamiss' `result_metan' 2E-8 "MAR `measure'"

// compare with metan under M=F
dicmd metamiss r1 f1 m1 r2 f2 m2, imor(0) nograph `measure'
local result_metamiss = r(ES)
dicmd metan r1 fm1 r2 fm2, nograph fixedi `measure'
local result_metan = r(ES)
mycompare `result_metamiss' `result_metan' 2E-8 "M=F `measure'"

// compare opposite ways round
dicmd metamiss r1 f1 m1 r2 f2 m2, imor(1 2) sdlogimor(.2 .4) nograph `measure'
local result_metamiss = r(ES)
dicmd metamiss r2 f2 m2 r1 f1 m1, imor(2 1) sdlogimor(.4 .2) nograph `measure'
if "`measure'" == "rd" local result_metan = -r(ES)
else local result_metan = 1/r(ES)
mycompare `result_metamiss' `result_metan' 1E-8 "Opp `measure'"

// check error message
dicmd cap noi metamiss r1 f1 m1 r2 f2 m2, imor(logimor) sdlogimor(.2 .4) nograph `measure'
assert _rc==498

// compare with metamiss2
dicmd metamiss r1 f1 m1 r2 f2 m2, logimor(1) nograph  `measure'
local result1 = r(ES)
dicmd metamiss2 r1 f1 m1 r2 f2 m2, impmean(1) fixed metanopt(nograph) `measure'
local result2 = r(ES)
mycompare `result1' `result2' 1E-8 "metamiss2 `measure'"

dicmd metamiss r1 f1 m1 r2 f2 m2, logimor(-1) sdlogimor(1) nograph method(Taylor) `measure'
local result1 = r(ES)
dicmd metamiss2 r1 f1 m1 r2 f2 m2, impmean(-1) impsd(1) fixed metanopt(nograph) `measure'
local result2 = r(ES)
mycompare `result1' `result2' 1E-8 "metamiss2 impsd `measure'"

}

// check help file
runhelpfile using metamiss.sthlp

// REPORT SUCCESS
di as result _n "*****************************************" ///
	_n "*** METAMISS HAS PASSED ALL ITS TESTS ***" ///
	_n "*****************************************"

log close
