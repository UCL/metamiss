/*
Test script for metamiss
Extended 27sep2018 to include simple tests under MAR and M=F
*/

cd c:\ado\ian\metamiss
set linesize 79
log using metamisstest.log, replace

cscript metamiss
which metamiss

clear
input study r1 f1 m1 r2 f2 m2
1 10 80 10 10 80 0
2 20 70 10 10 50 40
end
gen logimor=_n

// compare with metan under MAR
gen zero = 0
metamiss r1 f1 zero r2 f2 zero, logimor(0) nograph
local result_metamiss = r(ES)
metan r1 f1 r2 f2, nograph fixedi
local result_metan = r(ES)
assert float(`result_metamiss') == float(`result_metan')

// compare with metan under M=F
gen fm1=f1+m1
gen fm2=f2+m2
metamiss r1 f1 m1 r2 f2 m2, imor(0) nograph
local result_metamiss = r(ES)
metan r1 fm1 r2 fm2, nograph fixedi
local result_metan = r(ES)
assert float(`result_metamiss') == float(`result_metan')

// compare opposite ways round
metamiss r1 f1 m1 r2 f2 m2, imor(1 2) sdlogimor(.2 .4) nograph
local result_metamiss = r(ES)
metamiss r2 f2 m2 r1 f1 m1, imor(2 1) sdlogimor(.4 .2) nograph
local result_metan = 1/r(ES)
assert reldif(`result_metamiss', `result_metan') < 1E-8
 
// check error message
cap noi metamiss r1 f1 m1 r2 f2 m2, imor(logimor) sdlogimor(.2 .4) nograph
assert _rc==498

// compare with metamiss2
metamiss r1 f1 m1 r2 f2 m2, logimor(1) nograph 
local result1 = r(ES)
metamiss2 r1 f1 m1 r2 f2 m2, impmean(1) fixed metanopt(nograph)
local result2 = r(ES)
assert float(`result1') == float(`result2')

metamiss r1 f1 m1 r2 f2 m2, logimor(-1) sdlogimor(1) nograph method(Taylor)
local result1 = r(ES)
metamiss2 r1 f1 m1 r2 f2 m2, impmean(-1) impsd(1) fixed metanopt(nograph)
local result2 = r(ES)
assert float(`result1') == float(`result2')

// REPORT SUCCESS
di as result _n "*******************************************" ///
	_n "*** MVMETAMISS HAS PASSED ALL ITS TESTS ***" ///
	_n "*******************************************"

log close
