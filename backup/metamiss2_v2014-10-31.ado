/*Missing data in meta-analysis with continuous outcomes, 31.10.2014, Anna Chaimani */

program metamiss2
	syntax varlist(min=8 max=8) [if] [in], [IMDom(string) IMRom(string) SDIMDom(string) SDIMRom(string) RHO(string) MD SMD ROM Bootstrap Taylor REPS(string) NOKEEP FIXED COMPare(string asis) SENSitivity *]
	
	tokenize `varlist'
	
	/*input data: observed, unobserved, mean, sd in treatment and control arm - consistent with metan*/
	local nt `1'
	local mt `2'
	local yt `3'
	local sdt `4'
	local nc `5'
	local mc `6'
	local yc `7'
	local sdc `8'
	
	/*mean and sd for lambda in treatment and control arms, default is lambda(0 0 0 0)*/
	if "`imdom'"!="" & "`imrom'"!=""{	
		di as err "Option {it:imdom()} may not be combined with option {it:imrom()}"
		exit
	}
	if "`sdimdom'"!="" & "`sdimrom'"!=""{	
		di as err "Option {it:sdimdom()} may not be combined with option {it:sdimrom()}"
		exit
	}
	if "`sdimdom'"!="" & "`imrom'"!=""{	
		di as err "Option {it:sdimdom()} may not be combined with option {it:imrom()}"
		exit
	}	
	if "`imdom'"!="" & "`sdimrom'"!=""{	
		di as err "Option {it:imdom()} may not be combined with option {it:sdimrom()}"
		exit
	}
	
	/*mean and sd for lambda*/
	local mulambdat=0
	local mulambdac=0
	local sdlambdat=0
	local sdlambdac=0
	local model "imdom"
	
	if "`imdom'"!=""{
		local model "imdom" //model IMDoM, default
		tokenize `imdom'
		local mulambdat `1'
		if "`2'"==""{			
			local mulambdac `1'
		}	
		else{
			local mulambdac `2'
		}
		if "`3'"!=""{
			di as err "Option {it:imdom()} requires 1 or 2 input values/variables"
		}		
	}
	if "`sdimdom'"!=""{
		if "`imdom'"==""{
			local model "imdom"
		}
		tokenize `sdimdom'
		local sdlambdat `1'
		if "`2'"==""{			
			local sdlambdac `1'
		}	
		else{
			local sdlambdac `2'
		}
		if "`3'"!=""{
			di as err "Option {it:sdimdom()} requires 1 or 2 input values/variables"
		}		
	}
	
	if "`imrom'"!=""{
		local model "imrom" //model IMRoM
		cap{
			assert `yt'>=0 & `yc'>=0
		}
		if _rc!=0{
			cap{
				assert `yt'<=0 & `yc'<=0
			}
			if _rc!=0{
				di as err "Effect sizes with opposite signs cannot be incorporated in the {it:imrom} model"
				exit
			}
		}		
		tokenize `imrom'
			local mulambdat `1'
		if "`2'"==""{			
			local mulambdac `1'
		}	
		else{
			local mulambdac `2'
		}
		if "`3'"!=""{
			di as err "Option {it:imrom()} requires 1 or 2 input values/variables"
		}		
	}
	if "`sdimrom'"!=""{
		if "`imrom'"==""{
			local model "imrom"
		}
		tokenize `sdimrom'
		local sdlambdat `1'
		if "`2'"==""{			
			local sdlambdac `1'
		}	
		else{
			local sdlambdac `2'
		}
		if "`3'"!=""{
			di as err "Option {it:sdimrom()} requires 1 or 2 input values/variables"
		}		
	}	
	
	/*correlation of lambdaT and lambdaC, default is rho(0)*/
	if "`rho'"!=""{
		local rho `rho'
	}
	else{
		local rho=0
	}	

	tempvar mulc sdlc mult sdlt rhol
	qui gen `mulc'=`mulambdac'
	qui gen `sdlc'=`sdlambdac'
	qui gen `mult'=`mulambdat'
	qui gen `sdlt'=`sdlambdat'
	qui gen `rhol'=`rho'
	
	if "`bootstrap'"==""{
		local method "taylor" //estimate relative effects and variances using parametric bootstrap, default
	}
	if "`bootstrap'"!=""{
		local method "bootstrap" //estimate relative effects and variances using Taylor approximation
	}
	if "`taylor'"!="" & "`bootstrap'"!=""{
		di as err "Option {it:bootstrap} may not be combined with option {it:taylor}"
		exit
	}
	if "`sdlambdat'"=="0" & "`sdlambdac'"=="0" & "`sensitivity'"==""{
		local method "taylor"
		if "`bootstrap'"!=""{
			di _newline " Note: Option {it:bootstrap} is not allowed when sd(IMDOM)/sd(IMROM) is zero for both experimental and control arms - reset to default option ({it:taylor})"
		}
	}

	if "`smd'"=="" & "`rom'"==""{
		local measure "md" //estimate mean differences, default
	}
	if "`smd'"!="" & "`rom'"=="" & "`md'"==""{
		local measure "smd" //estimate standardized mean differences
	}
	if "`rom'"!="" & "`smd'"=="" & "`md'"==""{
		local measure "rom" //estimate ratios of means
		cap{
			assert `yt'>=0 & `yc'>=0
		}
		if _rc!=0{
			cap{
				assert `yt'<=0 & `yc'<=0
			}
			if _rc!=0{
				di as err "Effect sizes with opposite signs cannot be used when ratios of means are estimated"
				exit
			}
		}		
	}
	if ("`smd'"!="" & "`rom'"!="") | ("`md'"!="" & "`rom'"!="") | ("`smd'"!="" & "`md'"!=""){
		di as err "Two or more effect measures specified"
		exit
	}
	
	if "`reps'"!="" & "`bootstrap'"==""{
		di as err "Option {it:reps()} requires option {it:bootstrap}"
		exit
	}	
	
	/*number of random samlpes for the bootstrap method, default is reps(10000)*/
	if "`reps'"!=""{
		local reps `reps'
	}
	else{
		local reps=10000
	}
*di "Starting analysis"
*foreach mac in model method measure sdlambdac sdlambdat mulambdac mulambdat rho {
*    mac list _`mac'
*}
	
	tempvar _pc _pt _ytotc _ytott _covf _f1 _sumc _sumt _f2 _vartot _f3 _f4
	qui{
		/*probabilities of the observations in control and treatment arms*/
		gen `_pc'=`nc'/(`nc'+`mc') `if' `in'
		gen `_pt'=`nt'/(`nt'+`mt') `if' `in'
		cap drop _es _sees
		
		if "`method'"=="taylor"{ //estimation via Taylor approximation
			if "`model'"=="imrom"{
				gen `_ytotc'=`_pc'*`yc'+(1-`_pc')*(exp(`mulc'+(`sdlc'^2)/2)*`yc') `if' `in' //estimate of 'true' mean in control arm
				gen `_ytott'=`_pt'*`yt'+(1-`_pt')*(exp(`mult'+(`sdlt'^2)/2)*`yt') `if' `in' //estimate of 'true' mean in treatment arm
				gen `_covf'=`rhol'*`sdlt'*`sdlc'*(1-`_pt')*(1-`_pc')*`yt'*`yc'*exp(`mult'+0.5*`sdlt'^2)*exp(`mulc'+0.5*`sdlc'^2) `if' `in'
				gen `_sumc'=((`_pc'*(1-`_pc'))/(`nc'+`mc'))*(1-2*exp(`mulc'+0.5*`sdlc'^2))*`yc'^2+(`sdc'^2/`nc')*((1-`_pc')^2*exp(2*`mulc'+2*`sdlc'^2)+2*`_pc'*(1-`_pc')*exp(`mulc'+0.5*`sdlc'^2)+`_pc'^2)+exp(2*`mulc'+`sdlc'^2)*(exp(`sdlc'^2)-1)*((1-`_pc')*`yc')^2 `if' `in'
				gen `_sumt'=((`_pt'*(1-`_pt'))/(`nt'+`mt'))*(1-2*exp(`mult'+0.5*`sdlt'^2))*`yt'^2+(`sdt'^2/`nt')*((1-`_pt')^2*exp(2*`mult'+2*`sdlt'^2)+2*`_pt'*(1-`_pt')*exp(`mult'+0.5*`sdlt'^2)+`_pt'^2)+exp(2*`mult'+`sdlt'^2)*(exp(`sdlt'^2)-1)*((1-`_pt')*`yt')^2	`if' `in'		
			}
			if "`model'"=="imdom"{
				gen `_ytotc'=`_pc'*`yc'+(1-`_pc')*(`mulc'+`yc') `if' `in' //estimate of 'true' mean in control arm
				gen `_ytott'=`_pt'*`yt'+(1-`_pt')*(`mult'+`yt') `if' `in' //estimate of 'true' mean in treatment arm
				gen `_covf'=`rhol'*`sdlt'*`sdlc'*(1-`_pt')*(1-`_pc') `if' `in'
				gen `_sumc'=((`_pc'*(1-`_pc'))/(`nc'+`mc'))*(`mulc'^2+`sdlc'^2)+(`sdc'^2/`nc')+`sdlc'^2*(1-`_pc')^2 `if' `in'
				gen `_sumt'=((`_pt'*(1-`_pt'))/(`nt'+`mt'))*(`mult'^2+`sdlt'^2)+(`sdt'^2/`nt')+`sdlt'^2*(1-`_pt')^2 `if' `in'
			}
			if "`measure'"=="md"{
				gen `_f1'=1 `if' `in'
				gen `_f2'=1 `if' `in'
				gen `_f4'=1 `if' `in'
				gen `_f3'=1 `if' `in'
				gen _es=`_ytott'-`_ytotc' `if' `in' //estimate of 'true' MD
			}
			if "`measure'"=="smd"{
				gen `_f1'=1 `if' `in'
				gen `_f2'=(`nc'+`nt'-2)/((`nc'-1)*`sdc'^2+(`nt'-1)*`sdt'^2) `if' `in'
				gen `_f4'=1 `if' `in'
				gen `_f3'=1	`if' `in'			
				gen _es=(`_ytott'-`_ytotc')/(1/sqrt(`_f2')) `if' `in' //estimate of 'true' SMD
			}
			if "`measure'"=="rom"{
				gen `_f1'=1/(`_ytotc'*`_ytott') `if' `in'
				gen `_f2'=1 `if' `in'
				gen `_f3'=1/(`yc'^2) `if' `in'
				gen `_f4'=1/(`yt'^2) `if' `in'
				gen _es=log(abs(`_ytott')/abs(`_ytotc')) `if' `in' //estimate of 'true' RoM
			}
			gen `_vartot'=`_f2'*(`_sumc'*`_f3'+`_sumt'*`_f4')-2*`_f1'*`_covf' `if' `in' //estimate of 'true' variances depending on model and measure
			gen _sees=sqrt(`_vartot') `if' `in'
		}
	}	
	if "`method'"=="bootstrap"{	//estimation via parametric bootstrap	
		mkmat `nc' `if' `in',mat(Mnc) nomissing
		mkmat `mc' `if' `in',mat(Mmc) nomissing
		mkmat `nt' `if' `in',mat(Mnt) nomissing
		mkmat `mt' `if' `in',mat(Mmt) nomissing
		mkmat `yc' `if' `in',mat(Myc) nomissing
		mkmat `yt' `if' `in',mat(Myt) nomissing
		mkmat `sdc' `if' `in',mat(Msdc) nomissing
		mkmat `sdt' `if' `in',mat(Msdt) nomissing
		mkmat `mulc' `if' `in',mat(Mmulc) nomissing
		mkmat `sdlc' `if' `in',mat(Msdlc) nomissing
		mkmat `mult' `if' `in',mat(Mmult) nomissing
		mkmat `sdlt' `if' `in',mat(Msdlt) nomissing
		mkmat `rhol' `if' `in',mat(Mrho) nomissing
		mat Mnc=Mnc'
		mat Mmc=Mmc'
		mat Mnt=Mnt'
		mat Mmt=Mmt'
		mat Myc=Myc'
		mat Myt=Myt'
		mat Msdc=Msdc'
		mat Msdt=Msdt'
		mat Mmulc=Mmulc'
		mat Msdlc=Msdlc'
		mat Mmult=Mmult'
		mat Msdlt=Msdlt'
		mat Mrho=Mrho'
		scalar Mreps=real("`reps'")
		scalar Mns=colsof(Mnc)
		mat Mmulc=Mmulc[1,1..`=Mns']
		mat Mmult=Mmult[1,1..`=Mns']
		mat Msdlc=Msdlc[1,1..`=Mns']
		mat Msdlt=Msdlt[1,1..`=Mns']
		mat Mrho=Mrho[1,1..`=Mns']
		if "`model'"=="imdom"{
			scalar Mod=1
		}
		if "`model'"=="imrom"{
			scalar Mod=2
		}
		if "`measure'"=="md"{
			scalar Meas=1
		}
		if "`measure'"=="smd"{
			scalar Meas=2
		}		
		if "`measure'"=="rom"{
			scalar Meas=3
		}		
		
		/*get the estimated 'true' relative effects and variances from the function random() - end of code - and store them in the dataset*/
		mata random()
		
		cap drop _es _sees
		mat _estot=_estot'
		mat _vartotnew=_vartotnew'
		tempvar _rank _use
		qui gen `_rank'=_n
		qui gen `_use'=1 `if' `in'
		qui sort `_use'
		qui svmat _estot, name(_es)
		qui rename _es1 _es
		qui svmat _vartotnew
		qui gen _sees=sqrt(_vartotnew)
		cap drop _vartotnew
		qui sort `_rank'
			
	}
	
	/*compare the results of different assumptions for the missing parameters*/
	if `"`compare'"'!=""{
		global varlist `varlist'
		global if `if'
		global in `in'
		global compare `compare'
		global options `options'
		
		compare `varlist'
	}
	
	/*run a sensitivity analysis using a range of standard deviations for the missing parameter*/
	if "`sensitivity'"!=""{
		global model `model'
		global measure `measure'
		global reps `reps'
		global imdom `imdom'
		global imrom `imrom'
		global rho `rho'
		global method `method'
		global varlist `varlist'
		
		sensitivity `varlist'
	}
	
	else{
		/*run a meta-analysis on the estimated 'true' relative effects*/
		if "`fixed'"!=""{
			metan _es _sees, fixed `options'
		}
		else{
			metan _es _sees, random `options'
		}
		if "`nokeep'"==""{
			qui rename _es _ES
			qui rename _sees _seES
			if "`measure'"=="md"{
				label var _ES MD
				label var _seES "se(MD)"
			}
			if "`measure'"=="smd"{
				label var _ES SMD
				label var _seES "se(SMD)"
			}	
			if "`measure'"=="rom"{
				label var _ES lnRoM
				label var _seES "se(lnRoM)"
			}			
		}
		else{
			cap drop _es _sees
			cap drop _LCI _UCI _SS _WT
		}
	}
	cap scalar drop Mreps Mns Mod Meas
	cap mat drop Mnc Mmc Mnt Mmt Myc Myt Msdc Msdt Mmulc Msdlc Mmult Msdlt Mrho
end


/*generate random samples for the bootstrap method*/	
mata
void random()
{	
	
	/*transfer the observed values into mata*/
	mulambdac=st_matrix("Mmulc")
	sdlambdac=st_matrix("Msdlc")
	mulambdat=st_matrix("Mmult")
	sdlambdat=st_matrix("Msdlt")	
	reps=st_numscalar("Mreps")
	ns=st_numscalar("Mns")
	mod=st_numscalar("Mod")
	meas=st_numscalar("Meas")	
	rho=st_matrix("Mrho")
	yc=st_matrix("Myc")
	yt=st_matrix("Myt")	
	nc=st_matrix("Mnc")
	nt=st_matrix("Mnt")	
	mc=st_matrix("Mmc")
	mt=st_matrix("Mmt")		
	sdt=st_matrix("Msdt")
	sdc=st_matrix("Msdc")	
	ntott=nt+mt
	ntotc=nc+mc
	pc=nc
	pt=nt
	for (i=1; i<=cols(pc); i++){
		pc[1,i]=nc[1,i]/ntotc[1,i]
		pt[1,i]=nt[1,i]/ntott[1,i]
	}
	sigmac=sdc
	sigmat=sdt
	for (i=1; i<=cols(sigmac); i++){
		sigmac[1,i] = sdc[1,i]/sqrt(nc[1,i])
		sigmat[1,i] = sdt[1,i]/sqrt(nt[1,i])
	}	

	/*generate random samples*/
	
	LC=rnormal(reps,1,mulambdac,sdlambdac) //bivariate normal distribution on lambdaT and lambdaC
	Mlambdat=LC
	Slambdat=LC
	for (i=1; i<=rows(Mlambdat); i++){
		for (j=1; j<=cols(Mlambdat); j++){
			Mlambdat[i,j]=mulambdat[1,j]+rho[1,j]*(LC[i,j]-mulambdac[1,j])
			Slambdat[i,j]=sdlambdat[1,j]*sqrt(1-rho[1,j]*rho[1,j])
		}
	}
	LT=rnormal(1,1,Mlambdat,Slambdat)	

	PC=rbeta(reps,1,nc,mc) //beta distribution on the number of observations
	PT=rbeta(reps,1,nt,mt)	

	df=ns-1
	TC=rt(reps,ns,df) //t-distribution (with df=#studies-1) on means from observed participants
	TT=rt(reps,ns,df)
	
	YC=TC
	YT=TT
	for (i=1; i<=rows(TC); i++){
		for (j=1; j<=cols(TC); j++){	
			YC[i,j]=TC[i,j]*sqrt(sdc[1,j]^2/nc[1,j])+yc[1,j]
			YT[i,j]=TT[i,j]*sqrt(sdt[1,j]^2/nt[1,j])+yt[1,j]
		}
	}

	if (mod==2){
		YMC=YC
		YMT=YT
		for (i=1; i<=rows(YMC); i++){
			for (j=1; j<=cols(YMC); j++){
				YMC[i,j]=exp(LC[i,j])*YC[i,j] //estimate the 'true' means from missing participants for treatment and control arm from IMRoM model
				YMT[i,j]=exp(LT[i,j])*YT[i,j]
			}	
		}
	}
	if (mod==1){	
		YMC=LC+YC //estimate the 'true' means from missing participants for treatment and control arm from IMDoM model
		YMT=LT+YT	
	}
	
	YCtrue=YMC
	YTtrue=YMT
	for (i=1; i<=rows(YCtrue); i++){
		for (j=1; j<=cols(YCtrue); j++){
			YCtrue[i,j] = PC[i,j]*YC[i,j]+(1-PC[i,j])*YMC[i,j] //estimate the 'true' means from all participants for treatment and control arm
			YTtrue[i,j] = PT[i,j]*YT[i,j]+(1-PT[i,j])*YMT[i,j]
		}
	}	
	if (meas==1){
		Ytot=YTtrue-YCtrue //estimate of 'true' MD
	}
	if (meas==2){
		spooled=YTtrue
		Ytot=YTtrue-YCtrue
		for (i=1; i<=rows(spooled); i++){
			for (j=1; j<=cols(spooled); j++){			
				spooled[i,j]=sqrt(((NT[i,j]-1)*sdt[1,j]^2+(NC[i,j]-1)*sdc[1,j]^2)/(NT[i,j]+NC[i,j]-2)) //estimate of Spooled
			}
		}
		for (i=1; i<=rows(Ytot); i++){
			for (j=1; j<=cols(Ytot); j++){		
				Ytot[i,j]=(YTtrue[i,j]-YCtrue[i,j])/spooled[i,j] //estimate of 'true' SMD
			}
		}		
	}
	if (meas==3){
		Ytot=YTtrue
		for (i=1; i<=rows(Ytot); i++){
			for (j=1; j<=cols(Ytot); j++){			
				Ytot[i,j]=log(abs(YTtrue[i,j])/abs(YCtrue[i,j])) //estimate of 'true' RoM
			}
		}
	}
	/*estimate for each study the means and variances of the 'true' relative effects*/	
	estot=mean(Ytot)
	vartot=variance(Ytot)
	vartotnew=J(1,ns,0)	
	for (i=1; i<=cols(vartotnew); i++){ //keep the diagonals of the variance-covariance matrices
		vartotnew[1,i]=vartot[i,i]	
	}

	/*transfer the estimated means and variances to Stata*/
	st_matrix("_estot",estot)
	st_matrix("_vartotnew",vartotnew)
}
end
	
	
program def compare
	syntax varlist(min=8 max=8) [if] [in]
	
	tempvar double
	qui expand 2 $if $in,gen(`double')
	lab def _anal 0 "Primary analysis" 1 "Secondary analysis"
	lab val `double' _anal
	
	metamiss2 $varlist $if $in, notable nograph $compare by(`double') nooverall
	cap rename _ES _second
	cap rename _seES _sesecond
	
	metamiss2 $varlist $if $in, notable nograph `md' `smd' `rom' `taylor' `bootstrap' reps(`reps') imdom(`imdom') imrom(`imrom') rho(`rho') sdimdom(`sdimdom') sdimrom(`sdimrom') `options' by(`double') nooverall

	cap replace _ES=_second if `double'==1
	cap replace _seES=_sesecond if `double'==1
	cap rename _ES _es
	cap rename _seES _sees
	metan _es _sees, random `options' by(`double') nooverall
	cap rename _es _ES
	cap rename _sees _seES
	cap rename _second _ESsecond
	cap rename _sesecond _seESsecond
	cap drop if `double'==1
	cap lab drop _anal
	
end
	
	
program def sensitivity
	syntax varlist(min=8 max=8) [if] [in]
	
	local nobs=`=_N'
	cap set obs 100
	
	tempvar _sdlambda
	
	qui gen `_sdlambda' in 1=0

	if "$model"=="imdom"{
		forvalues i=2/100{
			qui replace `_sdlambda' in `i'=`_sdlambda'[`=`i'-1']+5/99
			local sdlambda`i'=`_sdlambda'[`i']
		}
	}
	if "$model"=="imrom"{
		forvalues i=2/100{
			qui replace `_sdlambda' in `i'=`_sdlambda'[`=`i'-1']+0.5/99
			local sdlambda`i'=`_sdlambda'[`i']
		}
	}
	qui metamiss2 $varlist $if $in, notable nokeep nograph "$measure" imdom("$imdom") imrom("$imrom")
	local ES1=r(ES)
	local LCI1=r(ci_low)
	local UCI1=r(ci_upp)	
	forvalues i=2/100{
		if "$model" =="imdom"{
			local sd `"sdimdom(`sdlambda`i'')"'
		}
		if "$model" =="imrom"{
			local sd `"sdimrom(`sdlambda`i'')"'
		}	
		if "$method" =="taylor"{
			qui metamiss2 $varlist $if $in, notable nokeep nograph "$measure" `sd' imdom("$imdom") imrom("$imrom") rho("$rho")
		}
		if "$method" =="bootstrap"{
			qui metamiss2 $varlist $if $in, notable nokeep nograph "$measure" b reps("$reps") `sd' imdom("$imdom") imrom("$imrom") rho("$rho")
		}
		local ES`i'=r(ES)
		local LCI`i'=r(ci_low)
		local UCI`i'=r(ci_upp)
	}
	
	tempvar _ESlambda _LCIlambda _UCIlambda
	qui{
		gen `_ESlambda' in 1=`ES1'
		gen `_LCIlambda' in 1=`LCI1'
		gen `_UCIlambda' in 1=`UCI1'

		forvalues i=2/100{
			replace `_ESlambda' in `i'=`ES`i''
			replace `_LCIlambda' in `i'=`LCI`i''
			replace `_UCIlambda' in `i'=`UCI`i''	
		}
	}

	if "`measure'"=="md"{
		local ytitle "Summary effect - Mean difference"
	}
	if "`measure'"=="smd"{
		local ytitle "Summary effect - Standardized mean difference"
	}
	if "`measure'"=="rom"{
		local ytitle "Summary effect - Ratio of means"
	}
	if "`model'"=="imdom"{
		local xtitle "IMDoM"
	}
	if "`model'"=="rom"{
		local xtitle "IMRoM"
	}		
	sc `_ESlambda' `_sdlambda', mcol(black)|| sc `_LCIlambda' `_sdlambda',mcol(black) msymb(plus) || sc `_UCIlambda' `_sdlambda', mcol(black) msymb(plus) legend(off) xtitle("Standard deviation of lambda - `xtitle' model") ytitle(`ytitle')
	cap drop if _n>`nobs'
	cap drop _LCI _UCI _WT
end	
	
	
	
	
	
