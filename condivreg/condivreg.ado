*! Version 2.0.4  17apr2006
* version 1 Marcelo Moreira and Brian Poi
* Estimated IV regression by 2sls or liml
* version 1.0.0 created 20021105 by Brian Poi
* version 2.0.0 created 20051001 by Anna Mikusheva
* version 2.0.1 20051205 code check and reformat for SJ inclusion by BPP
* version 2.0.2 changed final display format
* version 2.0.3 bug fixes
* version 2.0.4 bug fix regarding first-stage F statistic

program define condivreg, eclass byable(recall)

	version 8.2
	
	if replay() {
		if "`e(cmd)'" != "condivreg" {
			error 301
		}
		syntax [, level(integer -1) ]
		if `level' > -1 {
			di as error "cannot specify {cmd:level()} on replay"
			exit 198
		}
		myprint
		myprintregions
		exit
	}
		 
	/* First get the list of variables.	*/
	local n 0
	gettoken lhs 0 : 0, parse(" ,[") match(paren)
	IsStop `lhs'
	if `s(stop)' { 
		error 198 
	}
	while `s(stop)'==0 {
		if "`paren'"=="(" {
			local n = `n' + 1
			gettoken p lhs : lhs, parse(" =")
			while "`p'"!="=" {
				local end`n' `end`n'' `p'
				gettoken p lhs : lhs, parse(" =")
			}
			tsunab end`n' : `end`n''
			tsunab rinst : `lhs'
		}
		else {
			local exog `exog' `lhs'
		}
		gettoken lhs 0 : 0, parse(" ,[") match(paren)
		IsStop `lhs'
	}
	local 0 `"`lhs' `0'"'

	tsunab exog : `exog'
	tokenize `exog'
	local ry1 "`1'"
	local 1 " "
	local rexog `*'
	loc ry2 "`end1'"
	syntax [if] [in], [nocons noinstcons liml 2sls ar lm interval 	/*
		*/	   level(integer $S_level) test(real 0)]
	
	marksample touse
	markout `touse' `ry1' `ry2' `rexog' `rinst'

/* Some syntax checking.	*/
	if ("`liml'" != "" & "`2sls'" != "") {
		di as error "Only one of liml or 2sls is allowed."
		exit 198
	}
	if ("`liml'" != "") {
		local est_type "liml"
	}
	else {
		local est_type "2sls"
	}

	quietly {
/* check for multicollinearity */
		tempvar one  AAA  vv dd 
		loc k : word count `rinst'
		loc p : word count `rexog'
		gen `one' = 1
	 	loc XZ "`rexog' `rinst'"
		if ("`cons'" == "") {
			mat accum `AAA' = `XZ'
			loc k = `k' + 1
		}
		else{
			mat accum `AAA' = `XZ', nocons
		}
		mat symeigen `vv' `dd' = `AAA'
		if (`dd'[1,`k'+`p']<0.000000001) {
			noi di as error "Multicollinearity!"
			exit
		}

		tempname bbb 
		sca `bbb'=`test'
		if ("`ar'"!="") { 
			loc AR "ar" 
		}
		if ("`lm'"!="") { 
			loc LM "lm" 
		}
 
/* Compute the 2sls and liml results.	*/
		tempname beta var rmse rss mss rsquared adjrsq fstat 	/*
			*/	xpx limlbeta
		myliml `ry1' `ry2' "`rexog'" "`rinst'" "`cons'" 	/*
			*/	"`instcons'" `touse'
		mat `beta' = r(beta)
		sca `limlbeta' = `beta'[1,1]

		if "`est_type'" == "2sls" {
			my2sls `ry1' `ry2' "`rexog'" "`rinst'" "`cons'"	/*
				*/	"`instcons'" `touse'  
		}
		mat `beta' = r(beta)
		mat `var' = r(var)
		if "`instcons'" == "" {
			mat colnames `beta' = `ry2' `rexog' _cons
			mat rownames `var' = `ry2' `rexog' _cons
			mat colnames `var' = `ry2' `rexog' _cons
								
		}
		else {
			mat colnames `beta' = `ry2' `rexog'
			mat rownames `var' = `ry2' `rexog'
			mat colnames `var' = `ry2' `rexog'
		}
		sca `rss' = r(rss)
/* Compute the general stuff that goes along. */
		tempvar junk
		gen double `junk' = `ry1'^2 if `touse'
		su `junk', meanonly
		loc capn = r(N)
		sca `mss' = r(sum) - `rss'
		if ("`cons'" == "") {
			su `ry1' if `touse', meanonly
			sca `mss' = `mss' - r(sum)^2 / `capn'
		}
		loc dfm : word count `rexog' `ry2'
		loc dfr = `capn' - `dfm'
		if ("`cons'" == "") {
			loc c = 1
			loc dfr = `dfr' - 1
		}
		else {
			loc c = 0
		}
		sca `rmse' = sqrt(`rss'/`dfr')
		sca `rsquared' = `mss' / (`mss' + `rss')
		sca `adjrsq' = 1 - (1-`rsquared')*(`capn'-`c') /	/*
			*/		(`capn'-`c'-`dfm')
		/* F statistic for 2SLS/Wald for LIML. */
		if "`cons'" != "" {
			sca `fstat' = .
		}
		else {
			tempname c dif vinv
			mat `c' = J(1, (`dfm' + 1), 0)
			su `ry1' if `touse', meanonly
			mat `c'[1, (`dfm'+1)] = r(mean) 
			mat `dif' = `beta' - `c'
			mat `vinv' = inv(`var')
			mat `fstat' = `dif'*`vinv'*`dif''
			sca `fstat' = trace(`fstat')
			if "`est_type'" == "2sls" {
				sca `fstat' = `fstat' / `dfm'
			}
		}

/* Get the first-stage stats.	*/
		tempname f1 r21 r2_a1 df_m1 df_r1
		if "`instcons'" == "" {
			reg `ry2' `rexog' `rinst'
		}
		else {
			reg `ry2' `rexog' `rinst', nocons
		}
		test `rinst'
		sca `f1' = r(F)
		sca `r21' = e(r2)
		sca `r2_a1' = e(r2_a)
		sca `df_m1' = r(df)
		sca `df_r1' = r(df_r)
		
/* Post results.  Use two different syntaxes of 2sls and liml
   so that we get t-stat for 2sls and z-stat for liml. */
   		tempvar touse2
   		gen byte `touse2' = `touse'	/* gets destroyed by eret post */
		if "`est_type'" == "2sls" {
			eret post `beta' `var', esample(`touse2') 	/*
				*/	depname("`ry1'") obs(`capn') dof(`dfr')
		}
		else {
			eret post `beta' `var', esample(`touse2') 	/*
				*/	depname("`ry1'") obs(`capn')
		} 
		eret local depvar "`ry1'"
		eret scalar df_m = `dfm'
		eret scalar df_r = `dfr'
		if "`est_type'" == "2sls" {
			eret scalar F = `fstat'
		}
		else {
			eret scalar wald = `fstat'
		}
		eret scalar H0_b = `bbb'
		eret scalar r2 = `rsquared'
		eret scalar r2_a = `adjrsq'
		eret scalar rmse = `rmse'
		eret scalar mss = `mss'
		eret scalar rss = `rss'
		eret scalar F_first = `f1'
		eret scalar df_m_first = `df_m1'
		eret scalar df_r_first = `df_r1'
		eret scalar r2_first = `r21'
		eret scalar r2_a_first = `r2_a1'
		eret scalar b_instd_liml = `limlbeta'
		eret scalar level = `level'
		eret local exog "`rexog'"
		eret local inst "`rinst'"
		eret local insts "`rexog' `rinst'"
		eret local instd "`ry2'"
		if "`est_type'" == "2sls" {
			eret local model "2sls"
		}
		else {
			eret local model "liml"
		}
		if "`instcons'" == "" {
			eret local instcons "yes"
		}
		else {
			eret local instcons "no"
		}
		if "`cons'" == "" {
			eret local cons "yes"
		}
		else {
			eret local cons "no"
		}

/* Common block for testing and conf sets */
 		if "`cons'" == "" {
			loc exog "`rexog' `one'"
		}
		else {
			loc exog "`rexog'"
		}
/* preserve currently posted results; -reg- will clobber them otherwise */
		tempname myest
		_estimates hold `myest'
		tempvar y1 y2
		if ("`exog'" != "") {
			foreach v in y1 y2 {
				reg `r`v'' `exog' if `touse', nocons
				predict double ``v'' if `touse', residuals
			}
		}
		else {
			gen double `y1' = `ry1' if `touse'
			gen double `y2' = `ry2' if `touse'
		}
/* Regress instruments on exog.	*/
		tempname ehold
		loc inst = ""
		loc n = 1
		foreach v in `rinst' {
			tempvar inst`n'
			if ("`exog'" != "") {
				reg `v' `exog' if `touse', nocons
				predict double `inst`n'' if `touse', residuals
			}
			else {
				gen double `inst`n'' = `v' if `touse'
			}
			loc inst "`inst' `inst`n''"
		}
/* Compute Omega. */
		tempname mzy1 mzy2 n k omega
		if "`instcons'" == "" {
			reg `y1' `inst' if `touse'
		}
		else {
			reg `y1' `inst' if `touse', nocons
		}
		predict `mzy1' if `touse', residuals
		if "`instcons'" == "" {
			reg `y2' `inst' if `touse'
		}
		else {
			reg `y2' `inst' if `touse', nocons
		}
		predict `mzy2' if `touse', residuals
		mat accum `omega' = `mzy1' `mzy2', noconstant
/* Safe to restore my results */
		_estimates unhold `myest'
		count if `touse'
		loc n = r(N)
		loc k : word count `inst'
		loc p : word count `exog'
		
		if ("`instcons'" == "") {
			loc k = `k' + 1
		}
		mat `omega' = `omega' / (`n'-`k'-`p')
		if (("`instcons'" == "")& ("`cons'" == "")) {
			loc df = `k' - 1	
		}
		else {
			loc df = `k' 
		}

/*------Confidence sets -------------*/

/*
Test		Result type		Interval
-----------------------------------------------------------------------
CLR		1			Empty set
		2			[x1, x2]
		3			(-infty, +infty)
		4		    (-infty, x1] U [x2, infty)
		
AR		1			Empty set
		2			[x1, x2]
		3			(-infty, +infty)
		4		    (-infty, x1] U [x2, infty)
		
LM		1			Not used (not possible)
                2			[x1, x2]                                
		3			(-infty, +infty)
		4		(-infty, x1] U [x2, infty)
		5		(-infty, x1] U [x2, x3] U [x4, infty)
		6		    [x1, x2] U [x3, x4]

-----------------------------------------------------------------------
NB: BPP modified the type codes used by AR and LM in the prequel to 
    match the ones used by CLR.
*/
		
		tempname cross zpz sqrtzpzi zpy MM ypz sqrtomegai v 	/*
			*/	d M N alpha C A D aa x1 x2 g type
		if ("`instcons'" == "") {
			mat accum `cross' = `inst' `one' `y1' `y2' 	/*
				*/	if `touse', noconstant
		}
		else {
			mat accum `cross' = `inst' `y1' `y2' 		/*
				*/	if `touse', noconstant
		}
	  
		mat `zpz' = `cross'[1..`k', 1..`k']
		mat `zpy' = `cross'[1..`k', (`k'+1)...]
		mat_inv_sqrt `zpz' `sqrtzpzi'
		mat_inv_sqrt `omega' `sqrtomegai'
		mat `ypz'=`zpy''
		mat `MM' = `sqrtomegai'*`ypz'*inv(`zpz')*`zpy'*`sqrtomegai'
		loc `k'=`df'
		mat symeigen `v' `d' = `MM'
		sca `M' = `d'[1,1]
		sca `N' =`d'[1,2]
		sca `alpha' = 1-`level'/100

/* inversion of CLR*/
		inversefun `M' `df' `alpha' `C'

		mat `A' =inv(`omega')*`ypz'*inv(`zpz')*`zpy'*inv(`omega')- /*
			*/	`C'*inv(`omega')
		sca `D' = -det(`A')
		sca `aa' = `A'[1,1]
		if (`aa'<0) {
	 		if (`D' <0) {
				sca `type'=1
			}
	 		else{
				sca `type'=2
				sca `x1'= (-`A'[1,2] + sqrt(`D'))/`aa'
				sca `x2' = (-`A'[1,2] - sqrt(`D'))/`aa'
				mat `g'=(`x1'\ `x2')
 			}
 		}
		else{
			if (`D'<0) {
		  		sca `type'=3
			}
	 		else {
		  		sca `type'=4
		  		sca `x1'= (-`A'[1,2]-sqrt(`D'))/`aa'
		  		sca `x2'= (-`A'[1,2]+sqrt(`D'))/`aa'
		  		mat `g'=(`x1' \ `x2')
 			}
	 	}
	 	ereturn local LR_type `"`=`type''"'
	 	if `type' == 2 | `type' == 4 {
	 		eret scalar LR_x1 = `x1'
	 		eret scalar LR_x2 = `x2'
	 	}

/* inversion of LM */

if "`LM'" != ""	{		/* BPP add 20051220 */
		tempname lmcv q1 q2 A1 A2 D1 D2 y1 y2 y3 y4 type1
		sca `lmcv' = invchi2tail(1, (1-`level'/100))
		if (`df'==1) {
			sca `q1' = `M'-`lmcv'
			mat `A1' =inv(`omega')*`ypz'*inv(`zpz')*`zpy'*	/*
				*/	inv(`omega')-`q1'*inv(`omega')
			sca `D1' = -4*det(`A1')
			sca `y1'= (-2*`A1'[1,2]+ sqrt(`D1'))/2/`A1'[1,1]
			sca `y2'= (-2*`A1'[1,2]-sqrt(`D1'))/2/`A1'[1,1]
			if (`A1'[1,1]>0) { 
				if (`D1'>0) {
					sca `type1'=4 
						/* two infinite intervals*/
					eret scalar LM_x1 = `y2'
					eret scalar LM_x2 = `y1'
				}
				else {
					sca `type1'=3
				}
			}
			else{
				if (`D1'>0) {
					sca `type1'=2 /*one interval */
					eret scalar LM_x1 = `y1'
					eret scalar LM_x2 = `y2'

				}
				else {
					sca `type1'=3
				}
			}
		}
		else {
			if ((`M' +`N' - `lmcv')^2-4*`M'*`N'<0) { 
		 		sca `type1' = 3
		 	}
			else {
			    sca `q1' = (`M'+ `N' - `lmcv' - 	/*
	 			*/ sqrt((`M'+`N' - `lmcv')^2 -  /*
	 			*/	4*`M'*`N'))/2
	 		    sca `q2' = (`M'+`N' - `lmcv' + 		/*
	 			*/ sqrt((`M'+`N'-`lmcv')^2 - 	/*
	 			*/	4*`M'*`N'))/2
	 		    if ((`q1' < `N') | (`q2' > `M')) {
	 			sca `type1' = 3
	 		    }
	 		    else { 		
				mat `A1' = inv(`omega')*`ypz'*inv(`zpz')*   /*
					*/ `zpy'*inv(`omega')-`q1'*inv(`omega')
				mat `A2' = inv(`omega')*`ypz'*inv(`zpz')*   /*
					*/ `zpy'*inv(`omega')-`q2'*inv(`omega')
		 		sca `D1' = -4*det(`A1')
				sca `D2' = -4*det(`A2')
	 			if (`A1'[1,1]>0) { 
		  			if (`A2'[1,1]>0) { 
						sca `type1' = 5
						sca `y1' = (-2*`A1'[1,2] +  /*
						    */ sqrt(`D1'))/2/`A1'[1,1]
						sca `y2' = (-2*`A1'[1,2] -  /*
						    */ sqrt(`D1'))/2/`A1'[1,1]
						sca `y3' = (-2*`A2'[1,2] +  /*
						    */ sqrt(`D2'))/2/`A2'[1,1]
						sca `y4' = (-2*`A2'[1,2] -  /*
						    */ sqrt(`D2'))/2/`A2'[1,1]
						eret scalar LM_x1 = `y4'
						eret scalar LM_x2 = `y2'
						eret scalar LM_x3 = `y1'
						eret scalar LM_x4 = `y3'
					}
		  			else {
						sca `type1' = 6
						sca `y1' = (-2*`A1'[1,2] +  /*
						    */ sqrt(`D1'))/2/`A1'[1,1]
						sca `y2' = (-2*`A1'[1,2] -  /*
						    */ sqrt(`D1'))/2/`A1'[1,1]
						sca `y3' = (-2*`A2'[1,2] +  /*
						    */ sqrt(`D2'))/2/`A2'[1,1]
						sca `y4' = (-2*`A2'[1,2] -  /*
						    */ sqrt(`D2'))/2/`A2'[1,1]
						eret scalar LM_x1 = `y3'
						eret scalar LM_x2 = `y4'
						eret scalar LM_x3 = `y2'
						eret scalar LM_x4 = `y1'
			  		}
			  	}
				if (`A1'[1,1]<=0) {
					sca `type1' =5
			  		sca `y1' = (-2*`A1'[1,2] + 	/*
			  			*/  sqrt(`D1'))/2/`A1'[1,1]
					sca `y2' = (-2*`A1'[1,2] - 	/*
						*/  sqrt(`D1'))/2/`A1'[1,1]
					sca `y3' = (-2*`A2'[1,2] + 	/*
						*/  sqrt(`D2'))/2/`A2'[1,1]
					sca `y4' = (-2*`A2'[1,2] -	/*
						*/  sqrt(`D2'))/2/`A2'[1,1]
					eret scalar LM_x1 = `y1'
					eret scalar LM_x2 = `y3'
					eret scalar LM_x3 = `y4'
					eret scalar LM_x4 = `y2'
		  		}
			    }
			}
		}
		eret local LM_type `"`=`type1''"'
}	/* End of LM block */

/* inversion of AR */
if "`AR'" != ""	{		/* BPP add 20051220 */
		tempname lmcv1  AAA type2 xx1 xx2 DDD aaa
		sca `lmcv1' = invchi2tail(`df', (1-`level'/100))
		mat `AAA' =`ypz'*inv(`zpz')*`zpy'-`lmcv1'*`omega'
		sca `DDD' = -det(`AAA')
		sca `aaa' = `AAA'[2,2]
		if (`aaa'<0) {
	 		if (`DDD' <0) {
				sca `type2'=3
			}
	 		else{
				sca `type2'=4
		 		sca `xx1'= (`AAA'[1,2] + sqrt(`DDD'))/`aaa'
		 		sca `xx2' = (`AAA'[1,2] - sqrt(`DDD'))/`aaa'
				eret scalar AR_x1 = `xx1'
				eret scalar AR_x2 = `xx2'
 			 }
		}
		else {
			if (`DDD'<0) {
				sca `type2'=1
			}
	 		else {
		  		sca `type2'=2
		  		sca `xx1'= (`AAA'[1,2]-sqrt(`DDD'))/`aaa'
		  		sca `xx2'= (`AAA'[1,2]+sqrt(`DDD'))/`aaa'
				eret scalar AR_x1 = `xx1'
				eret scalar AR_x2 = `xx2'
 			}
	 	}
		eret local AR_type `"`=`type2''"'
}	/* End of AR block */

	} /*---end of quiet--*/

/*------Testing block-----*/
	tempname  a b aprime bprime	
	
	/* Set up a and b vectors.	*/
	mat `a' = (`bbb'\1)
	mat `b' = (1\(-1*`bbb'))
	mat `aprime' = `a''		/* Macro parser chokes 	*/
	mat `bprime' = `b''		/* otherwise.		*/

 	tempname oia
	mat `oia' = inv(`omega')*`a'
	/* Compute scalars b'*O*b and a'*O^-1*a	*/
	tempname bob aoia
	matrix `bob' = `bprime'*`omega'*`b'
	scalar `bob' = trace(`bob')
	matrix `aoia' = `aprime'*inv(`omega')*`a'
	scalar `aoia' = trace(`aoia')
	
	tempname  sbar tbar
	mat `sbar' = `sqrtzpzi'*`zpy'*`b'/sqrt(`bob')
	mat `tbar' = `sqrtzpzi'*`zpy'*`oia'/sqrt(`aoia')
 
	tempname ar lm lr wald qt pval_new
	calcstat1 `sbar' `tbar' `ar' `lm' `qt'
	calcstat2 `sbar' `tbar' `aoia' `bob' `bbb' `omega' `lr' 
	loc kk=`k'-1

	if (("`instcons'" == "")&("`cons'" == "")) {
		new_try `kk' `qt' `lr' `pval_new'
	}
	else {
		new_try `k' `qt' `lr' `pval_new'
	}
		
	if "`AR'" != "" {
		tempname arcv arpv 
		if (("`instcons'" == "")&("`cons'" == "")) {
			sca `arcv' = invchi2tail((`k'-1), (1-`level'/100))
			sca `arpv' = 1-chi2((`k'-1), `ar')
		}
		else {
			sca `arcv' = invchi2tail(`k', (1-`level'/100))
			sca `arpv' = 1-chi2(`k', `ar')	
		}
		eret scalar AR_p = `arpv'
	}
	if "`LM'" != "" {
		tempname lmpv
		sca `lmpv' = 1-chi2(1, `lm')
		eret scalar LM_p = `lmpv'
	}
	
	eret scalar LR_p = `pval_new'

	eret local cmd "condivreg"

	myprint		/* Level stored in e(level) */
	myprintregions "`interval'" "`AR'" "`LM'"
		
end
  
	

/*----------------------------------------------------------------------*/

/* Computes AR and LM statistics.	*/
prog def calcstat1

	args sbar tbar ar lm qt
	tempname sbarp tbarp
	
	mat `sbarp' = `sbar''
	mat `tbarp' = `tbar''	
	mat `ar' = `sbarp'*`sbar'
	sca `ar' = trace(`ar')
	mat `lm' = (trace(`sbarp'*`tbar')^2) / trace(`tbarp'*`tbar')
	sca `lm' = trace(`lm')
	mat `qt' = `tbarp'*`tbar'
	sca `qt' = trace(`qt')

	
end 

/* Computes LR and Wald statistics.	*/
prog def calcstat2

	args sbar tbar aoia bob beta omega lr 
	tempname sbarp tbarp
	
	mat `sbarp' = `sbar''
	mat `tbarp' = `tbar''		
	tempname ss st tt
	mat `ss' = `sbarp'*`sbar'
	sca `ss' = trace(`ss')
	mat `tt' = `tbarp'*`tbar'
	sca `tt' = trace(`tt')
	mat `st' = `sbarp'*`tbar'
	sca `st' = trace(`st')
	sca `lr' = 0.5*(`ss' - `tt' + sqrt((`ss' + `tt')^2 -	/*
			*/	4*(`ss'*`tt' - (`st')^2)))	
  
end



program define new_try

	args k qt lrstat pval_new

	tempname gamma pval  u s2 qs farg1 farg2 farg wt
	
	sca `gamma' = 2*exp(lngamma(`k'/2)) / sqrt(_pi) / exp(lngamma((`k'-1)/2))

	if("`k'" == "1") {
		sca `pval' = 1 - chi2(`k', `lrstat')
	}
	else if ("`k'"== "2") {
		local ni 20
		mat `u' = J(`ni'+1,1,0)
		mat `s2' = J(`ni'+1,1,0)
		mat `qs' = J(`ni'+1,1,0)
		mat `wt' = J(1,`ni'+1,2)
		mat `farg1' = J(`ni'+1,1,0)
		mat `qs'[1,1] = (`qt'+`lrstat')
		mat `farg1'[1,1] = `gamma'*chi2(`k',`qs'[1,1])
		forv i =1(1)`ni'{
			mat `u'[`i'+1,1] = `i'*_pi/2/`ni'
			mat `s2'[`i'+1,1] = sin(`u'[`i'+1,1])
			mat `qs'[`i'+1,1] = (`qt'+`lrstat') / 	/*
				*/ (1+(`qt'/`lrstat')*`s2'[`i'+1,1]*`s2'[`i'+1,1])
			mat `farg1'[`i'+1,1] = `gamma'*chi2(`k',`qs'[`i'+1,1])
		}
		mat `wt'[1,1] = 1
		mat `wt'[1,`ni'+1] = 1
		local ni = `ni'/2
		forv i =1(1)`ni'{
			mat `wt'[1,`i'*2] = 4
		}
		local ni = `ni'*2
		mat `wt' = `wt'*_pi/2/3/`ni'
		mat `pval' = `wt'*`farg1'
		sca `pval' = 1-trace(`pval')
	}
	else if ("`k'"== "3") {
		local ni 20
		mat `s2' = J(`ni'+1,1,0)
		mat `qs' = J(`ni'+1,1,0)
		mat `wt' = J(1,`ni'+1,2)
		mat `farg1' = J(`ni'+1,1,0)
		mat `qs'[1,1] = (`qt'+`lrstat')
		mat `farg1'[1,1] = `gamma'*chi2(`k',`qs'[1,1])
		forv i =1(1)`ni'{
			mat `s2'[`i'+1,1] = `i'/`ni'
			mat `qs'[`i'+1,1] = (`qt'+`lrstat') / 	/*
				*/ (1+(`qt'/`lrstat')*`s2'[`i'+1,1]*`s2'[`i'+1,1])
			mat `farg1'[`i'+1,1] = `gamma'*chi2(`k',`qs'[`i'+1,1])
		}
		mat `wt'[1,1] = 1
		mat `wt'[1,`ni'+1] = 1
		local ni = `ni'/2
		forv i =1(1)`ni'{
			mat `wt'[1,`i'*2] = 4
		}
		local ni = `ni'*2
		mat `wt' = `wt'/3/`ni'
		mat `pval' = `wt'*`farg1'
		sca `pval' = 1-trace(`pval')
	}
	else if ("`k'"== "4") {
		local eps .02
		local ni 50
		mat `s2' = J(`ni'+1,1,0)
		mat `qs' = J(`ni'+1,1,0)
		mat `wt' = J(1,`ni'+1,2)
		mat `farg' = J(`ni'+1,1,0)
		mat `farg1' = J(`ni'+1,1,0)
		mat `farg2' = J(`ni'+1,1,1)
		mat `qs'[1,1] = (`qt'+`lrstat')
		mat `farg1'[1,1] = `gamma'*chi2(`k',`qs'[1,1])
		mat `farg'[1,1] = `farg1'[1,1]*`farg2'[1,1]
		forv i = 1(1)`ni'{
			mat `s2'[`i'+1,1] = `i'/`ni'*(1-`eps')
			mat `qs'[`i'+1,1] = (`qt'+`lrstat') / 	/*
				*/ (1+(`qt'/`lrstat')*`s2'[`i'+1,1]*`s2'[`i'+1,1])
			mat `farg1'[`i'+1,1] = `gamma'*chi2(`k',`qs'[`i'+1,1])
			mat `farg2'[`i'+1,1] = sqrt(1-`s2'[`i'+1,1]*`s2'[`i'+1,1])
			mat `farg'[`i'+1,1] = `farg1'[`i'+1,1]*`farg2'[`i'+1,1]
		}
		mat `wt'[1,1] = 1
		mat `wt'[1,`ni'+1] = 1
		local ni = `ni'/2
		forv i = 1(1)`ni'{
			mat `wt'[1,`i'*2] = 4
		}
		local ni = `ni'*2
		mat `wt' = `wt'/3/`ni'*(1-`eps')
		mat `pval' = `wt'*`farg'
		sca `pval' = 1-trace(`pval')
		sca `s2' = 1-`eps'/2
		sca `qs' = (`qt'+`lrstat')/(1+(`qt'/`lrstat')*`s2'*`s2')
		sca `farg1' = `gamma'*chi2(`k',`qs')
		sca `farg2' = 0.5*(asin(1)-asin(1-`eps'))-(1-`eps') / /*
			*/	2*sqrt(1-(1-`eps')*(1-`eps'))
		sca `pval' = `pval'-`farg1'*`farg2'
	}
	else {
		local ni 20
		mat `s2' = J(`ni'+1,1,0)
		mat `qs' = J(`ni'+1,1,0)
		mat `wt' = J(1,`ni'+1,2)
		mat `farg' = J(`ni'+1,1,0)
		mat `farg1' = J(`ni'+1,1,0)
		mat `farg2' = J(`ni'+1,1,1)
		mat `qs'[1,1] = (`qt'+`lrstat')
		mat `farg1'[1,1] = `gamma'*chi2(`k',`qs'[1,1])
		mat `farg'[1,1] = `farg1'[1,1]*`farg2'[1,1]
		forv i =1(1)`ni'{
			mat `s2'[`i'+1,1] = `i'/`ni'
			mat `qs'[`i'+1,1] = (`qt'+`lrstat') / 	/*
				*/ (1+(`qt'/`lrstat')*`s2'[`i'+1,1]*`s2'[`i'+1,1])
			mat `farg1'[`i'+1,1] = `gamma'*chi2(`k',`qs'[`i'+1,1])
			if "`i'" == "`ni'"{
				mat `farg2'[`i'+1,1] = 0
			}
			else{
				mat `farg2'[`i'+1,1] = /*
				*/   (1-`s2'[`i'+1,1]*`s2'[`i'+1,1])^((`k'-3)/2)
			}
			mat `farg'[`i'+1,1] = `farg1'[`i'+1,1]*`farg2'[`i'+1,1]
		}
		mat `wt'[1,1] = 1
		mat `wt'[1,`ni'+1] = 1
		local ni = `ni'/2
		forv i = 1(1)`ni'{
			mat `wt'[1,`i'*2] = 4
		}
		local ni = `ni'*2
		mat `wt' = `wt'/3/`ni'
		mat `pval' = `wt'*`farg'
		sca `pval' = 1-trace(`pval')
	}
	
	sca `pval_new' = `pval'

end 


prog def inversefun

	args M k alpha C

	tempname eps a  b x fa fb lrstat fx

	sca `eps' = 0.000001
	sca `a' = `eps' 
	sca `b' = `M' - `eps'
	sca `lrstat'= `M' - `a'

	new_try `k' `a' `lrstat' `fa'
	sca `lrstat' = `M' - `b'
	new_try `k' `b' `lrstat' `fb'

	if(`fa' > `alpha') {
		sca `C' = `a'
	}
	else  if ( `fb' <`alpha') {
		sca `C' = `b'
	}
	else {
		while (`b'-`a'>`eps') {
			sca `x' = (`b'-`a')/2+`a'
			sca `lrstat'= `M'-`x'
			new_try `k' `x' `lrstat' `fx'
			if (`fx' >`alpha') {
			 	  sca `b' = `x'
			}
			else {
				sca `a' = `x'
			}				 
		}
		sca `C' = `x'
	}
	  
end



prog def mat_inv_sqrt

	args in out
	tempname v vpri lam srlam

	loc k = rowsof(`in')
	mat symeigen `v' `lam' = `in'
	mat `vpri' = `v''
	/* Get sqrt(lam)	  */
	mat `srlam' = diag(`lam')
	forv i = 1/`k' {
		mat `srlam'[`i', `i'] = 1/sqrt(`srlam'[`i', `i'])
	}
	mat `out' = `v'*`srlam'*`vpri'

end

program define my2sls, rclass

	args ry1 ry2 rexog rinst cons instcons touse

	quietly {
		tempname b2sls var2sls rss rssold
		tempvar y2hat
		if "`instcons'" == "" {
			reg `ry2' `rexog' `rinst' if `touse'
		}
		else {
			reg `ry2' `rexog' `rinst' if `touse', noconstant
		}
		predict double `y2hat' if `touse', xb
		reg `ry1' `y2hat' `rexog' if `touse', `cons'
		mat `b2sls' = e(b)
		mat `var2sls' = e(V)
		sca `rssold' = e(rss)
		tempvar resids hat
		replace `y2hat' = `ry2'
		predict double `resids' if `touse', residuals
		replace `resids' = `resids'^2
		su `resids', meanonly
		sca `rss' = r(sum)
		mat `var2sls' = `var2sls'*`rss'/`rssold'
		return scalar rss = `rss'
		return matrix beta `b2sls'
		return matrix var `var2sls'
	}

end

program define myliml, rclass

	args ry1 ry2 rexog rinst cons instcons touse

	quietly {
		/* LIML estimator.	*/
		/* Follows Davidson and MacKinnon (1993).	*/
		tempvar one
		gen double `one' = 1
		if "`instcons'" == "" {
			loc x2 `rinst' `one'
		}
		else {
			loc x2 `rinst'
		}
		if "`cons'" == "" {
			loc x1 `rexog' `one'
		}
		else {
			loc x1 `rexog'
		}
	
		tempname xpx x1px1 x1py1 y1py1
		mat accum `xpx' = `ry2' `x1' if `touse', noconstant
		mat `y1py1' = `xpx'[1, 1]

		/* Compute Y1'MxY1.	*/
		tempname y1mxy1
		loc k4 : word count "`ry2'"
		if `k4' != 1 {
			di as error /*
			*/ "Multiple endogenous RHS variables not implemented."
			exit 198
		}
		reg `ry2' `x1' `x2' if `touse', noconstant 
			/* May be two ones here, but reg will catch it. */
		tempvar y2resid
		predict double `y2resid' if `touse', residuals
		replace `y2resid' = `y2resid'^2
		summ `y2resid', meanonly
		sca `y1mxy1' = r(sum)
	
		/* Compute kappa.	*/
		/* First get Y'MxY.	*/
		tempvar y1h y2h y1hy2h
		tempname s11 s12 s22 ymxy ym1y
		reg `ry1' `x1' `x2' if `touse', noconstant
		predict double `y1h' if `touse', residuals
		reg `ry2' `x1' `x2' if `touse', noconstant
		predict double `y2h' if `touse', residuals
		gen double `y1hy2h' = `y1h'*`y2h' if `touse'
		replace `y1h' = `y1h'^2
		replace `y2h' = `y2h'^2
		summ `y1hy2h' if `touse', meanonly
		scalar `s12' = r(sum)
		summ `y1h' if `touse', meanonly
		scalar `s11' = r(sum)
		summ `y2h' if `touse', meanonly
		scalar `s22' = r(sum)
		mat `ymxy' = J(2,2,0)
		mat `ymxy'[1, 1] = `s11'
		mat `ymxy'[1, 2] = `s12'
		mat `ymxy'[2, 1] = `s12'
		mat `ymxy'[2, 2] = `s22'

		/* And Y'M1Y.	*/
		tempvar y1hh y2hh y1hhy2hh
		reg `ry1' `x1' if `touse', noconstant
		predict double `y1hh' if `touse', residuals
		reg `ry2' `x1' if `touse', noconstant
		predict double `y2hh' if `touse', residuals
		gen double `y1hhy2hh' = `y1hh'*`y2hh' if `touse'
		replace `y1hh' = `y1hh'^2
		replace `y2hh' = `y2hh'^2
		summ `y1hhy2hh' if `touse', meanonly
		scalar `s12' = r(sum)
		summ `y1hh' if `touse', meanonly
		scalar `s11' = r(sum)
		summ `y2hh' if `touse', meanonly
		scalar `s22' = r(sum)
		mat `ym1y' = J(2,2,0)
		mat `ym1y'[1, 1] = `s11'
		mat `ym1y'[1, 2] = `s12'
		mat `ym1y'[2, 1] = `s12'
		mat `ym1y'[2, 2] = `s22'

		/* Now form the matrix to get kappa from. */
		tempname ymxyih kmat vals vecs kappa bliml varliml
		mat_inv_sqrt `ymxy' `ymxyih'
		mat `kmat' = `ymxyih'*`ym1y'*`ymxyih'
		mat symeigen `vecs' `vals' = `kmat'
		loc junk = colsof(`vals')
		scalar `kappa' = `vals'[1, `junk'] 
		
		/* Now assemble the "X'X" matrix.	*/
		tempname xpxi xpy xpy1 xpy2 bliml varliml rss sighatsq
		mat `xpx'[1, 1] = `y1py1' - `kappa'*`y1mxy1'
		mat `xpxi' = syminv(`xpx')
		/* And the "X'y" matrix.	*/
		mat accum `xpy1' = `ry2' `ry1' if `touse', noconstant
		mat `xpy1' = `xpy1'[1, 2]
		mat `xpy1' = `xpy1' - `kappa'*`ymxy'[1, 2]
		mat accum `xpy2' = `x1' `ry1' if `touse', noconstant
		loc k1 : word count `x1'
		loc k2 = `k1' + 1
		mat `xpy2' = `xpy2'[1..`k1', `k2'...]
		mat `xpy' = `xpy1' \ `xpy2'
		mat `bliml' = (`xpxi'*`xpy')'
		
		/* Now compute sigma_squared to get the covariance matrix. */
		tempvar linpred residsq
		mat score `linpred' = `bliml' if `touse'
		gen double `residsq' = (`ry1' - `linpred')^2 if `touse'
		su `residsq' if `touse', meanonly
		sca `rss' = r(sum)
		sca `sighatsq' = r(sum) / r(N)
		mat `varliml' = `xpxi' * `sighatsq'

		return matrix beta `bliml'
		return matrix var `varliml'
		return scalar rss = `rss'
	}	
	
end	

prog def myprint

	local level `e(level)'
	di
	if "`e(model)'" == "2sls" {
		di as text "Instrumental variables (2SLS) regression"
	}
	else {
		di as text "Instrumental variables (LIML) regression"
	}
	di
	di as text "First-stage results" _col(56) "Number of obs =" 	/*
		*/ as result %8.0f `e(N)'
	if "`e(model)'" == "2sls" {
		di as text "{hline 23}" _col(56) "F(" %3.0f `e(df_m)' 	/*
			*/ "," %6.0f `e(df_r)' ") =" as result %8.2f `e(F)'
		di as text "F(" %3.0f `e(df_m_first)' "," %6.0f 	/*
			*/ `e(df_r_first)' ") =" as result %8.2f 	/*
			*/ `e(F_first)' _col(56) as text "Prob > F      =" /*
			*/ as result %8.4f Ftail(`e(df_m)', `e(df_r)', `e(F)')
	}
	else {
		di as text "{hline 23}" _col(56) "Wald chi2(" %2.0f 	/*
			*/ `e(df_m)' ") =" as result %8.2f `e(wald)'
		di as text "F(" %3.0f `e(df_m_first)' "," %6.0f 	/*
			*/ `e(df_r_first)' ") =" as result %8.2f 	/*
			*/ `e(F_first)' _col(56) as text 		/*
			*/ "Prob > w      =" as result %8.4f 	/*
			*/ chi2tail(`e(df_m)', `e(wald)')
	}	
	di as text "Prob > F      =" as result %8.4f 		/*
		*/ Ftail(`e(df_m_first)', `e(df_r_first)', `e(F_first)') /*
		*/ _col(56) as text "R-squared     =" as result %8.4f `e(r2)'
	di as text "R-squared	  =" as result %8.4f `e(r2_first)'	/*
		*/ _col(56) as text "Adj R-squared =" as result %8.4f `e(r2_a)'
	di as text "Adj R-squared =" as result %8.4f `e(r2_a_first)'	/*
		*/ _col(56) as text "Root MSE      =" as result 	/*
		*/ %8.3f `e(rmse)' 
	di
	eret display, level(`level')
	di as text "Instrumented:  " _c
	Disp `e(instd)'
	di as text "Instruments:   " _c
	if "`e(instcons)'" == "yes" {
		Disp `e(insts)'
	}
	else {
		Disp `e(insts)' (No constant included)
	}
	di as text "Confidence set and p-value for " 			/*
		*/ abbrev("`e(instd)'", 13) " "				/*
		*/ "are based on normal approximation"
	di as text "{hline 78}"
	
	
end

program define myprintregions

	args interval AR LM

	di _n _n
	di as text "{hline 78}"
	if ("`interval'"=="") {
		local what "confidence set"
	}
	else {
		local what "confidence interval"
	}
	local plural
	if "`AR'`LM'" != "" {
		local plural "s"
	}
	local title : di "Coverage-corrected `what'`plural' and p-value`plural'"
	local titcol `=int((78 - length("`title'"))/2 + 1)'
	di as text _col(`titcol') "`title'"
	
	local title : di "for Ho: _b[`e(instd)'] = " %-9.0g e(H0_b)
	local titcol `=int((78 - length("`title'"))/2 + 1)'
	di as text _col(`titcol') "`title'"
	
	local title : di "LIML estimate of _b[`e(instd)'] = " /*
		*/	%-9.0g `e(b_instd_liml)'
	local titcol `=int((78 - length("`title'"))/2 + 1)'
	di as text _col(`titcol') "`title'"
	
	di
	di as text "{hline 78}"
	local title : di proper("`what'")
	local titcol `=int((70 - 15 - length("`title'"))/2) + 15'
	di as text _col(2) "Test" _col(`titcol') "`title'" _col(71) "p-value"
	di as text "{hline 78}"	
	
	local lrpval  "as result _col(72) %6.4f e(LR_p)"
	if `e(LR_type)' == 1 {
		local lrint "empty"
	}
	else if `e(LR_type)' == 2 {
		local lrint : di "[" 					/*
			*/	 %9.0g e(LR_x1)				/*
			*/	 ", "					/*
			*/	 %9.0g e(LR_x2)				/*
			*/	 "]"
	}
	else if `e(LR_type)' == 3 {
		local lrint "(-inf, +inf)"
	}
	else {
		if "`interval'" == "" {
			local lrint : di "(-inf, " 			/*
				*/	 %9.0g e(LR_x1) 		/*
				*/       "] U [" 			/*
				*/	 %9.0g e(LR_x2) 		/*
				*/	 ", +inf)" 
		}
		else {
			local lrint "(-inf, +inf)"
		}
	}
	local len `=length(`"`lrint'"')'
	local len `=int((70 - 15 - `len')/2)'
	local len `=`len'+15'
	di _col(2) "Conditional LR" _col(`len') as result "`lrint'" `lrpval'

	if ("`AR'"!="") {
		local arpval "as result _col(72) %6.4f e(AR_p)"
		if `e(AR_type)' == 1 {
			local arint "empty"
		}
		else if `e(AR_type)' == 2 {
			local arint : di "[" 				/*
				*/	 %9.0g e(AR_x1)			/*
				*/	 ", "				/*
				*/	 %9.0g e(AR_x2)			/*
				*/	 "]"
		}
		else if `e(AR_type)' == 3 {
			local arint "(-inf, +inf)"
		}
		else {
			if "`interval'" == "" {
				local arint : di "(-inf, " 		/*
					*/	 %9.0g e(AR_x1) 	/*
					*/       "] U [" 		/*
					*/	 %9.0g e(AR_x2) 	/*
					*/	 ", +inf)" 
			}
			else {
				local arint "(-inf, +inf)"
			}
		}
		local len `=length(`"`arint'"')'
		local len `=int((70 - 15 - `len')/2)'
		local len `=`len'+15'
		di _col(2) as text "Anderson-Rubin" 			/*
			*/	_col(`len') as result "`arint'" `arpval'
	}
	if ("`LM'"!="") {
		local lmpval "as result _col(72) %6.4f e(LM_p)"
		if `e(LM_type)' == 2 {
			local lmint : di "["				/*
				*/	 %9.0g e(LM_x1)			/*
				*/	 ", "				/*
				*/	 %9.0g e(LM_x2)			/*
				*/	"]"
		}
		else if `e(LM_type)' == 3 {
			local lmint "(-inf, +inf)"
		}
		else if `e(LM_type)' == 4 {
			if "`interval'" == "" {
				local lmint : di "(-inf, " 		/*
					*/	 %9.0g e(LM_x1) 	/*
					*/       "] U [" 		/*
					*/	 %9.0g e(LM_x2) 	/*
					*/	 ", +inf)" 
			}
			else {
				local lmint "(-inf, +inf)"
			}
		}
		else if `e(LM_type)' == 5 {
			if "`interval'" == "" {
				local lmint : di "(-inf, "		/*
					*/	 %6.0g e(LM_x1)		/*
					*/	 "] U ["		/*
					*/	 %6.0g e(LM_x2)		/*
					*/	 ", "			/*
					*/	 %6.0g e(LM_x3)		/*
					*/	 "] U ["		/*
					*/	 %6.0g e(LM_x4)		/*
					*/	 ", +inf)"	
			}
			else {
				local lmint "(-inf, +inf)"
			}
		}
		else if `e(LM_type)' == 6 & e(LM_x1) < e(LM_x3) {
			if "`interval'" == "" {
				local lmint : di "["			/*
					*/	 %9.0g e(LM_x1)		/*
					*/	 ", "			/*
					*/	 %9.0g e(LM_x2)		/*
					*/	 "] U ["		/*
					*/	 %9.0g e(LM_x3)		/*
					*/	 ", "			/*
					*/	 %9.0g e(LM_x4)		/*
					*/	 "]"
			}
			else {
				local lmint : di "["			/*
					*/	 %9.0g e(LM_x1)		/*
					*/	 ", "			/*
					*/	 %9.0g e(LM_x4)		/*
					*/	 "]"
			}
		}
		else {
			if "`interval'" == "" {
				local lmint : di "["			/*
					*/	 %9.0g e(LM_x3)		/*
					*/	 ", "			/*
					*/	 %9.0g e(LM_x4)		/*
					*/	 "] U ["		/*
					*/	 %9.0g e(LM_x1)		/*
					*/	 ", "			/*
					*/	 %9.0g e(LM_x2)		/*
					*/	 "]"			
			}
			else {
				local lmint : di "["			/*
					*/	 %9.0g e(LM_x3)		/*
					*/	 ", "			/*
					*/	 %9.0g e(LM_x2)		/*
					*/	 "]"
			}
		}
		local len `=length(`"`lmint'"')'
		local len `=int((70 - 15 - `len')/2)'
		local len `=`len'+15'
		di _col(2) as text "Score (LM)" 			/*
			*/ _col(`len') as result "`lmint'" `lmpval'
	}
	di as text "{hline 78}"

end


/* Lifted straight out of ivreg.ado.	*/
program define IsStop, sclass

	/* sic, must do tests one-at-a-time, 0 may be very large */
	if `"`0'"' == "[" {
		sret local stop 1
		exit
	}
	if `"`0'"' == "," {
		sret local stop 1
		exit
	}
	if `"`0'"' == "if" {
		sret local stop 1
		exit
	}
	if `"`0'"' == "in" {
		sret local stop 1
		exit
	}
	if `"`0'"' == "" {
		 sret local stop 1
		 exit
	}
	else	 sret local stop 0

end

/* More code borrowed from ivreg.ado.	*/
program define Disp
	local first ""
	local piece : piece 1 64 of `"`0'"'
	local i 1
	while "`piece'" != "" {
		di in gr "`first'`piece'"
		local first "              "
		local i = `i' + 1
		local piece : piece `i' 64 of `"`0'"'
	}
	if `i'==1 { 
		di 
	}
end

exit

