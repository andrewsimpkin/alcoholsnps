cd "H:/My Papers/LDA Alcohol"
use "Data/Alc Polgene Long_final.dta", clear
* set up drinkers/non-drinkers
sort uniqid age
bysort uniqid: gen visit =_n
bysort uniqid: egen max_alc = max(alc_week)
bysort uniqid: gen all0 = max_alc == 0

* log transform units
gen logalc = log(alc+1) 

* split into 80-20 samples 5 times for ums and kids separately:
* generate a random number for each mum and child:
set seed 12345
gen mum_unif = runiform() if mum_kid == 0
	bysort uniqid: egen mum_rand = max(mum_unif) if mum_kid == 0
* sort mums and give each a number from 1-5, named "mum_drop". we will then model with drop==1, ==2 in turn etc.
sort mum_rand
egen mum_rgroup = group(mum_rand) if mum_kid == 0
 centile mum_rgroup if visit == 1, centile(20 40 60 80)
	gen mum_drop = 1 if mum_kid == 0
		replace mum_drop = 2 if mum_rgroup > r(c_1)
		replace mum_drop = 3 if mum_rgroup > r(c_2)
		replace mum_drop = 4 if mum_rgroup > r(c_3)
		replace mum_drop = 5 if mum_rgroup > r(c_4)
			replace mum_drop = . if mum_kid != 0 // make sure no kids included
 

* repeat for kids
set seed 12345
gen kid_unif = runiform() if mum_kid == 1
	bysort uniqid: egen kid_rand = max(kid_unif) if mum_kid == 1
* sort mums and give each a number from 1-5, named "mum_drop". we will then model with drop==1, ==2 in turn etc.
sort kid_rand
egen kid_rgroup = group(kid_rand) if mum_kid == 1
 centile kid_rgroup if visit == 1, centile(20 40 60 80)
	gen kid_drop = 1 if mum_kid == 1
		replace kid_drop = 2 if kid_rgroup > r(c_1)
		replace kid_drop = 3 if kid_rgroup > r(c_2)
		replace kid_drop = 4 if kid_rgroup > r(c_3)
		replace kid_drop = 5 if kid_rgroup > r(c_4)
			replace kid_drop = . if mum_kid != 1 // make sure no kids included
	
* save, will need these groups to make scores later...
save "H:\My Papers\LDA Alcohol\Results\dropgroup.dta", replace 

* run for i=1, 2, .. 5 (can do in separate Stata GUIs to save time)
use "H:\My Papers\LDA Alcohol\Results\dropgroup.dta", clear
cap postclose units
postfile units group str10 SNP b_m se_m t_m p_m r_m r_mSNP ///
						 b_k se_k t_k p_k r_k r_kSNP ///
						 using "H:/My Papers/LDA Alcohol/Results/logunits.dta", replace
set more off
	quietly{
			forvalues i = 1(1)5 {
				foreach var of varlist  rs* {
						* model log units for mums
						cap xtmixed logalc `var' age PC* if mum_kid == 0 & all0 == 0 & mum_drop != `i' || aln:
									local b_m 	= _b[`var']
									local se_m 	= _se[`var']
									local t_m 	= _b[`var']/_se[`var']
									local p_m 	= (2*(1 - normal(abs(_b[`var']/_se[`var']))))
								* get R2 for SNP + age, PC
								mltrsq
								local r_m  = e(sb_rsq_l1)
										* get R2 for just SNP
										cap xtmixed logalc `var' if mum_kid == 0 & all0 == 0 & mum_drop != `i' || aln:
											mltrsq
												local r_mSNP  = e(sb_rsq_l1)
						* model log units for kids 
						cap xtmixed logalc `var' age PC* if mum_kid == 1 & all0 == 0  & kid_drop != `i' || aln:
									local b_k 	= _b[`var']
									local se_k 	= _se[`var']
									local t_k	= _b[`var']/_se[`var']
									local p_k 	= (2*(1 - normal(abs(_b[`var']/_se[`var']))))
								* get R2 for SNP + age, PC
								mltrsq
								local r_k  = e(sb_rsq_l1)
										* get R2 for just SNP
										cap xtmixed logalc `var' if mum_kid == 1 & all0 == 0 & kid_drop != `i' || aln:
											mltrsq
												local r_kSNP  = e(sb_rsq_l1)
									
										post  units (`i') ("`var'") (`b_m') (`se_m') ///
																	(`t_m') (`p_m') (`r_m') (`r_mSNP') ///
																	(`b_k') (`se_k') ///
																	(`t_k') (`p_k') (`r_k') (`r_kSNP') 							
																			
								noisily di "*`count'" _c
							
							local count = `count' +1

						}  	// end loop over cg's
					} // end loop over drop group
				}  	// end loop over quietly

		postclose units

* set up coefficient by threshold datasets to enable score making	
cd "H:/My Papers/LDA Alcohol/Results/"
* loop over mums/kids
set more off
foreach i in m k {
		* loop over p-value thresholds
		foreach thresh in 0.5 0.4 0.3 0.2 0.1 0.05 0.01 0.001 {
					* load data
					use logunits, clear
						* drop SNPs which are not under threshold
						* keep if p_`i' < `thresh'
							* make local macro = to name of SNP
							levelsof SNP, local(levels) 
								* loop over SNPs to create variable containing beta coefficients (which will be used to make score)
								foreach snp of local levels {
									* loop over drop groups, i.e. label for which 80% were used as training data to estimate the beta coefficient
									forvalues gp = 1(1)5 {
														* generate a variable labelled
														gen b_`snp'_`gp' = b_`i' if SNP == "`snp'" & group == `gp' & p_`i' < `thresh'
															egen max = max(b_`snp'_`gp')
																replace b_`snp'_`gp' = max
																	drop max
															} // end loop over groups (same dataset of coefficients)
													} // end loop over SNPs (same dataset of coefficients)
													
											* save for each threshold and separately for mums and kids
											* don't need model info, only beta coefficients
											drop group SNP b_m se_m t_m p_m r_m r_mSNP b_k se_k t_k p_k r_k r_kSNP
											* only need one row of data
											gen n=_n
												keep if n==1
													drop n
											
											* save
											save beta_`i'_`thresh'.dta, replace
				
											} // end loop over thresholds (one dataset per threshold)
									} // end loop over mums/kids
		
		
		
						
* merge sample data with coefficients for each threshold, derive score, post R-sq to results file for mums 
cd "H:/My Papers/LDA Alcohol/Results/"
cap postclose units
postfile units p gp r r_a ///
						 using "snp_resM.dta", replace
set more off
	quietly{
		foreach thresh in 0.5 0.4 0.3 0.2 0.1 0.05 0.01 {
			forvalues gp = 1(1)5 {
				use dropgroup, clear
					append using beta_m_`thresh'.dta
							* repeat coefficient value for each individual
							foreach var of varlist b_rs* {
								egen max = max(`var')
									replace `var' = max
										drop max
								} // end loop over coefficients
						
						* make scores for each drop group for everyone
						* multiply coef by SNP if coef != . i.e. SNP under threshold
						foreach var of varlist rs* {
							if b_`var'_`gp'!=. {
												gen S_`var'_`gp' = b_`var'_`gp'*`var'
												} // end if loop
										} // end loop over snps
								* sum products to get score
								egen S_`gp' = rsum(S_rs*_`gp')
				* model repeated measures of log units on scores in discovery data, i.e. those dropped in original analysis
				xtmixed logalc S_`gp' age PC* if mum_kid == 0 & all0 == 0 & mum_drop == `gp' || aln:
											* get R2 for allelic score + age, PC
											mltrsq
												local r_a  = e(sb_rsq_l1)
				xtmixed logalc S_`gp' if mum_kid == 0 & all0 == 0 & mum_drop == `gp' || aln:
											* get R2 for allelic score only
											mltrsq
												local r  = e(sb_rsq_l1)	
							* post out R-squared, threshold and dropped group
							post  units (`thresh') (`gp') (`r') (`r_a')
						
					} // end loop over group
				noisily di "*`thresh'" _c
				} // end loop over thresholds
		}	// end quietly loop	
		postclose units
				
* merge sample data with coefficients for each threshold, post R-sq to results file for kids
cd "H:/My Papers/LDA Alcohol/Results/"
cap postclose units
postfile units p gp r r_a ///
						 using "snp_resK.dta", replace
set more off
	quietly{
		foreach thresh in 0.5 0.4 0.3 0.2 0.1 0.05 {
			forvalues gp = 1(1)5 {
				use dropgroup, clear
					append using beta_k_`thresh'.dta
							* repeat coefficient value for each individual
							foreach var of varlist b_rs* {
								egen max = max(`var')
									replace `var' = max
										drop max
								} // end loop over coefficients
						
						* make scores for each drop group for everyone
						* multiply coef by SNP if coef != . i.e. SNP under threshold
						foreach var of varlist rs* {
							if b_`var'_`gp'!=. {
												gen S_`var'_`gp' = b_`var'_`gp'*`var'
												} // end if loop
										} // end loop over snps
								* sum products to get score
								egen S_`gp' = rsum(S_rs*_`gp')
				
				* model repeated measures of log units on scores in discovery data, i.e. those dropped in original analysis
				cap xtmixed logalc S_`gp' age PC* if mum_kid == 1 & all0 == 0 & kid_drop == `gp' || aln:
											* get R2 for allelic score + age, PC
											mltrsq
											local r_a  = e(sb_rsq_l1)
											gen r_a  = e(sb_rsq_l1)
				cap xtmixed logalc S_`gp' if mum_kid == 1 & all0 == 0 & kid_drop == `gp' || aln:
											* get R2 for allelic score only
											mltrsq
											local r  = e(sb_rsq_l1)	
											gen r  = e(sb_rsq_l1)
							post  units (`thresh') (`gp') (`r') (`r_a')
						
					} // end loop over group
				noisily di "*`thresh'" _c
				} // end loop over thresholds
		}	// end quietly loop	
		postclose units
				
		
		
		
		
		
		
		
		
		
		
		
		
********** analysis of everyone, i.e. no dropping.		
		
* run model 
use "H:\My Papers\LDA Alcohol\Results\dropgroup.dta", clear
cap postclose units
postfile units str10 SNP b_m se_m t_m p_m r_m r_mSNP ///
						 b_k se_k t_k p_k r_k r_kSNP ///
						 using "H:/My Papers/LDA Alcohol/Results/logunits_all.dta", replace
set more off
	quietly{
			local count 1
				foreach var of varlist  rs* {
						* model log units for mums
						cap xtmixed logalc `var' age PC* if mum_kid == 0 & all0 == 0 || aln:
									local b_m 	= _b[`var']
									local se_m 	= _se[`var']
									local t_m 	= _b[`var']/_se[`var']
									local p_m 	= (2*(1 - normal(abs(_b[`var']/_se[`var']))))
								* get R2 for SNP + age, PC
								mltrsq
								local r_m  = e(sb_rsq_l1)
										* get R2 for just SNP
										cap xtmixed logalc `var' if mum_kid == 0 & all0 == 0  || aln:
											mltrsq
												local r_mSNP  = e(sb_rsq_l1)
						* model log units for kids 
						cap xtmixed logalc `var' age PC* if mum_kid == 1 & all0 == 0 || aln:
									local b_k 	= _b[`var']
									local se_k 	= _se[`var']
									local t_k	= _b[`var']/_se[`var']
									local p_k 	= (2*(1 - normal(abs(_b[`var']/_se[`var']))))
								* get R2 for SNP + age, PC
								mltrsq
								local r_k  = e(sb_rsq_l1)
										* get R2 for just SNP
										cap xtmixed logalc `var' if mum_kid == 1 & all0 == 0 || aln:
											mltrsq
												local r_kSNP  = e(sb_rsq_l1)
									
										post  units ("`var'") (`b_m') (`se_m') ///
																	(`t_m') (`p_m') (`r_m') (`r_mSNP') ///
																	(`b_k') (`se_k') ///
																	(`t_k') (`p_k') (`r_k') (`r_kSNP') 							
																			
								noisily di "*`count'" _c
							
							local count = `count' +1

						}  	// end loop over cg's
					
				}  	// end loop over quietly

		postclose units

* set up beta datasets to enable score making later		
cd "H:/My Papers/LDA Alcohol/Results/"
* loop over mums/kids
set more off
foreach i in m k {
		* loop over p-value thresholds
		foreach thresh in 0.5 0.4 0.3 0.2 0.1 0.05 0.01 0.001 {
					* load data
					use logunits, clear
						* drop SNPs which are not under threshold
						* keep if p_`i' < `thresh'
							* make local macro = to name of SNP
							levelsof SNP, local(levels) 
								* loop over SNPs to create variable containing beta coefficients (which will be used to make score)
								foreach snp of local levels {
														* generate a variable labelled
														gen b_`snp' = b_`i' if SNP == "`snp'" & p_`i' < `thresh'
															egen max = max(b_`snp')
																replace b_`snp' = max
																	drop max
													} // end loop over SNPs (same dataset of coefficients)
													
											* save for each threshold and separately for mums and kids
											* don't need model info, only beta coefficients
											drop SNP b_m se_m t_m p_m r_m r_mSNP b_k se_k t_k p_k r_k r_kSNP
											* only need one row of data
											gen n=_n
												keep if n==1
													drop n
											
											* save
											save ALLbeta_`i'_`thresh'.dta, replace
				
											} // end loop over thresholds (one dataset per threshold)
									} // end loop over mums/kids
		
		
		
						
* merge sample data with coefficients for each threshold, get Score for mums, kids, each threshold: 
cd "H:/My Papers/LDA Alcohol/Results/"
use dropgroup, clear
	* setup data to contain scores:
	save units_long_scores, replace

	quietly{
	foreach i in m k {
		foreach thresh in 0.5 0.4 0.3 0.2 0.1 0.05 0.01 {
		set more off
				use units_long_scores, clear
					append using ALLbeta_`i'_`thresh'.dta
							* repeat coefficient value for each individual
							foreach var of varlist b_rs* {
								egen max = max(`var')
									replace `var' = max
										drop max
								} // end loop over coefficients
						
						* multiply coef by SNP if coef != . i.e. SNP under threshold
						foreach var of varlist rs* {
							if b_`var'!=. {
												gen S_`var' = b_`var'*`var'
												} // end if loop
										} // end loop over snps
								* sum products to get score
								* can't have var name with decimal point...multiply thresh by 100 to get usable name
								local name = `thresh'*100
								egen S`i'_`name' = rsum(S_rs*)
								* set to missing if mum data for kids score and vice versa:
								cap replace Sm_`name' = . if mum_kid == 1
								cap replace Sk_`name' = . if mum_kid == 0
									drop S_rs* b_rs*
									save units_long_scores, replace
									}
								}
							}

							
							
* model repeated measures of log units on scores in whole dataset, calculate R-squared
cd "H:/My Papers/LDA Alcohol/Results/"
* load scores
use units_long_scores, clear
* optimal threshold scores are S50 for mums:
set more off
xtmixed logalc Sm_50 if mum_kid == 0 & all0 == 0 || aln:
							* get R2 for allelic score only
							mltrsq
								gen r_m  = e(sb_rsq_l1)	
* S5 for kids
set more off
xtmixed logalc Sk_5 if mum_kid == 1 & all0 == 0 || aln:
							* get R2 for allelic score only
							mltrsq
												gen r_k  = e(sb_rsq_l1)	

