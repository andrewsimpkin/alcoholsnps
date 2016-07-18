cd "H:/My Papers/LDA Alcohol"
use "Data/Alc Polgene Long_final.dta", clear
* set up drinkers/non-drinkers
sort uniqid age
bysort uniqid: gen visit =_n
bysort uniqid: egen max_alc = max(alc_week)
bysort uniqid: gen all0 = max_alc == 0

* merge confounders
merge m:m aln mum_kid using "Data/Conf_for_AS/confounders_kids"
	drop _merge
merge m:m aln mum_kid using "Data/Conf_for_AS/confounders_mums"
	drop _merge
	
	* keep only visit 1 for mums/kids
	bysort aln mum: replace visit  = _n
	keep if visit ==1
	
	* all conf binary except meduc (ord) and cpd (count). msocclass has 65 for army, drop:
	drop if m_soc==65
	
	* associations with confounders
	putexcel set "Results/Tables.xlsx", modify sheet("Confounders")
	
	* non-binary first
	local i 2
	quietly{
	foreach snp of varlist rs* {
	ologit m_educ `snp'
		lincom `snp', eform
		putexcel B`i' = ("m_educ") C`i' = ("`snp'") D`i' = (r(estimate)) ///
										 E`i' = (r(se)) F`i' =(2*(1-normal(abs(_b[`snp']/_se[`snp'])))) ///
										 G`i' =(r(estimate)-1.96*r(se)) H`i'=(r(estimate)+1.96*r(se))
						local i = `i' + 1
				}
		}
		
	local i 92
	quietly{
			foreach snp of varlist rs* {
					foreach cigs in cpdF cpdG cpdH cpdK {
					regress `cigs' `snp'
						lincom `snp'
						putexcel B`i' = ("`cigs'") C`i' = ("`snp'") D`i' = (r(estimate)) ///
														 E`i' = (r(se)) F`i' =(2*(1-normal(abs(r(estimate)/r(se))))) ///
														 G`i' =(r(estimate)-1.96*r(se)) H`i'=(r(estimate)+1.96*r(se))
					
					local i = `i' + 1
					
					}
			}
	}
	
local i 452  
	quietly {	
			foreach snp of varlist rs* {
				foreach conf in ethnicity_kid ethnicity_mum anxiety7 anxiety10 ///
								anxiety13 anxiety15s depression7 depression10 ///
								depression13 depression15s conductD7 conductD10 ///
								conductD13 conductD15 adhd7 adhd10 adhd13 adhd15 ///
								b_eating18 b_eating13 b_eating14 b_eating16 asb11 ///
								asb13 asb14 asb15 asb18 asb19 asb21 cannabisF cannabisG ///
								cannabisH cannabisK antidF antidG antidH antidK ///
								amphF amphG amphH amphK op_cocF op_cocG op_cocH op_cocK {
								
								cap logistic `conf' `snp'
								cap lincom `snp'
								cap putexcel B`i' = ("`conf'") C`i' = ("`snp'") D`i' = (r(estimate)) ///
										 E`i' = (r(se)) F`i' =(2*(1-normal(abs(_b[`snp']/_se[`snp'])))) ///
										 G`i' =(r(estimate)-1.96*r(se)) H`i'=(r(estimate)+1.96*r(se))	
								
								local i = `i' + 1

				}
			}
		
		}

** PGRS	
cd "H:\My Papers\LDA Alcohol\"
use "Results\units_long_scores.dta", clear

* merge confounders
merge m:m aln mum_kid using "Data/Conf_for_AS/confounders_kids"
	drop _merge
merge m:m aln mum_kid using "Data/Conf_for_AS/confounders_mums"
	drop _merge
	
	* keep only visit 1 for mums/kids
	bysort aln mum_kid: replace visit  = _n
	keep if visit ==1
	
	* all conf binary except meduc (ord) and cpd (count). msocclass has 65 for army, drop:
	drop if m_soc==65
	
	* associations with confounders
	putexcel set "Results/PGRS_confounders.xlsx", modify
	
	* non-binary first
	local i 2
	quietly{
	foreach snp of varlist Sm* {
	cap ologit m_educ `snp'
		cap lincom `snp', eform
		cap putexcel B`i' = ("m_educ") C`i' = ("`snp'") D`i' = (r(estimate)) ///
										 E`i' = (r(se)) F`i' =(2*(1-normal(abs(_b[`snp']/_se[`snp'])))) 
						local i = `i' + 1
				}
		}	
		
		
	local i 9
	quietly{
			foreach snp of varlist Sm* {
					foreach cigs in cpdF cpdG cpdH cpdK {
					cap regress `cigs' `snp'
						cap lincom `snp'
						cap putexcel B`i' = ("`cigs'") C`i' = ("`snp'") D`i' = (r(estimate)) ///
														 E`i' = (r(se)) F`i' =(2*(1-normal(abs(r(estimate)/r(se)))))
					
					local i = `i' + 1
					
					}
			}
	}
	
local i 37  
	quietly {	
			foreach snp of varlist Sm* {
				foreach conf in cannabisF cannabisG ///
								cannabisH cannabisK antidF antidG antidH antidK ///
								amphF amphG amphH amphK op_cocG op_cocH op_cocK {
								
								cap logistic `conf' `snp'
								cap lincom `snp'
								cap putexcel B`i' = ("`conf'") C`i' = ("`snp'") D`i' = (r(estimate)) ///
										 E`i' = (r(se)) F`i' =(2*(1-normal(abs(_b[`snp']/_se[`snp']))))
								
								local i = `i' + 1

				}
			}
		
		}
			
local i 142  
	quietly {	
			foreach snp of varlist Sk* {
				foreach conf in ethnicity_kid ethnicity_mum anxiety7 anxiety10 ///
								anxiety13 anxiety15s depression7 depression10 ///
								depression13 depression15s conductD7 conductD10 ///
								conductD13 conductD15 adhd7 adhd10 adhd13 adhd15 ///
								b_eating18 b_eating13 b_eating14 b_eating16 asb11 ///
								asb13 asb14 asb15 asb18 asb19 asb21  {
								
								cap logistic `conf' `snp'
								cap lincom `snp'
								cap putexcel B`i' = ("`conf'") C`i' = ("`snp'") D`i' = (r(estimate)) ///
										 E`i' = (r(se)) F`i' =(2*(1-normal(abs(_b[`snp']/_se[`snp'])))) ///
										 G`i' =(r(estimate)-1.96*r(se)) H`i'=(r(estimate)+1.96*r(se))	
								
								local i = `i' + 1

				}
			}
		
		}
