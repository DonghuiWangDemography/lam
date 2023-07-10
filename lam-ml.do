*ml program to estimate attitude toward China from 1974-2019 using pooled survey series 
*prepared by Donghui Wang on 05/29/2023
net install cleanplots, from("https://tdmize.github.io/data/cleanplots")
set scheme cleanplots, perm




cd "\\NAS\home\Drive\US_project\cleaned_data"


****************************
*STEP1 : prepare initial values to feed the program 
*Note :stata ml evalautes individual likelihood function. The data lam-china is prepared as individual data. 
****************************
use lam-china, clear //cleaned data, marginal response of favorability expand to individual observation

  qui: oprobit k i.s
  mat rs=r(table)
  mat s=rs[1, " k:"]


  levelsof q, local(survey)
  foreach x of local survey { 
  qui: oprobit k i.s if q==`x' 
  qui: mat c=r(table)
  qui: mat tau`x'=c[1, " /:"]
  	
  }
 
  	//rename tau
	forval i= 1/13{
	
    mat colname tau`i'=_cons
	local n= `= colsof(tau`i')'
	local ctau`i' ""
	forval j=1/`n'{	
	local ctau`i' "`ctau`i'' tau`i'_`j'"  
	
	}
	mat coleq tau`i'= `ctau`i''
	}
	
  *put together intial values 
  mat e0=s,tau1, tau2, tau3, tau4, tau5,tau6, tau7,tau8,tau9,tau10,tau11, tau12, tau13
  
******************************
* STEP 2 : ml program 
****************************
 capture program drop lam  
 
  program define lam  
   
    #delimit ; 
	args  lnf mu tau1_1   
	             tau2_1  tau2_2 tau2_3 tau2_4 tau2_5 tau2_6 tau2_7 tau2_8 tau2_9  
				 tau3_1  tau3_2 tau3_3  
				 tau4_1  tau4_2 tau4_3 
				 tau5_1  tau5_2 tau5_3 tau5_4
				 tau6_1  tau6_2
				 tau7_1  tau7_2 tau7_3 tau7_4 tau7_5 tau7_6 tau7_7 tau7_8 tau7_9
				 tau8_1  tau8_2 tau8_3
				 tau9_1  tau9_2 tau9_3 tau9_4	
                 tau10_1  tau10_2 tau10_3  
				 tau11_1  tau11_2 tau11_3 
                 tau12_1  tau12_2 tau12_3 
				 tau13_1  tau13_2 tau13_3 			 

	;
	#delimit cr
	

	quietly {
  
	  // 1. ABC : 2 level 
    replace `lnf' =ln(normal(`tau1_1' -`mu'))                             if $ML_y1 ==1 & q==1
    replace `lnf' =ln(1 - normal(`tau1_1'  -`mu' ))                       if $ML_y1 ==2 & q==1
	
	
 	// 2. GSS
	replace `lnf' =ln(normal( `tau2_1' -`mu'))                            if $ML_y1 ==1 & q==2
	replace `lnf' =ln(normal(`tau2_2' -`mu')  - normal( `tau2_1' -`mu'))  if $ML_y1 ==2 & q==2
    replace `lnf' =ln(normal(`tau2_3' -`mu')  - normal( `tau2_2' -`mu'))  if $ML_y1 ==3 & q==2
	replace `lnf' =ln(normal(`tau2_4' -`mu')  - normal( `tau2_3' -`mu'))  if $ML_y1 ==4 & q==2
	replace `lnf' =ln(normal(`tau2_5' -`mu')  - normal( `tau2_4' -`mu'))  if $ML_y1 ==5 & q==2
    replace `lnf' =ln(normal(`tau2_6' -`mu')  - normal( `tau2_5' -`mu'))  if $ML_y1 ==6 & q==2
	replace `lnf' =ln(normal(`tau2_7' -`mu')  - normal( `tau2_6' -`mu'))  if $ML_y1 ==7 & q==2
	replace `lnf' =ln(normal(`tau2_8' -`mu')  - normal( `tau2_7' -`mu'))  if $ML_y1 ==8 & q==2
    replace `lnf' =ln(normal(`tau2_9' -`mu')  - normal( `tau2_8' -`mu'))  if $ML_y1 ==9 & q==2
    replace `lnf' =ln(1 - normal(`tau2_9'  -`mu'))                        if $ML_y1 ==10 & q==2
			
	
	// 3. PEW 
	replace `lnf' =ln(normal(`tau3_1' -`mu'))                             if $ML_y1 ==1 & q==3
	replace `lnf' =ln(normal(`tau3_2' -`mu')  - normal( `tau3_1' -`mu'))  if $ML_y1 ==2 & q==3
    replace `lnf' =ln(normal(`tau3_3' -`mu')  - normal( `tau3_2' -`mu'))  if $ML_y1 ==3 & q==3
    replace `lnf' =ln(1 - normal(`tau3_3'  -`mu'))                        if $ML_y1 ==4 & q==3
	
	
 	//4. TRA_4 
	
	replace `lnf' =ln(normal(`tau4_1' -`mu'))                             if $ML_y1 ==1 & q==4
	replace `lnf' =ln(normal(`tau4_2' -`mu')  - normal( `tau4_1' -`mu'))  if $ML_y1 ==2 & q==4
    replace `lnf' =ln(normal(`tau4_3' -`mu')  - normal( `tau4_2' -`mu'))  if $ML_y1 ==3 & q==4
    replace `lnf' =ln(1 - normal(`tau4_3'  -`mu'))                        if $ML_y1 ==4 & q==4

	
	//5.  TRA 5 
	replace `lnf' =ln(normal(`tau5_1' -`mu'))                             if $ML_y1 ==1 & q==5
	replace `lnf' =ln(normal(`tau5_2' -`mu')  - normal( `tau5_1' -`mu'))  if $ML_y1 ==2 & q==5
    replace `lnf' =ln(normal(`tau5_3' -`mu')  - normal( `tau5_2' -`mu'))  if $ML_y1 ==3 & q==5
	replace `lnf' =ln(normal(`tau5_4' -`mu')  - normal( `tau5_3' -`mu'))  if $ML_y1 ==4 & q==5
    replace `lnf' =ln(1 - normal(`tau5_4'  -`mu')) 					      if $ML_y1 ==5 & q==5	
	
	
	//6. USCBS_3 		
	replace `lnf' =ln(normal(`tau6_1' -`mu'))                             if $ML_y1 ==1 & q==6
	replace `lnf' =ln(normal(`tau6_2' -`mu')  - normal( `tau6_1' -`mu'))  if $ML_y1 ==2 & q==6
    replace `lnf' =ln(1 - normal(`tau6_2'  -`mu'))                        if $ML_y1 ==3 & q==6
		
	
	// 7. USGALLUP_10
	replace `lnf' =ln(normal( `tau7_1' -`mu'))                            if $ML_y1 ==1 & q==7
	replace `lnf' =ln(normal(`tau7_2' -`mu')  - normal( `tau7_1' -`mu'))  if $ML_y1 ==2 & q==7
    replace `lnf' =ln(normal(`tau7_3' -`mu')  - normal( `tau7_2' -`mu'))  if $ML_y1 ==3 & q==7
	replace `lnf' =ln(normal(`tau7_4' -`mu')  - normal( `tau7_3' -`mu'))  if $ML_y1 ==4 & q==7
	replace `lnf' =ln(normal(`tau7_5' -`mu')  - normal( `tau7_4' -`mu'))  if $ML_y1 ==5 & q==7
    replace `lnf' =ln(normal(`tau7_6' -`mu')  - normal( `tau7_5' -`mu'))  if $ML_y1 ==6 & q==7
	replace `lnf' =ln(normal(`tau7_7' -`mu')  - normal( `tau7_6' -`mu'))  if $ML_y1 ==7 & q==7
	replace `lnf' =ln(normal(`tau7_8' -`mu')  - normal( `tau7_7' -`mu'))  if $ML_y1 ==8 & q==7
    replace `lnf' =ln(normal(`tau7_9' -`mu')  - normal( `tau7_8' -`mu'))  if $ML_y1 ==9 & q==7
    replace `lnf' =ln(1 - normal(`tau7_9'  -`mu'))                        if $ML_y1 ==10 & q==7
				
	
	// 8. gallup_4
	replace `lnf' =ln(normal(`tau8_1' -`mu'))                             if $ML_y1 ==1 & q==8
	replace `lnf' =ln(normal(`tau8_2' -`mu')  - normal( `tau8_1' -`mu'))  if $ML_y1 ==2 & q==8
    replace `lnf' =ln(normal(`tau8_3' -`mu')  - normal( `tau8_2' -`mu'))  if $ML_y1 ==3 & q==8
    replace `lnf' =ln(1 - normal(`tau8_3'  -`mu')) 	                      if $ML_y1 ==4 & q==8	
	
	
	// 9. USKN_5
	replace `lnf' =ln(normal(`tau9_1' -`mu'))                             if $ML_y1 ==1 & q==9
	replace `lnf' =ln(normal(`tau9_2' -`mu')  - normal( `tau9_1' -`mu'))  if $ML_y1 ==2 & q==9
    replace `lnf' =ln(normal(`tau9_3' -`mu')  - normal( `tau9_2' -`mu'))  if $ML_y1 ==3 & q==9
	replace `lnf' =ln(normal(`tau9_4' -`mu')  - normal( `tau9_3' -`mu'))  if $ML_y1 ==4 & q==9
    replace `lnf' =ln(1 - normal(`tau9_4'  -`mu')) 						  if $ML_y1 ==5 & q==9	
	
	// 10. USORC_4
	replace `lnf' =ln(normal(`tau10_1' -`mu'))                              if $ML_y1 ==1 & q==10
	replace `lnf' =ln(normal(`tau10_2' -`mu')  - normal( `tau10_1' -`mu'))  if $ML_y1 ==2 & q==10
    replace `lnf' =ln(normal(`tau10_3' -`mu')  - normal( `tau10_2' -`mu'))  if $ML_y1 ==3 & q==10
    replace `lnf' =ln(1 - normal(`tau10_3'  -`mu')) 	                    if $ML_y1 ==4 & q==10	
	
	//11. USPSRA_4 
	replace `lnf' =ln(normal(`tau11_1' -`mu'))                              if $ML_y1 ==1 & q==11
	replace `lnf' =ln(normal(`tau11_2' -`mu')  - normal( `tau11_1' -`mu'))  if $ML_y1 ==2 & q==11
    replace `lnf' =ln(normal(`tau11_3' -`mu')  - normal( `tau11_2' -`mu'))  if $ML_y1 ==3 & q==11
    replace `lnf' =ln(1 - normal(`tau11_3'  -`mu')) 	                    if $ML_y1 ==4 & q==11	
	
	
	//12 zoby4
	replace `lnf' =ln(normal(`tau12_1' -`mu'))                              if $ML_y1 ==1 & q==12
	replace `lnf' =ln(normal(`tau12_2' -`mu')  - normal( `tau12_1' -`mu'))  if $ML_y1 ==2 & q==12
    replace `lnf' =ln(normal(`tau12_3' -`mu')  - normal( `tau12_2' -`mu'))  if $ML_y1 ==3 & q==12
    replace `lnf' =ln(1 - normal(`tau12_3'  -`mu')) 	                    if $ML_y1 ==4 & q==12	
	
	// abcwp4
	replace `lnf' =ln(normal(`tau13_1' -`mu'))                              if $ML_y1 ==1 & q==13
	replace `lnf' =ln(normal(`tau13_2' -`mu')  - normal( `tau13_1' -`mu'))  if $ML_y1 ==2 & q==13
    replace `lnf' =ln(normal(`tau13_3' -`mu')  - normal( `tau13_2' -`mu'))  if $ML_y1 ==3 & q==13
    replace `lnf' =ln(1 - normal(`tau13_3'  -`mu')) 	                    if $ML_y1 ==4 & q==13	
		
  }
 end 

*initial model : constant cutpoint 
    #delimit ; 
    ml model lf lam (k: k= i.s, noconstant) 
	                 (tau1_1:)
					 (tau2_1:) (tau2_2:) (tau2_3:) (tau2_4:) (tau2_5:) (tau2_6:) (tau2_7:) (tau2_8:) (tau2_9:) 
					 (tau3_1:) (tau3_2:) (tau3_3:) 
					 (tau4_1:) (tau4_2:) (tau4_3:)
					 (tau5_1:) (tau5_2:) (tau5_3:) (tau5_4:)
					 (tau6_1:) (tau6_2:)
					 (tau7_1:) (tau7_2:) (tau7_3:) (tau7_4:) (tau7_5:) (tau7_6:) (tau7_7:) (tau7_8:) (tau7_9:) 
					 (tau8_1:) (tau8_2:) (tau8_3:)
					 (tau9_1:) (tau9_2:) (tau9_3:) (tau9_4:)
					 (tau10_1:) (tau10_2:) (tau10_3:)
					 (tau11_1:) (tau11_2:) (tau11_3:)
					 (tau12_1:) (tau12_2:) (tau12_3:) 
					 (tau13_1:) (tau13_2:) (tau13_3:) 

	;
    #delimit cr
		
ml init e0  // intial values obtained from ordered probit rountine 
ml maximize , difficult  
   estimates store m0 

   mat list r(table)
   mat c=r(table)
   mat e_mu = c[1, "k:"]'
   mat e_ci=c[5..6,"k:"]'


 
 
*graph results  
   mat e_mu = c[1, "k:"]'
   mat e_ci=c[5..6,"k:"]'

   mat e_tau1 = c[1, "tau1_1:"]
   mat e_tau2 = c[1, "tau2_1:" .. "tau2_9:"]'
   mat e_tau3 = c[1, "tau3_1:" .. "tau3_3:"]'
   mat e_tau4 = c[1, "tau4_1:" .. "tau4_3:"]'
   mat e_tau5 = c[1, "tau5_1:" .. "tau5_4:"]'
   mat e_tau5 = c[1, "tau5_1:" .. "tau5_4:"]'
   mat e_tau6 = c[1, "tau6_1:" .. "tau6_2:"]'
   mat e_tau7 = c[1, "tau7_1:" .. "tau7_9:"]'
   mat e_tau8 = c[1, "tau8_1:" .. "tau8_3:"]'
   mat e_tau9 = c[1, "tau9_1:" .. "tau9_4:"]'
   mat e_tau10 = c[1, "tau10_1:".. "tau10_3:"]'
   mat e_tau11 = c[1, "tau11_1:" .. "tau11_3:"]'
   mat e_tau12 = c[1, "tau12_1:" ..  "tau12_3:"]'
   mat e_tau13 = c[1, "tau13_1:" ..  "tau13_3:"]'

   forval i=1/13 {
   svmat e_tau`i', names(matcol)
   }     
   
   levelsof s, matrow(s)
   mat mu= s, e_mu
   
   svmat  mu
   rename mu1 year   
   rename mu2 Esimates 
   
   svmat e_ci
   rename e_ci1 ul
   rename e_ci2 ll        
   twoway (line Esimates year) ///
          (rscatter  ul ll year, recast(rarea) color(grey%20)) ,  /// 
           xlab(1974(5)2020) ylab(-1(0.2)1) legend(ring(0)  order (1 "Estimates" 2 "Confidence Interval"))
   graph save "Graph" "C:\Users\Donghui\SynologyDrive\US_project\image\fig6.gph"
   
   graph use  "C:\Users\Donghui\SynologyDrive\US_project\image\fig6.gph"

 graph export "fig6.jpg",replace  width(8000) height(6000)  
  
   
    
***************************************************************
* STEP 3 : Alternative model specifications (footnote 2)
***************************************************************

*m2 : different taus before and after 1989 

gen a = s>1989
*only to the surveys that conducted both before and after 1989 
*2.GSS 6.CBS  7.GALLUP10 8.Gallup4  
     #delimit ; 
    ml model lf ml (k: k= i.s, noconstant) 
	                 (tau1_1:  )
					 (tau2_1: a) (tau2_2: a) (tau2_3: a) (tau2_4: a) (tau2_5: a) (tau2_6: a) (tau2_7: a) (tau2_8: a) (tau2_9: a) 
					 (tau3_1: ) (tau3_2: ) (tau3_3: ) 
					 (tau4_1: ) (tau4_2: ) (tau4_3: )
					 (tau5_1: ) (tau5_2: ) (tau5_3: ) (tau5_4: )
					 (tau6_1: a) (tau6_2: a)
					 (tau7_1: a) (tau7_2: a) (tau7_3: a) (tau7_4: a) (tau7_5: a) (tau7_6: a) (tau7_7: a) (tau7_8: a) (tau7_9: a) 
					 (tau8_1: a) (tau8_2: a ) (tau8_3: a )
					 (tau9_1: ) (tau9_2: ) (tau9_3: ) (tau9_4: )
					 (tau10_1: ) (tau10_2: ) (tau10_3: )
					 (tau11_1: ) (tau11_2: ) (tau11_3: )
					 (tau12_1: ) (tau12_2: ) (tau12_3: ) 
					 (tau13_1: ) (tau13_2: ) (tau13_3: ) 

	;
    #delimit cr
   ml init e0  
   ml maximize , difficult	
estat ic
 
 
*m3 varying tau_2 (tau8_2) in gallup 4
    #delimit ; 
    ml model lf lam (k: k= i.s, noconstant) 
	                 (tau1_1:)
					 (tau2_1:) (tau2_2:) (tau2_3:) (tau2_4:) (tau2_5:) (tau2_6:) (tau2_7:) (tau2_8:) (tau2_9:) 
					 (tau3_1:) (tau3_2:) (tau3_3:) 
					 (tau4_1:) (tau4_2:) (tau4_3:)
					 (tau5_1:) (tau5_2:) (tau5_3:) (tau5_4:)
					 (tau6_1:) (tau6_2:)
					 (tau7_1:) (tau7_2:) (tau7_3:) (tau7_4:) (tau7_5:) (tau7_6:) (tau7_7:) (tau7_8:) (tau7_9:) 
					 (tau8_1: ) (tau8_2: i.s ) (tau8_3:)
					 (tau9_1:) (tau9_2:) (tau9_3:) (tau9_4:)
					 (tau10_1:) (tau10_2:) (tau10_3:)
					 (tau11_1:) (tau11_2:) (tau11_3:)
					 (tau12_1:) (tau12_2:) (tau12_3:) 
					 (tau13_1:) (tau13_2:) (tau13_3:) 

	;
    #delimit cr
	
	
   ml init e0 
   ml maximize , difficult 
  estimates store m3
estat ic 

*m4: taus to be a linear function of time 
    #delimit ; 
    ml model lf lam (k: k= i.s, noconstant) 
	                 (tau1_1: s)
					 (tau2_1: s) (tau2_2: s) (tau2_3: s) (tau2_4: s) (tau2_5: s) (tau2_6: s) (tau2_7:s) (tau2_8:s) (tau2_9:s) 
					 (tau3_1: s) (tau3_2: s) (tau3_3: s) 
					 (tau4_1: s) (tau4_2: s) (tau4_3: s)
					 (tau5_1: s) (tau5_2: s) (tau5_3: s) (tau5_4: s)
					 (tau6_1: s) (tau6_2: s)
					 (tau7_1: s) (tau7_2: s) (tau7_3: s) (tau7_4: s) (tau7_5: s) (tau7_6: s) (tau7_7: s) (tau7_8: s) (tau7_9: s) 
					 (tau8_1: s) (tau8_2: s) (tau8_3: s)
					 (tau9_1: s) (tau9_2: s) (tau9_3: s) (tau9_4: s)
					 (tau10_1: s) (tau10_2: s) (tau10_3: s)
					 (tau11_1: s) (tau11_2: s) (tau11_3: s)
					 (tau12_1: s) (tau12_2: s) (tau12_3: s) 
					 (tau13_1: s) (tau13_2: s) (tau13_3: s) 

	;
    #delimit cr
	
   ml init e0  // intial values obtained from ordered probit rountine 
   ml maximize , difficult  
estat ic 
