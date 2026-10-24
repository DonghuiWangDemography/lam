*! Original fixed-threshold likelihood from lam-ml.do
*! Evaluator renamed to keep the production command available
program define lam_reference

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
    replace `lnf' =ln(1 - normal(`tau5_4'  -`mu'))                           if $ML_y1 ==5 & q==5


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
    replace `lnf' =ln(1 - normal(`tau8_3'  -`mu'))                           if $ML_y1 ==4 & q==8


    // 9. USKN_5
    replace `lnf' =ln(normal(`tau9_1' -`mu'))                             if $ML_y1 ==1 & q==9
    replace `lnf' =ln(normal(`tau9_2' -`mu')  - normal( `tau9_1' -`mu'))  if $ML_y1 ==2 & q==9
    replace `lnf' =ln(normal(`tau9_3' -`mu')  - normal( `tau9_2' -`mu'))  if $ML_y1 ==3 & q==9
    replace `lnf' =ln(normal(`tau9_4' -`mu')  - normal( `tau9_3' -`mu'))  if $ML_y1 ==4 & q==9
    replace `lnf' =ln(1 - normal(`tau9_4'  -`mu'))                           if $ML_y1 ==5 & q==9

    // 10. USORC_4
    replace `lnf' =ln(normal(`tau10_1' -`mu'))                              if $ML_y1 ==1 & q==10
    replace `lnf' =ln(normal(`tau10_2' -`mu')  - normal( `tau10_1' -`mu'))  if $ML_y1 ==2 & q==10
    replace `lnf' =ln(normal(`tau10_3' -`mu')  - normal( `tau10_2' -`mu'))  if $ML_y1 ==3 & q==10
    replace `lnf' =ln(1 - normal(`tau10_3'  -`mu'))                         if $ML_y1 ==4 & q==10

    //11. USPSRA_4
    replace `lnf' =ln(normal(`tau11_1' -`mu'))                              if $ML_y1 ==1 & q==11
    replace `lnf' =ln(normal(`tau11_2' -`mu')  - normal( `tau11_1' -`mu'))  if $ML_y1 ==2 & q==11
    replace `lnf' =ln(normal(`tau11_3' -`mu')  - normal( `tau11_2' -`mu'))  if $ML_y1 ==3 & q==11
    replace `lnf' =ln(1 - normal(`tau11_3'  -`mu'))                         if $ML_y1 ==4 & q==11


    //12 zoby4
    replace `lnf' =ln(normal(`tau12_1' -`mu'))                              if $ML_y1 ==1 & q==12
    replace `lnf' =ln(normal(`tau12_2' -`mu')  - normal( `tau12_1' -`mu'))  if $ML_y1 ==2 & q==12
    replace `lnf' =ln(normal(`tau12_3' -`mu')  - normal( `tau12_2' -`mu'))  if $ML_y1 ==3 & q==12
    replace `lnf' =ln(1 - normal(`tau12_3'  -`mu'))                         if $ML_y1 ==4 & q==12

    // abcwp4
    replace `lnf' =ln(normal(`tau13_1' -`mu'))                              if $ML_y1 ==1 & q==13
    replace `lnf' =ln(normal(`tau13_2' -`mu')  - normal( `tau13_1' -`mu'))  if $ML_y1 ==2 & q==13
    replace `lnf' =ln(normal(`tau13_3' -`mu')  - normal( `tau13_2' -`mu'))  if $ML_y1 ==3 & q==13
    replace `lnf' =ln(1 - normal(`tau13_3'  -`mu'))                         if $ML_y1 ==4 & q==13

  }
 end
