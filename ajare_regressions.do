/*------------------------------------------------------------------------------
Replication code for discrete choice experiment regressions

1) run all regressions in paper
2) create all regression tables in paper
3) show regressions that did not converge that are not in paper
________________________________________________________________________________
File ~root\regressions.do
Last Modifed: 3-29-2021
Created by: Daniel Brent
------------------------------------------------------------------------------*/

/*** stata setup ***/
clear matrix 
clear mata
clear
set matsize 11000
set maxvar 30000
set more off
capture log close

/* set path */
global path "/Users/dab320/Dropbox/Monash/CRC/CRCshare/Salient/DAB/Empirics_AJARE"
cd ${path}

/* load data */
set more off, perm
use Data/data_ajare, clear

/* create negative cost for lognormal regression */
gen lncost = -1*cost

/* set seed for replication */
set seed 365247

/*********************************************************************
********             TABLE 3 - main regressions               ********
**********************************************************************/

/*** column (1) - base regression ***/
#delimit ;
mixlogit choice_water status_quo, 
group(caseid) id(id) cluster(id) nrep(5000) burn(500)
rand(
cost
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo);
estimates save Estimates/base_norm_nocorr5000, replace; /* save estimates */

/*** column (2) - cost*treatment ***/
#delimit ;
mixlogit choice_water status_quo
cost_salient, 
group(caseid) id(id) cluster(id) nrep(5000) burn(500)
rand(
cost
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo);
estimates save Estimates/treatint_norm_nocorr5000, replace; /* save estimates */

/*** column (3) - cost*treatment*income ***/
#delimit ;
mixlogit choice_water status_quo
inc_low_cost inc_high_cost
cost_salient_inc_low cost_salient_inc_med cost_salient_inc_high, 
group(caseid) id(id) cluster(id) nrep(5000) burn(500)
rand(
cost
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo);
estimates save Estimates/treatint_het_norm_nocorr5000, replace; /* save estimates */

/*** column (4) - attributes*treatment ***/
#delimit ;
mixlogit choice_water status_quo
alt_flood_both_salient alt_water_stage3_4_salient
alt_water_none_salient alt_stream_medium_salient
alt_stream_high_salient alt_req_both_salient 
alt_temp_minustwo_salient, 
group(caseid) id(id) cluster(id) nrep(5000) burn(500)
rand(
cost
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo);
estimates save Estimates/treatint_att_norm_nocorr5000, replace; /* save estimates */

/***********************
***  CREATE TABLE 3  ***
***********************/

/* load estimates */
#delimit ;
eststo clear;
estimates use Estimates/base_norm_nocorr5000;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;
estimates use Estimates/treatint_norm_nocorr5000;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;
estimates use Estimates/treatint_het_norm_nocorr5000;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;
estimates use Estimates/treatint_att_norm_nocorr5000;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;

/* view table */
#delimit ;
esttab, 
se(4) b(4) star(* 0.10 ** 0.05 *** 0.01) nonotes label nodiscrete  
scalars("bicn BIC/N" "aicn AIC/N" "N Observations" "N_clust Individuals") sfmt(%10.0fc)
noobs compress nogaps unstack
mtitles("Base" "Cost" "Cost*Income" "Attributes")
mgroups("None" "Salient Interactions", pattern(0 1 0 0) 
prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}));

/* export to tex */
#delimit ;
esttab using Output/main.tex, replace
se(4) b(4) star(* 0.10 ** 0.05 *** 0.01) nonotes label nodiscrete  
scalars("bicn BIC/N" "aicn AIC/N" "N Observations" "N_clust Individuals") sfmt(%10.0fc)
noobs compress nogaps unstack
mtitles("Base" "Cost" "Cost*Income" "Attributes")
mgroups("None" "Salient Interactions", pattern(0 1 0 0) 
prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}));







/*********************************************************************
********          TABLE A.4 - specificcation tests            ********
**********************************************************************/

/*** column (1) - normal cost, uncorrelated (replicates column (1) of Table (3) ***/
#delimit ;
mixlogit choice_water status_quo, 
group(caseid) id(id) cluster(id) nrep(5000) burn(500)
rand(
cost
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo);
estimates save Estimates/base_norm_nocorr, replace;


/*** column (2) - normal cost, correlated ***/
#delimit ;
mixlogit choice_water status_quo, 
group(caseid) id(id) cluster(id)  nrep(5000) burn(500) corr
rand(
cost
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo);
estimates save Estimates/base_norm_corr, replace;

/*** column (3) - fixed cost, uncorrelated ***/
#delimit ;
mixlogit choice_water status_quo 
cost, 
group(caseid) id(id) cluster(id) nrep(5000) burn(500)
rand(
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo);
estimates save Estimates/base_fixed_nocorr, replace;

/*** column (4) - fixed income*cost, uncorrelated ***/
#delimit ;
mixlogit choice_water status_quo 
inc_low_cost inc_med_cost inc_high_cost, 
group(caseid) id(id) cluster(id) nrep(5000) burn(500)
rand(
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo);
estimates save Estimates/base_fixedhet_nocorr, replace;

/*** column (5) - lognormal cost, uncorrelated ***/

#delimit ;
mixlogit choice_water status_quo,
group(caseid) id(id) cluster(id) nrep(5000) burn(500)
rand(alt_flood_both alt_water_stage3_4
alt_water_none alt_stream_medium alt_stream_high
alt_req_both alt_temp_minustwo lncost) ln(1);
estimates save Estimates/base_ln_nocorr, replace;
eststo base_ln_nocorr;

/* mean of cost distribution */
#delimit ;
nlcom (mean_cost: -1*exp([Mean]_b[lncost]+0.5*[SD]_b[lncost]^2))
(med_cost: -1*exp([Mean]_b[lncost]))
(sd_cost: exp([Mean]_b[lncost]+0.5*[SD]_b[lncost]^2)*
sqrt(exp([SD]_b[lncost]^2)-1));

/* median of cost distribution */
#delimit ;
nlcom (med_cost: -1*exp([Mean]_b[lncost]))
(water_none: [Mean]_b[alt_water_none])
(wtp_water_none:  [Mean]_b[alt_water_none]/exp([Mean]_b[lncost]));


/*************************
***  CREATE TABLE A.4  ***
*************************/

/*** load estimates ***/
eststo clear;
estimates use Estimates/base_norm_nocorr;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;
estimates use Estimates/base_norm_corr;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;
estimates use Estimates/base_fixed_nocorr;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;
estimates use Estimates/base_fixedhet_nocorr;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;
estimates use Estimates/base_ln_nocorr;
/* add mean and median of underlying cost distribution */
nlcom (mean_cost: -1*exp([Mean]_b[lncost]+0.5*[SD]_b[lncost]^2)) 
(med_cost: -1*exp([Mean]_b[lncost]));
estadd scalar mean_cost = r(b)[1,1];
estadd scalar med_cost = r(b)[1,2];
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;

/* label variables */
label var lncost "LN(Cost)";
label var inc_med_cost "Med_Income*Cost";

/* view table**/
#delimit ;
esttab, 
se(4) b(4) star(* 0.10 ** 0.05 *** 0.01) nonotes label nodiscrete  
scalars("mean_cost Mean Cost" "med_cost Median Cost" "bicn BIC/N" "aicn AIC/N"
"N Observations" "N_clust Individuals")  sfmt(%10.0gc %10.0gc %10.0fc)
 noobs compress nogaps unstack 
drop(l11: l21: l31: l41: l51: l61: l71: l81: l22: l32: l42: l52: l62: l72: l82: 
l33: l43: l53: l63: l73: l83: l44: l54: l64: l74: l84: l55: l65: l75: l85: 
l66: l76: l86: l77: l87: l88:)
mtitles("Normal" "Normal Corr" "Fixed" "Fixed*Income" "Lognormal");

/* export to tex */
#delimit ;
esttab using Output/spec.tex, replace
se(4) b(4) star(* 0.10 ** 0.05 *** 0.01) nonotes label nodiscrete  
scalars("mean_cost Mean Cost" "med_cost Median Cost" "bicn BIC/N" "aicn AIC/N"
"N Observations" "N_clust Individuals")  sfmt(%10.0gc %10.0gc %10.0fc)
 noobs compress nogaps unstack 
drop(l11: l21: l31: l41: l51: l61: l71: l81: l22: l32: l42: l52: l62: l72: l82: 
l33: l43: l53: l63: l73: l83: l44: l54: l64: l74: l84: l55: l65: l75: l85: 
l66: l76: l86: l77: l87: l88:)
mtitles("Normal" "Normal Corr" "Fixed" "Fixed*Income" "Lognormal");






/*********************************************************************
********            TABLE A.5 - robustness to seed            ********
**********************************************************************/
/* set new seed */
#delimit ;
set seed 3652472;

/*** now replicate Table 3 regressions ***/

/*** column (1) - base regression ***/
mixlogit choice_water status_quo, 
group(caseid) id(id) cluster(id) nrep(5000) burn(500)
rand(
cost
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo);
estimates save Estimates/base_norm_nocorr5000_seed2, replace;

/*** column (2) - cost*treatment ***/
#delimit ;
mixlogit choice_water status_quo
cost_salient, 
group(caseid) id(id) cluster(id) nrep(5000) burn(500)
rand(
cost
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo);
estimates save Estimates/treatint_norm_nocorr5000_seed2, replace;

/*** column (3) - cost*treatment*income ***/
#delimit ;
mixlogit choice_water status_quo
inc_low_cost inc_high_cost
cost_salient_inc_low cost_salient_inc_med cost_salient_inc_high, 
group(caseid) id(id) cluster(id) nrep(5000) burn(500)
rand(
cost
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo);
estimates save Estimates/treatint_het_norm_nocorr5000_seed2, replace;

/*** column (4) - attributes*treatment ***/
#delimit ;
mixlogit choice_water status_quo
alt_flood_both_salient alt_water_stage3_4_salient
alt_water_none_salient alt_stream_medium_salient
alt_stream_high_salient alt_req_both_salient 
alt_temp_minustwo_salient, 
group(caseid) id(id) cluster(id) nrep(5000) burn(500)
rand(
cost
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo);
estimates save Estimates/treatint_att_norm_nocorr5000_seed2, replace;

/*************************
***  CREATE TABLE A.5  ***
*************************/

/*** load estimates ***/
#delimit ;
eststo clear;
estimates use Estimates/base_norm_nocorr5000_seed2;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;
estimates use Estimates/treatint_norm_nocorr5000_seed2;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;
estimates use Estimates/treatint_het_norm_nocorr5000_seed2;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;
estimates use Estimates/treatint_att_norm_nocorr5000;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;

/* view table */
#delimit ;
esttab, 
se(4) b(4) star(* 0.10 ** 0.05 *** 0.01) nonotes label nodiscrete  
scalars("bicn BIC/N" "aicn AIC/N" "N Observations" "N_clust Individuals") sfmt(%10.0fc)
noobs compress nogaps unstack
mtitles("Base" "Cost" "Cost*Income" "Attributes")
mgroups("None" "Salient Interactions", pattern(0 1 0 0) 
prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}));

/* export to tex */
#delimit ;
esttab using Output/robust_seed.tex, replace
se(4) b(4) star(* 0.10 ** 0.05 *** 0.01) nonotes label nodiscrete  
scalars("bicn BIC/N" "aicn AIC/N" "N Observations" "N_clust Individuals") sfmt(%10.0fc)
noobs compress nogaps unstack
mtitles("Base" "Cost" "Cost*Income" "Attributes")
mgroups("None" "Salient Interactions", pattern(0 1 0 0) 
prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}));



/*************************************************************************
********         TABLE A.6 - robustness to starting values        ********
**************************************************************************/
/* reset original seed */
#delimit ;
set seed 365247;

/* load primary treatment effect regression (column (2) of Table 3) */
estimates use Estimates/treatint_norm_nocorr5000;

/* extract starting values */
matrix b = e(b);
matrix e = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

/* permute starting values */
#delimit ;
forvalues j = 1/4 {; /* 4 different permutations */
forvalues i = 1/17 {; /* each of the 17 different paramters */
matrix e[1,`i'] = runiform(-abs(b[1,`i'])/`j',abs(b[1,`i'])/`j'); /* generate disturbance */
};
matrix b`j' = b + e; /* create new starting value */
};

/* run regressions with new starting values */
#delimit ;
forvalues j = 1/4 {;
mixlogit choice_water status_quo
cost_salient, 
group(caseid) id(id) cluster(id) nrep(5000) burn(500)
from(b`j',copy)
rand(
cost
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo);
estimates save Estimates/treatint_norm_nocorr5000_start`j', replace;
};


/*************************
***  CREATE TABLE A.6  ***
*************************/

/* load estimates */
eststo clear;
estimates use Estimates/treatint_norm_nocorr5000_start4;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;
estimates use Estimates/treatint_norm_nocorr5000_start3;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;
estimates use Estimates/treatint_norm_nocorr5000_start2;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;
estimates use Estimates/treatint_norm_nocorr5000_start1;
estadd scalar aicn = (2*e(k) - 2*e(ll))/e(N_clust); /* adjust aic to aic/N */
estadd scalar bicn = (ln(e(N))*e(k) - 2*e(ll))/e(N_clust); /* adjust bic to bic/N */
estadd scalar N = e(N)/3, replace; /* correct # of observations */
eststo;

/* view table */
#delimit ;
esttab, 
se(4) b(4) star(* 0.10 ** 0.05 *** 0.01) nonotes label nodiscrete  
scalars("bicn BIC/N" "aicn AIC/N" "N Observations" "N_clust Individuals") sfmt(%10.0fc)
noobs compress nogaps unstack
mtitles("Start 0.25" "Start 0.33" "Start 0.5" "Start 1");

/* export to tex */
#delimit ;
esttab using Output/robust_start.tex, replace
se(4) b(4) star(* 0.10 ** 0.05 *** 0.01) nonotes label nodiscrete  
scalars("bicn BIC/N" "aicn AIC/N" "N Observations" "N_clust Individuals") sfmt(%10.0fc)
noobs compress nogaps unstack
mtitles("Start 0.25" "Start 0.33" "Start 0.5" "Start 1");




/**************************************************************/
************          SCALE HETEROGENEITY          *************
************   (NOT IN PAPER - DID NOT CONVERGE    *************
***************************************************************/

/*** SCALE HET ***/
#delimit ;
gmnl choice_water status_quo 
cost
alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo,
het(salient)
group(caseid) id(id) vce(cluster id) 
difficult nrep(500) burn(50);
estimates save Estimates/smnl_treat_fixed_nocorr, replace;

/*** GMNL ***/
# delimit ;
gmnl choice_water status_quo, gamma(1)
rand(cost alt_water_stage3_4 alt_water_none 
alt_flood_both 
alt_stream_high  alt_stream_medium
alt_req_both
alt_temp_minustwo)
het(salient)
group(caseid) id(id) vce(cluster id)
difficult nrep(500) burn(50);
estimates save Estimates/gmnl_treat_norm_nocorr, replace;




