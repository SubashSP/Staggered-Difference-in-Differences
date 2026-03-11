/*******************************************************************************
*
*  STAGGERED DIFFERENCE-IN-DIFFERENCES
*  Callaway & Sant'Anna (2021) Estimator — Full Implementation
*  ---------------------------------------------------------------
*  Dataset : minwage_panel.csv (simulated minimum wage panel)
*  Outcome : lemp    — log teen employment
*  Cohorts : States adopt MW increases in 2004, 2007, 2010, or 2013
*  Period  : 2000–2019 (balanced panel, 50 states × 20 years)
*
*  PACKAGES REQUIRED (install once):
*    ssc install csdid
*    ssc install drdid
*    ssc install bacondecomp
*    ssc install coefplot
*    ssc install pretrends
*    ssc install eventstudyinteract
*    ssc install reghdfe
*    ssc install ftools
*    ssc install estout
*
*  NOTE: csdid requires Stata 16+
*
*******************************************************************************/


/*******************************************************************************
*  PART 0: SETUP
*******************************************************************************/

clear all
set more off
capture log close

** ─── SET YOUR WORKING DIRECTORY ──────────────────────────────────────────────
** Replace with the path to the cloned repository folder
**   Windows:  cd "C:\Users\YourName\staggered-did"
**   Mac/Linux: cd "/Users/YourName/staggered-did"
** ─────────────────────────────────────────────────────────────────────────────
* cd "CHANGE_THIS_PATH"

log using "outputs/staggered_did_log.txt", text replace

** Install packages if not already installed
* ssc install csdid
* ssc install drdid
* ssc install bacondecomp
* ssc install coefplot
* ssc install reghdfe
* ssc install ftools
* ssc install estout


/*******************************************************************************
*  PART 1: LOAD & PREPARE DATA
*******************************************************************************/

** Import the simulated panel dataset
import delimited "data/minwage_panel.csv", clear

** Declare panel structure
xtset state_id year
xtdescribe

** Quick summary
describe
sum
tab gvar, missing


/*******************************************************************************
*  PART 2: DATA QUALITY CHECKS
*******************************************************************************/

** Check for duplicates
duplicates report state_id year

** Verify irreversibility: once treated, always treated
sort state_id year
by state_id: gen reversal = (mw_treat < L.mw_treat)
tab reversal    // Should be all zeros
drop reversal

** Balance check: compare treated vs. never-treated pre-treatment
ttest lpop    if year == 2000, by(ever_treated)
ttest lincome if year == 2000, by(ever_treated)

** Distribution of cohorts
tab gvar, missing

** Pre-treatment summary by cohort
tabstat lemp lpop lincome if year == 2000, by(gvar) stats(n mean sd)


/*******************************************************************************
*  PART 3: NAIVE TWFE AS BENCHMARK
*******************************************************************************/

** Standard Two-Way Fixed Effects DiD (BENCHMARK — potentially biased)
reghdfe lemp mw_treat, absorb(state_id year) cluster(state_id)
estimates store twfe
estadd local fe "State + Year"

di "TWFE ATT estimate: " _b[mw_treat]
di "NOTE: This may be biased due to heterogeneous treatment effects."


/*******************************************************************************
*  PART 4: GOODMAN-BACON DECOMPOSITION
*******************************************************************************/

** Decompose the TWFE estimate into its 2x2 DiD components
** Reveals how much weight falls on 'Already Treated vs Timing' comparisons
bacondecomp lemp mw_treat, ddetail msymbols(O T)
graph export "figures/bacon_decomposition.png", replace width(1200)

di "If large share of weight is from 'Already Treated vs Timing',"
di "  this is strong motivation to use CS-DiD instead of TWFE."


/*******************************************************************************
*  PART 5: TWFE EVENT STUDY (BENCHMARK)
*******************************************************************************/

** Construct event-time variable (relative time to treatment)
gen event_time = year - treat_year if treat_year != .
replace event_time = -99 if event_time == .    // code never-treated

** TWFE event study: use l=-1 as the base/reference period
reghdfe lemp ib(-1).event_time if event_time >= -4 & event_time <= 5, ///
    absorb(state_id year) cluster(state_id)
estimates store twfe_es

di "TWFE event study stored as: twfe_es"


/*******************************************************************************
*  PART 6: CALLAWAY & SANT'ANNA (2021) — MAIN ESTIMATION
*******************************************************************************/

** ── 6a. Unconditional Parallel Trends (no covariates) ───────────────────────
csdid lemp, ivar(state_id) tvar(year) gvar(gvar) ///
    method(dripw) notyet cluster(state_id)
estimates store cs_nocov

** ── 6b. Conditional Parallel Trends (with covariates) ───────────────────────
** lpop and lincome are time-varying but used to control for secular trends
csdid lemp lpop lincome, ivar(state_id) tvar(year) gvar(gvar) ///
    method(dripw) notyet cluster(state_id)
estimates store cs_cov

di "CS estimates stored: cs_nocov, cs_cov"


/*******************************************************************************
*  PART 7: AGGREGATIONS — INTERPRETING ATT(g,t)
*******************************************************************************/

** Restore the main CS model (with covariates)
estimates restore cs_cov

** ── 7a. Simple (overall) ATT ─────────────────────────────────────────────────
estat simple
** Interpretation: single weighted-average ATT across all groups and periods

** ── 7b. Group/Cohort ATTs ────────────────────────────────────────────────────
estat group
** Shows: average ATT for each treatment cohort
** Detects whether early vs. late adopters respond differently

** ── 7c. Calendar-Time ATTs ───────────────────────────────────────────────────
estat calendar
** Shows: ATT at each calendar year, averaged over treated cohorts

** ── 7d. Dynamic / Event-Study ATTs ──────────────────────────────────────────
estat event
** MOST INFORMATIVE:
**   l < 0 (Tm*) = pre-treatment periods → should be ~0 if PT holds
**   l = 0 (T0)  = effect at treatment onset
**   l > 0 (T1+) = dynamic post-treatment effects

** Save event-study estimates as e() scalars for plotting
estat event, post
estimates store cs_event


/*******************************************************************************
*  PART 8: PRE-TRENDS TESTING
*******************************************************************************/

** Restore and re-run to get event-study with post for testing
estimates restore cs_cov
estat event, post

** Test individual pre-treatment coefficients (Tm2, Tm1 = l=-2,-1)
** Note: l=-1 is often the reference period; Tm2 is typically the first testable pre-period
capture test _b[Tm2] = 0
capture test _b[Tm3] = 0

** Joint test of all pre-treatment periods
** (Update coefficient names based on your actual output — check with: matrix list r(b))
** Example: test _b[Tm3] = 0 & _b[Tm2] = 0

** Using pretrends package (if installed)
* pretrends, numpre(3) power(0.5)
** Tests: are pre-trends consistent with a parallel-trends violation?

di "PRE-TRENDS INTERPRETATION:"
di "  p > 0.05 for pre-period coefficients = no evidence against parallel trends"
di "  p < 0.05 = potential violation; re-examine design or use conditional PT"


/*******************************************************************************
*  PART 9: EVENT STUDY PLOTS
*******************************************************************************/

** ── 9a. Built-in csdid_plot ──────────────────────────────────────────────────
csdid lemp lpop lincome, ivar(state_id) tvar(year) gvar(gvar) ///
    method(dripw) notyet cluster(state_id) asinr

csdid_plot, ///
    title("Callaway-Sant'Anna Event Study") ///
    subtitle("Effect of Minimum Wage Increase on Log Teen Employment") ///
    xtitle("Years Relative to Treatment Onset") ///
    ytitle("ATT Estimate (log employment)") ///
    xline(0, lpattern(dash) lcolor(red)) ///
    yline(0, lpattern(dot) lcolor(black)) ///
    scheme(s2color)

graph export "figures/cs_event_study.png", replace width(1200)

** ── 9b. TWFE vs. CS comparison with coefplot ─────────────────────────────────
estimates restore cs_event    // CS event study
coefplot (twfe_es, label(TWFE) mcolor(red) ciopts(lcolor(red%40))) ///
         (cs_event, label(CS DiD) mcolor(navy) ciopts(lcolor(navy%40))), ///
    vertical yline(0, lpattern(dot)) ///
    xline(4.5, lwidth(medium) lpattern(dash) lcolor(gs8)) ///
    title("TWFE vs. Callaway-Sant'Anna Event Study") ///
    xtitle("Event Time (years since MW increase)") ///
    ytitle("Estimate") ///
    scheme(s2color)

graph export "figures/twfe_vs_cs_comparison.png", replace width(1200)


/*******************************************************************************
*  PART 10: SENSITIVITY ANALYSIS
*******************************************************************************/

** ── 10a. Only never-treated as controls (more conservative) ─────────────────
csdid lemp lpop lincome, ivar(state_id) tvar(year) gvar(gvar) ///
    method(dripw) notyetnever cluster(state_id)
estat simple
estimates store cs_nevertreated

** ── 10b. Regression Adjustment only (alternative method) ────────────────────
csdid lemp lpop lincome, ivar(state_id) tvar(year) gvar(gvar) ///
    method(reg) notyet cluster(state_id)
estat simple
estimates store cs_ra

** ── 10c. IPW only (alternative method) ──────────────────────────────────────
csdid lemp lpop lincome, ivar(state_id) tvar(year) gvar(gvar) ///
    method(ipw) notyet cluster(state_id)
estat simple
estimates store cs_ipw

** ── 10d. Wild bootstrap inference (recommended for small N clusters) ─────────
** (slower — uncomment when needed)
* csdid lemp lpop lincome, ivar(state_id) tvar(year) gvar(gvar) ///
*     method(dripw) notyet wboot reps(999) seed(42)
* estat simple


/*******************************************************************************
*  PART 11: RESULTS TABLE
*******************************************************************************/

** ── Restore each model and get simple ATT ────────────────────────────────────
estimates restore cs_nocov
estat simple, post
estimates store att_nocov

estimates restore cs_cov
estat simple, post
estimates store att_cov

estimates restore cs_nevertreated
estat simple, post
estimates store att_nevertreated

estimates restore cs_ra
estat simple, post
estimates store att_ra

** ── Export comparison table ──────────────────────────────────────────────────
esttab twfe att_nocov att_cov att_nevertreated att_ra, ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    mtitles("TWFE" "CS No Cov" "CS Cov (DR)" "CS Never-Treated" "CS Reg-Adj") ///
    title("Effect of Minimum Wage Increase on Log Teen Employment") ///
    note("Doubly-robust (dripw) unless noted. Clustered SEs at state level.") ///
    using "outputs/results_table.rtf", replace

esttab twfe att_nocov att_cov att_nevertreated att_ra, ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    mtitles("TWFE" "CS No Cov" "CS Cov (DR)" "CS Never-Treated" "CS Reg-Adj") ///
    title("Effect of Minimum Wage Increase on Log Teen Employment") ///
    note("Doubly-robust (dripw) unless noted. Clustered SEs at state level.")


/*******************************************************************************
*  PART 12: QUICK REFERENCE SUMMARY
*******************************************************************************/

di ""
di "============================================================"
di "  CSDID QUICK REFERENCE"
di "============================================================"
di ""
di "  BASIC:     csdid Y, ivar(id) tvar(t) gvar(g)"
di "  WITH COV:  csdid Y X1 X2, ivar(id) tvar(t) gvar(g) method(dripw)"
di ""
di "  AGGREGATIONS:"
di "    estat simple    // overall ATT"
di "    estat group     // ATT by cohort"
di "    estat calendar  // ATT by calendar time"
di "    estat event     // ATT by event time (event study)"
di ""
di "  PLOTTING:"
di "    csdid Y, ivar(id) tvar(t) gvar(g) asinr"
di "    csdid_plot"
di "============================================================"


/*******************************************************************************
*  CLOSE LOG
*******************************************************************************/

log close
