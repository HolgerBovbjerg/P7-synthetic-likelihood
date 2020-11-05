The code contained in this folder implements BSL and uBSL for the Ricker model example from Section 4.3 of the main paper.
For completeness, we have also included our ABC implementation.
The Ricker model was presented in Wood (2010), and the model simulation and summary statistics have been replicated here.

For more details on model, refer to the main paper or to Wood (2010).

The 'run.m' file contains everything required to run BSL, uBSL and ABC for this example.

The files contained here are:

simulate_ricker			-	simulates from the Ricker model
ricker_summstats		-	computes the same summary statistics as in Wood (2010)
sl_log_like_ghuryeolkin	-	gets the unbiased estimator for the synthetic likelihood (when the summary statistics are MVN)
bayes_sl_ricker_wood	-	performs MCMC BSL
bayes_sl_ricker_wood_go	-	performs MCMC uBSL
abc_ricker_wood			-	performs ABC
