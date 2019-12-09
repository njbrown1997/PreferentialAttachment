### NOTES DEC 9 MONDAY

Made boxplots of assortativity coefficients (Figure 2)

### NOTES DEC 8
Added all age functions for n repetitions, sweeping through parameters

Changed poisson distribution to normal distribution for that age function, as its a bit simpler to adjust and still has a curve.

Plots CCDFs for each function (w/different parameter settings)

Created build_network.m for reducing code repetition

Created ccdf.m for finding CCDF

Added a few power law files that can be used to quantify the CCDFs shown. May or may not be able to acquire much info, since most probably aren't closely fit with a power-law

To do next: more network statistics, now that networks are generated

### NOTES DEC 2

Driver broken into sections by age function. Each section sweeps relevant parameter and stores network. 

Plotting enabled to get an idea of how parameters affect the networks.

Network size held constant across all settings.

To be done: 
-implement normalized Poisson aging
-calculate network features
-add/adjust analysis as needed
-figure generation of network features
