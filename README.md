# AEC-Modelling
DTU Alkaline electrolyser modelling built during visitng PhD period @ DTU. <br>

The objective is to build a model that can describe the polarisation curve, which depends on the setting parameters such as: <br>
- Wt% of the electrolyte, in this specific case we are talking about KOH 
- Temperature
- Operating pressure, which directly influences electrolyte density

## Files description

There are basically 4 types of files in this repository <br>
- Main, is the main operating code 
- func_class, is the classs of the actual AEC model, containing every function built 
- dataframe and txt files, are polynomial or Arrenhius coefficients for $i_0$ and $\alpha$, obtained through the fitting of the experimental data from Kibria and Miles experiments
- .csv files, are data from bachelor thesis ( including both Comsol and experimental data) and the one from Henrik

### References
Below i report all literature used, divided by category:

