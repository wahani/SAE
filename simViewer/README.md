ShinyApp_simViewer
=======================

ShinyApp for visualizing simulation results. At this time only ment for the results in the context of small area statistics. Under development - nothing working yet:

Feature list
-----------------------
### ready:
* load data button

### wanted:
* standard plots for bias-like and mse-like measures
* optional log-scale for axis
* optional reference-line
* Filter for scenarios
* Filter for data - parameters-Domains
* maybe possibility to read-in multiple files
* export for plots and tables
* table of measures for specific area
* table of results
* plots and tables for parameter estimation - density/boxplot/table
* plots concordance analysis

Data Format
----------------------
The app will expect a specific data format. Essentially you can read in a *.RData* file which will contain an object called *simData*. All other objects in that file will be ignored. *simData* is a named list containing two *data.frames*. They are called *parData*, for the results of the parameter estimates, and *estData*, containing the predicted domain-level statistic.

### *parData*
```
something
```