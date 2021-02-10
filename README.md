# Model development for APRANK

Here you can find all you need to recreate the two models used by APRANK to predict antigenicity in proteins and peptides.

Inside the **code** folder you will find several scripts that should be run in the order indicated by their numeric prefixes.

The first of the scripts will run APRANK for all the 15 species, so APRANK is needed. Also, keep in mind that this step might take a long while depending on your PC (you can set the number of parallel processes in the **Config - Run Predictors** sections to increase the speed; you can also edit the **organisms** variable to calculate some species one day and some another).

You will have to edit the first lines of all the scripts so they match to the location of both these files and APRANK, such as:

```r
aprank_development_folder <- "/home/user/Desktop/APRANK/ModelDevelopment"
aprank_folder <- "/home/user/Desktop/APRANK"
```
