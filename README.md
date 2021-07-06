# Isotopic Niche Overlap : considering both intra- and interindividual variability in the identification of isotopic niches with standard ellipses
• **OBJECTIVES** : 
<br>1) **Identify isotopic niches at species level** with isotopic data including both **intra- and interindividual variability** (e.g. of study cases: serveral isotopic values measured along seal whiskers, bird feathers, tooth growth layers).</br>
<br>2) Measure **isotopic niche overlap** between species.</br>

• **METHOD** : estimate **isotopic standard ellipses** with a **novel hierarchical model developed in a Bayesian framework** to consider **two levels of isotopic variability**: an **intra-individual level** (characterised by several isotopic measurements for each individual, e.g. along a whisker) and an **interindividual** one. Our model therefore expanded on standard model with one level of isotopic variability (an interindividual one; *cf*. [Jackson et al., 2011](https://doi.org/10.1111/j.1365-2656.2011.01806.x)).

• **LANGUAGE** : R (version 4.0.2). Model implemented in software Stan (package rstan, version 2.21.0).

• **CASE STUDY** : isotopic niche overlap was quantified here using δ<sup>13</sup>C and δ<sup>15</sup>N isotopic values from the whiskers of 8 harbour and 10 grey seals. 

## Script developped as part of: 
Planque Y, Spitz J, Authier M, Guillou G, Vincent C, Caurant F. Trophic niche overlap between sympatric harbour seals (Phoca vitulina) and grey seals (Halichoerus grypus) at the southern limit of their European range (Eastern English Channel). Ecol Evol. 2021;00:1– 22. https://doi.org/10.1002/ece3.7739

## Data used in this script are freely available on SEANOE:
Planque Yann, Vincent Cécile, Guillou Gaël, Lebreton Benoit, Caurant Florence (2020). δ13C and δ15N stable isotope compositions of the whisker of 8 harbour seals (Phoca vitulina) and 10 grey seals (Halichoerus grypus) captured in the baie de Somme, France, in 2008 and 2012, for telemetry tracking. SEANOE. https://doi.org/10.17882/76528

Before running the script, please download the file repository in the ZIP file and place it on your desktop. Place the data previously dowloaded in the subfolder "Input".

First publication on GitHub : 2020-11-05

Last update : 2021-04-12 (Version 1.3)

## Authors : Yann Planque<sup>(1)</sup>, Matthieu Authier<sup>(2)(3)</sup>
 Affiliations : 
 
    (1) Centre d'Etudes Biologiques de Chizé (CEBC, UMR 7372 CNRS - La Rochelle Université), La Rochelle, France
    
    (2) Observatoire PELAGIS (UMS 3462 CNRS - La Rochelle Université), La Rochelle, France
    
    (3) ADERA, Pessac, France

## Contact : yann.planque@univ-lr.fr ; yann.planque@hotmail.fr
