---
author: Michael Matschiner
level: Intermediate
title: CladeAge Tutorial
subtitle: Fossil-based divergence-time estimation with CladeAge
beastversion: 2.5.0
---


# Background

In Bayesian divergence-time estimation, phylogenies are commonly time calibrated through the specification of calibration densities on nodes representing clades with known fossil occurrences. Unfortunately, the optimal shape of these calibration densities is usually unknown and they are therefore often chosen arbitrarily, which directly impacts the reliability of the resulting age estimates. CladeAge overcomes this limitation by calculating optimal calibration densities for clades with fossil records, based on estimates for diversification rates and the sampling rate for fossils. CladeAge thus shares similarities with the Fossilized Birth-Death (FBD) process; however, while the FBD model assumes that the fossil record is either completely or randomly sampled, the CladeAge model assumes that only information about the oldest fossil of each clade is available. CladeAge allows uncertainty in the diversification- and sampling-rate estimates, but unlike with the FBD model, these parameters must be known _a priori_ and can not be estimated as part of the analysis.


<!--Please start the tutorial by adding some background about the tutorial in this section, clearly explaining the question/problem and the type of analysis that the methods in the tutorial should be used for. In the next section please add a short description of all the programs or packages used in the tutorial. The tutorial exercise should follow this part. Please add a short explanation on the dataset used in the tutorial before starting with the exercise. Please also add a section after the exercise interpreting the results. End your tutorial with some useful links.-->

----

# Programs used in this Exercise 

### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

BEAST2 is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees {% cite Bouckaert2014 --file CladeAge-Tutorial/master-refs.bib %}. This tutorial uses the BEAST2 version 2.5.0.

### BEAUti2 - Bayesian Evolutionary Analysis Utility

BEAUti2 is a graphical user interface tool for generating BEAST2 XML configuration files.

Both BEAST2 and BEAUti2 are Java programs, which means that the exact same code runs on all platforms. For us it simply means that the interface will be the same on all platforms. The screenshots used in this tutorial are taken on a Mac OS X computer; however, both programs will have the same layout and functionality on both Windows and Linux. BEAUti2 is provided as a part of the BEAST2 package so you do not need to install it separately.

### TreeAnnotator

TreeAnnotator is used to summarise the posterior sample of trees to produce a maximum clade credibility tree. It can also be used to summarise and visualise the posterior estimates of other tree parameters (e.g. node height).

TreeAnnotator is provided as a part of the BEAST2 package so you do not need to install it separately.

### Tracer

Tracer ([http://tree.bio.ed.ac.uk/software/tracer](http://tree.bio.ed.ac.uk/software/tracer)) is used to summarize the posterior estimates of the various parameters sampled by the Markov Chain. This program can be used for visual inspection and to assess convergence. It helps to quickly view median estimates and 95% highest posterior density intervals of the parameters, and calculates the effective sample sizes (ESS) of parameters. It can also be used to investigate potential parameter correlations. We will be using Tracer v1.6.0.

### FigTree

FigTree ([http://tree.bio.ed.ac.uk/software/figtree](http://tree.bio.ed.ac.uk/software/figtree)) is a program for viewing trees and producing publication-quality figures. It can interpret the node-annotations created on the summary trees by TreeAnnotator, allowing the user to display node-based statistics (e.g. posterior probabilities). We will be using FigTree v1.4.2.


----

# Practical: Fossil-based divergence-time estimation with CladeAge

In this tutorial, we are going to use a multi-marker sequence dataset in combination with fossil information to estimate the divergence times of cichlid fishes with the CladeAge model.

The aim of this tutorial is to:

- Learn which fossil information is required to apply the CladeAge model.
- Learn how to estimate divergence times with CladeAge.


## The Data

The analyses in this tutorial will be based on the dataset of Near et al. {% cite Near2013 -A --file CladeAge-Tutorial/master-refs.bib %}, comprising alignments for ten nuclear genes. In their study, Near et al. {% cite Near2013 -A --file CladeAge-Tutorial/master-refs.bib %} used this dataset to estimate divergence times of spiny-rayed fishes (=Acanthomorphata), and to identify shifts in diversification rates among different groups of these fishes. While the dataset of Near et al. {% cite Near2013 -A --file CladeAge-Tutorial/master-refs.bib %} did not focus on cichlid diversification, it included nine cichlid species among the 520 species sampled for the extensive phylogeny. Thus, we can here use part of this dataset of Near et al. {% cite Near2013 -A --file CladeAge-Tutorial/master-refs.bib %} to estimate early divergences among cichlid fishes. To facilitate the analyses in this tutorial, we will reduce the dataset of Near et al. {% cite Near2013 -A --file CladeAge-Tutorial/master-refs.bib %} to sequences of 24 selected species. These species represent divergent cichlid lineages as well as the most ancestral groups of spiny-rayed fishes so that the fossil record of these lineages can be employed for calibration.

## The CladeAge model

The approach of CladeAge is similar to the more traditional node-dating approach in the sense that prior densities are defined for the ages of different clades, and the minimum ages of these prior densities are provided by the oldest fossils of these clades. However, important differences exist between the CladeAge approach and node dating: First, the shape of age-prior densities is informed by a model of diversification and fossil sampling in the CladeAge approach, whereas in node dating, parametric distributions (e.g. lognormal or gamma) with more or less arbitrarily chosen parameters were usually applied. Second, because of the quantitative model used in the CladeAge approach, the clades used for calibration should also not be chosen at will. Instead, strictly all clades included in the phylogeny that (i) have a fossil record, (ii) are morphologically recognizable, and (iii) have their sister lineage also included in the phylogeny should be constrained according to the age of their oldest fossil. A consequence of this is that clades are constrained even when their known sister lineage has an older fossil record, and that the same fossil may be used to constrain not just one clade, but multiple nested clades, if the more inclusive clades do not have an even older fossil record. More details about these criteria can be found in our paper on CladeAge {% cite Matschiner2017 --file CladeAge-Tutorial/master-refs.bib %}, and further information on using CladeAge is given in our [Rough Guide to CladeAge](http://evoinformatics.eu/cladeage.pdf).

## Divergence-time estimation with CladeAge

### Preparing the molecular dataset

As a first step, download the molecular dataset of Near et al. {% cite Near2013 -A --file CladeAge-Tutorial/master-refs.bib %} from the Dryad data repository connected to their publication. On the command-line of a Mac or Linux computer, this could be done using the `wget` utility.

> Execute the following command:

```
wget https://datadryad.org/bitstream/handle/10255/dryad.50839/Near_et_al.nex
```

Alternatively, you could use your browser to open the Dryad repository [https://datadryad.org/resource/doi:10.5061/dryad.d3mb4](https://datadryad.org/resource/doi:10.5061/dryad.d3mb4) and download the file `Near_et_al.nex`.

This file in Nexus format contains the sequence alignments for ten nuclear markers, sequenced for 608 species of spiny-rayed fishes. As this dataset is far too large to be analyzed in this tutorial, we'll extract the sequences of 24 species that represent major groups of spiny-rayed fishes {% cite Betancur2017 --file CladeAge-Tutorial/master-refs.bib %} as well as the most divergent groups of cichlid fishes. The 24 species that we will focus on are listed in the table below.

<center>

| ID                       | Species                   | Group                |
|--------------------------|---------------------------|----------------------|
| Oreochromis_niloticus    | _Oreochromis niloticus_   | African cichlids     |
| Heterochromis_multidensA | _Heterochromis multidens_ | African cichlids     |
| Cichla_temensisA         | _Cichla temensis_         | Neotropical cichlids |
| Heros_appendictulatusA   | _Heros appendictulatus_   | Neotropical cichlids |
| Etroplus_maculatusA      | _Etroplus maculatus_      | Indian cichlids      |
| Oryzias_latipes          | _Oryzias latipes_         | Atherinomorphae      |
| Trachinotus_carolinusA   | _Trachinotus carolinus_   | Carangaria           |
| Channa_striataA          | _Channa striata_          | Anabantiformes       |
| Monopterus_albusA        | _Monopterus albus_        | Synbranchiformes     |
| Gasterosteus_acuC        | _Gasterosteus aculeatus_  | Eupercaria           |
| Astrapogon_stellatusA    | _Astrapogon tellatus_     | Gobiaria             |
| Aulostomus_chinensisA    | _Aulostomus chinensis_    | Syngnatharia         |
| Thunnus_albacaresA       | _Thunnus albacares_       | Pelagaria            |
| Porichthys_notatusA      | _Porichthys notatus_      | Batrachoidiaria      |
| Diplacanthopoma_brunneaA | _Diplacanthopoma brunnea_ | Ophidiaria           |
| Sargocentron_cornutumA   | _Sargocentron cornutum_   | Holocentrimorphaceae |
| Rondeletia_loricataA     | _Rondeletia loricata_     | Beryciformes         |
| Monocentris_japonicaA    | _Monocentris japonica_    | Trachichthyiformes   |
| Polymixia_japonicaA      | _Polymixia japonica_      | Polymixiipterygii    |
| Regalecus_Glesne         | _Regalecus glesne_        | Lampripterygii       |
| Percopsis_omiscomaycusA  | _Percopsis omiscomaycus_  | Percopsaria          |
| Zenopsis_conchiferaB     | _Zenopsis conchifera_     | Zeiariae             |
| Stylephorus_chordatusB   | _Stylephorus chordatus_   | Stylephoriformes     |
| Gadus_morhua             | _Gadus morhua_            | Gadiformes           |

</center>

A list of the 24 species IDs in plain text format is also in file `Near_et_al_ids.txt`. On the command line, we can use that file to extract the sequences of these species from the full alignment, and write them to a new file in Nexus format named `Near_et_al_red.nex`.

> Execute the following commands on the command line:

```
head -n 7 Near_et_al.nex | sed 's/ntax=608/ntax=24/g'> Near_et_al_red.nex
grep -f Near_et_al_ids.txt -e "\[" Near_et_al.nex | sed -e $'s/\[/\\\n\[/g' >> Near_et_al_red.nex
tail -n 35 Near_et_al.nex | sed 's/paup/assumptions/g' >> Near_et_al_red.nex
```

If you should not be able to execute these commands on the command line, you could instead download the reduced alignment file `Near_et_al_red.nex` using the link in the left-hand column under "Data".


### Installing the CladeAge package

To use fossil constraints as calibrations points according to the CladeAge model, we'll first have to install the CladeAge add-on package for BEAST2.

> To do so, open BEAUti, and click on "Manage Packages" in the "File" menu.
 
This will open a window for the BEAST2 Package Manager. In this window, select "CA" and click "Install/Upgrade" as shown in the screenshot below.

<figure>
	<a id="fig:beauti1"></a>
	<img style="width:80%;" src="figures/beauti1.png" alt="BEAUti">
	<figcaption>Figure 1: Install the CladeAge package.</figcaption>
</figure>

> Close and reopen BEAUti.

You should then see that an additional tab has been added named "Clade Ages", as in the screenshot below.

<figure>
	<a id="fig:beauti2"></a>
	<img style="width:80%;" src="figures/beauti2.png" alt="BEAUti">
	<figcaption>Figure 2: The BEAUti interface after installing the CladeAge package.</figcaption>
</figure>


### Generating the analysis file with BEAUti

> Click on "Import Alignment" in BEAUti's "File" menu, and select the alignment file `Near_et_al_red.nex`.

BEAUti should then recognize 30 different partitions, one for each codon position of each of the ten markers. The BEAUti window should then look as shown in the screenshot below.

<figure>
	<a id="fig:beauti3"></a>
	<img style="width:80%;" src="figures/beauti3.png" alt="BEAUti">
	<figcaption>Figure 3: The BEAUti interface after importing the alignment file.</figcaption>
</figure>

> Select all partitions, and click "Link Trees" as well as "Link Clock Models", as shown below.

<figure>
	<a id="fig:beauti4"></a>
	<img style="width:80%;" src="figures/beauti4.png" alt="BEAUti">
	<figcaption>Figure 4: Linking trees and clock models.</figcaption>
</figure>

> Move on to the "Site Model" tab to select the site model for all partitions.

Instead of selecting a model such as HKY or GTR, I highly recommend the use of the model averaging implemented in the bModelTest package {% cite Bouckaert2017 --file CladeAge-Tutorial/master-refs.bib %}. If you did not already install this package, you can do so with the BEAST2 Package Manager from BEAUti as described above for the installation of the CladeAge package (don't forget to close and reopen BEAUti after installation to see changes to the interface). More information on model averaging with the bModelTest package is provided in the [Substitution Model Averaging](https://taming-the-beast.org/tutorials/Substitution-model-averaging/) tutorial. While recommended, model averaging with bModelTest is not required for this CladeAge tutorial, and you could also pick a model such as HKY or GTR instead. The description given here, however, assumes that the bModelTest package has been installed.

> Select "BEAST Model Test" from the drop-down menu at the top of the window, as shown below.

<figure>
	<a id="fig:beauti5"></a>
	<img style="width:80%;" src="figures/beauti5.png" alt="BEAUti">
	<figcaption>Figure 5: Selecting model averaging with the bModelTest model.</figcaption>
</figure>

> Select "namedSelected" from the drop-down menu that at first had "transitionTransversionSplit" selected.
> Leave the checkbox next to "Empirical" unticked to allow estimation of nucleotide frequencies.
> Then, set the tick to the right of "Mutation Rate" to specify that this rate should be estimated.

The window should then look as in the next screenshot.

<figure>
	<a id="fig:beauti6"></a>
	<img style="width:80%;" src="figures/beauti6.png" alt="BEAUti">
	<figcaption>Figure 6: Selecting settings for model averaging with the bModelTest model.</figcaption>
</figure>

> Select all partitions in the list at the left of the window, and click "OK" to clone the site model from the first partition to all other partitions.

<figure>
	<a id="fig:beauti7"></a>
	<img style="width:80%;" src="figures/beauti7.png" alt="BEAUti">
	<figcaption>Figure 7: Cloning the site model.</figcaption>
</figure>

> Continue to the "Clock Model" tab, and select the "Relaxed Clock Log Normal" clock model from the drop-down menu.

<figure>
	<a id="fig:beauti8"></a>
	<img style="width:80%;" src="figures/beauti8.png" alt="BEAUti">
	<figcaption>Figure 8: Selecting the clock model.</figcaption>
</figure>

Because we are going to calibrate the molecular clock with fossil constraints, we should allow the clock rate to be estimated. By default, however, the checkbox that we would need to tick, next to "estimate" at the bottom right, is disabled.

> To enable this checkbox, click on "Automatic set clock rate" in BEAUti's "Mode" menu.

<figure>
	<a id="fig:beauti9"></a>
	<img style="width:80%;" src="figures/beauti9.png" alt="BEAUti">
	<figcaption>Figure 9: Enabling the estimation of the clock rate.</figcaption>
</figure>

The checkbox for clock-rate estimation should then be enabled.

> Set a tick in the checkbox at the bottom right to estimate the clock rate.

The "Clock Model" tab should then look as shown in the screenshot below.

<figure>
	<a id="fig:beauti10"></a>
	<img style="width:80%;" src="figures/beauti10.png" alt="BEAUti">
	<figcaption>Figure 10: Selecting estimation of the clock rate.</figcaption>
</figure>

> Move on to the "Priors" tab. Select the "Birth Death Model" as the tree prior, from the drop-down menu at the very top of the window.

<figure>
	<a id="fig:beauti11"></a>
	<img style="width:80%;" src="figures/beauti11.png" alt="BEAUti">
	<figcaption>Figure 11: Selecting the tree prior.</figcaption>
</figure>

Instead of specifying constraints on monophyly and divergence times in the "Priors" tab, this is done in the separate tab named "Clade Ages" when the CladeAge package is installed.

> Open the "Clade Ages" tab and click on the "+ Add Prior" button

<figure>
	<a id="fig:beauti12"></a>
	<img style="width:80%;" src="figures/beauti12.png" alt="BEAUti">
	<figcaption>Figure 12: Selecting a first taxon set for calibration.</figcaption>
</figure>

We will now have to specify a rather long list of constraints to make the best possible use of the information provided by the fossil record and to obtain divergence time estimates that are as reliable as possible given our dataset. We'll start simple by specifying that the origin of African cichlid tribe *Heterochromini*, represented in our dataset by *Heterochromis multidens*, must have occurred at least 15.97-33.9 million years ago (Ma), since this is the age of the oldest fossil of Heterochromi, which was reported from the Baid Formation of Saudi Arabia by Lippitsch and Micklich {% cite Lippitsch1998 -A --file CladeAge-Tutorial/master-refs.bib %}.

> Specify "Heterochromini" in the field next to "Taxon set label", and select only "Heterochromis\_multidensA" as the ingroup of this taxon set.

<figure>
	<a id="fig:beauti13"></a>
	<img style="width:80%;" src="figures/beauti13.png" alt="BEAUti">
	<figcaption>Figure 13: Selecting a first taxon set for Heterochromini.</figcaption>
</figure>

After you click "OK", you'll see the window shown below, in which you can specify minimum and maximum values for the parameters net diversification rate, turnover rate, and sampling rate. Based on these values as well as the fossil age, the CladeAge model is going to automatically determine the optimal shape of the prior density used for each constraint.

<figure>
	<a id="fig:beauti14"></a>
	<img style="width:80%;" src="figures/beauti14.png" alt="BEAUti">
	<figcaption>Figure 14: Specifying parameters for CladeAge calibrations.</figcaption>
</figure>

While it might sometimes be preferable to estimate the net diversification (speciation minus extinction) and turnover (extinction divided by speciation) rates as part of the analysis, this is not implemented in the CladeAge model yet. Instead, the values for the two parameters have to be specified _a priori_, together with the "sampling rate", the rate at which fossils that are eventually sampled by scientists were once deposited in the fossil record. Fortunately, estimates for these three rates can be found in the literature, and the confidence intervals for these estimates can be accounted for by the CladeAge model. Here, as in Matschiner et al. {% cite Matschiner2017 -A --file CladeAge-Tutorial/master-refs.bib %}, we adopt teleost-specific estimates for net diversification and turnover from Santini et al. {% cite Santini2009 -A --file CladeAge-Tutorial/master-refs.bib %}, and a sampling-rate estimate for bony fishes from Foote and Miller {% cite Foote2007 -A --file CladeAge-Tutorial/master-refs.bib %}. These estimates are 0.041-0.081 for the net diversification rate, 0.0011-0.37 for the turnover, and 0.0066-0.01806 for the sampling rate. 

> Specify 0.041-0.081 for the net diversification rate, 0.0011-0.37 for the turnover, and 0.0066-0.01806 for the sampling rate. Don't change the selection of "empirical CladeAge" in the drop-down menu to the left.

<figure>
	<a id="fig:beauti15"></a>
	<img style="width:80%;" src="figures/beauti15.png" alt="BEAUti">
	<figcaption>Figure 15: Specifying parameters for CladeAge calibrations.</figcaption>
</figure>

While these three rate estimates are going to apply to all fossil constraints, information specific to the fossil constraint still must be added.

> To add more information specific to the Heterochromini constraint, click on the triangle to the left of "Heterochromini".

You'll see four more fields in which you can specify minimum and maximum values for the "First occurrence age" (the age of the oldest fossil of the clade) and a "Sampling gap". This sampling gap is supposed to represent a period of time right after the origin of the clade during which one might assume that either fossilization is unlikely (for example because the clade's population might still be too small) or that fossils would not be recognized as belonging to that clade (because it might not have evolved diagnostic characters yet). We are going to ignore the possibility of a sampling gap here.

> Specify the known age of the oldest fossil of Heterochromini in million of years, 15.97-33.9.

<figure>
	<a id="fig:beauti16"></a>
	<img style="width:80%;" src="figures/beauti16.png" alt="BEAUti">
	<figcaption>Figure 16: Specifying further parameters for the Heterochromini calibration.</figcaption>
</figure>

It is also possible to preview the shape of the prior densities as calculated by CladeAge based on the specified parameters.

> Click on the "Preview" button below the CladeAge icon on the right (you may have to increase the window size to see it).

A plot outlining the prior densities should then appear as in the screenshot below.

<figure>
	<a id="fig:beauti17"></a>
	<img style="width:80%;" src="figures/beauti17.png" alt="BEAUti">
	<figcaption>Figure 17: Previewing the prior density for the Heterochromini calibration.</figcaption>
</figure>

From this plot, you can see that under the assumption that all specified model parameters are correct, there's a good probability that Heterochromini originated some time between 25 and 80 Ma. While this range is rather wide, it is based on only one fossil; we will obtain more precise estimates when we run the phylogenetic analysis with multiple fossil constraints.

> Click the triangle to the left of "Heterochromini" again to close the section with details on this fossil constraint.
>
> Add further fossil constraints for each of the clades listed below (see the Supplementary Material of Matschiner et al. {% cite Matschiner2017 -A --file CladeAge-Tutorial/master-refs.bib %} if interested in details and references):

* **"Other African cichlid tribes"**<br>Ingroup: *Oreochromis niloticus* ("Oreochromis\_niloticus")<br>Oldest fossil species: *Mahengechromis* spp.<br>First occurrence age: 45.0-46.0 Ma
* **"African cichlids"**<br>Ingroup: *Heterochromis multidens* ("Heterochromis\_multidensA"), *Oreochromis niloticus* ("Oreochromis\_niloticus")<br>Oldest fossil species: *Mahengechromis* spp.<br>First occurrence age: 45.0-46.0 Ma
* **"Retroculini and Cichlini"**<br>Ingroup: *Cichla temensis* ("Cichla\_temensisA")<br>Oldest fossil species: *Palaeocichla longirostrum*<br>First occurrence age: 5.332-23.03 Ma
* **"Other Neotropical cichlid tribes"**<br>Ingroup: *Heros appendictulatus* ("Heros\_appendictulatusA")<br>Oldest fossil species: *Plesioheros chauliodus*<br>First occurrence age: 39.9-45.0 Ma
* **"Neotropical cichlids"**<br>Ingroup: *Cichla temensis* ("Cichla\_temensisA"), *Heros appendictulatus* ("Heros\_appendictulatusA")<br>Oldest fossil species: *Plesioheros chauliodus*<br>First occurrence age: 39.9-45.0 Ma
* **"Afro-American cichlids"**<br>Ingroup: *Cichla temensis* ("Cichla\_temensisA"), *Heros appendictulatus* ("Heros\_appendictulatusA"), *Heterochromis multidens* ("Heterochromis\_multidensA"), *Oreochromis niloticus* ("Oreochromis\_niloticus")<br>Oldest fossil species: *Mahengechromis* spp.<br>First occurrence age: 45.0-46.0 Ma
* **"Cichlids"**<br>*Cichla temensis* ("Cichla\_temensisA"), *Etroplus maculatus* ("Etroplus\_maculatusA"), *Heros appendictulatus* ("Heros\_appendictulatusA"), *Heterochromis multidens* ("Heterochromis\_multidensA"), *Oreochromis niloticus* ("Oreochromis\_niloticus")<br>Oldest fossil species: *Mahengechromis* spp.<br>First occurrence age: 45.0-46.0 Ma
* **"Atherinomorphae"**<br>Ingroup: *Oryzias latipes* ("Oryzias\_latipes")<br>Oldest fossil species: *Rhamphexocoetus volans*<br>First occurrence age: 49.1-49.4 Ma
* **"Ovalentaria"**<br>Ingroup: *Cichla temensis* ("Cichla\_temensisA"), *Etroplus maculatus* ("Etroplus\_maculatusA"), *Heros appendictulatus* ("Heros\_appendictulatusA"), *Heterochromis multidens* ("Heterochromis\_multidensA"), *Oreochromis niloticus* ("Oreochromis\_niloticus"), *Oryzias latipes* ("Oryzias\_latipes")<br>Oldest fossil species: *Rhamphexocoetus volans*<br>First occurrence age: 49.1-49.4 Ma
* **"Carangaria"**<br>Ingroup: *Trachinotus carolinus* ("Trachinotus\_carolinusA")<br>Oldest fossil species: *Trachicaranx tersus*<br>First occurrence age: 55.8-57.23 Ma
* **"Anabantiformes"**<br>Ingroup: *Channa striata* ("Channa\_striataA")<br>Oldest fossil species: *Osphronemus goramy*<br>First occurrence age: 45.5-50.7 Ma
* **"Anabantaria"**<br>Ingroup: *Channa striata* ("Channa\_striataA"), *Monopterus albus* ("Monopterus\_albusA")<br>Oldest fossil species: *Osphronemus goramy*<br>First occurrence age: 45.5-50.7 Ma
* **"Eupercaria"**<br>Ingroup: *Gasterosteus aculeatus* ("Gasterosteus_acuC")<br>Oldest fossil species: *Cretatriacanthus guidottii*<br>First occurrence age: 83.5-99.6 Ma
* **"Gobiaria"**<br>Ingroup: *Astrapogon tellatus* ("Astrapogon\_stellatusA")<br>Oldest fossil species: *"Gobius" gracilis*<br>First occurrence age: 30.7-33.9 Ma
* **"Syngnatharia"**<br>Ingroup: *Aulostomus chinensis* ("Aulostomus\_chinensisA")<br>Oldest fossil species: *Prosolenostomus lessinii*<br>First occurrence age: 49.1-49.4 Ma
* **"Pelagaria"**<br>Ingroup: *Thunnus albacares* ("Thunnus\_albacaresA")<br>Oldest fossil species: *Eutrichiurides opiensis*<br>First occurrence age: 56.6-66.043 Ma
* **"Batrachoidiaria"**<br>Ingroup: *Porichthys notatus* ("Porichthys\_notatusA")<br>Oldest fossil species: *Louckaichthys novosadi*<br>First occurrence age: 27.82-33.9 Ma
* **"Ophidiaria"**<br>Ingroup: *Diplacanthopoma brunnea* ("Diplacanthopoma\_brunneaA")<br>Oldest fossil species: *Eolamprogrammus senectus*<br>First occurrence age: 55.8-57.23 Ma
* **"Percomorphaceae"**<br>Ingroup: *Astrapogon tellatus* ("Astrapogon\_stellatusA"), *Aulostomus chinensis* ("Aulostomus\_chinensisA"), *Channa striata* ("Channa\_striataA"), *Cichla temensis* ("Cichla\_temensisA"), *Diplacanthopoma brunnea* ("Diplacanthopoma\_brunneaA"), *Etroplus maculatus* ("Etroplus\_maculatusA"), *Gasterosteus aculeatus* ("Gasterosteus_acuC"), *Heros appendictulatus* ("Heros\_appendictulatusA"), *Heterochromis multidens* ("Heterochromis\_multidensA"), *Monopterus albus* ("Monopterus\_albusA"), *Oreochromis niloticus* ("Oreochromis\_niloticus"), *Oryzias latipes* ("Oryzias\_latipes"), *Porichthys notatus* ("Porichthys\_notatusA"), *Thunnus albacares* ("Thunnus\_albacaresA"), *Trachinotus carolinus* ("Trachinotus\_carolinusA")<br>Oldest fossil species: *Cretatriacanthus guidottii*<br>First occurrence age: 83.5-99.6 Ma
* **"Holocentrimorphaceae"**<br>Ingroup: *Sargocentron cornutum* (Sargocentron\_cornutumA)<br>Oldest fossil species: *Caproberyx pharsus*<br>First occurrence age: 97.8-99.1 Ma
* **"Acanthopterygii"**<br>Ingroup: *Astrapogon tellatus* ("Astrapogon\_stellatusA"), *Aulostomus chinensis* ("Aulostomus\_chinensisA"), *Channa striata* ("Channa\_striataA"), *Cichla temensis* ("Cichla\_temensisA"), *Diplacanthopoma brunnea* ("Diplacanthopoma\_brunneaA"), *Etroplus maculatus* ("Etroplus\_maculatusA"), *Gasterosteus aculeatus* ("Gasterosteus_acuC"), *Heros appendictulatus* ("Heros\_appendictulatusA"), *Heterochromis multidens* ("Heterochromis\_multidensA"), *Monocentris japonica* ("Monocentris\_japonicaA"), *Monopterus albus* ("Monopterus\_albusA"), *Oreochromis niloticus* ("Oreochromis\_niloticus"), *Oryzias latipes* ("Oryzias\_latipes"), *Porichthys notatus* ("Porichthys\_notatusA"), *Rondeletia loricata* ("Rondeletia\_loricataA"), *Thunnus albacares* ("Thunnus\_albacaresA"), *Trachinotus carolinus* ("Trachinotus\_carolinusA"), *Sargocentron cornutum* (Sargocentron\_cornutumA)<br>Oldest fossil species: *Caproberyx pharsus*<br>First occurrence age: 97.8-99.1 Ma
* **"Polymixiipterygii"**<br>Ingroup: *Polymixia japonica* ("Polymixia\_japonicaA")<br>Oldest fossil species: *Homonotichthys rotundus*<br>First occurrence age: 93.5-96.0 Ma
* **"Percopsaria"**<br>Ingroup: *Percopsis omiscomaycus* ("Percopsis\_omiscomaycusA")<br>Oldest fossil species: *Mcconichthys longipinnis*<br>First occurrence age: 61.1-66.043 Ma
* **"Zeiariae"**<br>Ingroup: *Zenopsis conchifera* ("Zenopsis\_conchiferaB")<br>Oldest fossil species: *Cretazeus rinaldii*<br>First occurrence age: 69.2-76.4 Ma
* **"Gadiformes"**<br>Ingroup: *Gadus morhua* ("Gadus\_morhua")<br>Oldest fossil species: *Protacodus* sp.<br>First occurrence age: 59.7-62.8 Ma
* **"Paracanthopterygii"**<br>Ingroup: *Gadus morhua* ("Gadus\_morhua"), *Percopsis omiscomaycus* ("Percopsis\_omiscomaycusA"), *Stylephorus chordatus* ("Stylephorus\_chordatusB"), *Zenopsis conchifera* ("Zenopsis\_conchiferaB")<br>Oldest fossil species: *Cretazeus rinaldii*<br>First occurrence age: 69.2-76.4 Ma

Once all these constraints are added, the BEAUti window should look as shown below.

<figure>
	<a id="fig:beauti18"></a>
	<img style="width:80%;" src="figures/beauti18.png" alt="BEAUti">
	<figcaption>Figure 18: Adding all calibrations.</figcaption>
</figure>

> Move on to the "MCMC" tab. Specify again an MCMC chain length of 100 million iterations, name the log file `Near_et_al_red.log` and the tree file `Near_et_al_red.trees`.

The BEAUti window should then look as in the next screenshot.

<figure>
	<a id="fig:beauti19"></a>
	<img style="width:80%;" src="figures/beauti19.png" alt="BEAUti">
	<figcaption>Figure 19: Specifying chain length and output.</figcaption>
</figure>

> Save the analysis settings to a new file named `Near_et_al_red.xml` by clicking "Save As" in BEAUti's "File" menu.

Open BEAST2, load file `Near_et_al_red.xml`, and try running the MCMC analysis.

<figure>
	<a id="fig:beast1"></a>
	<img style="width:50%;" src="figures/beast1.png" alt="BEAST">
	<figcaption>Figure 20: Starting the BEAST analysis.</figcaption>
</figure>

Most likely, the MCMC analysis is going to crash right at the start with an error message as shown below.

<figure>
	<a id="fig:beast2"></a>
	<img style="width:50%;" src="figures/beast2.png" alt="BEAST">
	<figcaption>Figure 21: Error message due to failure to initialise the MCMC chain.</figcaption>
</figure>

This is a common problem when several fossil constraints are specified: According to the error message, BEAST2 could not find a proper state to initialise. This means that even after several attempts, no starting state of the MCMC chain could be found that had a non-zero probability. Most often, the issue is that the tree that BEAST2 randomly generates to start the chain is in conflict with one or more fossil constraints. Unfortunately, the only way to fix this issue is to manually edit the XML file and specify a starting tree that is in agreement with the specified fossil constraints. In particular, because all fossil constraints imposed hard minimum ages on the origin of the respective clades, this clades must at least be as old as this minimum age in the starting tree. In case of doubt, it is usually safer to make the starting tree too old rather than too young, the course of the MCMC chain should, at least after the burnin, not be influenced by the starting state anymore anyway. Some helpful advice on how to specify starting trees is provided on the [BEAST2](https://www.beast2.org/fix-starting-tree/) webpage. With trees of hundreds of taxa, generating a suitable starting tree can be a tricky task in itself, but with the small number of 24 species used here, writing a starting tree by hand is feasible.

> Copy and paste the below starting tree string into a new FigTree window.

```
((((((((((((((Oreochromis_niloticus:50,Heterochromis_multidensA:50):10,(Cichla_temensisA:50,Heros_appendictulatusA:50):10):10,Etroplus_maculatusA:70):30,Oryzias_latipes:100):10,(Trachinotus_carolinusA:70,(Channa_striataA:60,Monopterus_albusA:60):10):40):10,Gasterosteus_acuC:120):10,Astrapogon_stellatusA:130):10,(Aulostomus_chinensisA:80,Thunnus_albacaresA:80):60):10,Porichthys_notatusA:150):10,Diplacanthopoma_brunneaA:160):10,Sargocentron_cornutumA:170):10,(Rondeletia_loricataA:100,Monocentris_japonicaA:100):80):10,Polymixia_japonicaA:190):10,((((Gadus_morhua:70,Stylephorus_chordatusB:70):10,Zenopsis_conchiferaB:80):10,Percopsis_omiscomaycusA:90):10,Regalecus_Glesne:100):100)
```

As you'll see, I just arbitrarily specified for most branches a length of 10 million years, and I made sure that particularly the more recent divergence events agree with the respective fossil constraints (e.g. by placing the divergence of *Oreochromis niloticus* and *Heterochromis multidens* at 50 Ma because the origin of the clade "Other African cichlid tribes", represented by *Oreochromis niloticus* is constrained to be at least 45 Ma).

> Open file `Near_et_al_red.xml` in a text editor and find the block shown below on lines 335-340.

```
		<init id="RandomTree.t:ZIC_2nd" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:enc_1st">
			<taxa id="ZIC_2nd" spec="FilteredAlignment" data="@Near_et_al_red" filter="7697-8576\3"/>
			<populationModel id="ConstantPopulation0.t:ZIC_2nd" spec="ConstantPopulation">
				<parameter id="randomPopSize.t:ZIC_2nd" name="popSize">1.0</parameter>
			</populationModel>
		</init>
```

-------

# Tutorial style guide

## Text styling

This is how to write _italic text_.

This is how to write **bold text**.

This is how to write **_bold and italic text_**.

Do text superscripts like this 7^th, x^2y or  x^(2y + 3z).


## Lists

### Unnumbered lists

- Lorem ipsum dolor sit amet, consectetur adipiscing elit.
- Integer pharetra arcu ut nisl mollis ultricies.
	- Fusce nec tortor at enim cursus dictum.
	- Phasellus nec urna quis velit eleifend convallis sodales nec augue.
- In iaculis turpis in massa facilisis, quis ultricies nibh ultricies.
- Nam vitae turpis eu lacus imperdiet mollis id at augue.
- Sed sed turpis ac dolor mollis accumsan.


### Numbered lists

1. Lorem ipsum dolor sit amet, consectetur adipiscing elit.
2. Integer pharetra arcu ut nisl mollis ultricies.
	1. Fusce nec tortor at enim cursus dictum.
	2. Phasellus nec urna quis velit eleifend convallis sodales nec augue.
1. In iaculis turpis in massa facilisis, quis ultricies nibh ultricies.
1. Nam vitae turpis eu lacus imperdiet mollis id at augue.
1. Sed sed turpis ac dolor mollis accumsan.

### Mixed lists

1. Lorem ipsum dolor sit amet, consectetur adipiscing elit.
2. Integer pharetra arcu ut nisl mollis ultricies.
	* Fusce nec tortor at enim cursus dictum.
	* Phasellus nec urna quis velit eleifend convallis sodales nec augue.
1. In iaculis turpis in massa facilisis, quis ultricies nibh ultricies.
1. Nam vitae turpis eu lacus imperdiet mollis id at augue.
1. Sed sed turpis ac dolor mollis accumsan.


## Figures


<figure>
	<a id="fig:example1"></a>
	<img style="width:25%;" src="figures/Logo_bw.png" alt="">
	<figcaption>Figure 1: This figure is 25% of the page width.</figcaption>
</figure>


<figure>
	<a id="fig:example2"></a>
	<img style="width:10%;" src="figures/Logo_bw.png" alt="">
	<figcaption>Figure 2: This figure is only 10% of the page width.</figcaption>
</figure>



# Code

A bit of inline monospaced font can be made `like this`. Larger code blocks can be made by using the code environment:

Java:

```java
public class HelloWorld {

    public static void main(String[] args) {
        // Prints "Hello, World" to the terminal window.
        System.out.println("Hello, World");
    }

}
```

XML:

```xml
	<BirthDeathSkylineModel spec="BirthDeathSkylineModel" id="birthDeath" tree="@tree" contemp="true">
	      <parameter name="origin" id="origin" value ="100" lower="0."/>    
	      <parameter name="R0" id="R0" value="2" lower="0." dimension ="10"/>
	      <parameter name="becomeUninfectiousRate" id="becomeUninfectiousRate" value="1" lower="0." dimension ="10"/>
	      <parameter name="samplingProportion" id="samplingProportion" value="0."/>
	      <parameter name="rho" id="rho" value="1e-6" lower="0." upper="1."/>
	</BirthDeathSkylineModel>
```

R:

```R
	> myString <- "Hello, World!"
	> print (myString)
	[1] "Hello, World!"
```

# Equations

Inline equations: {% eqinline \dot{x} = \sigma(y-x) %}

Displayed equations: 
{% eq \left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right) %}



## Instruction boxes

Use block-quotes for step-by-step instruction that the user should perform (this will produce a framed box on the website):

> The data we have is not the data we want, and the data we need is not the data we have.
> 
> We can input **any** formatted text in here:
>
> - Even
> - Lists
>
> or equations:
>
> {% eq (x_1, \ldots, x_n) \left( \begin{array}{ccc}
      \phi(e_1, e_1) & \cdots & \phi(e_1, e_n) \\
      \vdots & \ddots & \vdots \\
      \phi(e_n, e_1) & \cdots & \phi(e_n, e_n)
    \end{array} \right)
  \left( \begin{array}{c}
      y_1 \\
      \vdots \\
      y_n
    \end{array} \right) %}






# Hyperlinks

Add links to figures like this: 

- [Figure 1](#fig:example1) is 25% of the page width.
- [Figure 2](#fig:example2) is 10% of the page width. 

Add links to external URLs like [this](http://www.google.com). 

Links to equations or different sections within the same document are a little buggy.


----

# Useful Links

- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file CladeAge-Tutorial/master-refs.bib %}
- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- BEAST 1 website and documentation: [http://beast.bio.ed.ac.uk](http://beast.bio.ed.ac.uk)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users) 

----

# Relevant References

{% bibliography --cited --file CladeAge-Tutorial/master-refs.bib %}

