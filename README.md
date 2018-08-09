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

As a first step, download the molecular dataset of Near et al. {% cite Near2013 -A --file CladeAge-Tutorial/master-refs.bib %} from the Dryad data repository connected to their publication. On the command-line of a Mac or Linux computer, this could be done using the `wget` utility:

```sh
wget https://datadryad.org/bitstream/handle/10255/dryad.50839/Near_et_al.nex
```

Alternatively, you could use your browser to open the Dryad repository [https://datadryad.org/resource/doi:10.5061/dryad.d3mb4](https://datadryad.org/resource/doi:10.5061/dryad.d3mb4) and download the file `Near_et_al.nex`.

This file in Nexus format contains the sequence alignments for ten nuclear markers, sequenced for 608 species of spiny-rayed fishes. As this dataset is far too large to be analyzed in this tutorial, we'll extract the sequences of 24 species that represent major groups of spiny-rayed fishes {% cite Betancur2017 -A --file CladeAge-Tutorial/master-refs.bib %} as well as the most divergent groups of cichlid fishes. The 24 species that we will focus on are listed in the table below.

<center>

| ID                       | Species                   | Group                |
|--------------------------|---------------------------|----------------------|
| Oreochromis_niloticus    | *Oreochromis niloticus*   | African cichlids     |
| Heterochromis_multidensA | *Heterochromis multidens* | African cichlids     |
| Cichla_temensisA         | *Cichla temensis*         | Neotropical cichlids |
| Heros_appendictulatusA   | *Heros appendictulatus*   | Neotropical cichlids |
| Etroplus_maculatusA      | *Etroplus maculatus*      | Indian cichlids      |
| Oryzias_latipes          | *Oryzias latipes*         | Atherinomorphae      |
| Trachinotus_carolinusA   | *Trachinotus carolinus*   | Carangaria           |
| Channa_striataA          | *Channa striata*          | Anabantiformes       |
| Monopterus_albusA        | *Monopterus albus*        | Synbranchiformes     |
| Gasterosteus_acuC        | *Gasterosteus aculeatus*  | Eupercaria           |
| Astrapogon_stellatusA    | *Astrapogon tellatus*     | Gobiaria             |
| Aulostomus_chinensisA    | *Aulostomus chinensis*    | Syngnatharia         |
| Thunnus_albacaresA       | *Thunnus albacares*       | Pelagaria            |
| Porichthys_notatusA      | *Porichthys notatus*      | Batrachoidiaria      |
| Diplacanthopoma_brunneaA | *Diplacanthopoma brunnea* | Ophidiaria           |
| Sargocentron_cornutumA   | *Sargocentron cornutum*   | Holocentrimorphaceae |
| Rondeletia_loricataA     | *Rondeletia loricata*     | Beryciformes         |
| Monocentris_japonicaA    | *Monocentris japonica*    | Trachichthyiformes   |
| Polymixia_japonicaA      | *Polymixia japonica*      | Polymixiipterygii    |
| Regalecus_Glesne         | *Regalecus glesne*        | Lampripterygii       |
| Percopsis_omiscomaycusA  | *Percopsis omiscomaycus*  | Percopsaria          |
| Zenopsis_conchiferaB     | *Zenopsis conchifera*     | Zeiariae             |
| Stylephorus_chordatusB   | *Stylephorus chordatus*   | Stylephoriformes     |
| Gadus_morhua             | *Gadus morhua*            | Gadiformes           |

</center>

A list of the 24 species IDs in plain text format is also in file `Near_et_al_ids.txt`. On the command line, we can use that file to extract the sequences of these species from the full alignment, and write them to a new file in Nexus format named `Near_et_al_red.nex`:

```
head -n 7 Near_et_al.nex | sed 's/ntax=608/ntax=24/g'> Near_et_al_red.nex
grep -f Near_et_al_ids.txt -e "\[" Near_et_al.nex | sed -e $'s/\[/\\\n\[/g' >> Near_et_al_red.nex
tail -n 35 Near_et_al.nex | sed 's/paup/assumptions/g' >> Near_et_al_red.nex
```

If you should not be able to execute these commands on the command-line, you could instead download the reduced alignment file `Near_et_al_red.nex` using the link in the left-hand column under "Data".


### Installing the CladeAge package

To specify fossil constraints as calibrations points in BEAUti according to the CladeAge model of Matschiner et al. {% cite Matschiner2017 -A --file CladeAge-Tutorial/master-refs.bib %}, we'll first have to install the CladeAge add-on package for BEAST2. To do so, open BEAUti, and click on "Manage Packages" in the "File" menu. This will open a window for the BEAST2 Package Manager. In this window, select "CA" and click "Intstall/Upgrade" as shown in the screenshot below.

<figure>
	<a id="fig:beauti1"></a>
	<img style="width:80%;" src="figures/beauti1.png" alt= BEAUti"">
	<figcaption>Figure 1: Install the CladeAge package.</figcaption>
</figure>



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

