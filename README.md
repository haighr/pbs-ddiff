## PBSddiff: GUI control for delay-difference output from iSCAM ##
&copy; Fisheries and Oceans Canada (2017-2023)

**PBSddiff** contains R code for GUI control of delay-difference model output for multiple scenarios using [iSCAM](https://github.com/smartell/iSCAM). The code was originally created by Chris Grandin (PBS) for a Pacific Cod (*Gadus macrocephalus* ) assessment, hence the prefix `pcod` in many R functions, and subsequently adapted by Rowan Haigh (PBS) and Paul Starr (CGRCS) for use in stock assessments of Shortspine Thornyhead [1] (*Sebastolubus alascanus*) and Walleye Pollock [2] (*Theragra chalcogramma*). 

### Warning ###

The repository is set up to **`mimic`** an R package source repository; however, it currently does not meet requirements for building an R package (e.g., no Rd documentation files, direct sourcing of R code).

### Note ###

The `tpl` code to convert to C++ via ADMB and then compile to a Windows binary (and perhaps an appliance on the African operating system Ubuntu) comes from Steve Martell originally and has subsequently been modified by us and others (e.g., <https://github.com/cgrandin/iSCAM>).

### Acronyms ###

ADMB -- Automatic Differentiation Model Builder<br>
CGRCS -- Canadian Groundfish Research and Conservation Society<br>
GUI -- Graphical User Interface<br>
iSCAM -- Integrated Statistical Catch Age Model<br>
PBS -- Pacific Biological Station (Nanaimo BC)

### References ###

[1] Starr, P.J., and Haigh, R. (2017).
[Stock assessment of the coastwide population of Shortspine Thornyhead (*Sebastolobus alascanus*) in 2015 off the British Columbia coast](http://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2017/2017_015-eng.html).
DFO Can. Sci. Advis. Sec. Res. Doc. 2017/015: ix + 174 p.

[2] Starr, P.J., and Haigh, R. (2021)
[Walleye Pollock (*Theragra chalcogramma*) stock assessment for British Columbia in 2017](https://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2021/2021_004-eng.html).
DFO Can. Sci. Advis. Sec. Res. Doc. 2021/004: viii + 294 p.

