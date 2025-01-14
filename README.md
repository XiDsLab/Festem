# Festem

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Last
Commit](https://badgen.net/github/last-commit/XiDsLab/Festem/main)
![Commits Since
Latest](https://img.shields.io/github/commits-since/XiDsLab/Festem/latest/main)
![GitHub Downloads](https://img.shields.io/github/downloads/XiDsLab/Festem/total)
<!-- badges: end -->

Festem is a statistical method for the direct selection of cell-type markers for downstream clustering. Festem distinguishes marker genes with heterogeneous distribution across cells that are cluster informative. 

<img src="https://github.com/XiDsLab/Festem/blob/main/img/graphical_abstract.png?raw=true">

For analysis codes for paper "Directly selecting differentially expressed genes for single-cell clustering analyses", see https://github.com/XiDsLab/Festem_paper.

## Installation

To install Festem, start R and enter:

    if (!require("pak", quietly = TRUE))
        install.packages("pak")
    pak::pkg_install("XiDsLab/Festem")

## Usage

For detailed usage of Festem, please refer to our protocol at [STAR Protocol](https://doi.org/10.1016/j.xpro.2024.103514).

## Citations

- Chen, Z., Wang, C., Huang, S., Shi, Y., & Xi, R. (2024). Directly selecting cell-type marker genes for single-cell clustering analyses. Cell Reports Methods, 4(7). <https://doi.org/10.1016/j.crmeth.2024.100810>
- Wang, C., Chen, Z., & Xi, R. (2023). Feature screening for clustering analysis. arXiv preprint arXiv:2306.12671. <https://arxiv.org/abs/2306.12671>
