# Pan-cancer quantification of neoantigen-mediated immunoediting in cancer evolution

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![version](https://img.shields.io/badge/version-dev-green.svg)](https://shields.io/)

<details>
<summary>Table of content</summary>

## Table of content
   * [Overview](#Overview)
   * [Contents](#Contents)
   * [Citation](#Citation)
   * [Acknowledgement](#Acknowledgement)
   * [LICENSE](#License)

</details>

----

## Overview

This repository provides the analysis reports, code and data for readers who are interest in this project and make it easier to reproduce the whole analysis procedure.

## Contents

* [code](./code) 
  * [R](./code/R/) : R functions and scripts that do not include in the analysis report.
  * [Julia](./code/julia) : Julia scripts for simulating tumor growth, modified from Eszter [Lakatos et.al](https://www.nature.com/articles/s41588-020-0687-1)
  * [Python](./code/python) : Snakemake scripts to do HLA typing.
  * [Shell](./code/shell) : Shell scripts to do mutation calling, RNA-seq quantification, Neoantigen prediction.
* [data](./data) The data used and produced by analysis report.
* [docs](./docs) Website pages and figures used for showing analysis reports.
* [report](./report) Rmarkdown files of analysis report and related html web page files.

## Citation

Wu T, Wang G, Wang X, et al. Quantification of neoantigen-mediated immunoediting in cancer evolution[J]. Cancer research, 2022: canres. 3717.2021.

## Acknowledgement

We thank ShanghaiTech University High Performance Computing Public Service Platform for computing services.This work was supported by Shanghai Science and Technology Commission (21ZR1442400), the National Natural Science Foundation of China (31771373), and startup funding from ShanghaiTech University.

## License

***

**Cancer Biology Group @ShanghaiTech**

**Research group led by Xue-Song Liu in ShanghaiTech University**