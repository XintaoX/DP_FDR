# DP_FDR

This repository contains code and supporting materials for the paper:

**Xintao Xia and Zhanrui Cai (2023).  
"Adaptive False Discovery Rate Control with Privacy Guarantee."  
Journal of Machine Learning Research, 24(252):1–35.**

The repository includes the original code used in the paper, data and figures from the numerical experiments, and correction materials associated with the erratum. The `correction/` folder contains the erratum, corrected code, and updated results after correction. The corrected results are almost identical to the original results.

## Overview

This repository is organized to support reproducibility for the paper **"Adaptive False Discovery Rate Control with Privacy Guarantee."** The paper studies adaptive false discovery rate control under differential privacy constraints.

The main purpose of this repository is to provide:

- The original code used for the paper.
- Data and figures associated with the original experiments.
- Code for the standard AdaPT procedure, based on Lihua Lei's AdaPT GitHub repository.
- The erratum of the paper.
- Corrected code and updated results after correction.

## Repository Structure

```text
DP_FDR/
├── AdaPT/
├── code/
├── correction/
├── LICENSE
└── README.md
```

## Folder Descriptions

### `AdaPT/`

This folder contains code for the standard AdaPT procedure, based on Lihua Lei's [GitHub repository](https://lihualei71.github.io/):

```text
https://github.com/lihualei71/AdaPT
```

AdaPT stands for **Adaptive P-value Thresholding**. It is a procedure for multiple testing with side information and serves as the standard non-private AdaPT reference implementation used for comparison.

### `code/`

This folder contains the original code used for the paper:

**"Adaptive False Discovery Rate Control with Privacy Guarantee."**

It includes the previous version of the code, together with the data and figures used in the original numerical experiments and empirical studies.

This folder is intended mainly for reproducing the original results reported in the paper.

### `correction/`

This folder contains the correction materials for the paper.

It includes:

- The erratum of the paper.
- The corrected code.
- The updated results after correction.

The results after correction are almost identical to the previous results and do not materially change the empirical conclusions of the paper.

For users who want to reproduce the corrected version of the experiments, this folder should be used as the primary reference.

## Paper

The paper is available from the *Journal of Machine Learning Research*:

```text
https://jmlr.org/papers/v24/23-0039.html
```

## Citation

If you use this repository, please cite the following paper:

```bibtex
@article{xia2023adaptive,
  title   = {Adaptive False Discovery Rate Control with Privacy Guarantee},
  author  = {Xia, Xintao and Cai, Zhanrui},
  journal = {Journal of Machine Learning Research},
  volume  = {24},
  number  = {252},
  pages   = {1--35},
  year    = {2023}
}
```

## Related Reference

The `AdaPT/` folder is based on the standard AdaPT implementation from Lihua Lei's GitHub repository:

```text
https://github.com/lihualei71/AdaPT
```

If you use the AdaPT code, please also cite the AdaPT paper:

```bibtex
@article{lei2018adapt,
  title   = {AdaPT: An Interactive Procedure for Multiple Testing with Side Information},
  author  = {Lei, Lihua and Fithian, William},
  journal = {Journal of the Royal Statistical Society: Series B},
  volume  = {80},
  number  = {4},
  pages   = {649--679},
  year    = {2018}
}
```

## License

This repository is released under the MIT License. See the `LICENSE` file for details.

## Notes

The original code and corrected code are both kept in this repository for transparency and reproducibility. Users interested in the final corrected version should refer to the materials in the `correction/` folder.
