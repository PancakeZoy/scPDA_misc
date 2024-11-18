# scPDA_misc
This Repo contains the code to reproduce the results and figures in manuscript of scPDA

## Folder Explanation

| Name | Content |
|-----------------|-------------|
| [data](data) | Four datasets (in the form of `.rds`) demonstrated in Results section of the manuscript|
| [code](code) | Code of applying protein counts denosing methods (`GMM`, `DSB`, `scAR`, `DecontPro`, `scPDA`) applied to each dataset in `data` folder. The corresponding results are saved in `results` folder|
| [results](results) | Denoised counts of each dataset resulted from each denoising method|
| [fig_reprod](fig_reprod) | Code of reproducing each figure in the manuscript|
| [scPDA](scPDA)| The developing version of `scPDA`|

## Downloading required folders
Please download the `data` and `results` folders from the figshare repository: [https://doi.org/10.6084/m9.figshare.27898152](https://doi.org/10.6084/m9.figshare.27898152), and place them in the root directory. Optionally, you can also generate the `results` folder by running scripts in the `code` folder (see below), if the `data` folder is present. Note: Running these scripts may take significant time due to `DecontPro`. 

**To generate the `results` folder (skip if you download it from figshare)**
1. Create a conda environment:
```
conda create -n denoise
```
2. Install each denoising methods:
```
# scPDA
cd scPDA
pip install -e .

# scAR
conda install bioconda::scar

# DSB
Rscript -e "install.packages('dsb')"

# DecontPro
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
Rscript -e "BiocManager::install('decontX')"

# GMM
Rscript -e "install.packages('mclust')"
```

3. Run the script of protein denoising for each dataset:
- `code/dsb/dsb.R`
- `code/wnn25/wnn25.R`
- `code/titr188/titr188.R`
- `code/tea/tea.R`


**The structure of the starting data is shown in the tree diagram below:**
```
├── LICENSE
├── README.md
├── code
│   ├── dsb
│   │   ├── dsb.R
│   │   ├── dsb_DecontPro.R
│   │   ├── dsb_scAR.py
│   │   └── dsb_scPDA.py
│   ├── tea
│   │   ├── tea.R
│   │   ├── tea_DecontPro.R
│   │   ├── tea_scAR.py
│   │   └── tea_scPDA.py
│   ├── titr188
│   │   ├── titr188.R
│   │   ├── titr188_scAR.py
│   │   └── titr188_scPDA.py
│   ├── tools.R
│   └── wnn25
│       ├── wnn25.R
│       ├── wnn25_DecontPro.R
│       ├── wnn25_scAR.py
│       └── wnn25_scPDA.py
├── data
│   ├── dsb.rds
│   ├── dsb_empty.rds
│   ├── supp
│   │   ├── 10k_neg_prot1.rds
│   │   ├── 10k_neg_prot2.rds
│   │   ├── 10k_neg_prot3.rds
│   │   ├── 10k_neg_prot4.rds
│   │   ├── 10k_pos_prot.rds
│   │   ├── 5`_neg_prot2.rds
│   │   ├── 5`_pos_prot.rds
│   │   ├── 5k_neg_prot2.rds
│   │   ├── 5k_pos_prot.rds
│   │   ├── dsb_cells.rds
│   │   ├── dsb_hash.rds
│   │   └── dsb_lib.rds
│   ├── teaseq.rds
│   ├── titr188.rds
│   └── wnn25.rds
├── fig_reprod
│   ├── Figure1_a-i.R
│   ├── Figure1_a-i.png
│   ├── Figure1_j-k.R
│   ├── Figure1_j-k.png
│   ├── Figure2.R
│   ├── Figure2.jpg
│   ├── Figure3.R
│   ├── Figure3.jpg
│   ├── Figure4.R
│   ├── Figure4.jpg
│   ├── Supplementary.pdf
│   └── supp
│       ├── Amb_Freq.R
│       ├── Hash_vs_Lib.jpeg
│       ├── Inconsistent_Definition.R
│       ├── Prop_10K.jpeg
│       ├── Prop_5K.jpeg
│       ├── Prop_5prime.jpeg
│       ├── Tr1_vs_Tr2.jpeg
│       ├── Tr2_vs_Tr3.jpeg
│       └── Tr3_vs_Tr4.jpeg
├── results
│   ├── dsb
│   │   ├── dsb_DSB.rds
│   │   ├── dsb_DecontPro.csv
│   │   ├── dsb_EmptyProfile.csv
│   │   ├── dsb_GMM.rds
│   │   ├── dsb_GMM_mu1.csv
│   │   ├── dsb_raw.csv
│   │   ├── dsb_scAR.h5
│   │   └── dsb_scPDA.h5
│   ├── tea
│   │   ├── tea_DecontPro.csv
│   │   ├── tea_EmptyProfile.csv
│   │   ├── tea_GMM_mu1.csv
│   │   ├── tea_dsb.rds
│   │   ├── tea_meta.csv
│   │   ├── tea_raw.csv
│   │   ├── tea_scAR.h5
│   │   └── tea_scPDA.h5
│   ├── titr188
│   │   ├── titr188_DecontPro.csv
│   │   ├── titr188_GMM_mu1.csv
│   │   ├── titr188_dsb.csv
│   │   ├── titr188_gmm.csv
│   │   ├── titr188_meta.csv
│   │   ├── titr188_raw.csv
│   │   ├── titr188_scAR.csv
│   │   ├── titr188_scAR.h5
│   │   ├── titr188_scPDA.csv
│   │   └── titr188_scPDA.h5
│   └── wnn25
│       ├── wnn25_DSB.rds
│       ├── wnn25_DecontPro.csv
│       ├── wnn25_GMM.rds
│       ├── wnn25_GMM_mu1.csv
│       ├── wnn25_meta.csv
│       ├── wnn25_raw.csv
│       ├── wnn25_scAR.h5
│       └── wnn25_scPDA.h5
└── scPDA
    ├── scPDA
    │   ├── __init__.py
    │   └── main
    │       ├── __init__.py
    │       ├── _loss.py
    │       ├── _network.py
    │       └── api.py
    └── setup.py

18 directories, 95 files
```

**To reproduce the figures in the manuscript**
Run the following scripts in the `fig_reprod`:
```
Figure1_a-i.R
Figure1_j-k.R
Figure2.R
Figure4.R
Figure3.R
```