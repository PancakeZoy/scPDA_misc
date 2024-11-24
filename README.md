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

### Downloading required folders
Please download the `data` and `results` folders from the figshare repository: [](), and place them in the root directory. Optionally, you can also generate the `results` folder by running scripts in the `code` folder, if the `data` folder is present. Note: Running these scripts may take significant time due to `DecontPro`. The structure of the starting data is shown in the tree diagram below:

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
│   └── Supplementary.pdf
└── results
    ├── dsb
    │   ├── dsb_DSB.rds
    │   ├── dsb_DecontPro.csv
    │   ├── dsb_EmptyProfile.csv
    │   ├── dsb_GMM.rds
    │   ├── dsb_GMM_mu1.csv
    │   ├── dsb_raw.csv
    │   ├── dsb_scAR.h5
    │   └── dsb_scPDA.h5
    ├── tea
    │   ├── tea_DecontPro.csv
    │   ├── tea_EmptyProfile.csv
    │   ├── tea_GMM_mu1.csv
    │   ├── tea_dsb.rds
    │   ├── tea_meta.csv
    │   ├── tea_raw.csv
    │   ├── tea_scAR.h5
    │   └── tea_scPDA.h5
    ├── titr188
    │   ├── titr188_DecontPro.csv
    │   ├── titr188_GMM_mu1.csv
    │   ├── titr188_dsb.csv
    │   ├── titr188_gmm.csv
    │   ├── titr188_meta.csv
    │   ├── titr188_raw.csv
    │   ├── titr188_scAR.csv
    │   ├── titr188_scAR.h5
    │   ├── titr188_scPDA.csv
    │   └── titr188_scPDA.h5
    └── wnn25
        ├── wnn25_DSB.rds
        ├── wnn25_DecontPro.csv
        ├── wnn25_GMM.rds
        ├── wnn25_GMM_mu1.csv
        ├── wnn25_meta.csv
        ├── wnn25_raw.csv
        ├── wnn25_scAR.h5
        └── wnn25_scPDA.h5

13 directories, 68 files
```

### To run the scripts in `code`
1. Create a conda environment:
```
conda create -n denoise
```

2. Install each denoising methods:
- `scPDA`
```
cd scPDA
pip install -e .
```

- `scAR`
```
conda install bioconda::scar
```

- `DSB`
```
install.packages('dsb')
```

- `DecontPro`
```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("decontX")
```
