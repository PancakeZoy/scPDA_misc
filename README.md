# scPDA_misc
This Repo contains the code to reproduce the results and figures in manuscript of scPDA

Please first download the `data` and `results` folders from this figshare repository: [][], and move them to the root directory. The structure of the starting data is shown in the tree diagram below.

```
|-dsb_normalization.Rproj
|-data
```

## Folder Explanation

| Name | Content |
|-----------------|-------------|
| [data](data) | Four datasets (in `.rds` form) demonstrated in Results section of the manuscript|
| [code](code) | Code of applying protein counts denosing methods (`GMM`, `DSB`, `scAR`, `DecontPro`, `scPDA`) applied to each dataset in `data` folder. The corresponding results are saved in `results` folder|
| [results](results) | Denoised counts of each dataset resulted from each denoising method|
| [fig_reprod](fig_reprod) | Code of reproducing each figure in the manuscript|
| [scPDA](scPDA)| The developing version of `scPDA`|





