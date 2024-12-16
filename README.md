# libra-run-template

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.INSERT-DOI.svg)](https://doi.org/10.5281/zenodo.INSERT-DOI)

This is a template repository for experimental runs of LIBRA.

This repository has the data for the run [**insert run name**].

## How to reproduce the results

### In Binder

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/LIBRA-project/INSERT-REPO-NAME/HEAD)

### Locally

1. Create a conda environment (requires conda):

```
conda env create -f environment.yml
```

2. Run the notebooks with the created environment `baby_1l_run_2`

## Todo list:
- [ ] [Link to Zenodo](https://zenodo.org/)
- [ ] Update Zenodo badge with new DOI
- [x] Change environment name in [`environment.yml`](environment.yml)
- [ ] Add general run data to [`data/general.json`](data/general.json)
- [ ] Add LSC data to [`data/tritium_detection`](data/tritium_detection)
- [ ] Add neutron detection data to [`data/neutron_detection`](data/neutron_detection)
- [ ] Add OpenMC model to [`analysis/neutron`](analysis/neutron)
- [ ] Add Tritium model to [`analysis/tritium`](analysis/tritium)
- [ ] Add the right version tags to [`environment.yml`](environment.yml)
- [ ] Add and update information in the README
- [ ] Modify [binder](https://mybinder.org/) badge by inserting the repo name
- [ ] Update [CI workflows](.github/workflows)
- [ ] Make first release on GitHub
- [ ] Link Zenodo record (created automatically) to the [LIBRA-project Zenodo community](https://zenodo.org/communities/libra-project/records)