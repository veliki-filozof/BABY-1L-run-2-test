# Neutron detection analysis

## Diamond detector

1. Download the raw data from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14810785.svg)](https://doi.org/10.5281/zenodo.14810785)
2. Unzip the `data.zip` archive and place its content in `data/neutron_detection/`
3. Navigate to the correct directory: `cd analysis/neutron/detection`
3. Run `python process_raw_data.py`

This will probably take a little while and you should see in the terminal:

```
Added file: ..\..\..\data\neutron_detection\20241210_part1\UNFILTERED\Data_CH4@DT5725_1360_20241210.CSV containing 202756 events
Added file: ..\..\..\data\neutron_detection\20241210_part1\UNFILTERED\Data_CH4@DT5725_1360_20241210_1.CSV containing 202758 events
Added file: ..\..\..\data\neutron_detection\20241210_part1\UNFILTERED\Data_CH4@DT5725_1360_20241210_10.CSV containing 202758 events
```

After it is done, you should see three `.npy` files in `data/neutron_detection/`.

4. Now run the `main.ipynb` notebook