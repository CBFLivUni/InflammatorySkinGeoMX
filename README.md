# Psoriasis spatial transcriptomics with Nanostring GeoMX

![image](https://github.com/CBFLivUni/InflammatorySkinGeoMX/assets/8311721/5a375f3e-5eec-4aa8-8307-d3a9a4b6d7e2)


## Running the analysis
### Reproducibly with singularity

1. Install [`Singularity`](https://docs.sylabs.io/guides/3.8/user-guide/)

2. Build the singularity container:
    ```console
     sudo singularity build runAnalysis.img Singularity
    ```
3. Run the analysis and render the html notebooks with a single command:
    ```console
     ./runAnalysis.img
    ```

### Alternatively using RScript

1.	Install the needed R packages
    ```console
     RScript install/install.R
    ```
2.	Run the analysis and render the html notebooks
    ```console
     RScript runAnalysis.R
    ```
