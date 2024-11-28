# 1000 Genomes PCA Visualization Pipeline

This repository contains a Nextflow pipeline for performing Principal Component Analysis (PCA) on 1000 Genomes Project data combined with ancient DNA samples, along with an interactive 3D visualization tool.

## Pipeline Overview

The pipeline:
1. Processes VCF files from the 1000 Genomes Project
2. Merges them with ancient DNA samples
3. Performs LD pruning
4. Runs PCA analysis
5. Generates visualization data

## Prerequisites

- Nextflow (21.04.0 or later)
- Docker or Singularity
- Python 3 (for serving the visualization)
- Modern web browser with WebGL support (for viewing the visualization)
- Access to 1000 Genomes Project data

## Running the Pipeline

1. Clone this repository:
```
bash
git clone <repository-url>
cd <repository-name>
```

```
bash
git clone <repository-url>
cd <repository-name>
```

2. Configure your input data in nextflow.config or provide parameters on the command line:

```
nextflow run workflow.nf \
  --thousand_genomes_vcf "/path/to/1kg/*.vcf.gz" \
  --ancient_vcfs "/path/to/ancient/*.vcf.gz" \
  --outdir results
```

The pipeline will generate several outputs in your results directory, including:
- `pca_data.json`: PCA coordinates and metadata for visualization
- Various QC reports and intermediate files
- PCA plots in PDF format

## Visualizing the Results

1. Navigate to the directory containing pca_visualization.html and pca_data.json:
```
cd results
```

2. Serve the visualization using Python's HTTP server:
```
python3 -m http.server 8080 --bind 0.0.0.0
```

3. Open your web browser and navigate to:
```
http://localhost:8080/pca_visualization.html
```

## Using the visualization
The interactive 3D visualization offers several features:

- Sample Types
Toggle between modern and ancient samples using checkboxes
- Ancient samples are shown as larger black diamonds with white borders
- Population Controls
  - Superpopulations: Filter major population groups (e.g., European, African)
  - Checking/unchecking a superpopulation automatically selects/deselects all its populations
  - Populations: Filter specific populations (e.g., CEU, YRI)
  - Individual populations can be toggled independently
  - Unchecking all populations in a superpopulation automatically unchecks the superpopulation
- Color Schemes
Toggle between two coloring modes using radio buttons:
  - Population: Each population gets a unique color within its superpopulation's color spectrum
  - Superpopulation: All populations within a superpopulation share the same color
  - AFR (African): Red spectrum
  - EUR (European): Blue spectrum
  - EAS (East Asian): Green spectrum
  - SAS (South Asian): Purple-Magenta spectrum
  - AMR (American): Orange-Yellow spectrum
- Ancient samples: Black
- 3D Controls
  - Rotate: Click and drag
  - Zoom: Mouse wheel or pinch gesture
  - Pan: Right-click and drag
  - Reset: Double-click to reset the view
- Interactive Features
  - Hover over points to see detailed sample information
  - Legend can be used to show/hide specific populations
  - All controls update the visualization in real-time

## Troubleshooting
- If the visualization doesn't load, check your browser's console for errors
- Ensure both pca_visualization.html and pca_data.json are in the same directory
- Make sure your browser supports WebGL for 3D visualization
- If you see CORS errors, ensure you're using a web server (like the Python command above) rather than opening the HTML file directly
- For performance issues, try:
  - Reducing the number of displayed populations
  - Using the superpopulation coloring mode
  - Using a more powerful computer/graphics card

## Citation
If you use this pipeline in your research, please cite:
- The 1000 Genomes Project
- Plotly.js for visualization
- This repository