# 1000 Genomes PCA and Ancient DNA Analysis Pipeline

This repository contains Nextflow pipelines for:
1. Performing Principal Component Analysis (PCA) on 1000 Genomes Project data combined with ancient DNA samples
2. Running Schmutzi contamination analysis on ancient DNA samples
3. Providing interactive 3D visualization tools for PCA results

## Pipeline Overview

The repository contains two main workflows:

### 1. PCA Pipeline (workflow.nf)
- Processes VCF files from the 1000 Genomes Project
- Merges them with ancient DNA samples
- Performs LD pruning
- Runs PCA analysis
- Generates visualization data

### 2. Schmutzi Pipeline (run_schmutzi.nf)
- Analyzes mitochondrial DNA contamination in ancient samples
- Prepares BAM files for analysis
- Estimates contamination using deamination patterns
- Generates contamination reports and consensus sequences

## Prerequisites

- Nextflow (21.04.0 or later)
- Docker or Singularity
- Python 3 (for serving the visualization)
- Modern web browser with WebGL support (for viewing the visualization)
- Access to 1000 Genomes Project data
- Schmutzi installation (for contamination analysis)

## Installation

1. Clone this repository:
```bash
git clone <repository-url>
cd <repository-name>
```

2. Install Schmutzi (if planning to use contamination analysis):
```bash
# Install to default location ($HOME/schmutzi)
./install_schmutzi.sh

# Or specify custom installation directory
./install_schmutzi.sh --prefix /path/to/install
```

## Running the Pipelines

### PCA Analysis Pipeline

Configure your input data in nextflow.config or provide parameters on the command line:

```bash
nextflow run workflow.nf \
  --thousand_genomes_vcf "/path/to/1kg/*.vcf.gz" \
  --ancient_vcfs "/path/to/ancient/*.vcf.gz" \
  --outdir results
```

### Schmutzi Contamination Analysis

Run the Schmutzi workflow on your BAM files:

```bash
nextflow run run_schmutzi.nf \
  --bams '/path/to/your/*.bam' \
  --schmutzi_path "$SCHMUTZI_PATH" \
  --outdir 'results_schmutzi' \
  --library_type 'double' \
  --length_deam 10
```

Key Schmutzi parameters:
- `--bams`: Path pattern to your BAM files (required)
- `--schmutzi_path`: Path to Schmutzi installation (required)
- `--library_type`: Library type ('single' or 'double')
- `--length_deam`: Length considered for deamination
- `--contamination_prior`: Prior contamination estimate (optional)

## Output Files

### PCA Pipeline Output
- `pca_data.json`: PCA coordinates and metadata for visualization
- Various QC reports and intermediate files
- PCA plots in PDF format

### Schmutzi Pipeline Output
For each sample, creates a directory containing:
- Initial contamination estimates
- Final contamination estimates
- Endogenous consensus sequence
- Contaminant consensus sequence
- Detailed analysis reports

## Visualizing PCA Results

1. Navigate to the directory containing pca_visualization.html and pca_data.json:
```bash
cd results
```

2. Serve the visualization using Python's HTTP server:
```bash
python3 -m http.server 8080 --bind 0.0.0.0
```

3. Open your web browser and navigate to:
```
http://localhost:8080/pca_visualization.html
```

## Using the PCA Visualization
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

### General Issues
- If the visualization doesn't load, check your browser's console for errors
- Ensure both pca_visualization.html and pca_data.json are in the same directory
- Make sure your browser supports WebGL for 3D visualization
- If you see CORS errors, ensure you're using a web server rather than opening the HTML file directly

### Schmutzi-specific Issues
- Ensure BAM files are aligned to the rCRS reference
- Check that BAM files are properly indexed
- Verify Schmutzi installation with `schmutzi contDeam.pl --help`
- For memory issues, try reducing the number of threads
- Check the log files in the work directory for detailed error messages

## Citation
If you use these pipelines in your research, please cite:
- The 1000 Genomes Project
- Schmutzi: Renaud G, Slon V, Duggan AT, Kelso J. Schmutzi: estimation of contamination and endogenous mitochondrial consensus calling for ancient DNA. Genome Biol. 2015
- Plotly.js for visualization
- This repository