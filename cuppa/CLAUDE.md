# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

CUPPA (Cancer of Unknown Primary Prediction Algorithm) is a tissue of origin classifier using whole genome sequencing (DNA) and whole transcriptome sequencing (RNA) data. It consists of two main components:

- **Java component**: Feature extraction from sequencing pipeline outputs (PURPLE, LINX, ISOFOX, VirusInterpreter)
- **Python component**: Classification using scikit-learn multinomial logistic regressions

## Build and Test Commands

### Java Component

Build the JAR file:
```bash
mvn clean package
```

Run tests:
```bash
mvn test
```

Run a specific test class:
```bash
mvn test -Dtest=CuppaDataPrepTest
```

### Python Component

Install the Python package in development mode:
```bash
cd src/main/python/pycuppa
pip install -e .
```

Install with development dependencies:
```bash
pip install -e ".[dev]"
```

Run Python tests:
```bash
cd src/main/python/pycuppa
pytest tests/
```

Run specific test file:
```bash
pytest tests/classifier/test_cuppa_prediction.py
```

## Architecture

### Data Flow

```
Raw Files (VCF, TSV, CSV from PURPLE/LINX/ISOFOX)
    ↓
Java: CuppaDataPrep extracts features
    ↓
TSV Files: cuppa_data.cohort.*.tsv.gz (feature matrix)
    ↓
Python: CuppaClassifier predicts cancer type
    ↓
Outputs: Probabilities TSV + Visualization PNG
```

### Java Architecture (Feature Extraction)

**Main Entry Point**: `com.hartwig.hmftools.cup.prep.CuppaDataPrep`

**Core Abstraction**: `CategoryPrep` interface - each feature category (SNV, SV, DRIVER, SAMPLE_TRAIT, GENE_EXP, ALT_SJ) has a concrete implementation that extracts features from pipeline output files.

**Key Data Structures**:
- `DataItem`: Single feature value with Index (source, type, key) and Value
- `CategoryType`: Enum of 6 feature categories (4 DNA, 2 RNA)
- `ItemType`: Enum of ~45K feature names (SNV96, GEN_POS, drivers, fusions, genes, etc.)
- `DataSource`: DNA or RNA

**Package Organization**:
- `cup.prep`: Orchestration (CuppaDataPrep, CategoryPrep, config, I/O, matrix operations)
- `cup.somatics`: SNV extraction - 96 trinucleotide contexts + 6,101 genomic position bins
- `cup.svs`: Structural variant features - LINE, simple DEL/DUP, telomeric events, chromothripsis proxy
- `cup.drivers`: Driver mutations, gene fusions, viral insertions from LINX/VirusInterpreter
- `cup.rna`: Gene expression (TPM) and alternative splice junctions from ISOFOX
- `cup.traits`: Sample-level features - sex, WGD, microsatellite indels from PURPLE
- `cup.liftover`: VCF genome version conversion utilities

**Configuration**: `PrepConfig` manages:
- Single or multi-sample mode (sample_id_file supports wildcards in directory paths)
- Input directories for each tool (PURPLE, LINX, ISOFOX, VirusInterpreter)
- Output mode: single file or split by category (memory efficient for large cohorts)
- Reference genome version (V37/V38)
- Threading for parallel sample processing

### Python Architecture (Classification)

**Main Entry Points**:
- Training: `cuppa/train.py` → `TrainingRunner`
- Prediction: `cuppa/predict.py` → `PredictionRunner`

**Core Classifier**: `CuppaClassifier` is a 3-layer pipeline:

1. **Sub-Classifiers** (5 independent classifiers per feature type):
   - `GenPosClassifier`: Cosine similarity to per-cancer-type genomic position profiles
   - `Snv96Classifier`: Trinucleotide context patterns
   - `EventClassifier`: Driver mutations, fusions, viruses, SV types, TMB, WGD, sex
   - `GeneExpClassifier`: Gene expression (RNA)
   - `AltSjClassifier`: Alternative splice junctions (RNA)

   Each uses LASSO (L1) regularized logistic regression.

2. **Meta-Classifiers**:
   - `DNA_COMBINED`: Weighs DNA sub-classifiers (GEN_POS, SNV96, EVENT)
   - `RNA_COMBINED`: Weighs RNA sub-classifiers (GENE_EXP, ALT_SJ)

   Each uses Ridge (L2) regularized logistic regression with probability calibration and overrides.

3. **Combined Layer**: Multiplies DNA and RNA probabilities (with 0.01 floor), normalizes to sum=1

**Key Components** (`cuppa/components/`):
- `ProfileSimilarityTransformer`: Builds consensus profiles per cancer type, calculates cosine similarity
- `NoiseProfileAdder`: Adds background noise to stabilize small cancer type cohorts
- `NonBestSimilarityScaler`: Exponentiates cosine similarities to spread values (exponent=5 for GEN_POS/GENE_EXP, 2 for ALT_SJ)
- `Chi2FeatureSelector`: Statistical feature selection (FDR < 0.001) for EVENT features
- `RollingAvgCalibration`: Isotonic regression-based probability calibration using rolling windows
- `FusionProbOverrider`: Forces high probability for pathognomonic fusions (64 fusion-cancer type mappings)
- `SexProbFilter`: Zeros probabilities for opposite-sex cancer types

**Module Organization**:
- `cuppa/classifier/`: Core classifier implementations
- `cuppa/components/`: Transformers, filters, calibration
- `cuppa/compose/`: Pipeline extensions with caching support
- `cuppa/runners/`: Training and prediction orchestration
- `cuppa/sample_data/`: Feature/metadata loading from TSV files
- `cuppa/visualization/`: Output plotting (PNG) and data export (TSV)
- `cuppa/performance/`: Cross-validation and metrics

## Important Design Patterns

1. **Strategy Pattern**: Each feature category uses `CategoryPrep` interface with pluggable implementations
2. **Pipeline Composition**: Python components chain transformers in sklearn pipelines
3. **Stratified K-Fold CV**: Training samples split evenly by cancer type and RNA read length to handle class imbalance (default 10 folds)
4. **Probability Calibration**: RollingAvgCalibration ensures "0.8 probability" means 80% empirical accuracy
5. **Biological Rule Overrides**: Pathognomonic fusions and sex force/suppress probabilities based on domain knowledge
6. **Caching**: Python training can cache expensive transformations for efficient reruns
7. **Class Balancing**: Logistic regressions use `class_weight="balanced"` to weight by inverse sample frequency

## Cancer Type Organization

- **Training Set**: 43 distinct cancer subtypes (e.g., "Lung: Non-small cell: LUAD", "Breast: Triple negative")
- **Grouping**: Subtypes grouped by prefix before colon (e.g., "Lung:", "Breast:") for visualization
- **Probabilities**: Sub-classifiers predict subtypes; group probabilities are sums of subtype probabilities

## Feature Categories and Counts

| Source | Category | Description | Count |
|--------|----------|-------------|-------|
| DNA | GEN_POS | SNVs per 500kb genomic bin | 6,101 |
| DNA | SNV96 | SNVs per 96 trinucleotide context | 96 |
| DNA | EVENT | Drivers, fusions, viruses, SVs, TMB, WGD, sex | ~777 |
| RNA | GENE_EXP | log(TPM+1) per gene | 37,985 |
| RNA | ALT_SJ | log(read count+1) per alt splice junction | 236,328 |

## Key Configuration Concepts

### Java (PrepConfig)
- **Wildcard Paths**: Multi-sample mode uses wildcards in directory paths (e.g., `"/data/datasets/*/purple/"`) which are replaced with sample IDs
- **Conditional Requirements**: Either `-sample` OR `-sample_id_file` required; input directories vary by category
- **Category Selection**: `-categories ALL|DNA|RNA` or specific categories (SNV, SV, DRIVER, SAMPLE_TRAIT, GENE_EXP, ALT_SJ)
- **Output Splitting**: `-write_by_category` creates separate TSV files per category (prevents memory issues on large cohorts)

### Python (Runners)
- **Cancer Type Filtering**: `--excl_classes` and `--min_samples_with_rna` filter training samples
- **Fusion Overrides**: `--fusion_overrides_path` specifies TSV with fusion-cancer type mappings for FusionProbOverrider
- **CV Predictions**: `--cv_predictions_path` allows reusing cross-validation predictions instead of recomputing
- **Caching**: `--no_cache_training` disables model caching during CV (default is cached for speed)

## Common Development Scenarios

### Adding a New Feature Type

1. **Java**: Create new `CategoryPrep` implementation in appropriate package
   - Define new `ItemType` enum values for feature names
   - Implement `extractSampleData()` to load input files and return `List<DataItem>`
   - Add to `CuppaDataPrep.createCategoryPrep()` switch statement

2. **Python**: Create new sub-classifier in `cuppa/classifier/`
   - Extend from sklearn Pipeline
   - Add transformers (scaling, feature selection, etc.)
   - Add logistic regression as final step
   - Register in `CuppaClassifier` constructor

### Modifying Feature Transformations

Edit components in `cuppa/components/` (e.g., change `NonBestSimilarityScaler` exponent, adjust `NoiseProfileAdder` noise counts, modify `Chi2FeatureSelector` FDR threshold).

### Debugging Predictions

- Check `cuppa.vis.tsv` output for feature contributions (`data_type=feat_contrib`)
- Examine SHAP-derived odds ratios for EVENT features
- Review sub-classifier probabilities (`clf_name` column) to identify which feature types drive predictions

## Testing Strategy

- **Java**: Unit tests in `src/test/java/` with sample pipeline outputs in `src/test/resources/`
- **Python**: Unit tests in `src/main/python/pycuppa/tests/` with mock data generators
- **Notebooks**: Example training/prediction workflows in `src/main/python/pycuppa/doc/notebooks/`

## Dependencies

### Java
- hmf-common (shared utilities)
- patient-db (database access for training set)
- Apache Commons (IO, Lang3)
- Immutables (code generation for value objects)
- JUnit (testing)

### Python
- pandas 2.0.* (fixed for compatibility)
- numpy <2.0.0 (v1 binary structure required)
- scikit-learn 1.3.0 (exact version required for pickle compatibility)
- pytest, ipython (dev dependencies)

## Resource Files

Pre-trained classifiers and reference files (e.g., `alt_sj.selected_loci.tsv.gz`) are available via the HMF Resource page (see README).
