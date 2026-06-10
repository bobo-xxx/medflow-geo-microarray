# Platform Detection Rules

This document defines detection rules for different microarray platform file types.

## File Pattern Detection

| File Pattern | Platform | Notes |
|--------------|----------|-------|
| `*.CEL(.gz)?` | Affymetrix | RMA normalization required |
| `*.GPR(.gz)?` | Agilent | GenePix/Agilent variants |
| `*.idat` + `*.BPM` | Methylation | SKIP! Not supported |
| `*.idat` only | Illumina | Expression array |
| `*.txt` only | TXT | Series matrix |
| `*.txt` (Agilent FE) | Agilent FE | Feature Extraction output in RAW.tar |
|
## GPR Header Detection

GPR files require header analysis to distinguish variants:

- `^ATF` → GenePix format
- `^TYPE` → Agilent format

## Agilent FE TXT Detection

TXT files in RAW.tar require header analysis:
- Contains `ProbeName`, `GeneName`, `gTotalGeneSignal`, etc. → Agilent FE format
- Use `limma::read.maimages()` to process
- Single or dual channel depends on column structure

## Gene Symbol Column Priorities

When extracting gene identifiers from processed data:

1. Primary: Gene Symbol column
2. Fallback: Gene ID / RefSeq columns
3. Last resort: Probe ID mapping
