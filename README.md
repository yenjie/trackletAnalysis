# trackletAnalysis
Tracklet analysis for dN/deta and multiplicity fluctuation analyses

## producing tracklet trees
in analysis/trackletProduction:
```bash
./configure.sh [5/8] [0/1/2]                     // 0: DATA, 1: EPOS, 2: HIJING
./produce_analysis_samples [5/8] [0/1/2] [label] // 0: BOTH, 1: DATA, 2: MC
```

## calculating final results
in analysis/dNdetaAnalysis:
```bash
./run-final-result.sh [5/8] [label]
```

## everything
in analysis/dNdetaAnalysis:
```bash
./run-all.sh
```
