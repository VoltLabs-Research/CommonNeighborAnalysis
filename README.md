# CommonNeighborAnalysis

`CommonNeighborAnalysis` classifies atoms by local crystal environment and exports the reconstructed state consumed by DXA-compatible downstream tools.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- CommonNeighborAnalysis
```

## CLI

Usage:

```bash
common-neighbor-analysis <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--crystalStructure <type>` | No | Input crystal structure: `FCC`, `BCC`, `HCP`, `CUBIC_DIAMOND`, `HEX_DIAMOND`. | `FCC` |
| `--dissolveSmallClusters` | No | Mark small clusters as `OTHER` after clustering. | `false` |
| `--help` | No | Print CLI help. | |
