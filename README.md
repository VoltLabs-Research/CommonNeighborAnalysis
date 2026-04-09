# CommonNeighborAnalysis

`CommonNeighborAnalysis` classifies atoms by local crystal environment and exports the reconstructed state consumed by DXA-compatible downstream tools.

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

## Build With CoreToolkit

```bash
cd /path/to/voltlabs-ecosystem/tools/CoreToolkit
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/StructureIdentification
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/CommonNeighborAnalysis
conan create . -nr
```
