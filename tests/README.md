# Existing tests

## tests/test_main.py (lightweight)
- `test_command` — Just runs `nanomonsv --version`

## tests/full_tests.py (E2E)
Requires GRCh38 reference genome (auto-downloaded on first run).
Tests run `nanomonsv parse/get/validate` against test BAM/CRAM files.

### parse
| Test | Input | Expected record counts (rearrangement/deletion/insertion) |
|------|-------|----------------------------------------------------------|
| `test_parse_tumor` | BAM | 96 / 179 / 165 |
| `test_parse_tumor_cram` | CRAM | 96 / 179 / 165 |

### get
| Test | Format | Control | Options | Expected (result / supporting_read) |
|------|--------|---------|---------|--------------------------------------|
| `test_get1_1` | CRAM | yes | (default: racon + single_bnd) | 12 / ~93 |
| `test_get1_2` | CRAM | yes | `--use_mafft --no_single_bnd` | 12 / ~92 |
| `test_get1_3` | CRAM | no | (default: racon + single_bnd) | 20 / ~168 |
| `test_get1_4` | CRAM | no | `--use_mafft --no_single_bnd` | ~19 / ~168 |
| `test_get1_5` | BAM | yes | (default: racon + single_bnd) | 12 / ~93 |
| `test_get1_6` | BAM | yes | `--use_mafft --no_single_bnd` | 12 / ~92 |
| `test_get1_7` | BAM | no | (default: racon + single_bnd) | 20 / ~168 |
| `test_get1_8` | BAM | no | `--use_mafft --no_single_bnd` | ~19 / ~168 |
| `test_get2` | BAM | yes | `--control_panel_prefix --use_mafft --no_single_bnd` | 11 / ~84 |
| `test_get3` | BAM | yes | `--no_single_bnd` | 12 / ~93 |
| `test_get4_1` | BAM | yes | `--processes 4 --use_mafft --no_single_bnd` | 12 / ~92 |
| `test_get4_2` | BAM | yes | `--processes 4 --use_mafft --no_single_bnd` (empty input) | 1 / 0 |
| `test_get4_3` | BAM | yes (swapped) | `--processes 4 --use_mafft --no_single_bnd` | 1 / 0 |

Note: `~` means the assertion allows ±5 tolerance.
None of the tests use `--simple_repeat_bed`.

### validate
| Test | Description |
|------|-------------|
| `test_validate` | Runs validate with control BAM, checks output matches input |
