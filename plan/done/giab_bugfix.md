# GIAB Workflow Test Bugfix

## Status: DONE
- Created: 2026-03-11
- Completed: 2026-03-16
- Remaining concerns documented in: `known_issues_generate_consensus.md`

Bugs found during v0.9.0b2 workflow test with GIAB HG008 samples.

## 1. csv.field_size_limit in sbnd VCF conversion

- **Symptom**: `nanomonsv get` crashes at final step (`sbnd2vcf_convert`) with `_csv.Error: field larger than field limit (131072)`
- **Affected sample**: HG008-T-ONT-UL-NE-GRCh38 (ONT ultra-long)
- **Cause**: ONT sbnd contig sequences exceed Python csv module's default field size limit (128KB). Max observed: 208KB.
- **Impact**: `nanomonsv get` fails entirely — no result.vcf, no sbnd.result.vcf output. Downstream filt.pass and insert_classify cannot run.
- **File**: `nanomonsv/vcf_convert.py` (`sbnd2vcf_convert`)

### Fix

Add `csv.field_size_limit(sys.maxsize)` at the top of `sbnd2vcf_convert`.

## 2. Excessively long sbnd consensus sequences

- **Symptom**: sbnd contig sequences up to 208KB in ONT data. These are too long to be meaningful single breakend events.
- **Affected sample**: HG008-T-ONT-UL-NE-GRCh38 (22 contigs > 100KB)
- **Cause**: No upper limit on consensus sequence length in `generate_consensus.py`
- **Impact**: Bloated output files, triggers issue #1, wastes computation in downstream steps
- **File**: `nanomonsv/generate_consensus.py` (`print_consensus_sbnd`)

### Fix

Skip sbnd consensus output when length exceeds 100KB (hardcoded constant). Add check before the final `print()` at line 351:

```python
if len(consensus) > 100000: return
```

## 3. minimap2 core dumps from RLIMIT_AS

- **Symptom**: 44 core dump files totaling 46GB in working directory
- **Affected samples**: HG008-T-ONT-UL-NE-GRCh38 (9 cores), HG008-T-ONT-UL-NE-CHM13 (35 cores)
- **Cause**: `preexec_fn` sets `RLIMIT_AS` to 2GB for minimap2 subprocess. When exceeded, minimap2 calls `abort()` -> SIGABRT -> core dump. nanomonsv catches the exception and continues, but core dumps are written.
- **Impact**: Disk space waste only. nanomonsv handles the error correctly (WARNING + skip).
- **File**: `nanomonsv/generate_consensus.py` (`preexec_fn`)

### Fix

Add `RLIMIT_CORE = 0` in `preexec_fn` to suppress core dump generation:

```python
def preexec_fn(self):
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    limit = self.max_memory_minimap2 * 1024 ** 3
    ...
```

## Priority

1 > 2 > 3. Issue #1 is a crash bug. Issue #2 prevents #1 and improves output quality. Issue #3 is cosmetic (disk waste).

Fixing #2 largely prevents #1 from occurring (contigs would be capped at 100KB < 128KB limit), but #1 should still be fixed as a safety net.
