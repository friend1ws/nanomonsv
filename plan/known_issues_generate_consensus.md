# Known issues: generate_consensus / vcf_convert

## Status: WATCH
- Created: 2026-03-16

These issues were found during GIAB workflow testing (v0.9.0b2).
Currently mitigated but not fully resolved.

## 1. Excessively long sbnd consensus sequences

- ONT ultra-long data can produce sbnd contigs up to 208KB
- No upper limit on consensus length in `generate_consensus.py`
- Current mitigation: csv.field_size_limit set to 10MB in vcf_convert.py, so no crash
- Remaining concern: bloated output, wasted computation in downstream steps (insert_classify etc.)
- Potential fix: skip sbnd consensus when length exceeds 100KB
  - `generate_consensus.py` L312 area, before final print

## 2. minimap2 core dumps from RLIMIT_AS on satellite repeats

- Satellite repeat sequences (CCATTC, ATGGAAT repeats, 100-600KB) cause minimap2 to consume excessive virtual address space
- `preexec_fn` sets RLIMIT_AS → minimap2 calls abort() → SIGABRT → core dump files
- nanomonsv catches the exception and continues (WARNING + skip), so results are correct
- Current mitigation: `-f 0.05` added to all minimap2 calls in generate_consensus.py (v0.9.0b2)
  - This filters repetitive minimizers, preventing memory explosion
  - Tested: 0 SIGABRTs, 0 core dumps across 8 samples
- Remaining concern: RLIMIT_CORE is not set to 0, so if a new pattern triggers the issue, core dumps will still be written
- Potential fix: add `resource.setrlimit(resource.RLIMIT_CORE, (0, 0))` in `preexec_fn`

## 3. RLIMIT_AS vs glibc malloc arenas

- glibc malloc creates per-thread arenas that consume virtual address space (VIRT) without using real memory (RSS)
- RLIMIT_AS limits virtual address space, not real memory
- With default malloc behavior, 2GB RLIMIT_AS can be hit even when RSS is well under 1GB
- Current setting: max_memory_minimap2 default changed from 2 to 4 (v0.9.0b2)
- The `-f 0.05` fix largely eliminates the problem, but the fundamental RLIMIT_AS vs malloc arena mismatch remains
- Alternative approaches considered but not implemented:
  - MALLOC_ARENA_MAX=1 (works but may degrade performance)
  - Switch to RLIMIT_RSS (not enforced on modern Linux)
