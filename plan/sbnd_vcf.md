# sbnd: VCF output

## Status: TODO

## Overview
Output sbnd results in VCF format, complementing the existing VCF conversion
for canonical SVs (`nanomonsv/vcf_convert.py`).

## VCF specification for single breakends

Per VCF 4.3 specification (Section 5.4 "Breakends"), a **single breakend** is a breakpoint
where only one side can be unambiguously placed on the genome.
The other end is represented by a period (`.`) in the ALT field.

### ALT field format

The ALT field can include the assembled contig sequence extending from the breakpoint:

| Breakpoint direction | ALT format | Meaning |
|---------------------|------------|---------|
| Forward (`+`) | `s[sequence].` | From reference base `s`, sequence extends rightward into unmapped region |
| Reverse (`-`) | `.[sequence]s` | Unmapped sequence joins from the left into reference base `s` |

Where `s` is the reference base at the breakpoint position, and `[sequence]` is the
assembled contig sequence at the breakpoint.

If no contig sequence is available, the ALT is simply `s.` or `.s`.

### Example VCF records

Forward single breakend at chr1:12345 (Dir=`+`) with contig `ATCGATCG`:
```
#CHROM  POS     ID      REF  ALT           QUAL  FILTER  INFO
chr1    12345   sbnd_0  A    AATCGATCG.    .     PASS    SVTYPE=BND
```

Reverse single breakend at chr3:67890 (Dir=`-`) with contig `GCTAGCTA`:
```
#CHROM  POS     ID      REF  ALT           QUAL  FILTER  INFO
chr3    67890   sbnd_1  G    .GCTAGCTAG    .     PASS    SVTYPE=BND
```

### Key differences from paired breakends (BND)
| Aspect | Paired breakend | Single breakend |
|--------|----------------|-----------------|
| ALT format | `s]chr:pos]`, `[chr:pos[s`, etc. | `s[seq].` or `.[seq]s` |
| MATEID | Required (links to mate record) | Not used |
| Number of VCF records | 2 per SV | 1 per sbnd |
| SVTYPE | BND | BND |

## Current state

### Canonical SV VCF conversion (`nanomonsv/vcf_convert.py`)
`genomesv2vcf_convert()` reads `*.nanomonsv.result.txt` and outputs VCF with:
- DEL, INS, DUP for intrachromosomal SVs
- Paired BND records for interchromosomal SVs (with MATEID)
- FORMAT fields: TR (total reads), VR (variant reads)

There is no sbnd VCF output yet.

### sbnd result columns
From `post_proc.py` L372:
```
Chr_1  Pos_1  Dir_1  Contig  SV_ID  Checked_Read_Num_Tumor  Supporting_Read_Num_Tumor
[Checked_Read_Num_Control  Supporting_Read_Num_Control]  (if control exists)
```

## Design

### Approach
Add a function (e.g., `sbnd2vcf_convert()`) in `vcf_convert.py` that converts
`*.nanomonsv.sbnd.result.txt` to VCF single breakend records.

### Mapping from sbnd result to VCF fields
| VCF field | Source |
|-----------|--------|
| CHROM | Chr_1 |
| POS | Pos_1 |
| ID | SV_ID |
| REF | Reference base at Pos_1 (from pysam.FastaFile) |
| ALT | `REF` + Contig + `.` if Dir_1=`+`; `.` + Contig + `REF` if Dir_1=`-` |
| QUAL | `.` |
| FILTER | Is_Filter (once implemented) or `PASS` |
| INFO | `SVTYPE=BND` |
| FORMAT | `TR:VR` with tumor (and optionally control) read counts |

### Open questions
- Should sbnd VCF records be appended to the same VCF file as canonical SVs,
  or written to a separate file?
- The Contig column contains the full assembled contig sequence. Should it be
  included entirely in ALT, or truncated/summarized?

## Related files
- `nanomonsv/vcf_convert.py` — Existing canonical SV VCF conversion
- `nanomonsv/post_proc.py` L371-374 — sbnd output header definition
- `nanomonsv/run.py` — Call site for VCF conversion
