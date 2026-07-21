# RP-01 Parity Architecture Selection

Status: frozen before calibration and validation. No architecture was executed
before this ranking.

## Candidates

| Architecture | Locality | Distance | Memory | Ambiguity | HAOS fit | Deletion resistance | Identity fit | Transparency | Negative value | Total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Local `[3,1,3]` repetition-coset blocks | 5 | 5 | 4 | 5 | 5 | 3 | 5 | 5 | 5 | 42 |
| Periodic cycle code | 4 | 4 | 5 | 3 | 5 | 4 | 4 | 4 | 5 | 38 |
| Random local LDPC | 4 | 5 | 4 | 3 | 4 | 5 | 4 | 2 | 5 | 36 |
| Hierarchical local/mesoscopic checks | 3 | 5 | 3 | 3 | 5 | 5 | 5 | 2 | 5 | 36 |

Scores are `1` to `5`; higher memory score means lower memory cost.

## Selected Architecture

Sixteen independent local blocks use:

```text
H_local = [[1, 1, 0],
           [0, 1, 1]]
```

The global `H` is block diagonal with shape `32 x 48`.

- stored parity bits: 32;
- check degree: 2;
- variable degree: 1 or 2;
- decoder radius: one coarse cell;
- maximum decoder visibility: 3 current symbols + 2 stored parity bits;
- code dimension: 16;
- minimum distance: 3;
- declared correction radius: 1 symbol error per local block.

Syndromes `10`, `11`, and `01` identify the first, second, and third local bit.
Two errors in one block can be miscorrected to a wrong codeword; three errors
are undetected. This is a falsifiable boundary, not a flexible decoder.

Parity bits are stored as two bits at each coarse-cell decoder. No node or
decoder has the full 32-bit memory. No global fallback is allowed in the
admissible mechanism.
