# RP-01 Relational Alphabet

Status: frozen before calibration and validation.

## Definition

The continuous system is an `8 x 8` periodic scalar grid. The grid is divided
into sixteen non-overlapping `2 x 2` coarse cells. Each cell contributes three
binary local relations at its upper-left anchor `(i,j)`:

```text
z_0 = 1[x_(i+1,j) - x_(i,j) >= 0]       horizontal order
z_1 = 1[x_(i,j+1) - x_(i,j) >= 0]       vertical order
z_2 = 1[x_(i+1,j+1) - x_(i,j) >= 0]     diagonal order
```

The resulting alphabet has `16 * 3 = 48` primary symbols.

## Symbol Audit

| Family | Local radius | Source | Quantization | Identity role | Function relation | Circularity | Ambiguity | Bits |
| --- | ---: | --- | --- | --- | --- | --- | --- | ---: |
| Horizontal order | 1 edge | current endpoint values | sign at zero | local orientation | constrains low-mode response indirectly | low | no magnitude | 16 |
| Vertical order | 1 edge | current endpoint values | sign at zero | local orientation | constrains low-mode response indirectly | low | no magnitude | 16 |
| Diagonal order | `sqrt(2)` | current diagonal endpoints | sign at zero | local sector ordering | constrains mixed-axis response indirectly | low | no magnitude | 16 |

The symbols do not encode the functional label, final recovery decision,
future state, node magnitude, or continuous trajectory. Multiple continuous
states and `2^16` distinct symbol vectors share each complete stored parity
vector under the selected code.

Unperturbed stability is measured on construction fixtures before any recovery
run. Ambiguous zero differences are resolved deterministically as bit `1` and
reported; they are not fitted against outcomes.

## Rejected Alphabets

- Direct quantization of the four functional probe responses was rejected
  because it would encode the evaluation target.
- Multi-bit edge magnitudes were rejected because their memory approaches a
  compressed continuous checkpoint.
- Future transition labels were rejected as target leakage.
