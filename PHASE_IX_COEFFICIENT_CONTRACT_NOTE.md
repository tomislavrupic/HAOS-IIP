"Coefficient ratios and rescaled descriptors defined in Phase IX are part of the frozen spectral-trace contract and must be reused unchanged by later phases."

Phase IX extends the frozen spectral-trace contract by locking the following branch-local descriptors:

- the scaled coefficient family `b0(h)`, `b1(h)`, `b2(h)` extracted on the frozen Phase VIII short-time window
- the ratio descriptors `b1(h) / b0(h)` and `b2(h) / b0(h)`
- the rescaled descriptors `I0`, `I1`, and `I2` defined in `phase9-invariants/phase9_manifest.json`
- the deterministic generic-control comparison rule recorded in `phase9-invariants/phase9_manifest.json`

Later phases may compare against these objects, extend them explicitly, or supersede them only through a new recorded authority step.

Later phases may not silently:

- redefine the coefficient ratios while keeping the same Phase IX contract language
- replace the rescaled descriptors with new normalization rules without declaring the change
- treat later coefficient behavior as comparable to Phase IX if the frozen descriptors were changed in place

Any later phase that requires different ratio rules, different rescalings, or different descriptor families must declare that break explicitly and record the superseding authority artifact.
