# HAOS-IIP Clean-Baseline Identity Resolution 51.2

## Purpose

Freeze the local `51.2` state after the clean-baseline follow-up program to `51.1`.

The repository no longer needs to treat the vector branch's clean-baseline identity question as a general unresolved branch-survival problem. The bounded closure now supports a narrower and stronger read:

- the clean torus carries the same active coexact branch;
- exact zero-clean symmetry obscures that branch at the raw mode-label level;
- tiny symmetry splitting plus family-level labeling resolves the branch as a stable first-shell family.

## What changed after 51.1

`51.1` closed the main `V1-V6b` validation line and isolated one open vector question:

> does the clean periodic baseline carry the same active coexact family as the same physical object under refinement, or only a symmetry-degenerate subspace whose identity is not well captured by raw mode-level transport?

`51.2` executes the clean-baseline closure path in bounded form:

1. `CB1-lite`: contiguous low-window subspace scanning
2. `CB2`: infinitesimal pinning continuation
3. `CB4`: tiny cycle-holonomy continuation
4. `CB3`: family-level symmetry labeling on the resolved holonomy window

## Bounded execution summary

### `CB1-lite`

Simple low-window sliding did **not** close the clean torus on its own.

- baseline `12 -> 16`: best projector affinity `0.320`
- baseline `16 -> 20`: best projector affinity `0.323`
- mild-disorder control remained stable under the same metric

This showed the clean problem was not just “pick a different 4-mode window”.

### `CB2`

Infinitesimal pinning gave the first positive resolution signal.

- exact clean pinning `0.000`: open, worst best affinity `0.366`
- tiny pinning `0.015`: resolved on `12 -> 16`, best affinity `0.956`

This supported the reading that the clean issue is degeneracy-sensitive, not branch death.

### `CB4`

Tiny cycle holonomy resolved the same clean family without leaving periodic topology.

- `12 -> 16`: zero-holonomy best affinity `0.467`, resolved at phase `0.050` with best affinity `0.974`
- `16 -> 20`: zero-holonomy best affinity `0.342`, resolved at phase `0.050` with best affinity `0.961`

This was the decisive branch-existence closure step.

### `CB3`

The resolved holonomy window was then labeled directly in bounded family form on sizes `n = 12, 16, 20`.

Resolved-window summary:

- first-shell momentum support stayed at least `0.999964`
- coarse field-axis support stayed at least `0.999981` on `A_x`
- window-level `k_y` fraction ranged from `0.456913` to `0.508369`
- window-level `k_z` fraction ranged from `0.491613` to `0.543052`
- the family label set stayed stable across sizes
- the raw mode ordering still shuffled inside that stable family

That is exactly the bounded closure we needed. The branch survives as the same labeled first-shell family even though the raw clean basis can still rotate inside it.

## What 51.2 now supports

The repository can now say:

> the clean periodic baseline carries the same active coexact vector branch as a symmetry-degenerate first-shell family. Exact zero-clean symmetry can reshuffle raw basis vectors, but the branch itself is no longer open in the substantive sense.

## What 51.2 does not claim

This note does **not** claim:

- a full exact-zero point-group irrep package;
- perfect one-to-one raw mode identity on the untwisted clean torus;
- a broader continuum claim than the bounded branch already validated through `51.1`;
- that no cleaner exact-zero labeling could still be developed later.

## Practical meaning

The vector branch's identity problem is no longer the main open question.

The remaining work, if desired, is packaging:

- cleaner exact-zero symmetry language;
- stronger group-theoretic labeling;
- larger-size confirmation on the already-resolved family.

Those are refinements, not the original identity obstruction.
