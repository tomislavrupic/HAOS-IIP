# Disciplined Memo

Question:

Does the fixed candidate become consistently more special under `anisotropic_local_disruptions`, the chosen HAOS-story perturbation proxy?

Answer:

Yes, under the explicit ranking rule used in this suite.

Decision rule:

Answer yes only if the story-family perturbation gives the fixed candidate the best or tied-best exhaustive rank among all `84` possible 6-position subsets in both baseline and bounded local repair, while also beating matched random and contiguous controls in both modes.

What supports the answer:

- Under `anisotropic_local_disruptions`, the fixed candidate ranks `7/84` at baseline and `3/84` under bounded local repair.
- The next-best family ranks are weaker: `10/84` at baseline and `5/84` under bounded local repair.
- Under `anisotropic_local_disruptions`, the fixed candidate beats sampled random controls by `+0.027819` at baseline and `+0.022395` under bounded local repair.
- Under `anisotropic_local_disruptions`, the fixed candidate beats sampled contiguous controls by `+0.068755` at baseline and `+0.054574` under bounded local repair.

What this does not show:

- It does not show that `anisotropic_local_disruptions` is a canonical HAOS perturbation family. That is an explicit proxy choice in the absence of a workspace-defined rule.
- It does not show that the fixed candidate is uniquely optimal. It remains below the top-ranked subset even in the story-family condition.
- It does not show that the story family maximizes every summary gap. Some larger control gaps appear under `burst_errors` or `corridor_limited_swaps`.
- It does not validate a HAOS-QATC bridge.

Conclusion:

Under the chosen directional-local perturbation proxy, the fixed candidate does become more special in exhaustive rank, but the result remains a bounded toy signal, not a theory validation.
