# Graph Input Contract

The external validation runner consumes graph artifacts from JSON.

## Required shape

```json
{
  "name": "graph_name",
  "nodes": [
    {"id": "n1", "sector": "alpha"}
  ],
  "edges": [
    {
      "source": "n1",
      "target": "n2",
      "strength": 0.9,
      "distance": 0.4,
      "curvature_penalty": 0.85
    }
  ]
}
```

## Rules

- `nodes[].id` must be unique
- `nodes[].sector` is optional in theory, but the current runner expects a sector string for sector-order output
- edges are undirected
- `source` and `target` must refer to existing node ids
- `strength`, `distance`, and `curvature_penalty` must be numeric
- the runner does not mutate the input artifact

## Current default artifact

- `external_validation/data/toy_graph.json`

## CLI

- default graph: `python3 external_validation/runners/run_validation.py`
- explicit graph: `python3 external_validation/runners/run_validation.py --graph-file /abs/path/to/graph.json`
