# HAOS-IIP Explorer

Static public-facing interface for the local HAOS-IIP repository.

## What it includes

- landing page with the minimal HAOS-IIP definition
- interactive stability playground tied to the public oracle demo vocabulary
- oracle architecture viewer with repo references
- curated paper cards with conceptual and technical modes
- experiment atlas fed by real repository artifacts
- roadmap and developer entry sections

## Run locally

Serve the repository root, then open `haos_explorer/`.

```bash
cd /Volumes/Samsung\ T5/2026/HAOS/HAOS\ DOCS/HAOS-IIP
python3 -m http.server 8000
```

Open:

- `http://localhost:8000/haos_explorer/`

## Notes

- The preset scenario cards use the same JSON scenarios that drive `python3 -m haos_iip.demo stability`.
- The slider surface is a client-side public interaction model layered on top of that frozen vocabulary.
- The site stays backend-free for now so it can be served from the repository immediately.
