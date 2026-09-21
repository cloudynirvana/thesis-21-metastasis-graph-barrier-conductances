# Anatomical metastasis graphs with edge-wise desmoplastic conductances: lumped burden and edge-rate identifiability

**Thesis #21.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Depends on:** Thesis #5 (metastasis as a stochastic process on anatomical graphs) and Thesis #11 (desmoplastic transport versus a lumped burden equation).

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

When each organ-graph edge carries an explicit desmoplastic conductance, do metastatic edge rates remain unidentified from lumped systemic burden, or do barrier observables restore an identifiable subset?

On a four-node toy, each directed edge carries a shedding rate and a desmoplastic conductance. Node burdens see only the series flux. A barrier reading, when scheduled, is the stalled fraction on that edge. The lumped sum leaves the six parameters at practical rank 2 of 6, with a three-dimensional structural kernel along the series level set. Site-resolved nodes, which identify the three fluxes when each edge has one rate, remain at rank 3 of 6 once the conductance is free. Barrier readings restore a subset: all six parameters when the nodes and the three stalled fractions are both recorded, and a thinner subset under thinner schedules. A hold-up generator is full structural rank under the lump and practical rank 3 of 6. Chapter Four is that laboratory. It is not the five-node single-rate experiment of Thesis #5, and it is not the three-shell lesion of Thesis #11.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/edge_conductance_identifiability.py` | Seeded graph, Fisher ranks, and profiles (seed 20260921) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Trajectories, level set, spectra, and profiles |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/edge_conductance_identifiability.py
python3 build_pdf.py
```

NumPy, SciPy, and Matplotlib are required for the toy. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Anatomical metastasis graphs with edge-wise desmoplastic conductances: lumped burden and edge-rate identifiability [Internet]. Thesis #21 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-21-metastasis-graph-barrier-conductances

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
