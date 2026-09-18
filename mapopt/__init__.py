"""mapopt: autoresearch hill-climbing for world map projections.

Pure-Python (stdlib only) toolkit that:
  * defines a flexible, symmetry-constrained *family* of projections that is
    not any named historical projection (mapopt.family),
  * measures distortion rigorously via Tissot's indicatrix (mapopt.metrics),
  * implements correctly-formulated classic projections as baselines
    (mapopt.classics),
  * builds an approximate human-population weighting (mapopt.population),
  * hill-climbs the family's coefficients along a multi-objective Pareto
    frontier of shape vs. area error (mapopt.optimize).

No third-party dependencies: everything runs on a stock Python 3 install so the
research is fully reproducible from the committed source.
"""

__all__ = ["family", "metrics", "classics", "population", "optimize"]
