"""mapopt: autoresearch hill-climbing for world map projections.

Pure-Python (stdlib only) toolkit that:
  * defines a flexible, symmetry-constrained *family* of projections that is
    not any named historical projection (mapopt.family),
  * measures distortion rigorously via Tissot's indicatrix (mapopt.metrics),
  * implements correctly-formulated classic projections as baselines
    (mapopt.classics),
  * hill-climbs the family's coefficients across the 3-objective Pareto surface
    of shape vs. area vs. distance error (mapopt.optimize).

No third-party dependencies: everything runs on a stock Python 3 install so the
research is fully reproducible from the committed source.
"""

__all__ = ["family", "metrics", "classics", "optimize"]
