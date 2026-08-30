Phylogeny-colored de Bruijn graph workflow
==========================================

This workflow builds a compacted colored de Bruijn graph with Cuttlefish and
recolors its unitigs by the parsimonious presence/absence reconstruction of each
unitig over a supplied phylogeny.

Outputs include the recolored GFA, compression metrics comparing raw color
vectors against the parsimony-assigned labels, and optional concordance
breakpoint, boundary severity, evolutionary persistence, and evolutionary
stratification analyses.
