"""The assembly route: each annotated element of genome A genotyped in another assembly B.

Stages, per ordered pair: anchors (anchor: flank pieces mapped onto B's minimap2 index, index, or a whole-genome
PAF; chained and netted, chain), windows around the elements (window), a LASTZ alignment of each window (lastz),
the state of each element (classify), written as the pair's call file (pair; calls.records). run.genotype drives
the pairs, after the elements stage (io.annotation.stage).
"""
