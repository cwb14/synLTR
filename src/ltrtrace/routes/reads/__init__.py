"""The reads route: the reference's annotated elements genotyped from long reads aligned to it (ltrkin's engine).

Modules: reads (per-read allele evidence at one element: full, solo LTR, empty site, other), phase (how many
distinct haplotypes carry each allele), model (the read set's depth and read-span model, likelihoods, priors and
calls), genotype (one read set: checkpointed evidence over every element, the model's passes, genotypes.tsv) and run
(the route over a samples sheet). Inside the route coordinates are 0-based half-open, from run.load_targets on.
"""
