"""Genome simulators with known LTR-RT truth (`synltr sim`): populations with gene flow (pop), and the scorer of
calls against that truth (score). Nothing is imported here, so building the command line loads no numpy; the
population model imports msprime only inside the functions that run it, so the package works without the extra
`sim`."""
