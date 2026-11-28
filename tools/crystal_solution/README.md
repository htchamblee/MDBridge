# Crystal against solution helix parameters (docs/reports/2026-11-28-crystal-vs-solution.md)

Measures base-pair-step helical twist, rise, inclination and x-displacement (Olson et al. 2001
standard reference frames, written from the published definitions) of duplexes solved both by X-ray
crystallography and by solution NMR, to estimate how far the crystal-derived geometry of the
non-Watson-Crick builder (DIFFERENCES CH-20) may be from solution.

    mkdir cif && cd cif
    for id in 404D 1RRD 355D 1BNA 1NAJ 1DUF 2G3S 1EKD 2BJ6 1EJZ 1AXP; do
        curl -sO https://files.rcsb.org/download/$id.cif; done
    sha256sum -c ../tools/crystal_solution/SHA256SUMS      # the files measured on 2026-11-28
    python3 -I -B ../tools/crystal_solution/bpframes.py params.json *.cif
    python3 -I -B ../tools/crystal_solution/summarize.py params.json

The structures are not committed (wwPDB files, CC0, fetched by id). Later wwPDB remediation can
change the files; SHA256SUMS records the ones measured. Needs NumPy.
