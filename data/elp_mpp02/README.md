# ELP/MPP02 source data and attribution

The 14 `elp_*` files are the full, untruncated ELP/MPP02 coefficient tables
distributed by Yuk Tung Liu's C++ implementation:

- Repository: https://github.com/ytliu0/ElpMpp02
- Pinned revision: `e25e08f6082a3f3211b3d14e062d935558ad3ea6`
- Theory: J. Chapront and G. Francou, “The lunar theory ELP revisited.
  Introduction of new planetary perturbations,” A&A 404, 735–742 (2003),
  https://doi.org/10.1051/0004-6361:20030529
- Implementation notes: https://ytliu0.github.io/ElpMpp02/docs/ElpMpp02.pdf
- Original theory distribution, as identified upstream:
  `ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02/`

These files match the upstream Git blobs. `manifest.json` records their SHA-256
hashes and term counts (35,901 total). The tables contain analytical model
coefficients, not a sampled lunar trajectory. They retain their upstream units
and ordering. No terms have been truncated, re-fitted or adjusted in the files.

The upstream repository distributes its implementation and data under GNU GPL
version 3; its license is preserved in `LICENSE`. Attribution: Yuk Tung Liu,
with the underlying lunar theory and original data by Chapront and Francou.
Spp's argument/parameter equations and J2000 rotation were adapted from Liu's
`ElpMpp02.cpp`. Spp replaces runtime file loading and pointer-owned arrays with
compiled immutable tables, adds typed interfaces and validation, and implements
series summation with analytical derivatives using its own storage layout.

`tools/elp_gen.cpp` validates row widths, counts and numbers and emits exact
hexadecimal binary64 literals. It preserves all seven main-series coefficients;
the fitted amplitude adjustment uses A and B1–B5, as in the reference code.
The generated include is a build artifact, never downloaded during builds.

Both upstream fitted sets are supported: lunar laser ranging (LLR), and
DE405/DE406. The initial default is DE405/DE406 to make the long-range comparison
configuration explicit. This is not a fit to DE440/DE441, and selecting it does
not establish uniform astronomical accuracy throughout the computational range.
