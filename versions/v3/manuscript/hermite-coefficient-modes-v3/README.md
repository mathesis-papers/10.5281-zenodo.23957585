# Hermite coefficient maxima

This package resolves the difference conjecture in OEIS A381423, Comments,
revision 8, and gives complete proofs of uniqueness, an explicit formula,
difference-block lengths and indices, and square-root growth. It also proves
the absolute-coefficient difference assertion in A381524 and the three formulas
for equal absolute maxima in A381737.

## Files and reproduction

- `note.tex`: the main source, with the bibliography included.
- `note.pdf`: the readable mathematical note.
- `verify.py`: exact-integer finite comparisons, requiring Python 3.10 or later
  and no third-party Python packages.
- `checks.json`: the deterministic finite-check summary.
- `metadata.json`: authors, license, mathematical coverage, and source locator.

Build the PDF with Tectonic 0.17.0:

```sh
tectonic --untrusted note.tex
```

Alternatively, use a standard LaTeX distribution:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is inside `note.tex`, so BibTeX is unnecessary. The source
uses standard AMS, Latin Modern, geometry, booktabs, and hyperref packages.

Reproduce the finite comparisons:

```sh
python3 verify.py --output checks-reproduced.json
cmp checks.json checks-reproduced.json
```

The checker constructs integer Hermite polynomials through degree 600 and
compares every coefficient with the finite expansion and each maximizing
exponent with the formulas, including mode interlacing and the equality
criterion for signed and absolute maximum values. It checks exact
successive-value ratio identities through degree 600, both selectors and
both block structures through degree one million, and the tie generating
function through coefficient 10,000. All comparisons use integers; no
square-root rounding is required.
Since `8*n+21` is 5 modulo 8, it is never a square. Consequently the floor in
the signed formula is exactly `(2*n+6-isqrt(8*n+21))//8`.
For the absolute formula, put `D=8*n+17`, `root=isqrt(D)` and
`correction=int(root*root != D)`; its mode index is exactly
`(2*n+5-root-correction)//4`. This includes the square-discriminant
boundaries associated with tied absolute maxima.

## Mathematical scope

The [pinned Lean declaration](https://github.com/the-omega-institute/trureturing/blob/a2d277ea0980f1230eff2c7e27838372effccb30/D5/S3/Analytic/Polynomial/HermiteSignedMode.lean)
`D5.S3.Analytic.Polynomial.HermiteSignedMode.result` proves the
least-signed-maximizer characterization and the difference assertion for every
natural degree. Its polynomial is defined by the physicists' finite expansion.
Uniqueness, both closed forms, both sets of block indices, mode interlacing,
maximal-value equivalence and ratio limits, and the absolute-coefficient
conclusions are established by the written proofs in the note.
The numerical comparisons have precisely the finite ranges stated above.

The exact posing locator is A381423, Comments, revision 8: proposed on
23 February 2025, with revision 8 dated 6 March 2025. The dated history is available at
https://oeis.org/history?seq=A381423 . A381524, Comments, and A381737,
Formula, contain the absolute-coefficient statements treated in the note.

The problem identifier for the principal conjecture is
`oeis-a381423-hermite-differences`. The other OEIS items are covered explicitly
in the note's coverage table.

## Authors and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

Copyright 2026 Haobo Ma and Wenlin Zhang. Licensed under Creative Commons
Attribution 4.0 International (CC BY 4.0),
https://creativecommons.org/licenses/by/4.0/.
