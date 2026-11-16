# GB-02 — E5: the two sealed forecasts of the twist energy, looked up (2026-11-16)

Module `mtft.maass.glueballs` (`e5`, `identify_column`, `transcription_check`); tests `tests/test_gb02_e5_sealed_forecasts.py`. This is the
experiment that GB-01 left sealed (`GB01_GLUEBALL_COMPARISON.md`, section 5). Class: PHENO.

**The answer.** One of the two forecasts could be tested, and it passed at 1.4σ. The other has no lattice value. By the seal's own words,
written beforehand, a pass is weak support and no more; this one is weaker than that sounds, for the reasons in section 4. H24's verdict
does not change: it fails as written (GB-01). The check of the lab's transcription found the six numbers right and one spin assignment
missing, which does not help H24.

## 1. What was sealed, and when

The lab's record MP20 registered two numbers at 08:21 New York time on 14 November 2026, before any lookup: with a constant energy μ on every
η-twisted form, a second odd C = −1 state at 2.87(7) and a second even one at 3.24(7), in units of the 0⁺⁺. The criterion was fixed at 06:54:50
on 16 November and hash-sealed as GB-S1 (sha256 `e75106afbfc3e3bd92c29ebcfeea81db48a0c67ab715c50fdc29dad1735500e7`): the second odd state is the first excited
state of the lattice's A₁⁻⁻ channel, the second even one that of A₁⁺⁻, in the continuum table of Athenodorou and Teper (JHEP 11 (2020) 172,
arXiv:2007.06422); σ = |forecast − lattice ratio|/√(0.0715² + e²); PASS to 2σ, TENSION to 3σ, FAIL beyond; no excited state listed: NO DATA; a fail
in either channel closes the description; a pass is weak support.

How a table is read was written down at 09:34:41 on 16 November, before the paper was opened (the build record of mtft 1.1.0): which table,
what "first excited" means, that a value at one lattice spacing only is NO DATA, how errors and flags are treated, two readings that must
agree digit for digit, and that the six transcribed values are compared at the same time. The seal still verifies; its outcome fields,
which are outside its hash, are now filled.

## 2. What the paper's table says

Read from arXiv:2007.06422v2, the PDF's text, by two separate automated extractions that agree digit for digit. The text of the tables is
flattened: the numbers of a row in order, with nothing to show the columns of a short row. **I did not see the PDF's layout.**

Table 14, "Continuum limit of glueball masses in units of the string tension, for all representations, R …", columns P=+,C=+; P=−,C=+;
P=+,C=−; P=−,C=−. The A₁ block:

    3.405(21)  5.276(45)  9.32(28)  9.78(46)
    5.855(41)  7.29(13)   10.39(37)
    7.515(50)  9.18(26)
    7.38(11)   9.37(22)

- **The six values of the lab's transcription are right**, digit for digit, table numbers included (`transcription_check()`): all six are in
  Table 14, and Table 17 has the four J = 0 values.
- **One excited C = −1 state of A₁ is listed, 10.39(37).** Its column is not in the flattened text. It is identified as P=−,C=−, the second
  A₁⁻⁻ state, by two routes (`identify_column()`), and that is an identification, not a reading:
  1. The paper's continuum extrapolation redone from its own Tables 1 and 8 to 13 (a straight line in a²σ). Controls: 3.403(21), χ² = 4.2 on 4 for the
     paper's 3.405(21), and 9.318(279), χ² = 1.6 on 3 for its 9.32(28). The P=−,C=− states listed first at the four finest spacings give
     **10.389(364), χ² = 1.5 on 2**; those listed second give 9.934(457), χ² = 3.6 on 2, beside the paper's 9.78(46). Under P=+,C=− the second row gives
     11.473(209), χ² = 21.3 on 2 and the second-lightest states 10.407(311), χ² = 6.3 on 2: fits the paper would have starred, and 10.39(37) has no star.
  2. The same authors' later SU(3) table (arXiv:2106.00364, Table 24), read from HTML that keeps the cells: the A₁ block repeats the three
     other columns digit for digit and has a second entry only under P=−,C=−, where the pair is now 9.93(49) and 10.03(47).
- A₁⁺⁻ has excited states at single spacings and none in the continuum table. By the reading rules that is NO DATA.

## 3. The sealed rule

`evaluate_open(excited_odd=(10.39, 0.37))`:

| channel | forecast | lattice | σ | verdict |
|---|---|---|---|---|
| A₁⁻⁻ first excited (second odd twisted state) | 2.8658 ± 0.0715 | 10.39(37)/3.405(21) = 3.051(110) | 1.41 | **PASS** |
| A₁⁺⁻ first excited (second even twisted state) | 3.2399 ± 0.0715 | none in the continuum table | | **NO DATA** |
| overall | | | | PASS (one channel) |

If 10.39(37) stood in the other column, the even forecast would pass at 1.43σ and the odd one would have no data: the overall verdict
is the same for either column.

## 4. What the pass is worth

- **Little power.** At the lattice's error any first excited A₁⁻⁻ state from 8.87 to 10.65 passes, and it cannot lie below the ground state,
  9.78(46). Only a state above 10.65 would have failed.
- **A nearly degenerate, noisy pair.** 9.78(46) and 10.39(37) in the sealed source; 9.93(49) and 10.03(47) in the authors' later table, where the
  rule would give 0.51σ. That is outside the seal, which names the 2020 paper.
- **One fitted number.** μ was fitted to the two C = −1 ground states; the excited state is the one point it did not see.
- **The source gives the A₁⁺⁻ ground state spin 4.** Its Table 16: "4⁺⁻ gs ∼** A₁⁺⁻gs + E⁺⁻ex1 + T₁⁺⁻ex3 + T₂⁺⁻ex2", two stars for "significant
  uncertainty"; its Table 17 lists that state at 9.02(10)**, which is the mean of the four members (9.02); its text doubts a 0⁺⁻ there. No A₁⁻⁻
  state is assigned a spin. The lab's record MP15 wrote "spins not identified" for both, which for A₁⁺⁻ left out the paper's own assignment.
  Under the dictionary's D2 (weight 0 is J = 0 only) the A₁⁺⁻ ground state is then not a state of the comparison. Two consequences, neither in
  H24's favour: the even twisted form at 1.92 times the edge (GB-01, section 3, fifth row) has no lattice partner at all; and μ was fitted to
  two ground states of which one probably should not have been used, so "one number for two measurements" is one number for one.

So the constant twist energy is not closed, and it is not supported beyond a 1.4σ agreement that almost any value would have given.

## 5. Scored

Written at 06:54:50, blind: the source lists an excited state in both channels, 0.5: it lists one. Given a value, the odd forecast passes, 0.15:
it passed, so that estimate was wrong, and the reason is that the gap to the first excited state was taken from the light channels (2.0 to 2.5
in √σ), where at 9 to 10 the paper's levels lie 0.3 to 1.5 apart. The even forecast passes, 0.35: not scored.

## 6. What is still to do

- **A look at the PDF.** Under which heading 10.39(37) stands in Table 14 can be read off the page in a moment; it was not read here.
- A dictionary that takes the source's spins as they are (A₁⁺⁻ ground state as 4⁺⁻) has not been written or tested. It would have to be
  registered first (`register_forecast`).

## 7. Ledger

GB-6 (the outcome; MODEL, caveat) and GB-7 (the spin of the A₁⁺⁻ ground state; STRUCTURAL, active), wave GB-02 in piece IV. GB-5 is kept as
written, with a pointer. The negative N24 is unchanged: a pass, however weak, is not a negative, and H24's failure does not rest on it.

## 8. Added at 11:08 on 16 November 2026: the column, read

Roger Tano read the PDF's page 43 (a screenshot of Table 14): in the A₁ block 9.32(28) stands alone under P=+,C=−, and 10.39(37) stands
beneath 9.78(46) under P=−,C=−. The same file's text with its layout kept shows the same (arXiv:2007.06422v2, sha256 f0eb6406…de23d0). So the
identification of section 2 is confirmed by reading, and it was the odd forecast that was tested. Everything the package holds from the
paper for this experiment was then compared with the layout, column by column, and agrees digit for digit; one A₁⁺⁻ state at a single
spacing, 1.77(4) at β = 6.235, was added to a list that no fit uses. Sections 1 to 7 stand as written before this reading.

