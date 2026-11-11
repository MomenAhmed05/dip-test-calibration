# Flat shoulders break the uniform calibration of the dip test

Code, deposited outputs, Lean files and paper for

> M. Ahmed. *Flat shoulders break the uniform calibration of the dip test, and a
> constant cutoff repairs it.* 2026. The paper is [`paper/main.pdf`](paper/main.pdf).

The dip test of Hartigan and Hartigan (1985) takes its critical values from the
uniform distribution, on the grounds that among unimodal laws the uniform gives
the largest dip. The paper shows that this is false. A density that is flat over
most of its range and then rises to a narrow peak at one end makes the nominal
5 percent test reject about one sample in four, however large the sample. Over
atomless unimodal laws with bounded support and a finite density at the mode, the
worst limiting rejection rate of the standard 5 percent test lies between 28.4 and
33.2 percent. Rejecting when the dip exceeds `c_alpha / sqrt(n)`, with
`c_0.05 = 0.7314` and `c_0.01 = 0.8599`, restores the nominal level on that class
in the limit. The paper also gives certified cutoffs for small and moderate samples.

## What is here

| directory | what it holds |
| --- | --- |
| [`code/`](code/) | the Python and R scripts that produce every number in the paper |
| [`deposits/`](deposits/) | the outputs of those scripts, as JSON, CSV and text files |
| [`lean/`](lean/) | the Lean 4 files that check the deterministic steps of several proofs |
| [`paper/`](paper/) | the paper (PDF and LaTeX source) and the script that draws its figures |

Every number in the paper is read from a deposit. The figures are drawn by
[`paper/figs/make_figs.py`](paper/figs/make_figs.py), which reads deposits only and
estimates nothing itself.

### Where each result comes from

| paper location | deposits | main scripts |
| --- | --- | --- |
| `F_M` and the cutoffs `c_alpha` (Sections 2 and 4.2, Table 1) | `dt8/`, `dt4/` | `fm_cert.py`, `fm_cert_check.py`, `fm_cert_check2.py`, `fm_table.py` |
| brackets for the uniform quantile `qbar_alpha` (Section 4.2, Table 1) | `dt8/`, `qbar/`, `dt13/` | `qbar_cert.py`, `dt4_qbar.py`, `dt13_qbar_analysis.py` |
| mechanism and block formula (Sections 3 and 4.4, Figures 1 to 3) | `dt1a/`, `dt1b/`, `dt2/`, `dt2b/`, `dt3c/` | `dt1a_mechanism.py`, `dt1b_confirm.py`, `dt2_farm.py`, `dt2b_pinning.py`, `dt3c_*.py` |
| edge statement at `t'` (Theorem 7, Section 5.4, Appendix B) | `dt16/` | `dt16_neprime_arb.py`, `dt16_hgeom_cert.py`, `dt16_pair_cert.py`, `dt16_pair_cert_summary.py`, `dt16_rebuild.py` |
| edge statement at `c_alpha` (Section 5.4) | `dt10/`, `dt11/` | `dt10_*.py`, `dt11_scalar_check.py` |
| size map for `4 <= n <= 1000` (Section 6.1, Figure 5) | `dt13/` | `dt13_confirm.py`, `dt13_analysis.py` |
| `n = 4` (Section 6.3) | `dt13/`, `dt14/` | `dt13_n4_*.py`, `dt14_n4_*.py`, `dt14_he2.py` |
| certified cutoffs for small and moderate `n` (Section 6.2, Table 2, Figure 6) | `dt14/`, `dt9/`, `dt8blk/` | `dt14_r1_n4.py`, `dt14_r1_small.py`, `dt14_r1_small_analysis.py`, `dt9_*.py`, `dt8_block_*.py` |
| size and power study (Section 7.1, Figure 7) | `dt5/` | `dt5_power.py`, `dt5_analysis.py`, `dt5_r5.py` |
| comparison with the multiscale tests of Dümbgen and Walther and Rufibach and Walther (Section 7.1) | `dt14/mh_power*` | `code/r/mh_*.R` |
| p-value mixtures and curved pile-ups (Section 7.2) | `dt0c/`, `dt0d/` | `dt0c_screen.py`, `dt0d_screen.py` |
| UniDip, SkinnyDip and DipMeans (Section 7.2) | `o5ml/`, `o5gap/` | `o5_ml.py`, `o5_gap.py`, `o5_analysis.py`, `o5_hx.py` |
| Cheng and Hall calibration (Section 7.2) | `dt5b/` | `dt5b_ch98.py`, `ch98.py`, `ch98_fidelity.py` |
| learned calibrator (Section 7.2) | `o8/`, `o8_redteam/` | `o8_calibrator.py`, `o8_audit.py`, `o8_redteam.py` |
| defects found after freezing (Appendix B) | `defects/` | |

Appendix B of the paper lists the confirmatory runs. Each one followed a protocol
that was fixed before it started, used fresh seeds, and was run once. The seed
roots are written inside each script and each deposit.

### What is not here

- **Raw per-shard simulation arrays.** The `.npy` and `.npz` files that the
  simulations write, about 2.9 GB in total, are not in this repository. The
  scripts regenerate them from the recorded seed roots. The summaries that the
  paper reads are included. The one raw array that a figure reads,
  `deposits/dt1b/dips_f0_n12800.npy`, is included.
- **Third-party texts.** `fm_table.py --build` parses the table of Balabdaoui and
  Pitman (2011) from a local copy of their paper, which is not redistributed. Its
  output, `deposits/dt4/bf_table_full.json`, is included.
- **Design notes.** Some docstrings refer to design notes (files such as
  `DT5-PREREG.md`). These are working records and are not part of this repository.

## Reproducing

Run every script from the repository root, for example `python code/dt13_analysis.py`.
The scripts find their deposits relative to their own location. Many scripts accept
`--quick` for a small smoke run that writes to a separate `_quick` folder.

Python 3.12 with

```
numpy 2.3.5   scipy 1.17.1   diptest 0.12.0   python-flint 0.9.0   mpmath 1.3.0
sympy 1.14.0  numba 0.63.1   matplotlib 3.10.9   torch 2.12.1 (CPU)
```

The UniDip, SkinnyDip and DipMeans scripts (`o5_*.py`) need ClustPy 0.0.3, which
was installed in a separate virtual environment with its own torch (2.15.0, CPU)
and scikit-learn 1.8.0.

The R scripts in `code/r/` need R 4.6 with `diptest` 0.77-2 and `modehunt` 1.0.8.
Run them from the repository root. If the packages are in a separate library, set
`DIPTEST_RLIB` to its path.

To redraw the figures:

```bash
python -I paper/figs/make_figs.py
```

## Lean

The files in [`lean/DipTest/`](lean/DipTest/) check the deterministic steps of
several arguments in Lean 4 (v4.35.0-rc3) with Mathlib. The probability arguments,
the published laws used as inputs and the Monte Carlo certificates lie outside Lean;
where a published result is used, it enters as an explicit hypothesis. Each file
is self-contained and imports only Mathlib.

```bash
cd lean
lake exe cache get
lake build
```

The Mathlib commit is pinned in `lakefile.toml` and `lake-manifest.json`. Each file
can also be checked alone with `lake env lean DipTest/<File>.lean`.

## Licence

MIT. See [`LICENSE`](LICENSE).
