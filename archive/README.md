# Archive and repository history

There is one active manuscript: [`../manuscript/Magenta_Gap.tex`](../manuscript/Magenta_Gap.tex). This directory holds historical material so it cannot be mistaken for that working text.

| Archive | What is preserved |
| --- | --- |
| [`july-2026/`](july-2026/) | Original 51-page PDF and recovered source. The source and PDF are byte-identical to their former `manuscript/jmp/` files. |
| [`legacy-pipeline/`](legacy-pipeline/) | Nine superseded calculation scripts, unchanged, plus the historical root README text. |

The corrected `recompute/` package, exact certificate, active algebra checks, and dated research proofs remain in their existing locations. This cleanup does not alter scientific content or erase Git history.

## Versions and branch history

- [September import at `cb997e8`](https://github.com/AlphaOmega-Institute/magenta-gap/tree/cb997e8b03e7e5c3955da07a641297eecbff1190/manuscript).
- [Reviewed Optica-format revision at `51b37d5`](https://github.com/AlphaOmega-Institute/magenta-gap/tree/51b37d506dee66de9228a1ee2c0e81aa7a193488/manuscript).
- [Reviewed venue-neutral manuscript and pre-cleanup layout at `44ee356`](https://github.com/AlphaOmega-Institute/magenta-gap/tree/44ee356a315fdd20b3106f595b7ebfb5656e3d01).

The following work-branch tips are all ancestors of `main` at `44ee356` (or equal to it). Removing these completed branch names does not remove their commits or files from history.

| Completed work branch | Preserved tip |
| --- | --- |
| `claude/corrected-bound-calculation` | [`a0282ae23d11c102550c88a9f44a8e2a99f99211`](https://github.com/AlphaOmega-Institute/magenta-gap/commit/a0282ae23d11c102550c88a9f44a8e2a99f99211) |
| `claude/figure-rebuild` | [`e07eb71eac24a7773838246590b395cc7245e505`](https://github.com/AlphaOmega-Institute/magenta-gap/commit/e07eb71eac24a7773838246590b395cc7245e505) |
| `codex/restore-conceptual-argument` | [`44ee356a315fdd20b3106f595b7ebfb5656e3d01`](https://github.com/AlphaOmega-Institute/magenta-gap/commit/44ee356a315fdd20b3106f595b7ebfb5656e3d01) |
| `codex/september-consolidation` | [`51b37d506dee66de9228a1ee2c0e81aa7a193488`](https://github.com/AlphaOmega-Institute/magenta-gap/commit/51b37d506dee66de9228a1ee2c0e81aa7a193488) |

A local branch can be recreated after fetching `main`, for example:

```sh
git fetch origin main
git branch recovered-september 51b37d506dee66de9228a1ee2c0e81aa7a193488
```

**Still unmerged:** `claude/independent-colour-checks`, inspected at [`a3a947c0c27d76dcb4b5e57d7b1329791fc771e7`](https://github.com/AlphaOmega-Institute/magenta-gap/tree/a3a947c0c27d76dcb4b5e57d7b1329791fc771e7). It contains unique reproduction scripts and exploratory analyses. Its branch is retained; preservation is not a claim that those analyses establish new manuscript results.
