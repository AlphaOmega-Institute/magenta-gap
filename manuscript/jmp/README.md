# The Magenta Gap — frozen July 2026 archive

This directory preserves the July 2026 manuscript as a historical baseline. It is frozen and is not a second active paper. The sole working manuscript is [`../Magenta_Gap.tex`](../Magenta_Gap.tex), with its [current PDF](../Magenta_Gap.pdf) and [build notes](../README.md).

- [Magenta_Gap_JMP.tex](Magenta_Gap_JMP.tex): recovered editable source, matched to the original 51-page PDF.
- [Magenta_Gap_JMP.pdf](Magenta_Gap_JMP.pdf): the original manuscript PDF, preserved byte for byte.

Its July claims and numerical results are preserved unchanged, including claims corrected in the later revision. They should not be read as the current conclusions. The archive may inform later framing, but ongoing corrections belong only in the canonical working manuscript. The recent Hilbert-geometry discussion is not incorporated.

## Source recovery

The named-author source was recovered from a saved TeX paste dated July 16, 2026. The PDF's embedded creation timestamp is approximately five minutes later. Two differences were found and restored in the recovered source:

1. The affiliation's postal address: 3724 Laurel Street, Suite 108, New Orleans, Louisiana, USA.
2. Parentheses around the reference to Proposition 5.1 in the sentence about the equatorial ridge (PDF page 35).

No other source content was changed. This is a verified recovery, not an export of the current Overleaf project.

## Verification

The recovered source, with those two edits, compiled successfully to 51 pages. The complete extracted text matched the original PDF after normalization of whitespace, line-end hyphenation, and typographic characters. All references resolved. The three embedded TikZ diagrams matched the original diagrams; the section structure, equations, numerical table, and all 21 bibliography entries were checked.

The PDF stored here is the original, not the verification build. Verification establishes correspondence with that July PDF, not mathematical correctness or agreement with any later Overleaf edits.

## Build

The source contains all three figures as TikZ code and the full bibliography. It requires no separate figure files or BibTeX database. With a complete TeX installation providing `elsarticle` and the packages loaded by the source, run from this directory:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error Magenta_Gap_JMP.tex
```

The verification build used pdfLaTeX, `elsarticle` 3.5, and CM-Super fonts.

## Checksums at import

Imported November 11, 2026.

| File | SHA-256 |
| --- | --- |
| Original saved TeX paste, before the two edits | `036a4ccdb5618f194372a148aa7d147faccfb512cc791aad62f460697e90a28a` |
| Recovered TeX in this directory | `bb85dba3c42039732efe2971919476d66af92460abbe5f82b752f0bcabbd778d` |
| Original PDF in this directory | `4ce1bcf9e9f3c798604b24d401954f43aab7faf12f2de03da305eaccd063e8b6` |

