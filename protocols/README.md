# Protocols

Each `DT*-PREREG.md` file is the written protocol of one computation, frozen before that
computation drew its confirmatory random numbers. A protocol fixes the code, the random
roots, the sample sizes, the decision rule and, where stated, predictions made before the run.

Tiers:

- **CONFIRM**: a confirmatory run on fresh roots. The paper's certified numbers come from
  these runs (Appendix B of the paper lists which run gives which number).
- **DEV**: an exploratory or development run. Its numbers are never used as evidence.
- **EXACT**: a computation without random numbers.

## What was removed from the working copies

These files are cleaned copies of the protocols as frozen. We removed identifiers of review
rounds, gate and receipt identifiers, references to project-management notes, quotes from
private messages, and three sections of `DT13-PREREG.md` that record process only (the
request that started the step, the plan for updating notes after the result, and the table
of review dispositions). The header of each file was rewritten to keep only the freeze time,
the tier, the roots and the code hashes. No design element, rule, threshold, seed, sample
size, code hash or prediction was changed. The original working copies remain with the author.

Some protocols refer to working notes that are not published, such as `CONVENTIONS.md`
(conventions CV-1 to CV-10, including the list of random roots in CV-6), `DT7-UNIF.md` and
`DT10-EDGE.md`. The results that the paper uses from those notes are proved in the paper.

## Clock times

The times typed into the protocols for DT0c, DT0d, DT1a, DT1b, DT2, DT2b, DT2ML and DT2ML2
were estimates that ran 1 to 5 hours ahead of the real clock. The headers of those files
give the true times instead: the last modification time of the protocol file and of the
result file of its run. In every case the protocol was last modified before its result file.
From DT3c on, every time was read from the system clock (`date -u`).
