# Project-Nonogram-Simulated-Annealing

This repository provides a C++ Nonogram solver using simulated annealing. The main source file is `src/nonogram_solver.cpp`.
This group project was carried out in 2024 during my master's degree at the Haute Ecole en Hainaut as part of the combinatorial optimisation course.

## Overview
The solver reads a `.pti` constraints file (row and column clues), applies a small logical pre-processing pass, and then performs simulated annealing to shift black cell groups and minimize a deviation score against the clues. Sequence comparison uses a modified Needleman–Wunsch scoring.

## Files and layout
- `src/nonogram_solver.cpp` — main C++ source with all core logic.
- `include/` — optional header directory if you split code later.
- `data/pti/` — input constraint files (e.g., `10x10.pti`).
- `output/pto/` — exported solutions in `.pto` format.

Key functions/symbols in `nonogram_solver.cpp`:
- `needlemanWunsch` — computes row/column sequence similarity scores.
- `movingRight`, `movingLeft` — neighborhood generation by shifting groups.
- `fitModif` — feasibility check for a proposed row modification.
- `MoveGroupsFirst` — utility for initial placement of groups.
- `createPTOFile` — writes the final solution as `.pto`.
- `main` — orchestrates I/O, preprocessing, annealing, and output.
- Important globals: `RowSize`, `ColSize`, `filename`, `ScoreTotal`, `ScoresCols`, `tableau2D`, `tableau2DReplica`.

## Build and run
1. Ensure the source file is saved as `src/nonogram_solver.cpp`.
2. Compile with g++ (C++17 recommended):
g++ -std=c++17 -O2 -Isrc -Iinclude -o build/nonogram_solver src/nonogram_solver.cpp

text
3. Place `.pti` inputs under `data/pti/`. Set the input path via the `filename` variable (default example: `data/pti/10x10.pti`).
4. Run the solver:
./build/nonogram_solver

text
The program will print score evolution and write a `.pto` output under `output/pto/<name>.pto` via `createPTOFile`.

## Parameters
Core penalties and annealing controls are defined near the top of the source:
- `MISMATCHPENALTY`, `GAPPENALTY`, `ERRORACCEPTANCEREDUCTOR`, `errorAcceptance`.

Tuning these parameters adjusts exploration vs. exploitation and the acceptance behavior during annealing.

## Algorithm outline
- Preprocessing: deterministically place mandatory black cells based on block-sum logic.
- Initial solution: place groups flush-left by default.
- Annealing loop:
  - Randomly select a row, a group, and a direction.
  - Propose a shift using `movingRight`/`movingLeft` on `tableau2DReplica`.
  - Validate feasibility with `fitModif`.
  - Recompute costs for affected columns using `needlemanWunsch`.
  - Acceptance: always accept improvements; otherwise accept with probability
    \( P = \exp\!\left(-\frac{\Delta E}{T}\right) \)
    as controlled by the acceptance variables/temperature schedule in code.

## Example settings and results
Tested on 10×10, 15×15, and 30×30 boards. Example parameters:
- `mismatch = -1`, `gap = -1`, `threads = 12`, initial acceptance ≈ 0.6, cooling factor ≈ 0.99.

See the accompanying report for detailed tables and figures.

## Limitations and future work
- Stronger logical preprocessing to reduce the search space.
- Smarter move heuristics to avoid unmovable groups.
- Alternative neighborhoods (bit flips) and in-place updates to avoid copying `tableau2DReplica`.
- Automated tuning for the temperature and cooling schedules.

## Notes
- To change the input file, edit the `filename` variable in `src/nonogram_solver.cpp`.
- Consider adding a small example `.pti` and the corresponding `.pto` output for quick validation by users.

## Sources

- [1] “Theory — nonogram solver,” Stevocity.me.uk, 2025. [Online]. Available : https://stevocity. me.uk/nonogram/theory
- [2] C. aux projets Wikimedia, “méthode d’optimisation,” Wikipedia.org, 08 2004. [Online]. Available : https://fr.wikipedia.org/wiki/Recuit_simul%C3%A9
- [3] lb@laurentbloch.org and W. de, “Algorithme de needleman et wunsch - [site www de laurent bloch],” Laurentbloch.net, 2019. [Online]. Available : https://www.laurentbloch.net/ MySpip3/Algorithme-de-Needleman-et-Wunsch
- [4] “Needleman-wunsch algorithm in c++ - javatpoint,” www.javatpoint.com, 2025. [Online]. Available : https://www.javatpoint.com/needleman-wunsch-algorithm-in-cpp
