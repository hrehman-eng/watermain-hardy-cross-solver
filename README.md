# Water Distribution Network Solver — Hardy-Cross Method

A MATLAB + Excel tool that solves flow distribution and residual pressure in a looped municipal watermain network using the classic Hardy-Cross iterative method and the Hazen-Williams head loss equation. Built on a synthetic 7-node, 8-pipe subdivision network ("Blackwood Meadows," a fictional Waterloo, ON layout) and checked against Ontario MECP design criteria.

## What it does

- Solves pipe flows and node pressures for a looped water distribution network
- Checks results against Ontario MECP design criteria:
  - Minimum residual pressure: 275 kPa (average-day demand), 140 kPa (fire flow)
  - Maximum design velocity: 2.5 m/s
- Runs three load cases:
  1. **Base case** — average-day demand, new pipe (Hazen-Williams C = 130)
  2. **Fire flow case** — average-day demand + FUS minimum required fire flow (66.7 L/s) at the hydraulically most remote node
  3. **Sensitivity case** — average-day demand, aged pipe (C = 100), to show how pressure margin erodes over time
- Outputs convergence history and pressure-margin comparison plots

## Key finding

Under fire-flow conditions, the cul-de-sac branch pipe (150 mm, serving the network's most remote node) exceeds the 2.5 m/s velocity limit at 3.94 m/s — a real design flag showing the branch would need to be upsized for fire protection, even though it comfortably handles average-day demand. Pressure at that node stays well above the fire-flow minimum, so the constraint is pipe velocity, not pressure.

Aged pipe (C = 100 vs. 130) barely affects pressure at these flow rates (~1 kPa difference) — a useful illustration that fire flow, not long-term pipe aging, governs this network's design.

## Files

| File | Description |
|---|---|
| `Watermain_Network_Input1.xlsx` | Input/output workbook — network geometry, demands, design criteria, and results (populated by the MATLAB script) |
| `hardy_cross_solver1.m` | MATLAB solver — reads the workbook, runs all three load cases, writes results back, and generates plots |

## How to run

1. Download both files into the same folder.
2. Open `hardy_cross_solver1.m` in MATLAB (base MATLAB only — no additional toolboxes required) and run it.
3. Results are written to the `Results`, `Results_FireFlow`, and `Results_Sensitivity` sheets of the workbook.
4. Two plots are saved as PNGs: `convergence_history.png` and `pressure_margin_comparison.png`.

To change the network (topology, demands, pipe sizing), edit the yellow-shaded cells in `Node_Data` and `Pipe_Data` — see the `Instructions` sheet for the network diagram and design-criteria reference. Note: the solver's initial flow guess (in `hardy_cross_solver1.m`) must satisfy mass balance at every node before the loop-correction iteration starts; if you change demands or topology, re-derive it by hand.

## Method

The Hardy-Cross method solves looped pipe networks by:
1. Assuming an initial flow distribution that satisfies continuity (mass balance) at every node
2. Computing head loss around each loop using Hazen-Williams
3. Applying a flow correction to every pipe in a loop until head loss around that loop converges to zero
4. Propagating node pressures outward from a fixed-head source via breadth-first search

## Limitations

This is a portfolio/teaching model, not a substitute for professional network modelling software (e.g., EPANET, Bentley WaterGEMS), which solve the full nodal system simultaneously via the global gradient algorithm rather than loop-by-loop. The network is synthetic (not tied to real GIS/SCADA data), demands are static snapshots rather than a diurnal demand pattern, and minor losses (fittings, valves) are neglected.

## Author

Hamza Rehman — [github.com/hrehman-eng](https://github.com/hrehman-eng) · [LinkedIn](https://www.linkedin.com/in/hrehman07)
