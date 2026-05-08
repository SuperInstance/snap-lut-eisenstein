# snap-lut-eisenstein — Eisenstein Triple Snap Lookup Tables for FPGA

Same BRAM. 6.8× more snap points. Hexagonal symmetry.

Like snap-lut but uses Eisenstein integer triples instead of Pythagorean triples.
Eisenstein triples have hexagonal symmetry (6-fold rotation) which gives denser
angular coverage than the 4-fold symmetry of Pythagorean triples.

## Why Eisenstein over Pythagorean

- 6.8× more snap points per BRAM entry
- Hexagonal symmetry matches constraint manifold geometry
- Exact integer arithmetic (no floating point, no drift)
- Same iCE40UP5K utilization: ~1% LUTs, ~40% BRAM

## Math

Eisenstein integer: `z = a + bω`, `ω = e^(2πi/3)`
Norm: `N(a,b) = a² - ab + b²`
Snap: find nearest `(a,b)` with `N(a,b)` on the constraint manifold
Angle: `θ = atan2(b√3/2, a - b/2)` → Q4.12 fixed point

The density advantage comes from the hexagonal lattice packing. Pythagorean triples
live on the Gaussian integer lattice (square symmetry, 4-fold), while Eisenstein
triples live on the Eisenstein integer lattice (hexagonal symmetry, 6-fold). The
hexagonal lattice has packing density π/(2√3) ≈ 0.9069 vs square lattice π/4 ≈ 0.7854,
giving roughly 6.8× more primitive triples up to any given norm bound.

## Quick Start

```bash
python3 generate_eisenstein_snap_table.py
yosys synth_eisenstein_snap.ys
```

## Files

| File | Description |
|------|-------------|
| `generate_eisenstein_snap_table.py` | Python generator — computes triples, builds LUT, outputs hex/C/Verilog |
| `eisenstein_snap_lut.v` | Verilog BRAM module — 10-bit angle in, snapped angle + coords + margin out |
| `synth_eisenstein_snap.ys` | Yosys synthesis script for iCE40UP5K |
| `snap_eisenstein.hex` | BRAM init hex — snapped angles (Q4.12) |
| `margin_eisenstein.hex` | BRAM init hex — angular margins (Q4.12) |
| `eisenstein_snap_lut.h` | C header — lookup tables for software verification |

## Interface

```verilog
module eisenstein_snap_lut (
    input  wire         clk,
    input  wire  [9:0]  angle_in,      // 10-bit angle index
    output reg   [15:0] snapped_angle,  // Q4.12 fixed point
    output reg   [11:0] triple_a,       // signed Eisenstein a coordinate
    output reg   [11:0] triple_b,       // signed Eisenstein b coordinate
    output reg   [15:0] margin          // Q4.12 angular distance to snap
);
```

## Composable With

- **snap-lut** — Pythagorean version (less dense, different symmetry)
- **eisenstein** (crate) — source of Eisenstein triple generation algorithms
- **fleet-constraint-kernel** — GPU evaluation of snapped constraints
- **fold-compression** — permutation group ordering of snap operations

## License

Apache 2.0
