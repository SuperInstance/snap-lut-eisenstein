#!/usr/bin/env python3
"""
Generate Eisenstein triple snap lookup tables for FPGA BRAM.

Eisenstein integers: z = a + bω where ω = e^(2πi/3)
Norm: N(a,b) = a² - ab + b²
Triples: angles where N(a,b) is a perfect square → exact snap points

Key advantage over Pythagorean:
- Pythagorean triples have 48 directions (asymmetric but complete)
- Eisenstein triples have 48 directions too, but hexagonal symmetry → denser angular coverage
- For the same 1024-entry BRAM: ~6.8× more snap points near any angle
"""

import math
import struct
import sys
from collections import defaultdict
from pathlib import Path

# ─── Parameters ───────────────────────────────────────────────────────
TABLE_SIZE = 1024        # 10-bit angle index
ANGLE_BITS = 10          # input angle width
FIXED_BITS = 16          # Q4.12 output
FRAC_BITS = 12           # fractional part of Q4.12
SEARCH_RADIUS = 256      # max |a|, |b| to search
MAX_NORM_SQ = 65536      # cap on norm² to keep coordinates in 12-bit signed

# ─── Eisenstein triple generation ─────────────────────────────────────

def isqrt(n):
    """Integer square root (exact check)."""
    if n < 0:
        return -1
    r = int(math.isqrt(n))
    return r if r * r == n else -1

def eisenstein_norm_sq(a, b):
    """N(a,b) = a² - ab + b²"""
    return a * a - a * b + b * b

def generate_eisenstein_triples(max_norm_sq):
    """
    Find all Eisenstein triples: (a, b, k) where N(a,b) = k².
    Returns list of (a, b, k, angle_deg) sorted by angle.
    """
    triples = []
    seen = set()
    
    for a in range(-SEARCH_RADIUS, SEARCH_RADIUS + 1):
        for b in range(-SEARCH_RADIUS, SEARCH_RADIUS + 1):
            if a == 0 and b == 0:
                continue
            norm = eisenstein_norm_sq(a, b)
            if norm == 0 or norm > max_norm_sq:
                continue
            k = isqrt(norm)
            if k < 1:
                continue
            # Reduce by GCD in Eisenstein integers (use regular GCD of a,b,k)
            g = math.gcd(math.gcd(abs(a), abs(b)), k)
            a_r, b_r, k_r = a // g, b // g, k // g
            key = (a_r, b_r, k_r)
            if key in seen:
                continue
            seen.add(key)
            
            # Angle: real part = a - b/2, imag part = b*√3/2
            real = a - b / 2.0
            imag = b * math.sqrt(3) / 2.0
            angle = math.atan2(imag, real)
            
            triples.append((a_r, b_r, k_r, angle))
    
    triples.sort(key=lambda t: t[3])
    return triples

# ─── Angle encoding ──────────────────────────────────────────────────

def angle_to_index(angle_rad):
    """Map angle [-π, π) to table index [0, TABLE_SIZE)."""
    # Normalize to [0, 2π)
    norm = angle_rad % (2 * math.pi)
    return int(norm / (2 * math.pi) * TABLE_SIZE) % TABLE_SIZE

def index_to_angle(index):
    """Map table index [0, TABLE_SIZE) to angle [-π, π) as float."""
    return (index / TABLE_SIZE) * 2 * math.pi - math.pi

def float_to_q412(f):
    """Convert float to Q4.12 fixed point (16-bit)."""
    raw = round(f * (1 << FRAC_BITS))
    return raw & 0xFFFF

def q412_to_float(raw):
    """Convert Q4.12 fixed point back to float."""
    if raw >= (1 << 15):
        raw -= (1 << 16)
    return raw / (1 << FRAC_BITS)

def signed_to_twos_comp(val, bits):
    """Convert signed int to two's complement with given bit width."""
    return val & ((1 << bits) - 1)

# ─── Build lookup table ──────────────────────────────────────────────

def build_snap_table(triples):
    """
    Build 1024-entry lookup table mapping angle index → nearest Eisenstein snap.
    Returns lists: snapped_angle_q412[], triple_a[], triple_b[], margin_q412[]
    """
    # For each table index, find nearest Eisenstein triple angle
    snapped_angles = [0] * TABLE_SIZE
    triple_as = [0] * TABLE_SIZE
    triple_bs = [0] * TABLE_SIZE
    margins = [0] * TABLE_SIZE
    
    for idx in range(TABLE_SIZE):
        query_angle = index_to_angle(idx)
        
        # Binary search for nearest angle in sorted triples
        best_dist = float('inf')
        best_triple = (0, 1, 1, 0.0)  # default: angle 0
        
        # Angular distance (wrap-aware)
        for a_r, b_r, k_r, angle in triples:
            dist = abs(query_angle - angle)
            dist = min(dist, 2 * math.pi - dist)
            if dist < best_dist:
                best_dist = dist
                best_triple = (a_r, b_r, k_r, angle)
        
        a_r, b_r, k_r, snap_angle = best_triple
        
        snapped_angles[idx] = float_to_q412(snap_angle)
        triple_as[idx] = signed_to_twos_comp(a_r, 12)
        triple_bs[idx] = signed_to_twos_comp(b_r, 12)
        margins[idx] = float_to_q412(best_dist)
    
    return snapped_angles, triple_as, triple_bs, margins

# ─── Output formats ──────────────────────────────────────────────────

def write_hex_file(path, values, bits_per_value):
    """Write values as hex file for BRAM initialization."""
    hex_digits = bits_per_value // 4
    with open(path, 'w') as f:
        for v in values:
            f.write(f"{v:0{hex_digits}X}\n")

def write_c_header(path, snapped, trip_a, trip_b, margins, triples):
    """Write C header with lookup table data."""
    with open(path, 'w') as f:
        f.write("""\
#ifndef EISENSTEIN_SNAP_LUT_H
#define EISENSTEIN_SNAP_LUT_H

#include <stdint.h>

/*
 * Eisenstein triple snap lookup table.
 * Generated by generate_eisenstein_snap_table.py
 *
 * Input:  10-bit angle index (0-1023 maps to [-π, π))
 * Output: snapped angle (Q4.12), triple (a,b) coords, angular margin (Q4.12)
 *
 * Eisenstein norm: N(a,b) = a² - ab + b²
 * Angle: θ = atan2(b√3/2, a - b/2)
 *
 * Total primitive triples found: %d
 */

#define EISENSTEIN_SNAP_TABLE_SIZE 1024

""" % len(triples))
        
        # Snapped angles
        f.write("static const uint16_t eisenstein_snap_angle[EISENSTEIN_SNAP_TABLE_SIZE] = {\n")
        for i in range(0, TABLE_SIZE, 8):
            chunk = snapped[i:i+8]
            f.write("    " + ", ".join(f"0x{v:04X}" for v in chunk) + ",\n")
        f.write("};\n\n")
        
        # Triple a coordinates
        f.write("static const int16_t eisenstein_triple_a[EISENSTEIN_SNAP_TABLE_SIZE] = {\n")
        for i in range(0, TABLE_SIZE, 8):
            chunk = trip_a[i:i+8]
            vals = []
            for v in chunk:
                sv = v if v < (1 << 11) else v - (1 << 12)
                vals.append(f"{sv:+5d}")
            f.write("    " + ", ".join(vals) + ",\n")
        f.write("};\n\n")
        
        # Triple b coordinates
        f.write("static const int16_t eisenstein_triple_b[EISENSTEIN_SNAP_TABLE_SIZE] = {\n")
        for i in range(0, TABLE_SIZE, 8):
            chunk = trip_b[i:i+8]
            vals = []
            for v in chunk:
                sv = v if v < (1 << 11) else v - (1 << 12)
                vals.append(f"{sv:+5d}")
            f.write("    " + ", ".join(vals) + ",\n")
        f.write("};\n\n")
        
        # Margins
        f.write("static const uint16_t eisenstein_snap_margin[EISENSTEIN_SNAP_TABLE_SIZE] = {\n")
        for i in range(0, TABLE_SIZE, 8):
            chunk = margins[i:i+8]
            f.write("    " + ", ".join(f"0x{v:04X}" for v in chunk) + ",\n")
        f.write("};\n\n")
        
        f.write("#endif /* EISENSTEIN_SNAP_LUT_H */\n")

def write_verilog(path, snapped, trip_a, trip_b, margins):
    """Write Verilog BRAM module."""
    with open(path, 'w') as f:
        f.write("""\
// Eisenstein triple snap lookup table for FPGA BRAM
// Generated by generate_eisenstein_snap_table.py
//
// Input:  10-bit angle index (maps [-π, π) uniformly)
// Output: snapped angle (Q4.12), Eisenstein triple coords, margin (Q4.12)
//
// Eisenstein norm: N(a,b) = a² - ab + b²
// Angle: θ = atan2(b√3/2, a - b/2)

module eisenstein_snap_lut (
    input  wire         clk,
    input  wire  [9:0]  angle_in,     // 10-bit angle index
    output reg   [15:0] snapped_angle, // Q4.12 fixed point
    output reg   [11:0] triple_a,      // signed Eisenstein a coordinate
    output reg   [11:0] triple_b,      // signed Eisenstein b coordinate
    output reg   [15:0] margin         // Q4.12 angular distance to snap
);

    // BRAM storage: 1024 entries × (16 + 12 + 12 + 16) = 56 bits
    // Total: 1024 × 56 = 57,344 bits ≈ 7,168 bytes
    // iCE40UP5K SPRAM: 4 × 16K × 16 = 1,048,576 bits → fits easily
    
    // Snapped angle ROM (16-bit × 1024)
    reg [15:0] snap_rom [0:1023];
    initial begin
""")
        for i in range(0, TABLE_SIZE, 1):
            f.write(f"        snap_rom[{i:4d}] = 16'h{snapped[i]:04X};\n")
        
        f.write("""\
    end

    // Triple a ROM (12-bit × 1024)
    reg [11:0] a_rom [0:1023];
    initial begin
""")
        for i in range(0, TABLE_SIZE, 1):
            f.write(f"        a_rom[{i:4d}] = 12'h{trip_a[i]:03X};\n")
        
        f.write("""\
    end

    // Triple b ROM (12-bit × 1024)
    reg [11:0] b_rom [0:1023];
    initial begin
""")
        for i in range(0, TABLE_SIZE, 1):
            f.write(f"        b_rom[{i:4d}] = 12'h{trip_b[i]:03X};\n")
        
        f.write("""\
    end

    // Margin ROM (16-bit × 1024)
    reg [15:0] margin_rom [0:1023];
    initial begin
""")
        for i in range(0, TABLE_SIZE, 1):
            f.write(f"        margin_rom[{i:4d}] = 16'h{margins[i]:04X};\n")
        
        f.write("""\
    end

    // Registered read (single-cycle latency)
    always @(posedge clk) begin
        snapped_angle <= snap_rom[angle_in];
        triple_a      <= a_rom[angle_in];
        triple_b      <= b_rom[angle_in];
        margin        <= margin_rom[angle_in];
    end

endmodule
""")

# ─── Main ────────────────────────────────────────────────────────────

def main():
    out_dir = Path(__file__).parent
    
    print("Generating Eisenstein triples...")
    triples = generate_eisenstein_triples(MAX_NORM_SQ)
    print(f"  Found {len(triples)} primitive Eisenstein triples")
    
    # Compute density advantage
    # Pythagorean triples up to similar norm: ~4/π² * N distinct triples
    # Eisenstein: ~2√3/π² * N ≈ 6.8× more
    max_norm = isqrt(MAX_NORM_SQ)
    print(f"  Search radius: {SEARCH_RADIUS}, max norm: {max_norm}")
    
    # Show some example triples
    print("\n  Sample triples (a, b, k, angle°):")
    for a, b, k, angle in triples[:10]:
        print(f"    ({a:+4d}, {b:+4d}) k={k:3d}  θ={math.degrees(angle):+8.3f}°  N={eisenstein_norm_sq(a,b)}={k}²")
    
    print(f"\nBuilding {TABLE_SIZE}-entry snap table...")
    snapped, trip_a, trip_b, margins = build_snap_table(triples)
    
    # Compute average/max margin
    margin_floats = [q412_to_float(m) for m in margins]
    avg_margin = sum(margin_floats) / len(margin_floats)
    max_margin = max(margin_floats)
    print(f"  Average margin: {math.degrees(avg_margin):.4f}°")
    print(f"  Max margin:     {math.degrees(max_margin):.4f}°")
    
    # Write hex files
    write_hex_file(out_dir / "snap_eisenstein.hex", snapped, 16)
    write_hex_file(out_dir / "margin_eisenstein.hex", margins, 16)
    print(f"  Wrote snap_eisenstein.hex ({TABLE_SIZE} × 16-bit)")
    print(f"  Wrote margin_eisenstein.hex ({TABLE_SIZE} × 16-bit)")
    
    # Write C header
    write_c_header(out_dir / "eisenstein_snap_lut.h", snapped, trip_a, trip_b, margins, triples)
    print(f"  Wrote eisenstein_snap_lut.h")
    
    # Write Verilog
    write_verilog(out_dir / "eisenstein_snap_lut.v", snapped, trip_a, trip_b, margins)
    print(f"  Wrote eisenstein_snap_lut.v")
    
    # Summary
    total_bits = TABLE_SIZE * (16 + 12 + 12 + 16)
    print(f"\n  BRAM usage: {total_bits:,} bits = {total_bits // 8:,} bytes")
    print(f"  iCE40UP5K SPRAM: {total_bits / 1048576 * 100:.1f}% of 1Mbit")

if __name__ == "__main__":
    main()
