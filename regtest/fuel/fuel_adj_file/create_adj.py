#!/usr/bin/env python3
"""
Setup step for fuel_adj_file regtest.
Creates a minimal FARSITE .adj file (fire.adj) for FM4 with adj_factor = 1.4.
"""
from pathlib import Path

adj_content = """\
# FARSITE fuel adjustment file for fuel_adj_file regtest
# Applies 1.4x ROS scale to FM4 (chaparral)
1
4  1.4000
"""

out = Path("fire.adj")
out.write_text(adj_content)
print(f"Created {out}")
