#!/usr/bin/env python3
"""
Setup step for fmd_moisture regtest.
Creates a minimal FARSITE .fmd file (moisture.fmd) with diurnal variation
for FM4: dry morning, moist evening.
"""
from pathlib import Path

# Simple 2-entry diurnal cycle (morning / afternoon)
fmd_content = """\
# FARSITE fuel moisture schedule for fmd_moisture regtest
# Format: MONTH DAY HOUR PRECIP NMODELS
#         MODEL 1hr% 10hr% 100hr% LHrb% LWdy%
7 1  800 0 1
4  8 10 15 90 120
7 1 1400 0 1
4 12 14 18 95 125
7 1 2000 0 1
4  6  9 13 88 118
"""

out = Path("moisture.fmd")
out.write_text(fmd_content)
print(f"Created {out}")
