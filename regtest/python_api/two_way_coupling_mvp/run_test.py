#!/usr/bin/env python3
"""
run_regression_test.sh wrapper for two-way coupling MVP test

This script is used by CMake ctest to run the regression test.
"""

import sys
import os

# Change to test directory
test_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(test_dir)

# Add parent directories to path
sys.path.insert(0, os.path.join(test_dir, '..', '..', '..', 'python'))

# Import and run test
try:
    from test_two_way_coupling_mvp import main
    sys.exit(main())
except ImportError as e:
    print(f"ERROR: Could not import test module: {e}")
    sys.exit(1)
