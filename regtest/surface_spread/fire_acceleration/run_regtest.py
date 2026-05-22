#!/usr/bin/env python3
"""
Fire Acceleration Regression Test
Validates size-based and FARSITE temporal acceleration models
"""

import sys
import os
import subprocess

def run_test(inputs_file):
    """Run wildfire simulation with given inputs file"""
    cmd = ['./wildfire_levelset', inputs_file]
    print(f"\n{'='*70}")
    print(f"Running: {inputs_file}")
    print(f"{'='*70}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("FAILED:", inputs_file)
        print(result.stdout)
        print(result.stderr)
        return False
    print("SUCCESS:", inputs_file)
    return True

def main():
    """Main test runner"""
    tests = [
        'inputs.size_based',
        'inputs.temporal_point',
        'inputs.wind_lag'
    ]
    
    # Check if wildfire_levelset executable exists
    if not os.path.exists('./wildfire_levelset'):
        print("ERROR: wildfire_levelset executable not found")
        print("Please build the code first: mkdir build && cd build && cmake .. && make")
        return 1
    
    passed = 0
    failed = 0
    
    for test in tests:
        if run_test(test):
            passed += 1
        else:
            failed += 1
    
    print(f"\n{'='*70}")
    print(f"Fire Acceleration Regression Test Results")
    print(f"{'='*70}")
    print(f"Passed: {passed}/{len(tests)}")
    print(f"Failed: {failed}/{len(tests)}")
    print(f"{'='*70}\n")
    
    return 0 if failed == 0 else 1

if __name__ == '__main__':
    sys.exit(main())
