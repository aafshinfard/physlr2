#!/usr/bin/env python3
"""Find the first mode after minimum in ntCard histogram output.

Faithful reimplementation of physlr find-ntcard-mode.
Usage: find-ntcard-mode.py <histogram_file>
Prints the mode count value to stdout.
"""
import sys

def find_mode(filename):
    freq_count = [int(line.rstrip().split("\t")[2]) for line in open(filename)
                  if line[0] != "k"]
    min_idx = 0
    min_val = freq_count[0]
    for idx, freq in enumerate(freq_count):
        if freq > min_val:
            min_idx = idx - 1
            break
        min_val = freq
    freq_count = freq_count[min_idx:]
    print(freq_count.index(max(freq_count)) + 1 + min_idx)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <histogram_file>", file=sys.stderr)
        sys.exit(1)
    find_mode(sys.argv[1])
