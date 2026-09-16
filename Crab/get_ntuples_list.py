#!/usr/bin/env python3

import subprocess
from pathlib import Path


# List of CRAB task directories
crab_dirs = [
    "./crab_muMuGamma_11Sep26_ParkingDoubleMuonLowMass0Run2025C-PromptReco-v1",
    # add more here...
]


output_file = "muMuGammaTree_ntuples_25_test.txt"

ntuples = []

for crab_dir in crab_dirs:

    print(f"Processing: {crab_dir}")

    result = subprocess.run(
        ["crab", "getoutput", "-d", crab_dir, "--dump"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True,
    )

    # Keep only ROOT files from the CRAB output
    files = [
        line.strip("1) PFN: ").split(" ")[0]
        for line in result.stdout.splitlines()
        if "PFN" in line
    ]

    print(f"  Found {len(files)} ROOT files")

    ntuples.extend(files)


# Remove duplicates while preserving order
ntuples = list(dict.fromkeys(ntuples))

print()
print(f"Total ROOT files: {len(ntuples)}")

with open(output_file, "w") as f:
    f.write("\n".join(ntuples) + "\n")

print(f"Written to: {output_file}")
