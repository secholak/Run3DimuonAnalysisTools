#!/usr/bin/env python3

import subprocess
from pathlib import Path


# List of CRAB task directories
crab_dirs = [


    # 26C is a low PU, skipping
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass0Run2026A-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass0Run2026B-PromptReco-v1",
    # "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass0Run2026C-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass0Run2026D-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass1Run2026A-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass1Run2026B-PromptReco-v1",
    # "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass1Run2026C-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass1Run2026D-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass2Run2026A-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass2Run2026B-PromptReco-v1",
    # "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass2Run2026C-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass2Run2026D-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass3Run2026A-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass3Run2026B-PromptReco-v1",
    # "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass3Run2026C-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass3Run2026D-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass4Run2026A-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass4Run2026B-PromptReco-v1",
    # "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass4Run2026C-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass4Run2026D-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass5Run2026A-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass5Run2026B-PromptReco-v1",
    # "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass5Run2026C-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass5Run2026D-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass6Run2026A-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass6Run2026B-PromptReco-v1",
    # "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass6Run2026C-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass6Run2026D-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass7Run2026A-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass7Run2026B-PromptReco-v1",
    # "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass7Run2026C-PromptReco-v1",
    "./crab_muMuGamma_16Sep26_ParkingDoubleMuonLowMass7Run2026D-PromptReco-v1",


    # add more here...
]


output_file = "muMuGammaTree_ntuples_26.txt"

ntuples = []

for crab_dir in crab_dirs:

    print(f"Processing: {crab_dir}")

    result = subprocess.run(
        ["crab", "getoutput", "-d", crab_dir, 
         "--dump",
         "--jobids=1,2,3,4,5,6,7,8,9,10",
         
         ],
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
