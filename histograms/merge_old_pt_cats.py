import os

import ROOT
import re

def run_mass_cut(input, output_dir, m_key, m_val, resolution = 0.1):
    INPUT_FILE = input

    MASS_MIN = m_val - m_val * resolution
    MASS_MAX = m_val + m_val * resolution
    # MASS_MIN = m_val - 0.01
    # MASS_MAX = m_val + 0.01

    pattern = re.compile(r"Cat_(\d+)_(bb|ee)_(SR|CR)$")

    f = ROOT.TFile.Open(INPUT_FILE)

    merged = {}

    for key in f.GetListOfKeys():

        name = key.GetName()

        m = pattern.match(name)
        if not m:
            continue

        n = int(m.group(1))
        xx = m.group(2)
        yy = m.group(3)

        if xx == "bb" and not (1 <= n <= 10):
            continue

        if xx == "ee" and not (1 <= n <= 5):
            continue

        h = f.Get(name)

        # find mass window
        bin_lo = h.FindBin(MASS_MIN)
        bin_hi = h.FindBin(MASS_MAX)

        nbins_new = bin_hi - bin_lo + 1

        xmin_new = h.GetBinLowEdge(bin_lo)
        xmax_new = h.GetBinLowEdge(bin_hi + 1)

        # create reduced histogram
        hwin = ROOT.TH1D(
            f"{name}_window",
            h.GetTitle(),
            nbins_new,
            xmin_new,
            xmax_new,
        )

        hwin.Sumw2()

        for new_bin, old_bin in enumerate(range(bin_lo, bin_hi + 1), start=1):
            hwin.SetBinContent(new_bin, h.GetBinContent(old_bin))
            hwin.SetBinError(new_bin, h.GetBinError(old_bin))

        merge_key = (xx, yy)

        if merge_key not in merged:

            merged_hist = hwin.Clone(f"mass_{xx}_{yy}")
            merged_hist.SetDirectory(0)

            merged[merge_key] = merged_hist

        else:
            merged[merge_key].Add(hwin)

    f.Close()

    # write one file per region

    # make results directory if it doesn't exist
    results_dir = f"./results/{output_dir}"
    if os.path.exists(results_dir) == False:
        os.makedirs(results_dir)

    for (xx, yy), hist in merged.items():

        fout = ROOT.TFile(f"{results_dir}/mass_{m_key}_{xx}_{yy}.root", "RECREATE")
        hist.Write()
        fout.Close()

        print(
            f"Saved mass_{m_key}_{xx}_{yy}"
            f" ({hist.GetNbinsX()} bins, "
            f"{hist.GetXaxis().GetXmin():.2f} - "
            f"{hist.GetXaxis().GetXmax():.2f})"
        )


def main():
    input_file = "/afs/cern.ch/work/s/secholak/cms/workdir/eta2mumug/histogram_production/Run3DimuonAnalysisTools/histograms/histosApr22.root"
    for m in [0.2850, 0.2875, 0.2900, 0.2925, 0.2950, 0.2975, 0.3000, 0.3025, 0.3050, 0.3075, 0.3100]:
        output_dir = "excess_0p2975"
        run_mass_cut(input_file, output_dir, f"{m:.4f}".replace(".", "p"), m, resolution=0.1)
    
    # excess at 0.4575
    for m in [0.4450, 0.4475, 0.4500, 0.4525, 0.4550, 0.4575, 0.4600, 0.4625, 0.4650, 0.4675, 0.4700]:
        output_dir = "excess_0p4575"
        run_mass_cut(input_file, output_dir, f"{m:.4f}".replace(".", "p"), m, resolution=0.07)



if __name__ == "__main__":
    main()
