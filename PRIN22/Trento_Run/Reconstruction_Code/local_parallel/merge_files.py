import ROOT
import pickle

def merge_files():
    path = "./"
    chain = ROOT.TChain("mtracks")

    with open("view_info.pkl", "rb") as file:
        data = pickle.load(file)

    Nviews = len(data["actual_views"])

    for n in range(0, Nviews, 200):
        filename = f"{path}mt.merged.temp.{n}.root"

        if not ROOT.TFile.Open(filename):
            print(f"File does not exist: {filename}, moving to the next file")
            continue

        chain.Add(filename)

    #outputfile = ROOT.TFile(f"{path}mt.merged.all.root", "RECREATE")

    chain.Merge(f"{path}mt.merged.all.root")

if __name__ == "__main__":
    merge_files()
