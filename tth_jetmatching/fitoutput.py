import ROOT
infile="test.root"
outfile="out.root"

f=ROOT.TFile(infile)
names=["Higgs","TopHad","TopLep","WHad"]
widths=[15,15,15,15]
suffixes=["best","reco","recoall"]
histos=[]

out=ROOT.TFile(outfile,"recreate")
for n in names:
    for s in suffixes:
        name="M_"+n+"_"+s
        print name
        h=f.Get(name)
        histos.append(h)
out.cd()
fits=[]
for h in histos:
    fn=h.GetName()+"_fit"
    imax=h.GetMaximumBin()
    xmax=h.GetBinCenter(imax)
    maxcontent=h.GetBinContent(imax)
    mx=0.
    for i in range(imax,h.GetNbinsX()):
        if h.GetBinContent(i) < 0.5*maxcontent:
            mx=h.GetBinCenter(i)
            break
    mn=0.
    for i in range(imax,0,-1):
        if h.GetBinContent(i) < 0.5*maxcontent:
            mn=h.GetBinCenter(i)
            break

    fit=ROOT.TF1(fn,"gaus",mn,mx)
    h.Fit(fn,"q","",mn,mx)
    h.Write()
    fits.append(fit)
    fit.Write()
    print h.GetName(),fit.GetParameter(1),fit.GetParameter(2)
