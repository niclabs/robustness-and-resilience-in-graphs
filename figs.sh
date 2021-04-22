target='/Users/elisa/Dropbox/Research/Topics/Resilience/survey/actual/ComNet/'
python3 experiments.py
cp nxcomplete.eps ${target}/fig1a.eps
cp nxtri.eps ${target}/fig1b.eps
cp nxsquare.eps ${target}/fig1c.eps
cp nxhex.eps ${target}/fig1d.eps
cp nxws.eps ${target}/fig1e.eps
cp nxregular.eps ${target}/fig1f.eps
cp nxerdos.eps ${target}/fig1g.eps
cp nxcirc.eps ${target}/fig1h.eps
cp nxladder.eps ${target}/fig1i.eps
cp nxba.eps ${target}/fig1j.eps
cp nxtree.eps ${target}/fig1k.eps
cp nxstar.eps ${target}/fig1l.eps
cp nxbarbell.eps ${target}/fig1m.eps
cp nxwheel.eps ${target}/fig1n.eps
cp nxpath.eps ${target}/fig1o.eps
cd ${target}
pdflatex schaeffer
