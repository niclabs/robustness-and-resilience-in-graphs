target='/Users/elisa/Dropbox/Research/Topics/Resilience/survey/actual/ComNet/'
# python3 experiments.py
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
Rscript analyze.R
cp single_scalar_10sec.eps ${target}/fig2a.eps
cp double_scalar_10sec.eps ${target}/fig2b.eps
cp corr_single_scalar.eps ${target}/fig3.eps
cp clust_single_scalar.eps ${target}/fig4.eps
cp poscor_g1.eps ${target}/fig5a.eps
cp poscor_g2.eps ${target}/fig5b.eps
cp poscor_g3.eps ${target}/fig5c.eps
cp poscor_g4.eps ${target}/fig5d.eps
cp values_splittingNumber_10sec.eps ${target}/fig6a.eps
cp values_connectivityRobustnessFunction_10sec.eps ${target}/fig6b.eps
cp values_electricalNodalRobustness_10sec.eps ${target}/fig6c.eps
cp values_percolatedPath_10sec.eps ${target}/fig6d.eps
cp values_RCB_10sec.eps ${target}/fig6e.eps
cp values_hubDensity_10sec.eps ${target}/fig6f.eps
# the journal wants the panels combined
echo "Combined"
cp examples.eps ${target}/final/fig1.eps
cp fig2.eps ${target}/final
cp corr_single_scalar.eps ${target}/final/fig3.eps
cp clust_single_scalar.eps ${target}/final/fig4.eps
cp fig5.eps ${target}/final
cp fig6.eps ${target}/final
cd ${target}/final
pdflatex schaeffer
