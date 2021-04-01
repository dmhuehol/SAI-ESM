# Use CDO to prepare files for later usage

import process_glens_fun as pgf

inPath = '/Users/dhueholt/Documents/GLENS_data/CompleteTest/'
cntrlCard = 'b.e15.B5505C5WCCML45BGCR.f09_g16.control.001.cam.h0.T.*.nc'
fdbckCard = 'b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.001.cam.h0.T.*.nc'

pgf.prep_raw_data(inPath,cntrlCard,'control')
print("Prepared control files")
pgf.prep_raw_data(inPath,fdbckCard,'feedback')
print("Prepared feedback files")
