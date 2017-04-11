import ROOT

# LOAD 
c=ROOT.KalmanMuonCalibrator("MC_80X_13TeV")

print 'correct the pt of a 40 GeV muon + at eta=phi=0.0'
print c.getCorrectedPt(40,0.0,0.0,1)

print 'correct the pt of a 40 GeV muon - at eta=phi=0.0'
print c.getCorrectedPt(40,0.0,0.0,-1)


print 'correct the dpt/pt  a 40 geV muon at eta =1 that has dpt/pt=0.01'
print 'Recall! This correction is only valid after smearing'
print c.getCorrectedError(40,1.0,0.01)


print 'smear a 40 geV muon at eta=1.0 to match the data'
print c.smear(40,1.0)




print 'propagate the statistical error of the calibratio to a + 40 GeV muon at eta=phi=0'
print 'first get number of parameters'
N=c.getN()
print N,'parameters'
for i in range(0,N):
    c.vary(i,+1)
    print 'variation',i,'ptUp', c.getCorrectedPt(40,0.0,0.0,1)
    c.vary(i,-1)
    print 'variation',i,'ptDwn', c.getCorrectedPt(40,0.0,0.0,1)


c.reset()

print 'propagate the closure error to a + 40 GeV muon at eta=phi=0'

c.varyClosure(+1)
print 'After closure shift pt=', c.getCorrectedPt(40,0.0,0.0,1)




 

