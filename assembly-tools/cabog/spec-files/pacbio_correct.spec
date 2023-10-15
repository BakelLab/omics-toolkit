# Use blasr for alignments
blasr = -bestn 14 -nCandidates 50 -noRefineAlign -ignoreQuality -minMatch 10 -maxLCPLength 15 -minPctIdentity 70.0 -maxAnchorsPerPosition 500

# original asm settings
utgErrorRate = 0.25
utgErrorLimit = 6.5

cnsErrorRate = 0.25
cgwErrorRate = 0.25
ovlErrorRate = 0.25

merSize      = 14
merThreshold = 500

merylMemory = 32000
merylThreads = 12

ovlStoreMemory = 48000

# grid info
useGrid = 0
scriptOnGrid = 0
frgCorrOnGrid = 0
ovlCorrOnGrid = 0

ovlHashBits = 24
ovlThreads = 4
ovlHashBlockLength =  1000000000
ovlRefBlockLength  =  6700000000

frgCorrThreads = 32
frgCorrBatchSize = 48000000
ovlCorrBatchSize = 2500000

mbtConcurrency = 16
ovlConcurrency = 10
cnsConcurrency = 32
frgCorrConcurrency = 16
ovlCorrConcurrency = 32 
