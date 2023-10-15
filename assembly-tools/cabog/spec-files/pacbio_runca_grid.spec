useGrid                = 1
scriptOnGrid           = 1
frgCorrOnGrid          = 1
ovlCorrOnGrid          = 1
cnsOnGrid              = 1

gridEngine             = LSF
sge                    = -q alloc -oo /hpc/users/vanbah01/lsf-output/%J.OU
sgeName                = gridASM
sgeScript              = -n 1   -W 144:00 -R span[ptile=1]  -R rusage[mem=2560]
sgeMerTrim             = -n 12  -W 24:00  -R span[ptile=12] -R rusage[mem=4437]
sgeOverlap             = -n 12  -W 96:00  -R span[ptile=12] -R rusage[mem=4437]
sgeMerOverlapSeed      = -n 12  -W 24:00  -R span[ptile=12] -R rusage[mem=4437]
sgeMerOverlapExtend    = -n 12  -W 24:00  -R span[ptile=12] -R rusage[mem=4437]
sgeConsensus           = -n 2   -W 96:00  -R span[ptile=2]  -R rusage[mem=4096]
sgeFragmentCorrection  = -n 12  -W 96:00  -R span[ptile=12] -R rusage[mem=4437]
sgeOverlapCorrection   = -n 12  -W 96:00  -R span[ptile=12] -R rusage[mem=4437]

utgGenomeSize          = 820000000
doOverlapBasedTrimming = 0
unitigger              = bogart

merSize                = 14
merylMemory            = 32000
merylThreads           = 12

ovlMinLen              = 400
ovlHashBits            = 24
ovlThreads             = 12
ovlHashBlockLength     = 800000000
ovlRefBlockSize        =  100000000
ovlStoreMemory         = 48000

frgCorrThreads         = 12
frgCorrBatchSize       = 500000
ovlCorrBatchSize       = 500000

utgErrorRate           = 0.07
ovlErrorRate           = 0.10
cgwErrorRate           = 0.10
cnsErrorRate           = 0.10

utgGraphErrorRate      = 0.07
utgGraphErrorLimit     = 3.25
utgMergeErrorRate      = 0.085 
utgMergeErrorLimit     = 5.25 
