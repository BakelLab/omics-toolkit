# -------------------------------------
#  GRID CONFIGURATION
# -------------------------------------
#
#  By default, do NOT use SGE; use manual override with command-line options.
#
useGrid                = 0
scriptOnGrid           = 0
frgCorrOnGrid          = 1
ovlCorrOnGrid          = 1
maxGridJobSize         = 10
sge                    = -m a -r n
sgeScript              = -l nodes=1:ppn=8,walltime=48:00:00,mem=16G
sgeOverlap             = -l nodes=1:ppn=8,walltime=48:00:00,mem=16G
sgeMerOverlapSeed      = -l nodes=1:ppn=8,walltime=48:00:00,mem=16G
sgeMerOverlapExtend    = -l nodes=1:ppn=8,walltime=48:00:00,mem=16G
sgeConsensus           = -l nodes=1:ppn=8,walltime=48:00:00,mem=16G
sgeFragmentCorrection  = -l nodes=1:ppn=8,walltime=48:00:00,mem=16G
sgeOverlapCorrection   = -l nodes=1:ppn=8,walltime=48:00:00,mem=16G
#
# -------------------------------------
#  ERROR RATES
# -------------------------------------
#
#  Expected rate of sequencing error. Allow pairwise alignments up to this rate.
#  
utgErrorRate        = 0.03  # Sanger reads can use values less than one. Titanium reads require 3% in unitig.
utgErrorLimit       = 2.5   # Allow mismatches over and above the utgErrorRate. This helps with Illumina reads.
ovlErrorRate        = 0.06  # Larger than utg to allow for correction.
cnsErrorRate        = 0.10  # Larger than utg to avoid occasional consensus failures
cgwErrorRate        = 0.10  # Larger than utg to allow contig merges across high-error ends
#
# -------------------------------------
# STAGE CONFIGURATION
# -------------------------------------
#
# MERYL
merylMemory         = 120000
merylThreads        = 16
# 
# OVERLAPPER
merSize             = 24    # default=22; use lower to combine across heterozygosity, higher to separate near-identical repeat copies 
overlapper          = ovl   # the mer overlapper for 454-like data is insensitive to homopolymer problems but requires more RAM and disk
ovlMemory           = 8GB --hashload 0.8 --hashstrings 400000 --hashdatalen 400000000
ovlThreads          = 8
ovlHashBlockSize    = 8000000
ovlRefBlockSize     = 32000000
ovlConcurrency      = 1          # Not needed when using grid
#
# OVERLAP STORE
ovlStoreMemory      = 90000      # Must be in Mbp; keep it a bit below the max available memory to avoid malloc errors
# 
# ERROR CORRECTION
frgCorrThreads      = 16
frgCorrBatchSize    = 48000000   # We want this to be as big as possible. A batch of 50M reads will use around 100 Gb of memory.
frgCorrConcurrency  = 1          # not needed when using grid
ovlCorrBatchSize    = 2500000    # 2.5M reads can run on a single GPC node (4M should work too); submitting 80 at a time works nicely
ovlCorrConcurrency  = 1          # not needed when using grid
#
# UNITIGGER
unitigger           = bog        # Use bog when dealing with Illumina reads
utgBubblePopping    = 1
# utgGenomeSize     =
#
# CONSENSUS
cnsConcurrency      = 16
#
# SCAFFOLDING
doExtendClearRanges = 0          # Let's skip extended clear ranges during scaffolding for now
#
# -------------------------------------
# FRG FILES
# -------------------------------------
#
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/100814_stitch_nobact_s_1.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/100814_stitch_nobact_s_2.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/100814_stitch_nobact_s_3.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/100814_stitch_nobact_s_5.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/100814_stitch_nobact_s_6.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/100814_stitch_nobact_s_7.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/100814_stitch_nobact_s_8.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/101020_stitch_nobact_s_1.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/101020_stitch_nobact_s_2.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/101020_stitch_nobact_s_3.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/101020_stitch_nobact_s_6.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/101020_stitch_nobact_s_7.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/101020_stitch_nobact_s_8.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/101111_stitch_nobact_s_2.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/101111_stitch_nobact_s_3.frg
/home/vanbakel/scratch/assemblies/cabog_short-inserts-stitchonly/101111_stitch_nobact_s_6.frg
