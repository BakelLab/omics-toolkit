map-array-probes_bpmap2tab.pl
-----------------------------

  Extracts a tab-delimited file of probe information from a bpmap file for a specified
  sequence group. The output simply consists of one column with the probe identifier and
  another column with the probe sequence. The difference with a standard blat input mapping
  file is that the probe identifier derived from the bpmap file is specially formatted to
  allow for specific remapping of the bpmap file after BLAT analysis of the probe file.
  The format of the probe identifier is '|' delimited and consists of:
  
  sequence group name | PM array X-coordinate | PM array Y-coordinate
  
  The sequence group name is later used to determine which sequence groups need to be remapped.
  Other sequence groups, such as controls etc. will not be remapped, but simply transferred
  to the new bpmap file. The last fields in the identifier together comprise a string that uniquely
  identifies each array probes and will be used to hash-lookup the blat hits for each probe during 
  the actual remapping of the bpmap file.



map-array-probes_run-blat.pl
----------------------------
  
  Takes a tab-delimited probe sequence file (probe_identifier <tab> probe_sequence) and perfoms
  a blat search against a series of fasta files of the reference genome. Each blat search against
  a genome sequence file will be submitted as separate jobs to the pbs queue. The output of all
  blat jobs will be collected in an output dir for further processing. To minimize the size of the
  output files, only hits with >18 matches will be included in the output.



map-array-probes_blat2gff.pl
----------------------------

  Processes a folder with blat output files generated by the map-array-probes_run-blat.pl command
  and produces a single GFF file with all genomic matches according to specified cutoffs.
  


map-array-probes_blat2bpmap.pl
------------------------------

  Processes a folder with blat output files that were generated by sequentially running
  map-array-probes_bpmap2tab.pl and then map-array-probes_run-blat.pl. It will produce a single
  BPMAP file with all the new mappings. 


