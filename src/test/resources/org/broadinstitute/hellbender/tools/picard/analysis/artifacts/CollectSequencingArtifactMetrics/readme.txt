Description of reads (may be out of order):

perfect_fr/1 - read 1+ with no artifacts
perfect_fr/2 - read 2- with no artifacts
nbases_fr/1 - read 1+ with no artifacts, but some Ns in REFERENCE (contexts with an N _anywhere_ will be skipped)
nbases_fr/2 - read 2- with no artifacts, but some Ns in READ (N sites will be skipped, but otherwise do not affect context)
pro_oxog_fr/1 - read 1+ consistent with OxoG (G>T)
pro_oxog_fr/2 - read 2- consistent with OxoG (G>T)
pro_oxog_rf/1 - read 1- consistent with OxoG (C>A)
pro_oxog_rf/2 - read 2+ consistent with OxoG (C>A)
pro_ffpe_fr/1 - read 1+ consistent with FFPE (C>T)
pro_ffpe_fr/2 - read 2- consistent with FFPE (C>T)
con_ffpe_rf/1 - read 1- INconsistent with FFPE (C>T)
con_ffpe_rf/2 - read 2+ INconsistent with FFPE (C>T)
indel_oxog_fr/1 - read 1+ with insertion, deletion and OxoG artifact (indel sites will be skipped, but otherwise do not affect context)
indel_oxog_fr/2 - read 2- with insertion, deletion and OxoG artifact (indel sites will be skipped, but otherwise do not affect context)
unmapped_mate/1 - mapped read 1 with unmapped mate (may be filtered out)
unmapped_mate/2 - unmapped read 2 with mapped mate (mate may be filtered out)
duplicate/1 - duplicate (should always be filtered out)
duplicate/2 - duplicate (should always be filtered out)
non_primary/1 - non-primary alignment (should always be filtered out)
non_primary/2 - non-primary alignment (should always be filtered out)
low_quality/1 - has some low BASE qualities (those sites may be skipped) 
low_quality/2 - has a low MAPPING quality (may be filtered out)