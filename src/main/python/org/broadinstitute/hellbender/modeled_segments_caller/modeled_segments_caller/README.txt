==============================
GATK CNV ModeledSegmentsCaller
==============================

`modeled_segments_caller` is a python module that identifies segments with cancerous copy
number variations.

This module makes use of both copy ratio and allele fraction data from the CNV pipeline,
as output by the java class ModeledSegments.

The data characterizing the copy ratio and allele fraction posterior distributions over
each segment in the genome are read from the output .seg files of ModeledSegments.
The data used for finding the cancerous segments are sampled from these distributions,
and clustered.

We use the copy ratio data to identify the range of copy ratios corresponding to normal
segments as follows: the one dimensional copy ratio data is clustered, and the normal
cluster is identified either as the one with the smallest or the second smallest copy
ratio values. To decide which one of these two is the normal segment, we check the
corresponding allele fraction values. If the cluster of the smallest copy ratio has
considerable weight with allele fraction values close to 0.5, then it is considered to
be the normal cluster. If not, then the cluster of the second smallest copy ratio is chosen
to be the normal one.

Finally, we fit Gaussians to the data in the two dimensional (copy ratio, allele fraction)
space, and chose those peaks to be normal which fall in the normal copy ratio region and
their mean allele fraction value is close to 0.5. We determine the responsibility of each
segment being normal by determining the responsibility of them being covered by the normal
Gaussian peak versus any other peak.

We plot the results showing the normal segments in black and the called segments in
red color. We also show the PHRED score corresponding to the probability of each segment
being normal.

