package org.broadinstitute.hellbender.tools.walkers.rnaseq;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3BaseData;
import htsjdk.tribble.gff.Gff3Feature;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.AlignmentAgreesWithHeaderReadFilter;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import org.broadinstitute.hellbender.utils.read.GATKRead;

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;


import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Evaluate gene expression from RNA-seq reads aligned to genome.
 *
 * <p>This tool evaluates gene expression from RNA-seq reads aligned to genome.  Features to evaluate expression over are defined in an input annotation file in gff3 fomat
 * (https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).  Output is a tsv listing sense and antisense expression for all stranded grouping features,
 * and expression (labeled as sense) for all unstranded grouping features.
 * </p>
 * 
 * <p>
 *     <h3>Input</h3>
 *     <ul>
 *     <li>BAM file of RNA-seq reads</li>
 *     <li>Gff3 file of feature annotations</li>
 *     </ul>
 * </p>
 * <p>
 *     <h3>Output</h3>
 *     TSV file of gene expression
 * </p>
 *
 * <p>
 *     <h3>Usage Examples</h3>
 *     <p>
 *         gatk GeneExpressionEvaluation
 *         -I input.bam
 *         -G geneAnnotations.gff3
 *         -O output.tsv
 *     </p>
 * </p>
 * 
 *<p>Reads are assumed to be paired-end.  Reads which are in a "good pair" are counted once together as a single fragment.  Reads which are not in a "good pair" are each counted separately.
 * A "good pair" is defined as:
 *  <li>Both reads are mapped</li>
 *  <li>Properly paired flag is set</li>
 *  <li>Both reads are on same contig</li>
 *  <li>Reads are on opposite strands</li>
 *  <li>Mapping quality of both reads is at least minimum-mapping-quality</li>
 *  <li>Reads are inward-facing</li>
 *</p>
 *
 * <p>Reads can be either from spliced are unspliced RNA.  If from spliced RNA, alignment blocks of reads are taken as their coverage.  
 * If from unspliced RNA, the entire region from the start of the earliest read to the end of the latest read is taken as the fragment
 * coverage.  Splice status is set through the command line.  By default, splice status is taken to be "spliced." </p>
 * <p>
 *     <h3>For spliced RNA</h3>
 *     <p>
 *         gatk GeneExpressionEvaluation
 *         -I input.bam
 *         -G geneAnnotations.gff3
 *         -O output.tsv
 *         --spliced true
 *     </p>
 * </p>
 * <p>
 *     <h3>For unspliced RNA</h3>
 *     <p>
 *         gatk GeneExpressionEvaluation
 *         -I input.bam
 *         -G geneAnnotations.gff3
 *         -O output.tsv
 *         --spliced true
 *     </p>
 * </p>
 * 
 * <p>Gene expression is aggregated over "groupingType" features, and overlap of fragments with features is determined by "overlapType" features.  By default,
 * "groupingType" is genes, and "overlapType" is exons.  Additional grouping_types and overlap_types can be used as well.</p>
 * <p>
 *     <h3>Use gene and pseudogene as grouping types</h3>
 *     <p>
 *     gatk GeneExpressionEvaluation
 *     -I input.bam
 *     -G geneAnnotations.gff3
 *     -O output.tsv
 *     --grouping-type gene
 *     --grouping-type pseudogene
 *     </p>
 * </p>
 * 
 * <p>The transcription read is by default read1, however this can be changed to read2.  The transcription read determines whether a fragment is sense or antisense.  
 * If the transcription read is on the same strand as the feature, the fragment is sense, if on the opposite strand, the fragment is antisense.
 * </p>
 * <p>
 *     <h3>Set transcription read to read2</h3>
 *     <p>
 *         gatk GeneExpressionEvaluation
 *         -I input.bam
 *         -G geneAnnotations.gff3
 *         -O output.tsv
 *         --transcription-read R2
 *         
 * </p>
 *
 * <p>Multi-overlapping fragments (fragment alignments which overlap multiple grouping features) can be handled in two ways.  Equals weight can be given to each grouping feature,
 * in which case each is given weight 1/N, for N overlapping grouping features.
 * </p>
 * <p>
 *     <h3>Equal weight for multi-overlapping fragments</h3>
 *     <p>
 *         gatk GeneExpressionEvaluation
 *         -I input.bam
 *         -G geneAnnotations.gff3
 *         -O output.tsv
 *         --multi-overlap-method EQUAL
 *     </p>
 * </p>
 * <p>Multi-overlapping fragments can also have weight distributed according to how much of the fragment overlaps each feature.  In this case, each grouping feature is given an unnormalized weight corresponding
 * to the fraction of the fragment that overlaps its overlapping features.  A "non-overlapping" option is also given an unnormalized weight corresponding the the fraction of the fragment that overlaps no features.
 * Final weights for each feature are found by normalizing by the sum of unnormalized weights.  This is the default behavior.
 * </p>
 * <p>
 *     <h3>Proportional weight for multi-overlapping fragments</h3>
 *     <p>
 *         gatk GeneExpressionEvaluation
 *         -I input.bam
 *         -G geneAnnotations.gff3
 *         -O output.tsv
 *         --multi-overlap-method PROPORTIONAL
 *     </p>
 * </p>
 *
 * <p>Multi-mapping fragments (fragments whose reads map to multiple locations) can also be handled in two ways.  They can be ignored, in which case only fragments with a single mapping are counted.
 * This is default behavior.
 * </p>
 * <p>
 *     <h3>Ignore multi-mapping fragments</h3>
 *     <p>
 *         gatk GeneExpressionEvaluation
 *         -I input.bam
 *         -G geneAnnotations.gff3
 *         -O output.tsv
 *         --multi-map-method IGNORE
 *     </p>
 * </p>
 * <p>Multi-mapping fragments can alternatively have their weight distributed equally between the different alignments.  If run in this setting, minimum-mapping-quality will be set to 0,
 * regardless of its value on the command line.
 * </p>
 * <p>
 *     <h3>Equally weight each alignment of multi-mapping fragments</h3>
 *     <p>
 *         gatk GeneExpressionEvaluation
 *         -I input.bam
 *         -G geneAnnotations.gff3
 *         -O output.tsv
 *         --multi-map-method EQUAL
 *     </p>
 * </p>
 * <p>The number of mapping for a particular fragment is extracted via the NH tag.</p>
 * <p>Currently this tool only works on a single sample BAM.  If a readgroup header line indicates multiple samples in the same BAM, the tool will throw an exception.</p>
 */
@CommandLineProgramProperties(
        summary = "This tool evaluates gene expression from RNA-seq reads aligned to genome.  Features to evaluate expression over are defined in an input annotation file in gff3 fomat " +
                "(https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).  Output is a tsv listing sense and antisense expression for all stranded grouping features, " +
                "and expression (labeled as sense) for all unstranded grouping features.",
        oneLineSummary = "Evaluate gene expression from RNA-seq reads aligned to genome.",
        programGroup = CoverageAnalysisProgramGroup.class
)
@DocumentedFeature
public final class GeneExpressionEvaluation extends ReadWalker {

    @Argument(
            doc = "Output file for gene expression.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputCountsFile = null;

    @Argument(doc="Gff3 file containing feature annotations", shortName = "G", fullName = "gff-file")
    private FeatureInput<Gff3Feature> gffFile;

    @Argument(doc="Feature types to group by", fullName = "grouping-type")
    private Set<String> groupingType = new HashSet<>(Collections.singleton("gene"));

    @Argument(doc="Feature overlap types", fullName = "overlap-type")
    private Set<String> overlapType = new HashSet<>(Collections.singleton("exon"));

    @Argument(doc="Attribute key to label features by in output", fullName = "feature-label-key")
    private String featureLabelKey = "Name";

    @Argument(doc = "Which read corresponds to the transcription strand", fullName = "transcription-read")
    private TrancriptionRead trancriptionRead = TrancriptionRead.R1;

    @Argument(doc = "How to distribute weight of alignments which overlap multiple features", fullName = "multi-overlap-method")
    private MultiOverlapMethod multiOverlapMethod = MultiOverlapMethod.PROPORTIONAL;

    @Argument(doc = "How to distribute weight of reads with multiple alignments", fullName = "multi-map-method")
    private MultiMapMethod multiMapMethod = MultiMapMethod.IGNORE;

    @Argument(doc = "Whether the rna is spliced.  If spliced, alignments must be from a splice aware aligner (such as star).  If unspliced, alignments must be from a non-splicing aligner (such as bwa). ")
    private boolean spliced = true;

    final private Map<Gff3BaseData, Coverage> featureCounts = new LinkedHashMap<>();

    final private OverlapDetector<Pair<Gff3BaseData, Interval>> featureOverlapDetector = new OverlapDetector<>(0,0);

    private String sampleName = null;

    final MappingQualityReadFilter mappingQualityFilter = new MappingQualityReadFilter();

    enum MultiOverlapMethod {

        EQUAL {
            @Override
            Map<Gff3BaseData, Float> getWeights(final List<Interval> alignmentIntervals, final OverlapDetector<Pair<Gff3BaseData, Interval>> featureOverlapDetector) {
                final Set<Gff3BaseData> overlappingFeatures = alignmentIntervals.stream().flatMap(i -> featureOverlapDetector.getOverlaps(i).stream().map(Pair::getLeft)).collect(Collectors.toSet());
                final int nOverlappingFeatures = overlappingFeatures.size();
                final Map<Gff3BaseData, Float> weights = new LinkedHashMap<>();
                for (final Gff3BaseData feature : overlappingFeatures) {
                    weights.put(feature, (float)1.0/nOverlappingFeatures);
                }
                return weights;
            }
        },

        PROPORTIONAL {
            @Override
            Map<Gff3BaseData, Float> getWeights(final List<Interval> alignmentIntervals, final OverlapDetector<Pair<Gff3BaseData, Interval>> featureOverlapDetector) {
                final List<Interval> mergedAlignmentIntervals = getMergedIntervals(alignmentIntervals);

                final int basesOnReference = mergedAlignmentIntervals.stream().map(Interval::getLengthOnReference).reduce(0, Integer::sum);
                int totalCoveredBases = 0;
                float summedUnNormalizedWeights = 0;
                final Map<Gff3BaseData, Float> weights = new LinkedHashMap<>();
                for (final Interval alignmentInterval : mergedAlignmentIntervals) {
                    final Set<Pair<Gff3BaseData,Interval>> overlaps = featureOverlapDetector.getOverlaps(alignmentInterval);
                    final Map<Gff3BaseData, List<Interval>> overlappingIntervalsByFeature = new LinkedHashMap<>();
                    final List<Interval> allOverlappingIntervals = new ArrayList<>();
                    for (Pair<Gff3BaseData, Interval> overlap : overlaps) {
                        final List<Interval> overlappingIntervals = overlappingIntervalsByFeature.computeIfAbsent(overlap.getLeft(), f -> new ArrayList<>());
                        overlappingIntervals.add(overlap.getRight());
                        allOverlappingIntervals.add(overlap.getRight());
                    }

                    final List<Interval> allOverlappingIntervalsMerged = getMergedIntervals(allOverlappingIntervals);
                    totalCoveredBases += allOverlappingIntervalsMerged.stream().map(alignmentInterval::getIntersectionLength).reduce(0,Integer::sum);
                    for (Map.Entry<Gff3BaseData, List<Interval>> overlap : overlappingIntervalsByFeature.entrySet()) {
                        final List<Interval> mergedOverlapIntervals = getMergedIntervals(overlap.getValue());
                        final float weight = (float)mergedOverlapIntervals.stream().map(alignmentInterval::getIntersectionLength).reduce(0,Integer::sum)/basesOnReference;
                        weights.compute(overlap.getKey(), (k,v) -> v == null? weight : v + weight);
                        summedUnNormalizedWeights += weight;
                    }
                }

                summedUnNormalizedWeights += 1.0 - (float)totalCoveredBases/basesOnReference;

                final float normalizationFactor = (float)1.0/summedUnNormalizedWeights;

                for (final Gff3BaseData feature : weights.keySet()) {
                    weights.compute(feature, (k,v) -> v*normalizationFactor);
                }

                return weights;
            }
        };

        abstract Map<Gff3BaseData, Float> getWeights(final List<Interval> alignmentIntervals, final OverlapDetector<Pair<Gff3BaseData, Interval>> featureOverlapDetector);
    }

    enum MultiMapMethod {
        IGNORE {
            @Override
            protected Map<Gff3BaseData, Float> getWeightsForMethod(final int nHits, final Map<Gff3BaseData, Float> previousWeights) {
                if (nHits == 1) {
                    return previousWeights;
                } else {
                    return Collections.emptyMap();
                }
            }
        },
        EQUAL {
            @Override
            protected Map<Gff3BaseData, Float> getWeightsForMethod(final int nHits, final Map<Gff3BaseData, Float> previousWeights) {
                if (nHits == 1) {
                    return previousWeights;
                } else {
                    final Map<Gff3BaseData, Float> newWeights = new HashMap<>(previousWeights.size());
                    for (final Map.Entry<Gff3BaseData, Float> entry : previousWeights.entrySet()) {
                        newWeights.put(entry.getKey(), entry.getValue()/(float)nHits);
                    }
                    return newWeights;
                }
            }
        };

        Map<Gff3BaseData, Float> getWeights(final int nHits, final Map<Gff3BaseData, Float> previousWeights) {
            if (nHits < 1) {
                throw new GATKException("nHits = " + nHits + ", cannot be less than 1");
            }

            return getWeightsForMethod(nHits, previousWeights);
        }

        abstract protected Map<Gff3BaseData, Float> getWeightsForMethod(final int nHits, final Map<Gff3BaseData, Float> previousWeights);

    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> readFilters = new ArrayList<>();
        readFilters.add(ReadFilterLibrary.VALID_ALIGNMENT_START);
        readFilters.add(ReadFilterLibrary.VALID_ALIGNMENT_END);
        readFilters.add(new AlignmentAgreesWithHeaderReadFilter(getHeaderForReads()));
        readFilters.add(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS);
        readFilters.add(ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH);
        readFilters.add(ReadFilterLibrary.SEQ_IS_STORED);
        readFilters.add(ReadFilterLibrary.MAPPED);
        readFilters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        readFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        readFilters.add(mappingQualityFilter);
        return readFilters;
    }

    @Override
    public void onTraversalStart() {
        validateOutputFile(outputCountsFile);
        if(multiMapMethod == MultiMapMethod.EQUAL) {
            mappingQualityFilter.minMappingQualityScore = 0;
        }
        final SAMSequenceDictionary dict = getBestAvailableSequenceDictionary();
        final SAMFileHeader header = getHeaderForReads();
        for (final SAMReadGroupRecord readGroupRecord : header.getReadGroups()) {
            if (sampleName == null) {
                sampleName = readGroupRecord.getSample();
            } else {
                if (!sampleName.equals(readGroupRecord.getSample())) {
                    throw new GATKException("Cannot run GeneExpressionEvaluation on multi-sample bam.");
                }
            }
        }
        if (dict == null) {
            throw new GATKException("sequence dictionary must be specified (" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME + ").");
        }

        logger.info("collecting list of features");
        final List<SimpleInterval> allIntervals = hasUserSuppliedIntervals()? getTraversalIntervals() : IntervalUtils.getAllIntervalsForReference(dict);
        for (final SimpleInterval interval : allIntervals) {
            final List<Gff3Feature> contigFeatures = features.getFeatures(gffFile, interval);
            logger.info("collecting features in " + interval.getContig() + ":" + interval.getStart() + "-" + interval.getEnd());
            for (final Gff3Feature feature : contigFeatures) {
                if (groupingType.contains(feature.getType())) {
                    final List<Interval> overlappingFeatures = feature.getDescendents().stream().filter(f -> overlapType.contains(f.getType())).map(f -> new Interval(f.getContig(), f.getStart(), f.getEnd())).collect(Collectors.toList());
                    final Gff3BaseData shrunkGroupingBaseData = shrinkBaseData(feature.getBaseData());
                    addGroupingFeature(shrunkGroupingBaseData, overlappingFeatures);
                }
            }
        }

        logger.info("Collecting read counts...");
    }

    private void validateOutputFile(final File file) {
        if (file.exists()) {
            if (!Files.isWritable(file.toPath())) {
                throw new UserException.CouldNotCreateOutputFile(file, " is not writable");
            }
        } else {
            if (!Files.isWritable(file.getParentFile().toPath())) {
                throw new UserException.CouldNotCreateOutputFile(file, " is not writable");
            }
        }
    }

    private Gff3BaseData shrinkBaseData(final Gff3BaseData baseData) {
        //remove all but featureLabelKey attributes
        final Map<String, String> shrunkAttributes = baseData.getAttributes().entrySet().stream().filter(e -> e.getKey().equals(featureLabelKey)).collect(Collectors.toMap(Map.Entry::getKey,Map.Entry::getValue));
        return new Gff3BaseData(baseData.getContig(), baseData.getSource(), baseData.getType(), baseData.getStart(), baseData.getEnd(), baseData.getStrand(), baseData.getPhase(), shrunkAttributes);
    }

    private void addGroupingFeature(final Gff3BaseData groupingBaseData, final List<Interval> overlappingFeatures) {
        final String geneID = groupingBaseData.getAttributes().get(featureLabelKey);
        if (geneID == null) {
            throw new UserException("no geneid field " + featureLabelKey + " found in feature at " + groupingBaseData.getContig() + ":" + groupingBaseData.getStart() + "-" + groupingBaseData.getEnd());
        }

        featureCounts.put(groupingBaseData, new Coverage(0, 0));
        for (final Interval overlappingFeature : overlappingFeatures) {
            featureOverlapDetector.addLhs(Pair.of(groupingBaseData, overlappingFeature), overlappingFeature);
        }
    }


    static boolean inGoodPair(final GATKRead read, int minimumMappingQuality) {

        boolean ret = !read.mateIsUnmapped() && read.isProperlyPaired() && read.getContig().equals(read.getMateContig()) &&
                read.isReverseStrand() != read.mateIsReverseStrand();

        if (ret) {
            if (!read.hasAttribute(SAMTag.MQ.toString())) {
                throw new GATKException("Mate quality must be included.  Consider running FixMateInformation.");
            }
            ret = ret && read.getAttributeAsInteger(SAMTag.MQ.toString()) >= minimumMappingQuality;
        }

        if (ret) {
            if (read.isReverseStrand()) {
                ret = ret && read.getEnd() >= read.getMateStart();
            } else {
                if (!read.hasAttribute(SAMTag.MC.toString())) {
                    throw new GATKException("Mate cigar must be present.  Consider running FixMateInformation.");
                }
                final Cigar mateCigar = TextCigarCodec.decode(read.getAttributeAsString(SAMTag.MC.toString()));
                ret = ret && read.getStart() <= read.getMateStart() + mateCigar.getReferenceLength();
            }
        }
        return ret;
    }

    static List<Interval> getAlignmentIntervals(final GATKRead read, final boolean spliced, final int minimumMappingQuality) {

        if (spliced) {
            final List<Interval> alignmentIntervals = new ArrayList<>();

            final List<AlignmentBlock> readAlignmentBlocks = SAMUtils.getAlignmentBlocks(read.getCigar(), read.getStart(), "read cigar");

            for( final AlignmentBlock block : readAlignmentBlocks) {
                alignmentIntervals.add(new Interval(read.getContig(), block.getReferenceStart(), block.getReferenceStart()+block.getLength() - 1));
            }

            boolean overlapsMate = false;
            if (inGoodPair(read, minimumMappingQuality)) {
                final String mateCigarString = read.getAttributeAsString(SAMTag.MC.toString());
                if(mateCigarString == null) {
                    throw new GATKException("Mate cigar must be present if using spliced reads");
                }
                final List<AlignmentBlock> mateAlignmentBlocks = SAMUtils.getAlignmentBlocks(TextCigarCodec.decode(mateCigarString), read.getMateStart(), "mate cigar");
                for( final AlignmentBlock block : mateAlignmentBlocks) {
                    final Interval alignmentBlockInterval = new Interval(read.getMateContig(), block.getReferenceStart(), block.getReferenceStart()+block.getLength() - 1);
                    alignmentIntervals.add(alignmentBlockInterval);

                    if (!overlapsMate && read.overlaps(alignmentBlockInterval)) {
                        overlapsMate = true;
                    }
                }
            }
            return getMergedIntervals(alignmentIntervals);

        } else {
            if (read.isUnmapped()) {
                return Collections.emptyList();
            }

            final boolean inGoodPair = inGoodPair(read, minimumMappingQuality);

            final int start = inGoodPair? Math.min(read.getStart(), read.getMateStart()) : read.getStart();
            final int end = inGoodPair? start + Math.abs(read.getFragmentLength()) - 1 : read.getEnd();
            return Collections.singletonList(new Interval(read.getContig(), start, end));
        }


    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        if ((!read.isReverseStrand() || !inGoodPair(read, mappingQualityFilter.minMappingQualityScore))) {
            final List<Interval> alignmentIntervals = getAlignmentIntervals(read, spliced, mappingQualityFilter.minMappingQualityScore);

            final Map<Gff3BaseData, Float> initalWeights = multiOverlapMethod.getWeights(alignmentIntervals, featureOverlapDetector);
            final Map<Gff3BaseData, Float> finalWeights = multiMapMethod.getWeights(read.getAttributeAsInteger(SAMTag.NH.toString()), initalWeights);

            final Strand fragmentStrand = getFragmentStrand(read, trancriptionRead);

            for (final Map.Entry<Gff3BaseData, Float> weight : finalWeights.entrySet()) {
                final Gff3BaseData feature = weight.getKey();
                final boolean isSense = (feature.getStrand() == Strand.NONE) || feature.getStrand() == fragmentStrand;
                if (isSense) {
                    featureCounts.get(feature).addSenseCount(weight.getValue());
                } else {
                    featureCounts.get(feature).addAntiSenseCount(weight.getValue());
                }
            }
        }
    }



    @Override
    public Object onTraversalSuccess() {
        logger.info(String.format("Writing read counts to %s...", outputCountsFile.getAbsolutePath()));
        try (final FragmentCountWriter writer = new FragmentCountWriter(outputCountsFile.toPath(), sampleName, featureLabelKey)) {
            int i=0;
            for (final File input_bam: readArguments.getReadFiles()) {
                writer.writeMetadata("input_bam_"+i, input_bam.toString());
                i++;
            }
            writer.writeMetadata("annotation_file", gffFile.toString());
            for (final Map.Entry<Gff3BaseData, Coverage> featureCoverage : featureCounts.entrySet()) {
                final Gff3BaseData feature = featureCoverage.getKey();
                final Coverage coverage = featureCoverage.getValue();
                writer.writeRecord(new SingleStrandFeatureCoverage(feature, coverage.sense_count, true));
                if (feature.getStrand() != Strand.NONE) {
                    writer.writeRecord(new SingleStrandFeatureCoverage(feature, coverage.antisense_count, false));
                }
            }
        } catch (IOException ex) {
            throw new UserException(ex.getMessage());
        }

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    public static class FragmentCountWriter extends TableWriter<SingleStrandFeatureCoverage> {
        final String name;
        final String gene_id_key;
        public FragmentCountWriter(final Path file , final String name, final String gene_id_key) throws IOException {
            super(file, new TableColumnCollection("gene_id", "contig", "start", "stop", "strand", "sense_antisense", name+"_counts"));
            this.name = name;
            this.gene_id_key = gene_id_key;
        }

        protected void composeLine(final SingleStrandFeatureCoverage fragmentCount, final DataLine dataLine) {
            final String gene_label = fragmentCount.baseData.getAttributes().get(gene_id_key);
            dataLine.set("contig", fragmentCount.baseData.getContig())
                    .set("start", fragmentCount.baseData.getStart())
                    .set("stop", fragmentCount.baseData.getEnd())
                    .set("strand", fragmentCount.baseData.getStrand().encode())
                    .set("sense_antisense", fragmentCount.sense? "sense" : "antisense")
                    .set(name != null? name+"_counts" : "counts", fragmentCount.count, 2)
                    .set("gene_id", gene_label == null? "" : gene_label);
        }
    }

    public static class FragmentCountReader extends TableReader<SingleStrandFeatureCoverage> {

        public FragmentCountReader(final Path file) throws IOException {
            super(file);
        }
        @Override
        protected GeneExpressionEvaluation.SingleStrandFeatureCoverage createRecord(final DataLine dataLine) {
            return new GeneExpressionEvaluation.SingleStrandFeatureCoverage(new Gff3BaseData(dataLine.get("contig"),".", ".",
                    dataLine.getInt("start"), dataLine.getInt("stop"), Strand.decode(dataLine.get("strand")), 0,
                    Collections.singletonMap("ID", dataLine.get("gene_id"))), (float)dataLine.getDouble(6), dataLine.get("sense_antisense").equals("sense"));
        }
    }

    public static class Coverage {
        private float sense_count;
        private float antisense_count;

        public Coverage(final float sense_count, final float antisense_count) {
            this.sense_count = sense_count;
            this.antisense_count = antisense_count;
        }

        public void addSenseCount(float count) {
            sense_count += count;
        }

        public void addAntiSenseCount(float count) {
            antisense_count += count;
        }

        public float getSenseCount() {
            return sense_count;
        }

        public float getAntisenseCount() {
            return antisense_count;
        }
    }

    public static class SingleStrandFeatureCoverage {
        public final Gff3BaseData baseData;
        public final float count;
        public final boolean sense;

        SingleStrandFeatureCoverage(final Gff3BaseData baseData, final float count, final boolean sense) {
            this.baseData = baseData;
            this.count = count;
            this.sense = sense;
        }
    }

    static Strand getFragmentStrand(final GATKRead read, final TrancriptionRead trancriptionRead) {
        return (read.isFirstOfPair() == (trancriptionRead==TrancriptionRead.R1))? (read.isReverseStrand() ? Strand.NEGATIVE : Strand.POSITIVE) :
                                                                                    (read.isReverseStrand() ? Strand.POSITIVE : Strand.NEGATIVE);
    }

    enum TrancriptionRead {
        R1,
        R2
    }

    private static List<Interval> getMergedIntervals(final List<Interval> intervals) {
        final List<Interval> intervalsCopy = new ArrayList<>(intervals);
        Collections.sort(intervalsCopy);
        final IntervalList.IntervalMergerIterator mergeringIterator = new IntervalList.IntervalMergerIterator(intervalsCopy.iterator(), true, true, false);
        final List<Interval> merged = new ArrayList<>();
        while (mergeringIterator.hasNext()) {
            merged.add(mergeringIterator.next());
        }
        return merged;
    }
}
