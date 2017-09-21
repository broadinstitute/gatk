package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredVariantInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.SparkSharder;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.FilterLongReadAlignmentsSAMSpark;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.ShardPartitioner;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.iterators.ArrayUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.reference.FastaReferenceWriter;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Created by valentin on 4/20/17.
 */
@CommandLineProgramProperties(summary = "genotype SV variant call files",
        oneLineSummary = "genotype SV variant call files",
        programGroup = StructuralVariationSparkProgramGroup.class)
public class GenotypeStructuralVariantsSpark extends GATKSparkTool {

    public static final String FASTQ_FILE_DIR_SHORT_NAME = "fastqDir";
    public static final String FASTQ_FILE_DIR_FULL_NAME = "fastqAssemblyDirectory";
    public static final String HAP_AND_CTG_FILE_SHORT_NAME = "assemblies";
    public static final String HAP_AND_CTG_FILE_FULL_NAME = "haplotypesAndContigsFile";
    public static final String SHARD_SIZE_SHORT_NAME = "shard";
    public static final String SHARD_SIZE_FULL_NAME = "shardSize";
    public static final String INSERT_SIZE_DISTR_SHORT_NAME = "insSize";
    public static final String INSERT_SIZE_DISTR_FULL_NAME = "insertSizeDistribution";

    private static final long serialVersionUID = 1L;

    private static final Comparator<List<AlignmentInterval>> LEFT_RIGHT_ALIGNMENT_COMPARATOR =
            Comparator.comparingInt(GenotypeStructuralVariantsSpark::unclippedStart)
                    .thenComparingInt(GenotypeStructuralVariantsSpark::unclippedEnd)
                    .thenComparingInt(GenotypeStructuralVariantsSpark::clippedStart)
                    .thenComparingInt(GenotypeStructuralVariantsSpark::clippedEnd);

    @ArgumentCollection
    private RequiredVariantInputArgumentCollection variantArguments = new RequiredVariantInputArgumentCollection();

    @Argument(doc = "fastq files location",
            shortName = FASTQ_FILE_DIR_SHORT_NAME,
            fullName = FASTQ_FILE_DIR_FULL_NAME)
    private String fastqDir = null;

    @Argument(doc = "assemblies SAM/BAM file location",
            shortName = HAP_AND_CTG_FILE_SHORT_NAME,
            fullName = HAP_AND_CTG_FILE_FULL_NAME)
    private String haplotypesAndContigsFile = null;

    @Argument(doc = "output VCF file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    private String outputFile = null;

    @Argument(doc = "shard size",
            shortName = SHARD_SIZE_SHORT_NAME,
            fullName = SHARD_SIZE_FULL_NAME,
            optional = true)
    private int shardSize = 1000;

    @Argument(doc = "insert size distribution",
            shortName = INSERT_SIZE_DISTR_SHORT_NAME,
            fullName = INSERT_SIZE_DISTR_FULL_NAME,
            optional = true)
    private InsertSizeDistribution insertSizeDistribution = new InsertSizeDistribution("N(309,149)");

    @ArgumentCollection
    private AlignmentPenalties penalties = new AlignmentPenalties();

    private VariantsSparkSource variantsSource;
    private ReadsSparkSource haplotypesAndContigsSource;

    @Override
    public boolean requiresReads() {
        return false;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    private void setUp(final JavaSparkContext ctx) {

        variantsSource = new VariantsSparkSource(ctx);
        haplotypesAndContigsSource = new ReadsSparkSource(ctx);
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        setUp(ctx);
        final List<SimpleInterval> intervals = hasIntervals() ? getIntervals()
                : IntervalUtils.getAllIntervalsForReference(getReferenceSequenceDictionary());
        final TraversalParameters traversalParameters = new TraversalParameters(intervals, false);
        final JavaRDD<GATKRead> haplotypeAndContigs = haplotypesAndContigsSource.getParallelReads(haplotypesAndContigsFile, referenceArguments.getReferenceFileName(), traversalParameters, 0);
        final JavaRDD<StructuralVariantContext> variants = variantsSource.getParallelVariantContexts(
                variantArguments.variantFiles.get(0).getFeaturePath(), getIntervals())
                .map(StructuralVariantContext::new).filter(GenotypeStructuralVariantsSpark::structuralVariantAlleleIsSupported);
        final SparkSharder sharder = new SparkSharder(ctx, getReferenceSequenceDictionary(), intervals, shardSize, 0, 10000);
        final JavaRDD<Shard<StructuralVariantContext>> variantSharded = sharder.shard(variants, StructuralVariantContext.class);
        final JavaRDD<Shard<GATKRead>> haplotypesSharded = sharder.shard(haplotypeAndContigs, GATKRead.class);
        final JavaPairRDD<Shard<StructuralVariantContext>, Shard<GATKRead>> variantAndHaplotypesSharded = sharder.cogroup(variantSharded, haplotypesSharded);
        final JavaPairRDD<StructuralVariantContext, Iterable<GATKRead>> variantAndHaplotypes = sharder.matchLeftByKey(variantAndHaplotypesSharded,
                StructuralVariantContext::getUniqueID,
                r -> r.getAttributeAsString("VC"));
        final ShardPartitioner<StructuralVariantContext> partitioner = sharder.partitioner(StructuralVariantContext.class, variantAndHaplotypes.getNumPartitions());
        final String fastqDir = this.fastqDir;
        final String fastqFileFormat = "asm%05d.fastq";
        final Pattern ASSEMBLY_NAME_ALPHAS = Pattern.compile("[a-zA-Z_]+");
        final JavaPairRDD<StructuralVariantContext, Tuple2<Iterable<GATKRead>, Iterable<Template>>> variantHaplotypesAndTemplates =
                variantAndHaplotypes.partitionBy(partitioner).mapPartitionsToPair(it -> {
                    final AssemblyCollection assemblyCollection = new AssemblyCollection(fastqDir, fastqFileFormat);
                    return Utils.stream(it)
                            .map(t -> {
                                final IntStream assemblyNumbers = Utils.stream(t._2())
                                        .filter(c -> c.getReadGroup().equals(ComposeStructuralVariantHaplotypesSpark.CONTIG_READ_GROUP))
                                        .map(GATKRead::getName)
                                        .map(n -> n.substring(0, n.indexOf(":")))
                                        .map(a -> ASSEMBLY_NAME_ALPHAS.matcher(a).replaceAll("0"))
                                        .mapToInt(Integer::parseInt)
                                        .sorted()
                                        .distinct();
                                final Stream<Template> allTemplates = assemblyNumbers
                                        .boxed()
                                        .flatMap(i -> assemblyCollection.templates(i).stream())
                                        .distinct();
                                return new Tuple2<>(t._1(), new Tuple2<>(t._2(), (Iterable<Template>) allTemplates.collect(Collectors.toList())));
                            }).iterator();
                });

        final JavaRDD<StructuralVariantContext> calls = processVariants(variantHaplotypesAndTemplates, ctx);
        final VCFHeader header = composeOutputHeader();
        SVVCFWriter.writeVCF(getAuthenticatedGCSOptions(), outputFile,
                referenceArguments.getReferenceFile().getAbsolutePath(), calls, header, logger);
        tearDown(ctx);
    }

    private JavaPairRDD<StructuralVariantContext, Tuple2<Iterable<GATKRead>, Iterable<Template>>> joinFragments(JavaPairRDD<StructuralVariantContext, Iterable<GATKRead>> variantAndHaplotypes) {
        return variantAndHaplotypes
                .mapToPair(tuple -> new Tuple2<>(tuple._1(), new Tuple2<>(tuple._2(), tuple._1().assemblyIDs())))
                .mapToPair(tuple -> new Tuple2<>(tuple._1(), new Tuple2<>(tuple._2()._1(), Utils.stream(tuple._2()._1())
                        .map(id -> String.format("%s/%s.fastq", fastqDir, id))
                        .flatMap(file -> SVFastqUtils.readFastqFile(file).stream())
                        .collect(Collectors.groupingBy(r -> r.getName()))
                        .entrySet()
                        .stream()
                        .map(entry -> Template.create(entry.getKey(),
                                removeRepeatedReads(entry.getValue()), read -> new Template.Fragment(read.getName(), read.getFragmentOrdinal(), read.getBases(), ArrayUtils.toInts(read.getQuals(), false)))
                        ).collect(Collectors.toList()))));
    }

    private VCFHeader composeOutputHeader() {
        final SAMFileHeader readHeader = getHeaderForReads();
        final List<String> samples = readHeader.getReadGroups().stream()
                .map(SAMReadGroupRecord::getSample)
                .distinct()
                .sorted()
                .collect(Collectors.toList());
        final VCFHeader result = new VCFHeader(Collections.emptySet(), samples);
        result.setSequenceDictionary(getReferenceSequenceDictionary());
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1, VCFHeaderLineType.Integer, "last base position of the variant"));
        return result;
    }


    private static File createTransientImageFile(final String name, final byte[] bases) {
        try {
            final File fasta = File.createTempFile(name, ".fasta");
            final File image = File.createTempFile(name, ".img");
            FastaReferenceWriter.writeSingleSequenceReference(fasta.toPath(), false, false, name, "", bases);
            BwaMemIndex.createIndexImageFromFastaFile(fasta.toString(), image.toString());
            fasta.delete();
            image.deleteOnExit();
            return image;
        } catch (final IOException ex) {
            throw new UncheckedIOException(ex);
        }
    }

    private JavaRDD<StructuralVariantContext> processVariants(final JavaPairRDD<StructuralVariantContext, Tuple2<Iterable<GATKRead>, Iterable<Template>>> input, final JavaSparkContext ctx) {
        return input.mapPartitions(it -> {
            final Stream<Tuple2<StructuralVariantContext, Tuple2<Iterable<GATKRead>, Iterable<Template>>>> variants = Utils.stream(it);
            final Map<String, File> imagesByName = new HashMap<>();
            final SAMFileHeader header = new SAMFileHeader();
            header.setSequenceDictionary(getReferenceSequenceDictionary());
            final GenotypeLikelihoodCalculator genotypeCalculator = new GenotypeLikelihoodCalculators().getInstance(2, 2);

            return variants.map(variant -> {
                final List<Template> templates = Utils.stream(variant._2()._2()).collect(Collectors.toList());
                final List<byte[]> sequences = templates.stream()
                        .flatMap(t -> t.fragments().stream())
                        .map(Template.Fragment::bases)
                        .collect(Collectors.toList());
                final List<SVHaplotype> haplotypes = Utils.stream(variant._2()._1())
                        .map(SVHaplotype::of).collect(Collectors.toList());
                final List<GATKRead> reads = templates.stream().map(t -> {
                    final SAMRecord record = new SAMRecord(header);
                    record.setReadName(t.name());
                    record.setReadUnmappedFlag(true);
                    record.setReferenceName(variant._1().getContig());
                    record.setAlignmentStart(variant._1().getStart());
                    return new SAMRecordToGATKReadAdapter(record);
                }).collect(Collectors.toList());
                final SVHaplotype ref = haplotypes.stream().filter(h -> h.getName().equals("ref")).findFirst().get();
                final SVHaplotype alt = haplotypes.stream().filter(h -> h.getName().equals("alt")).findFirst().get();
                final int refHaplotypeIndex = haplotypes.indexOf(ref);
                final int altHaplotypeIndex = haplotypes.indexOf(alt);

                final TemplateHaplotypeScoreTable scoreTable = new TemplateHaplotypeScoreTable(templates, haplotypes);
                for (int h = 0; h < haplotypes.size(); h++) {
                    final SVHaplotype haplotype = haplotypes.get(h);
                    final boolean isContig = haplotype.isContig();
                    final String imageName = isContig
                            ? haplotype.getName()
                            : haplotype.getVariantId() + "/" + haplotype.getName();
                    final File imageFile = imagesByName.computeIfAbsent(imageName, (in) ->
                            createTransientImageFile(in, haplotype.getBases()));
                    final BwaMemIndex index = BwaMemIndexCache.getInstance(imageFile.toString());
                    final BwaMemAligner aligner = new BwaMemAligner(index);
                    final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(sequences);
                    for (int i = 0; i < templates.size(); i++) {
                        final List<BwaMemAlignment> firstAlignment = alignments.get(i * 2);
                        final List<BwaMemAlignment> secondAlignment = alignments.get(i * 2 + 1);
                        final TemplateMappingInformation mappingInformation = calculateMappingInfo(haplotype,
                                templates.get(i), firstAlignment, secondAlignment);
                        scoreTable.setMappingInfo(h, i, mappingInformation);
                    }
                }
                // resolve the missing mapping scores to the worst seen + a penalty.
                for (int t = 0; t < templates.size(); t++) {
                    final OptionalDouble bestFirstAlignmentScore = scoreTable.getBestAlignmentScore(t, 0);
                    if (bestFirstAlignmentScore.isPresent()) {
                        final double missingAlignmentScore = bestFirstAlignmentScore.getAsDouble() - 0.1 * penalties.unmappedFragmentPenalty;
                        scoreTable.applyMissingAlignmentScore(t, 0, missingAlignmentScore);
                    }
                    final OptionalDouble bestSecondAlignmentScore = scoreTable.getBestAlignmentScore(t, 1);
                    if (bestSecondAlignmentScore.isPresent()) {
                        final double missingAlignmentScore = bestSecondAlignmentScore.getAsDouble() - 0.1 * penalties.unmappedFragmentPenalty;
                        scoreTable.applyMissingAlignmentScore(t, 1, missingAlignmentScore);
                    }
                }
                final ReadLikelihoods<SVHaplotype> likelihoods = new ReadLikelihoods<>(SampleList.singletonSampleList("sample"),
                        new IndexedAlleleList<>(ref, alt), Collections.singletonMap("sample", reads));
                final LikelihoodMatrix<SVHaplotype> sampleLikelihoods = likelihoods.sampleMatrix(0);
                sampleLikelihoods.fill(Double.NEGATIVE_INFINITY);
                final int refIdx = likelihoods.indexOfAllele(ref);
                final int altIdx = likelihoods.indexOfAllele(alt);
                for (final SVHaplotype haplotype : haplotypes) {
                    final double haplotypeRefScore = haplotype.getReferenceScore();
                    final double haplotypeAltScore = haplotype.getAlternativeScore();
                    for (int t = 0; t < templates.size(); t++) {
                        for (int h = 0; h < haplotypes.size(); h++) {
                            final double firstMappingScore = scoreTable.getMappingInfo(h, t).firstAlignmentScore.getAsDouble();
                            final double secondMappingScore = scoreTable.getMappingInfo(h, t).secondAlignmentScore.getAsDouble();
                            sampleLikelihoods.set(refIdx, t,
                                    Math.max(firstMappingScore + secondMappingScore + haplotypeRefScore, sampleLikelihoods.get(refIdx, t)));
                            sampleLikelihoods.set(altIdx, t * 2,
                                    Math.max(firstMappingScore + secondMappingScore + haplotypeAltScore, sampleLikelihoods.get(altIdx, t * 2)));
                        }
                    }
                }
                for (int t = 0; t < templates.size(); t++) {
                    final TemplateMappingInformation refMapping = scoreTable.getMappingInfo(refHaplotypeIndex, t);
                    final TemplateMappingInformation altMapping = scoreTable.getMappingInfo(altHaplotypeIndex, t);
                    if (refMapping.pairOrientation.isProper() == altMapping.pairOrientation.isProper()) {
                        if (refMapping.pairOrientation.isProper()) {
                            sampleLikelihoods.set(refIdx, t, sampleLikelihoods.get(refIdx, t) + insertSizeDistribution.logProbability(refMapping.insertSize.getAsInt()) / Math.log(10));
                            sampleLikelihoods.set(altIdx, t, sampleLikelihoods.get(altIdx, t) + insertSizeDistribution.logProbability(refMapping.insertSize.getAsInt()) / Math.log(10));

                        }
                    } else if (refMapping.pairOrientation.isProper()) {
                        sampleLikelihoods.set(altIdx, t, sampleLikelihoods.get(altIdx, t) - 0.1 * penalties.improperPairPenalty);
                    } else {
                        sampleLikelihoods.set(refIdx, t, sampleLikelihoods.get(refIdx, t) - 0.1 * penalties.improperPairPenalty);
                    }
                }
                likelihoods.normalizeLikelihoods(true, -0.1 * penalties.maximumLikelihoodDiffernenceCap);
                final GenotypeLikelihoods likelihoods1 = genotypeCalculator.genotypeLikelihoods(sampleLikelihoods);
                final VariantContextBuilder newVariantBuilder = new VariantContextBuilder(variant._1());
                newVariantBuilder.genotypes(
                        new GenotypeBuilder().name("sample")
                                .PL(likelihoods1.getAsPLs()).make());
                return new StructuralVariantContext(newVariantBuilder.make());
            }).iterator();
        });

        //final ReferenceMultiSource reference = getReference();
        // return null;// input
        //.partitionBy(new );
        //map(t -> processVariant(t._1(), t._2()._1(), t._2()._2()));
    }

    private static StructuralVariantContext processVariant(final StructuralVariantContext variant, final Iterable<GATKRead> haplotypesAndContigs, final Iterable<Template> templates) {
        Haplotype refHaplotype = null;
        Haplotype altHaplotype = null;
        final List<GenotypingContig> alignedContigs = new ArrayList<>();
        for (final GATKRead haplotypeOrContig : haplotypesAndContigs) {
            final String call = haplotypeOrContig.getAttributeAsString("HP");
            if (call == null) {
                throw new IllegalArgumentException("missing call in record: " + haplotypeOrContig);
            }
            if (haplotypeOrContig.getReadGroup().equals("HAP")) {
                if (call == ".") {
                    throw new IllegalArgumentException("the call cannot be . for a haplotype record: " + haplotypeOrContig);
                } else if (call == "ref") {
                    refHaplotype = Haplotype.fromGATKRead(haplotypeOrContig, true);
                } else if (call == "alt") {
                    altHaplotype = Haplotype.fromGATKRead(haplotypeOrContig, false);
                } else {
                    throw new IllegalArgumentException("illegal call string: " + call);
                }
            } else if (haplotypeOrContig.getReadGroup().equals("CTG")) {
                final GenotypingContig genotypingContig = new GenotypingContig(haplotypeOrContig);
                alignedContigs.add(genotypingContig);
            }

        }
        //final BwaVariantTemplateScoreCalculator calculator = new BwaVariantTemplateScoreCalculator();

     /*   input.collect().forEach(vt -> {
            logger.debug("Doing  " + vt._1().getContig() + ":" + vt._1().getStart());
            if (structuralVariantAlleleIsSupported(vt._1().getStructuralVariantAllele())) {
                System.err.println("" + vt._1().getContig() + ":" + vt._1().getStart() + " " + vt._2().size() + " templates ");
                final List<Haplotype> haplotypes = new ArrayList<>(2);
                haplotypes.add(vt._1().composeHaplotype(0, padding, reference));
                haplotypes.add(vt._1().composeHaplotype(1, padding, reference));
                final TemplateHaplotypeScoreTable table = new TemplateHaplotypeScoreTable(vt._2(), haplotypes);
                calculator.calculate(table);
                table.dropUninformativeTemplates();
                System.err.println("table.0 is " + Arrays.toString(table.getRow(0)));
                System.err.println("table.1 is " + Arrays.toString(table.getRow(1)));
                System.err.println("table.n is " + Arrays.toString(table.templates().stream().map(Template::name).toArray()));
                final GenotypeLikelihoods likelihoods = table.calculateGenotypeLikelihoods(2);
                System.err.println("likelihoods = " + likelihoods.getAsString());
            } else {
                System.err.println("" + vt._1().getContig() + ":" + vt._1().getStart() + " not supported ");
            }
        });
        return variants;*/
        return null;
    }

    private TemplateMappingInformation calculateMappingInfo(final SVHaplotype haplotype, final Template template, final List<BwaMemAlignment> first, final List<BwaMemAlignment> second) {

        final List<AlignmentInterval> firstIntervals = composeAlignmentIntervals(haplotype, template, 0, first);
        final List<AlignmentInterval> secondIntervals = composeAlignmentIntervals(haplotype, template, 1, second);
        final double firstScore = score(template.fragments().get(0).length(), firstIntervals);
        final double secondScore = score(template.fragments().get(1).length(), secondIntervals);

        if (firstIntervals.isEmpty() && secondIntervals.isEmpty()) {
            return new TemplateMappingInformation();
        } else if (secondIntervals.isEmpty()) {
            return new TemplateMappingInformation(firstScore, true);
        } else if (firstIntervals.isEmpty()) {
            return new TemplateMappingInformation(secondScore, false);
        } else {
            final Pair<List<AlignmentInterval>, List<AlignmentInterval>> sortedAlignments
                    = sortLeftRightAlignments(firstIntervals, secondIntervals);
            final Pair<SVFastqUtils.Strand, SVFastqUtils.Strand> strands = new ImmutablePair<>(
                    strand(sortedAlignments.getLeft()), strand(sortedAlignments.getRight()));
            final ReadPairOrientation orientation = ReadPairOrientation.fromStrands(strands.getLeft(), strands.getRight());
            if (orientation.isProper()) {
                return new TemplateMappingInformation(firstScore, secondScore,
                        unclippedEnd(sortedAlignments.getRight()) - unclippedStart(sortedAlignments.getLeft()));
            } else {
                return new TemplateMappingInformation(firstScore, secondScore, orientation);
            }
        }
    }

    private static double score(final int length, final List<AlignmentInterval> intervals) {
        return intervals.isEmpty() ? Double.NaN : AlignmentScore.calculate(length, intervals).getValue();

    }

    private static Pair<List<AlignmentInterval>, List<AlignmentInterval>> sortLeftRightAlignments(
            final List<AlignmentInterval> first, final List<AlignmentInterval> second) {
        final int cmp = LEFT_RIGHT_ALIGNMENT_COMPARATOR.compare(first, second);
        return (cmp <= 0) ? new ImmutablePair<>(first, second)
                : new ImmutablePair<>(second, first);
    }

    private static SVFastqUtils.Strand strand(final List<AlignmentInterval> firstIntervals) {
        final int mappedBasesOrientation = firstIntervals.stream()
                .mapToInt(ai -> (ai.forwardStrand ? 1 : -1) * CigarUtils.countAlignedBases(ai.cigarAlong5to3DirectionOfContig))
                .sum();
        if (mappedBasesOrientation != 0) {
            return mappedBasesOrientation < 0 ? SVFastqUtils.Strand.NEGATIVE : SVFastqUtils.Strand.POSITIVE;
        } else { // tie-break:
            final Comparator<AlignmentInterval> comparator0 = Comparator.comparingInt(a -> CigarUtils.countAlignedBases(a.cigarAlong5to3DirectionOfContig));
            final Comparator<AlignmentInterval> comparator1 = comparator0.thenComparingInt(a -> a.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                    .filter(ce -> ce.getOperator() == CigarOperator.I)
                    .mapToInt(CigarElement::getLength)
                    .sum());
            final Comparator<AlignmentInterval> comparator = comparator1.thenComparingInt(a -> a.startInAssembledContig).reversed();

            final boolean forwardStrand = firstIntervals.stream().sorted(comparator).findFirst().get().forwardStrand;
            return forwardStrand ? SVFastqUtils.Strand.POSITIVE : SVFastqUtils.Strand.NEGATIVE;
        }
    }

    private static int unclippedStart(final List<AlignmentInterval> intervals) {
        return intervals.stream()
                .mapToInt(ai -> ai.referenceSpan.getStart() - (ai.forwardStrand ? CigarUtils.countLeftClippedBases(ai.cigarAlong5to3DirectionOfContig)
                                                                                : CigarUtils.countRightClippedBases(ai.cigarAlong5to3DirectionOfContig)))
                .min().getAsInt();
    }

    private static int clippedStart(final List<AlignmentInterval> intervals) {
        return intervals.stream()
                .mapToInt(ai -> ai.referenceSpan.getStart())
                .min().getAsInt();
    }

    private static int unclippedEnd(final List<AlignmentInterval> intervals) {
        return intervals.stream()
                .mapToInt(ai -> ai.referenceSpan.getEnd() + (ai.forwardStrand ? CigarUtils.countRightClippedBases(ai.cigarAlong5to3DirectionOfContig)
                                                                              : CigarUtils.countLeftClippedBases(ai.cigarAlong5to3DirectionOfContig)))
                .max().getAsInt();
    }

    private static int clippedEnd(final List<AlignmentInterval> intervals) {
        return intervals.stream()
                .mapToInt(ai -> ai.referenceSpan.getEnd())
                .max().getAsInt();
    }

    private static List<AlignmentInterval> composeAlignmentIntervals(final SVHaplotype haplotype, final Template template, final int fragmentIndex, final List<BwaMemAlignment> alignerOutput) {
        if (alignerOutput.isEmpty() || alignerOutput.size() == 1 && SAMFlag.READ_UNMAPPED.isSet(alignerOutput.get(0).getSamFlag())) {
            return Collections.emptyList();
        } else {
            final List<String> haplotypeNameList = Collections.singletonList(haplotype.getName());
            final Set<String> haplotypeNameSet = Collections.singleton(haplotype.getName());
            final int readLength = template.fragments().get(fragmentIndex).length();
            final List<AlignmentInterval> allIntervals = alignerOutput.stream().map(bwa -> new AlignmentInterval(bwa, haplotypeNameList, readLength)).collect(Collectors.toList());
            final List<AlignmentInterval> bestIntervals = FilterLongReadAlignmentsSAMSpark.pickBestConfigurations(new AlignedContig(template.name(), template.fragments().get(fragmentIndex).bases(), allIntervals, false), haplotypeNameSet).get(0);
            return bestIntervals;
        }
    }

    private static List<SVFastqUtils.FastqRead> removeRepeatedReads(final List<SVFastqUtils.FastqRead> in) {
        if (in.size() <= 2) {
            return in;
        } else {
            boolean[] repeats = new boolean[in.size()];
            for (int i = 0; i < repeats.length; i++) {
                if (repeats[i]) continue;
                final SVFastqUtils.FastqRead iread = in.get(i);
                for (int j = i + 1; j < repeats.length; j++) {
                    if (repeats[j]) continue;
                    final SVFastqUtils.FastqRead jread = in.get(j);
                    if (Arrays.equals(iread.getBases(), jread.getBases()) &&
                            Arrays.equals(iread.getQuals(), jread.getQuals())) {
                        repeats[j] = true;
                    }
                }
            }
            return IntStream.range(0, repeats.length)
                    .filter(i -> !repeats[i])
                    .mapToObj(in::get)
                    .collect(Collectors.toList());
        }
    }

    private static boolean structuralVariantAlleleIsSupported(final StructuralVariantContext ctx) {
        switch (ctx.getStructuralVariantType()) {
            case INS:
            case DEL:
                return true;
            default:
                return false;
        }
    }



    private void tearDown(final JavaSparkContext ctx) {
        variantsSource = null;
    }
}
