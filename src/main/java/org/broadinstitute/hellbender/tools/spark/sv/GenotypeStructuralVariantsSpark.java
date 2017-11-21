package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.base.Functions;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.avro.generic.GenericData;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
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
import org.broadinstitute.hellbender.tools.spark.sv.utils.ArraySVHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalLocator;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.ShardPartitioner;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVContext;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SerializableBiFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignmentUtils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.iterators.ArrayUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.reference.FastaReferenceWriter;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import scala.Tuple2;
import scala.Tuple3;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.io.UncheckedIOException;
import java.util.*;
import java.util.function.IntFunction;
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


    private final boolean ignoreReadsThatDontOverlapBreakingPoint = true;

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
    private int shardSize = 100000;

    @Argument(doc = "parallelism factor", shortName = "pfactor", fullName = "parallelismFactor", optional = true)
    private int parallelismFactor = 4;

    @Argument(doc = "insert size distribution",
            shortName = INSERT_SIZE_DISTR_SHORT_NAME,
            fullName = INSERT_SIZE_DISTR_FULL_NAME,
            optional = true)
    private InsertSizeDistribution insertSizeDistribution = new InsertSizeDistribution("N(309,149)");

    @ArgumentCollection
    private AlignmentPenalties penalties = new AlignmentPenalties();

    private VariantsSparkSource variantsSource;
    private ReadsSparkSource haplotypesAndContigsSource;
    public static final Pattern ASSEMBLY_NAME_ALPHAS = Pattern.compile("[a-zA-Z_]+");
    private int parallelism = 0;

    @Override
    public boolean requiresReads() {
        return false;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    private void setUp(final JavaSparkContext ctx) {

        parallelism = ctx.defaultParallelism() * parallelismFactor;
        variantsSource = new VariantsSparkSource(ctx);
        haplotypesAndContigsSource = new ReadsSparkSource(ctx);
    }

    private static class Localized<E> implements Locatable, Serializable {

        private static final long serialVersionUID = 1L;
        private final Locatable location;
        private final E element;

        public Localized(final E element, final Locatable location) {
            this.location = location;
            this.element = element;
        }

        public String getContig() {
            return location.getContig();
        }

        public int getStart() {
            return location.getStart();
        }

        public int getEnd() {
            return location.getEnd();
        }

        public E get() {
            return element;
        }

        @SuppressWarnings("unchecked")
        public static <E> Class<Localized<E>> getSubClass(final Class<E> elementClass) {
            return (Class<Localized<E>>) (Class<?>) Localized.class;
        }
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        setUp(ctx);
        final List<SimpleInterval> intervals = hasIntervals() ? getIntervals()
                : IntervalUtils.getAllIntervalsForReference(getReferenceSequenceDictionary());
        final TraversalParameters traversalParameters = new TraversalParameters(intervals, false);
        final JavaRDD<SVHaplotype> haplotypeAndContigs = haplotypesAndContigsSource
                .getParallelReads(haplotypesAndContigsFile, referenceArguments.getReferenceFileName(), traversalParameters)
                .repartition(parallelism)
                .map(r -> r.getReadGroup().equals("CTG") ? SVContig.of(r) : ArraySVHaplotype.of(r));
        final JavaRDD<SVContext> variants = variantsSource.getParallelVariantContexts(
                variantArguments.variantFiles.get(0).getFeaturePath(), getIntervals())
                .repartition(parallelism)
                .map(SVContext::of).filter(GenotypeStructuralVariantsSpark::structuralVariantAlleleIsSupported);
        final SparkSharder sharder = new SparkSharder(ctx, getReferenceSequenceDictionary(), intervals, shardSize, 0);

        final JavaRDD<Shard<Localized<SVContext>>> variantSharded = sharder.shard(variants.map(v -> new Localized<>(v, new SimpleInterval(v.getContig(), v.getStart(), v.getStart()))));
        final JavaRDD<Shard<Localized<SVHaplotype>>> haplotypesSharded = sharder.shard(haplotypeAndContigs.map(h -> new Localized<>(h, h.getVariantLocation())));
        final JavaPairRDD<Shard<Localized<SVContext>>, Shard<Localized<SVHaplotype>>> variantAndHaplotypesSharded = sharder.cogroup(variantSharded, haplotypesSharded);
        final JavaPairRDD<Localized<SVContext>, Iterable<Localized<SVHaplotype>>> variantAndHaplotypesLocalized = sharder
                .matchLeftByKey(variantAndHaplotypesSharded, x -> x.get().getUniqueID(), x -> x.get().getVariantId());
        final JavaPairRDD<SVContext, Iterable<SVHaplotype>> variantAndHaplotypes = variantAndHaplotypesLocalized
                .mapToPair(tuple -> new Tuple2<>(tuple._1().get(), Utils.stream(tuple._2().iterator()).map(Localized::get).collect(Collectors.toList())))
                .mapValues(haplotypes -> {
                    if (haplotypes.size() <= 2) {
                        return haplotypes;
                    } else {
                        final Set<SVHaplotype> result = new LinkedHashSet<>(haplotypes.size());
                        hapLoop: for (final SVHaplotype haplotype : haplotypes) {
                            if (haplotype.getName().equals("ref") || haplotype.getName().equals("alt")) {
                                result.add(haplotype);
                            } else {
                                final SVContig contig = (SVContig) haplotype;
                                if (contig.isPerfectAlternativeMap() && contig.isPerfectReferenceMap()) {
                                    continue;
                                }
                                for (final SVHaplotype added : result) {
                                    if (Arrays.equals(contig.getBases(), added.getBases())) {
                                        continue hapLoop;
                                    }
                                }
                                result.add(contig);
                            }
                        }
                        return result;
                    }
                });

        final ShardPartitioner<SVContext> partitioner = sharder.partitioner(SVContext.class, parallelism);
        final String fastqDir = this.fastqDir;
        final String fastqFileFormat = "asm%06d.fastq";
        final Broadcast<SVIntervalLocator> locatorBroadcast = ctx.broadcast(SVIntervalLocator.of(getReferenceSequenceDictionary()));
        final Broadcast<InsertSizeDistribution> insertSizeDistributionBroadcast = ctx.broadcast(insertSizeDistribution);

        final JavaPairRDD<SVContext, Tuple3<Iterable<SVHaplotype>, Iterable<Template>, Iterable<int[]>>> variantHaplotypesAndTemplates =
                variantAndHaplotypes.mapPartitionsToPair(it -> {
                    final AssemblyCollection assemblyCollection = new AssemblyCollection(fastqDir, fastqFileFormat);
                    return Utils.stream(it)
                            .map(t -> {
                                final SVIntervalLocator locator = locatorBroadcast.getValue();
                                final InsertSizeDistribution insertSizeDistribution = insertSizeDistributionBroadcast.getValue();
                                final SVIntervalTree<SimpleInterval> coveredReference = Utils.stream(t._2())
                                        .filter(h -> !h.isContig())
                                        .flatMap(h -> h.getReferenceAlignmentIntervals().stream())
                                        .flatMap(ai -> ai.referenceCoveredIntervals().stream())
                                        .collect(locator.toTreeCollector(si -> si));
                                final IntStream assemblyNumbers = Utils.stream(t._2())
                                        .filter(SVHaplotype::isContig)
                                        .map(SVHaplotype::getName)
                                        .map(n -> n.substring(0, n.indexOf(":")))
                                        .map(a -> ASSEMBLY_NAME_ALPHAS.matcher(a).replaceAll("0"))
                                        .mapToInt(Integer::parseInt)
                                        .sorted()
                                        .distinct();
                                final List<Template> allTemplates = assemblyNumbers
                                        .boxed()
                                        .flatMap(i -> assemblyCollection.templates(i).stream())
                                        .distinct()
                                        .collect(Collectors.toList());
                             //           .filter(tt -> tt.fragments().stream().map(f -> AlignmentInterval.encode(f.alignmentIntervals())).filter(s -> s.contains("chr10,2480")).count() > 0);
                                final Stream<int[]> allTemplatesMaxMappingQualities = allTemplates.stream()
                                        .map(tt -> tt.maximumMappingQuality(coveredReference, locator, insertSizeDistribution));
                                return new Tuple2<>(t._1(), new Tuple3<>(t._2(),
                                        (Iterable<Template>) allTemplates,
                                        (Iterable<int[]>) allTemplatesMaxMappingQualities.collect(Collectors.toList())));
                            }).iterator();
                });

        final JavaRDD<SVContext> calls = processVariants(variantHaplotypesAndTemplates, ctx);
        final VCFHeader header = composeOutputHeader();
        header.addMetaDataLine(new VCFInfoHeaderLine("READ_COUNT", 1, VCFHeaderLineType.Integer, "number of reads"));
        header.addMetaDataLine(new VCFInfoHeaderLine("CONTIG_COUNT", 1, VCFHeaderLineType.Integer, "number of contigs"));
        header.addMetaDataLine(new VCFInfoHeaderLine("REF_COVERED_RANGE", 2, VCFHeaderLineType.Integer, "ref haplotype offset covered range"));
        header.addMetaDataLine(new VCFInfoHeaderLine("REF_COVERED_RATIO", 1, VCFHeaderLineType.Float, "ref covered effective length with total length ratio"));
        header.addMetaDataLine(new VCFInfoHeaderLine("ALT_REF_COVERED_SIZE_RATIO", 1, VCFHeaderLineType.Float, "alt/ref covered length ratio"));
        SVVCFWriter.writeVCF(getAuthenticatedGCSOptions(), outputFile,
                referenceArguments.getReferenceFileName(), calls, header, logger);
        tearDown(ctx);
    }

    private JavaPairRDD<SVContext, Tuple2<Iterable<GATKRead>, Iterable<Template>>> joinFragments(final JavaPairRDD<SVContext, Iterable<GATKRead>> variantAndHaplotypes) {
        return variantAndHaplotypes
                .mapToPair(tuple -> new Tuple2<>(tuple._1(), new Tuple2<>(tuple._2(), tuple._1().getUniqueID())))
                .mapToPair(tuple -> new Tuple2<>(tuple._1(), new Tuple2<>(tuple._2()._1(), Utils.stream(tuple._2()._1())
                        .map(id -> String.format("%s/%s.fastq", fastqDir, id))
                        .flatMap(file -> SVFastqUtils.readFastqFile(file).stream())
                        .collect(Collectors.groupingBy(SVFastqUtils.FastqRead::getName))
                        .entrySet()
                        .stream()
                        .map(entry -> Template.create(entry.getKey(),
                                removeRepeatedReads(entry.getValue()), Template.Fragment::of)
                        ).collect(Collectors.toList()))));
    }

    private VCFHeader composeOutputHeader() {
        final List<String> samples = Collections.singletonList("sample");
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

    private JavaRDD<SVContext> processVariants(final JavaPairRDD<SVContext, Tuple3<Iterable<SVHaplotype>, Iterable<Template>, Iterable<int[]>>> input, final JavaSparkContext ctx) {
        final Broadcast<SAMSequenceDictionary> broadCastDictionary = ctx.broadcast(getReferenceSequenceDictionary());
        final SerializableBiFunction<String, byte[], File> imageCreator =  GenotypeStructuralVariantsSpark::createTransientImageFile;
        final AlignmentPenalties penalties = this.penalties;
        final boolean ignoreReadsThatDontOverlapBreakingPoint = this.ignoreReadsThatDontOverlapBreakingPoint;
        final InsertSizeDistribution insertSizeDistribution = this.insertSizeDistribution;
        return input.mapPartitions(it -> {
            final SAMSequenceDictionary dictionary = broadCastDictionary.getValue();
            final Stream<Tuple2<SVContext, Tuple3<Iterable<SVHaplotype>, Iterable<Template>, Iterable<int[]>>>> variants =
                    Utils.stream(it);
            final Map<String, File> imagesByName = new HashMap<>();
            final SAMFileHeader header = new SAMFileHeader();
            header.setSequenceDictionary(dictionary);
            final GenotypeLikelihoodCalculator genotypeCalculator =
                    new GenotypeLikelihoodCalculators().getInstance(2, 2);

            return variants
                    .filter(variant -> !Utils.isEmpty(variant._2()._1()))
                    .map(variant -> {
                final List<Template> allTemplates = Utils.stream(variant._2()._2())
                        .collect(Collectors.toList());
                final List<int[]> allMapQuals = Utils.stream(variant._2()._3())
                        .collect(Collectors.toList());

                final List<Template> allInformativeTemplates = new ArrayList<>(allTemplates.size());
                final List<int[]> allInformativeMapQuals = new ArrayList<>(allTemplates.size());
                for (int i = 0; i < allTemplates.size(); i++) {
                    for (final int mq : allMapQuals.get(i)) {
                        if (mq > 0) {
                            allInformativeMapQuals.add(allMapQuals.get(i));
                            allInformativeTemplates.add(allTemplates.get(i));
                            break;
                        }
                    }
                }


                final List<SVHaplotype> haplotypes = Utils.stream(variant._2()._1())
                        .collect(Collectors.toList());

                final List<Template> templates;
                final List<int[]> mapQuals;
                if (allInformativeTemplates.size() <= 5000) {
                   templates = allInformativeTemplates;
                   mapQuals = allInformativeMapQuals;
                } else {
                    templates = new ArrayList<>(5000);
                    mapQuals = new ArrayList<>(5000);
                    final Random rdn = new Random(variant._1().getUniqueID().hashCode());
                    for (int i = 0; i < 5000; i++) {
                        final int idx = rdn.nextInt(allInformativeTemplates.size());
                        templates.add(allInformativeTemplates.get(idx));
                        mapQuals.add(allInformativeMapQuals.get(idx));
                    }
                }
                final List<byte[]> sequences = templates.stream()
                                .flatMap(t -> t.fragments().stream())
                                .map(Template.Fragment::bases)
                                .collect(Collectors.toList());
                final List<GATKRead> reads = templates.stream().map(t -> {
                            final SAMRecord record = new SAMRecord(header);
                            record.setReadName(t.name());
                            record.setReadUnmappedFlag(true);
                            record.setReferenceName(variant._1().getContig());
                            record.setAlignmentStart(variant._1().getStart());
                            return new SAMRecordToGATKReadAdapter(record);
                        }).collect(Collectors.toList());
                //if (true) {
                //    final VariantContextBuilder newVariantBuilder = new VariantContextBuilder(variant._1());
                //    newVariantBuilder.attribute("READ_COUNT", sequences.size());
                //    newVariantBuilder.attribute("CONTIG_COUNT", haplotypes.size());
                //    return SVContext.of(newVariantBuilder.make());
                //}
                final SVHaplotype ref = haplotypes.stream().filter(h -> h.getName().equals("ref")).findFirst().get();
                final int[] refBreakPoints = calculateBreakPoints(ref, variant._1(), dictionary);
                final SVHaplotype alt = haplotypes.stream().filter(h -> h.getName().equals("alt")).findFirst().get();
                final int[] altBreakPoints = calculateBreakPoints(alt, variant._1(), dictionary);
                final GenotypingAllele refAllele = GenotypingAllele.of(ref, variant._1());
                final GenotypingAllele altAllele = GenotypingAllele.of(alt, variant._1());
                final int refHaplotypeIndex = haplotypes.indexOf(ref);
                final int altHaplotypeIndex = haplotypes.indexOf(alt);

                final TemplateHaplotypeScoreTable scoreTable =
                        new TemplateHaplotypeScoreTable(templates, haplotypes);
                  for (int h = 0; h < haplotypes.size(); h++) {
                    final SVHaplotype haplotype = haplotypes.get(h);

                    final boolean isContig = haplotype.isContig();
                    final String imageName = isContig
                            ? haplotype.getName()
                            : haplotype.getVariantId() + "/" + haplotype.getName();
                    final File imageFile =
                            imagesByName.computeIfAbsent(imageName, (in) -> imageCreator.apply(in, haplotype.getBases()));
                    final BwaMemIndex index = BwaMemIndexCache.getInstance(imageFile.toString());
                    final BwaMemAligner aligner = new BwaMemAligner(index);
                    final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(sequences);
                    final IntFunction<String> haplotypeName = i -> i == 0 ? haplotype.getName() : null;
                    for (int i = 0; i < templates.size(); i++) {
                        final Template template = templates.get(i);
                        final List<BwaMemAlignment> firstAlignment = alignments.get(i * 2);
                        final List<BwaMemAlignment> secondAlignment = alignments.get(i * 2 + 1);
                        final List<AlignmentInterval> firstIntervals = BwaMemAlignmentUtils.toAlignmentIntervals(firstAlignment, haplotypeName, template.fragments().get(0).length());
                        final List<AlignmentInterval> secondIntervals = BwaMemAlignmentUtils.toAlignmentIntervals(secondAlignment, haplotypeName, template.fragments().get(1).length());

                        final TemplateMappingInformation mappingInformation = TemplateMappingInformation.fromAlignments(
                                firstIntervals, template.fragments().get(0).length(),
                                secondIntervals, template.fragments().get(1).length());
                        scoreTable.setMappingInfo(h, i, mappingInformation);
                    }
                }
             //    resolve the missing mapping scores to the worst seen + a penalty.
                for (int t = 0; t < templates.size(); t++) {

                    final OptionalDouble bestFirstAlignmentScore = scoreTable.getWorstAlignmentScore(t, 0);
                    if (bestFirstAlignmentScore.isPresent()) {
                        final double missingAlignmentScore = bestFirstAlignmentScore.getAsDouble() - 0.1 * penalties.unmappedFragmentPenalty;
                        scoreTable.applyMissingAlignmentScore(t, 0, missingAlignmentScore);
                    }
                    final OptionalDouble bestSecondAlignmentScore = scoreTable.getWorstAlignmentScore(t, 1);
                    if (bestSecondAlignmentScore.isPresent()) {
                        final double missingAlignmentScore = bestSecondAlignmentScore.getAsDouble() - 0.1 * penalties.unmappedFragmentPenalty;
                        scoreTable.applyMissingAlignmentScore(t, 1, missingAlignmentScore);
                    }
                }
                scoreTable.calculateBestMappingScores();
                final ReadLikelihoods<GenotypingAllele> likelihoods = new ReadLikelihoods<>(SampleList.singletonSampleList("sample"),
                        new IndexedAlleleList<>(refAllele, altAllele), Collections.singletonMap("sample", reads));
                final ReadLikelihoods<GenotypingAllele> likelihoodsFirst = new ReadLikelihoods<>(SampleList.singletonSampleList("sample"),
                                new IndexedAlleleList<>(refAllele, altAllele), Collections.singletonMap("sample", reads));
                final ReadLikelihoods<GenotypingAllele> likelihoodsSecond = new ReadLikelihoods<>(SampleList.singletonSampleList("sample"),
                                new IndexedAlleleList<>(refAllele, altAllele), Collections.singletonMap("sample", reads));

                final LikelihoodMatrix<GenotypingAllele> sampleLikelihoods = likelihoods.sampleMatrix(0);
                        final LikelihoodMatrix<GenotypingAllele> sampleLikelihoodsFirst = likelihoodsFirst.sampleMatrix(0);
                        final LikelihoodMatrix<GenotypingAllele> sampleLikelihoodsSecond = likelihoodsSecond.sampleMatrix(0);

                        sampleLikelihoods.fill(Double.NEGATIVE_INFINITY);
                final int refIdx = likelihoods.indexOfAllele(refAllele);
                final int altIdx = likelihoods.indexOfAllele(altAllele);
                final List<SVContig> contigs = haplotypes.stream().filter(SVHaplotype::isContig).map(SVContig.class::cast)
                        .collect(Collectors.toList());
                for (int t = 0; t < templates.size(); t++) {
                    sampleLikelihoodsFirst.set(refIdx, t,
                            scoreTable.getMappingInfo(refHaplotypeIndex, t).firstAlignmentScore.orElse(0));
                    sampleLikelihoodsSecond.set(refIdx, t,
                            scoreTable.getMappingInfo(refHaplotypeIndex, t).secondAlignmentScore.orElse(0));
                    sampleLikelihoodsFirst.set(altIdx, t,
                            scoreTable.getMappingInfo(altHaplotypeIndex, t).firstAlignmentScore.orElse(0));
                    sampleLikelihoodsSecond.set(altIdx, t,
                            scoreTable.getMappingInfo(altHaplotypeIndex, t).secondAlignmentScore.orElse(0));
                }
                for (int h = 0; h < contigs.size(); h++) {
                    final SVContig contig = contigs.get(h);
                    final int mappingInfoIndex = haplotypes.indexOf(contig);
                    double haplotypeAltScore = AlignmentScore.calculate(contig.getLength(), contig.getAlternativeAlignment()).getValue();
                    double haplotypeRefScore = AlignmentScore.calculate(contig.getLength(), contig.getReferenceAlignment()).getValue();
                    final double base = Math.max(haplotypeAltScore, haplotypeRefScore);
                    haplotypeAltScore -= base;
                    haplotypeRefScore -= base;
                    for (int t = 0; t < templates.size(); t++) {
                        final boolean noAlignment = !scoreTable.getMappingInfo(mappingInfoIndex, t).firstAlignmentScore.isPresent()
                                && !scoreTable.getMappingInfo(mappingInfoIndex, t).secondAlignmentScore.isPresent();
                        if (noAlignment) continue;
                        final double firstMappingScore = scoreTable.getMappingInfo(mappingInfoIndex, t).firstAlignmentScore.orElse(0.0);
                        final double secondMappingScore = scoreTable.getMappingInfo(mappingInfoIndex, t).secondAlignmentScore.orElse(0.0);
                        if (firstMappingScore == scoreTable.bestMappingScorePerFragment[t][0]) {
                            sampleLikelihoodsFirst.set(refIdx, t,
                                    Math.max(firstMappingScore + haplotypeRefScore, sampleLikelihoodsFirst.get(refIdx, t)));
                            sampleLikelihoodsFirst.set(altIdx, t,
                                    Math.max(firstMappingScore + haplotypeAltScore, sampleLikelihoodsFirst.get(altIdx, t)));
                        }
                        if (secondMappingScore == scoreTable.bestMappingScorePerFragment[t][1]) {
                            sampleLikelihoodsSecond.set(refIdx, t,
                                    Math.max(secondMappingScore + haplotypeRefScore, sampleLikelihoodsSecond.get(refIdx, t)));
                            sampleLikelihoodsSecond.set(altIdx, t,
                                    Math.max(secondMappingScore + haplotypeAltScore, sampleLikelihoodsSecond.get(altIdx, t)));
                        }
                    }
                }
                for (int t = 0; t < templates.size(); t++) {
                    for (int f = 0; f < 2; f++) {
                        final int maxMq = mapQuals.get(t)[f];
                        final LikelihoodMatrix<GenotypingAllele> matrix = f == 0 ? sampleLikelihoodsFirst : sampleLikelihoodsSecond;
                        final double base = Math.max(matrix.get(refIdx, t), matrix.get(altIdx, t));
                        final int maxIndex = matrix.get(refIdx, t) == base ? refIdx : altIdx;
                        final int minIndex = maxIndex == refIdx ? altIdx : refIdx;
                        matrix.set(minIndex, t, Math.max(matrix.get(maxIndex, t) - 0.1 * maxMq, matrix.get(minIndex, t)));
                       //         Math.min(matrix.get(maxIndex, t), matrix.get(minIndex, t) - base)));
                    //    matrix.set(minIndex, t,  Math.min(matrix.get(maxIndex, t), matrix.get(minIndex, t) - base));
                    }
                }
                for (int j = 0; j < sampleLikelihoods.numberOfReads(); j++) {
                    final boolean considerFirstFragment;
                    final boolean considerSecondFragment;
                    if (ignoreReadsThatDontOverlapBreakingPoint) {
                        considerFirstFragment = scoreTable.getMappingInfo(refHaplotypeIndex, j).crossesBreakPoint(refBreakPoints, 0)
                                || scoreTable.getMappingInfo(altHaplotypeIndex, j).crossesBreakPoint(altBreakPoints, 0);
                        considerSecondFragment = scoreTable.getMappingInfo(refHaplotypeIndex, j).crossesBreakPoint(refBreakPoints, 1)
                                || scoreTable.getMappingInfo(altHaplotypeIndex, j).crossesBreakPoint(altBreakPoints, 1);
                    } else {
                        considerFirstFragment = true;
                        considerSecondFragment = true;
                    }
                    for (int i = 0; i < sampleLikelihoods.numberOfAlleles(); i++) {
                        sampleLikelihoods.set(i, j, (considerFirstFragment ? sampleLikelihoodsFirst.get(i, j) : 0) +
                                ((considerSecondFragment) ? sampleLikelihoodsSecond.get(i, j) : 0));

                    }
                }
                        for (int t = 0; t < templates.size(); t++) {
                            final double base = Math.max(sampleLikelihoods.get(refIdx, t), sampleLikelihoods.get(altIdx, t));
                            final int maxIndex = sampleLikelihoods.get(refIdx, t) == base ? refIdx : altIdx;
                            final int minIndex = maxIndex == refIdx ? altIdx : refIdx;
                            sampleLikelihoods.set(minIndex, t, Math.min(sampleLikelihoods.get(maxIndex, t), sampleLikelihoods.get(minIndex, t) - base));
                        }

                        for (int t = 0; t < templates.size(); t++) {
                    final TemplateMappingInformation refMapping = scoreTable.getMappingInfo(refHaplotypeIndex, t);
                    final TemplateMappingInformation altMapping = scoreTable.getMappingInfo(altHaplotypeIndex, t);
                    if (refMapping.pairOrientation.isProper() == altMapping.pairOrientation.isProper()) {
                        if (refMapping.pairOrientation.isProper() && (scoreTable.getMappingInfo(refHaplotypeIndex, t).crossesBreakPoint(refBreakPoints) || scoreTable.getMappingInfo(altHaplotypeIndex, t).crossesBreakPoint(altBreakPoints))) {
                            sampleLikelihoods.set(refIdx, t, sampleLikelihoods.get(refIdx, t) + insertSizeDistribution.logProbability(refMapping.insertSize.getAsInt()) / Math.log(10));
                            sampleLikelihoods.set(altIdx, t, sampleLikelihoods.get(altIdx, t) + insertSizeDistribution.logProbability(altMapping.insertSize.getAsInt()) / Math.log(10));
                        }
                    } else if (refMapping.pairOrientation.isProper()) {
                    //    sampleLikelihoods.set(altIdx, t, sampleLikelihoods.get(altIdx, t) - 0.1 * penalties.improperPairPenalty);
                    } else {
                    //    sampleLikelihoods.set(refIdx, t, sampleLikelihoods.get(refIdx, t) - 0.1 * penalties.improperPairPenalty);
                    }
                }
                int minRefPos = ref.getLength();
                int maxRefPos = 0;
                for (int t = 0; t < scoreTable.numberOfTemplates(); t++) {
                    final TemplateMappingInformation mappingInfo = scoreTable.getMappingInfo(refHaplotypeIndex, t);
                    if (mappingInfo.minCoordinate < minRefPos) {
                        minRefPos = mappingInfo.minCoordinate;
                    }
                    if (mappingInfo.maxCoordinate > maxRefPos) {
                        maxRefPos = mappingInfo.maxCoordinate;
                    }
                }
                minRefPos = Math.max(0, minRefPos - 1);
                maxRefPos = Math.min(ref.getLength(), maxRefPos + 1);
                final int dp = likelihoods.readCount();
                likelihoods.removeUniformativeReads(0.0);
                likelihoods.normalizeLikelihoods(true, -0.1 * penalties.maximumTemplateScoreDifference);
                final int[] ad = new int[2];
                for (int t = 0; t < scoreTable.numberOfTemplates(); t++) {
                    final TemplateMappingInformation refMapping  = scoreTable.getMappingInfo(refHaplotypeIndex, t);
                    final TemplateMappingInformation altMapping  = scoreTable.getMappingInfo(altHaplotypeIndex, t);
                    final double refScore = refMapping.firstAlignmentScore.orElse(0) + refMapping.secondAlignmentScore.orElse(0);
                    final double altScore = altMapping.firstAlignmentScore.orElse(0) + altMapping.secondAlignmentScore.orElse(0);
                    if (refScore > altScore) {
                        ad[0]++;
                    } else if (altScore > refScore) {
                        ad[1]++;
                    }

                }
                ad[0] = 0; ad[1] = 0;
                for (int r = 0; r < likelihoods.readCount(); r++) {
                    if (likelihoods.sampleMatrix(0).get(0, r) > likelihoods.sampleMatrix(0).get(1, r)) {
                        ad[0]++;
                    } else {
                        ad[1]++;
                    }
                }
                final GenotypeLikelihoods likelihoods1 = genotypeCalculator.genotypeLikelihoods(sampleLikelihoods);
                final int pl[] = likelihoods1.getAsPLs();
                final int bestGenotypeIndex = MathUtils.maxElementIndex(likelihoods1.getAsVector());
                final int gq = GATKVariantContextUtils.calculateGQFromPLs(pl);
                final List<Allele> genotypeAlleles = gq == 0
                        ? Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)
                        : (bestGenotypeIndex == 0
                           ? Arrays.asList(variant._1().getReference(), variant._1().getReference())
                           : ((bestGenotypeIndex == 1)
                              ? Arrays.asList(variant._1().getReference(), variant._1().getAlternateAllele(0))
                              : Arrays.asList(variant._1().getAlternateAllele(0), variant._1().getAlternateAllele(0))));
                final VariantContextBuilder newVariantBuilder = new VariantContextBuilder(variant._1());
                newVariantBuilder.attribute("READ_COUNT", allTemplates.size() * 2);
                newVariantBuilder.attribute("CONTIG_COUNT", haplotypes.size());
                if (minRefPos <= maxRefPos) {
                    newVariantBuilder.attribute("REF_COVERED_RANGE", new int[] { minRefPos + 1, maxRefPos - 1 });
                    newVariantBuilder.attribute("REF_COVERED_RATIO", ((double) maxRefPos - minRefPos) / ref.getLength());
                    newVariantBuilder.attribute("ALT_REF_COVERED_SIZE_RATIO", (alt.getLength() + (maxRefPos - minRefPos )- ref.getLength()) / ((double) maxRefPos - minRefPos));
                }
                newVariantBuilder.genotypes(
                        new GenotypeBuilder().name("sample")
                                .PL(pl)
                                .GQ(gq)
                                .DP(dp)
                                .AD(ad)
                                .alleles(genotypeAlleles).make());
                return SVContext.of(newVariantBuilder.make());
           }).iterator();
        });
    }

    private static int[] calculateBreakPoints(final SVHaplotype haplotype, final SVContext context, final SAMSequenceDictionary dictionary) {
        final List<AlignmentInterval> intervals = haplotype.getReferenceAlignmentIntervals();
        final List<SimpleInterval> breakPoints = context.getBreakPointIntervals(0, dictionary, false);
        final List<Integer> result = new ArrayList<>(breakPoints.size());
        for (final SimpleInterval breakPoint : breakPoints) {
            for (final AlignmentInterval interval : intervals) {
                final Cigar cigar = interval.cigarAlongReference();
                if (interval.referenceSpan.overlaps(breakPoint)) {
                    int refPos = interval.referenceSpan.getStart();
                    int hapPos = interval.startInAssembledContig;
                    for (final CigarElement element : cigar) {
                        final CigarOperator operator = element.getOperator();
                        final int length = element.getLength();
                        if (operator.consumesReferenceBases() && breakPoint.getStart() >= refPos && breakPoint.getStart() < refPos + length) {
                            if (operator.isAlignment()) {
                                result.add(hapPos + breakPoint.getStart() - refPos);
                            } else {
                                result.add(hapPos);
                            }
                        }
                        if (operator.consumesReferenceBases()) {
                            refPos += length;
                        }
                        if (operator.consumesReadBases() || operator.isClipping()) {
                            hapPos += length;
                        }
                    }
                }
            }
        }
        Collections.sort(result);
        return result.stream().mapToInt(i -> i).toArray();
    }

    private static List<AlignmentInterval> composeAlignmentIntervals(final SVContig haplotype, final Template template, final int fragmentIndex, final List<BwaMemAlignment> alignerOutput) {
        if (alignerOutput.isEmpty() || alignerOutput.size() == 1 && SAMFlag.READ_UNMAPPED.isSet(alignerOutput.get(0).getSamFlag())) {
            return Collections.emptyList();
        } else {
            final List<String> haplotypeNameList = Collections.singletonList(haplotype.getName());
            final Set<String> haplotypeNameSet = Collections.singleton(haplotype.getName());
            final int readLength = template.fragments().get(fragmentIndex).length();
            final List<AlignmentInterval> allIntervals = alignerOutput.stream().map(bwa -> new AlignmentInterval(bwa, haplotypeNameList, readLength)).collect(Collectors.toList());
            final List<AlignmentInterval> bestIntervals = FilterLongReadAlignmentsSAMSpark.pickBestConfigurations(new AlignedContig(template.name(), template.fragments().get(fragmentIndex).bases(), allIntervals, false), haplotypeNameSet, 0.0).get(0);
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

    private static boolean structuralVariantAlleleIsSupported(final SVContext ctx) {
        switch (ctx.getStructuralVariantType()) {
            case INS:
            case DEL:
                return true;
            default:
                return false;
        }
    }

    private static class GenotypingAllele extends Allele {

        private static final long serialVersionUID = 1L;

        private final SVHaplotype haplotype;

        private static GenotypingAllele of(final SVHaplotype haplotype, final SVContext context) {
            if (haplotype.isContig()) {
                return new GenotypingAllele(haplotype, "<" + haplotype.getName() + ">", false);
            } else if (haplotype.getName().equals(SVHaplotype.REF_HAPLOTYPE_NAME)) {
                return new GenotypingAllele(haplotype, context.getReference().getBaseString(), true);
            } else { // assume is "alt".
                return new GenotypingAllele(haplotype, context.getAlternateAlleles().get(0).getDisplayString(), false);
            }
        }

        protected GenotypingAllele(final SVHaplotype haplotype, final String basesString, final boolean isRef) {
            super(basesString, isRef);
            this.haplotype = haplotype;
        }

        private boolean isReference(final SVHaplotype haplotype) {
            return haplotype.getName().equals("ref");
        }
    }


    private void tearDown(final JavaSparkContext ctx) {
        variantsSource = null;
    }
}
