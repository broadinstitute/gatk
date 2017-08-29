package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredVariantInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.SparkSharder;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.tools.spark.utils.ShardedPairRDD;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.iterators.ArrayUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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
    public static final String PADDING_SHORT_NAME = "padding";
    public static final String PADDING_FULL_NAME = "padding";
    public static final String INSERT_SIZE_DISTR_SHORT_NAME = "insSize";
    public static final String INSERT_SIZE_DISTR_FULL_NAME = "insertSizeDistribution";

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private RequiredVariantInputArgumentCollection variantArguments = new RequiredVariantInputArgumentCollection();

    @Argument(doc = "fastq files location",
            shortName = FASTQ_FILE_DIR_SHORT_NAME,
            fullName  = FASTQ_FILE_DIR_FULL_NAME)
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

    @Argument(doc = "padding",
            shortName = PADDING_SHORT_NAME,
            fullName = PADDING_FULL_NAME,
            optional = true)
    private int padding = 300;

    @Argument(doc = "insert size distribution",
              shortName = INSERT_SIZE_DISTR_SHORT_NAME,
              fullName = INSERT_SIZE_DISTR_FULL_NAME,
              optional = true)
    private InsertSizeDistribution dist = new InsertSizeDistribution("N(309,149)");

    @Argument(doc = "pair-hmm implementation")
    private StructuralVariantPairHMMImplementation pairHmm = new StructuralVariantPairHMMImplementation("Affine(45,10)");

    private VariantsSparkSource variantsSource;
    private ReadsSparkSource haplotypesAndContigsSource;

    @Override
    public boolean requiresReads() {
        return true;
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
        final JavaRDD<GATKRead> haplotypeAndContigs = haplotypesAndContigsSource.getParallelReads(haplotypesAndContigsFile, referenceArguments.getReferenceFileName(), getIntervals(), 0);
        final JavaRDD<StructuralVariantContext> variants = variantsSource.getParallelVariantContexts(
                variantArguments.variantFiles.get(0).getFeaturePath(), getIntervals())
                .map(StructuralVariantContext::new);

//        final SparkSharder sharder = new SparkSharder(ctx, getReference(), getIntervals(), 10000);

//        final SparkSharder.ShardedRDD<GATKRead> sharedReads = sharder.shard(haplotypeAndContigs, null);
//        final SparkSharder.ShardedRDD<StructuralVariantContext> sharedVCs = sharder.shard(variants, null);
//        final ShardedPairRDD<StructuralVariantContext, GATKRead> variantsAndHaplotypes = sharder.join(sharedVCs, sharedReads);
//        final JavaPairRDD<StructuralVariantContext, List<GATKRead>> matchedVariantsAndHaplotypes = variantsAndHaplotypes.matchLeft(
//                (a, r) -> { r.hasAttribute("VC") && r.getAttributeAsString("VC").equals(a.getUniqueID())}
//        final JavaPairRDD<StructuralVariantGenotypingContext, Tuple2<List<GenotypingContig> contigs>
  ///      final JavaPairRDD<StructuralVariantContext, List<GATKRead>> sharedReads.groupByRight((r,l) -> r.overlaps(l), );

        final JavaPairRDD<SimpleInterval, GATKRead> readsPerShard;
        final JavaPairRDD<SimpleInterval, StructuralVariantContext> variantsPerShard;
        final JavaPairRDD<SimpleInterval, Tuple2<List<StructuralVariantContext>, List<GATKRead>>> variantsHaplotypesAndContigsPerShard;
        final JavaPairRDD<StructuralVariantContext, List<GATKRead>> variantsHaplotypeAndContigs;
        final JavaRDD<StructuralVariantContext> outputVariants = processVariants(variants, ctx);
        final VCFHeader header = composeOutputHeader();
        SVVCFWriter.writeVCF(getAuthenticatedGCSOptions(), outputFile,
                referenceArguments.getReferenceFile().getAbsolutePath(), outputVariants, header, logger);
        tearDown(ctx);
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

    private JavaRDD<StructuralVariantContext> processVariants(final JavaRDD<StructuralVariantContext> variants, final JavaSparkContext ctx) {
        final JavaPairRDD<StructuralVariantContext, List<String>> variantAndAssemblyIds = variants.mapToPair(v -> new Tuple2<>(v, v.assemblyIDs()));
        final String fastqDir = this.fastqDir;

        final JavaPairRDD<StructuralVariantContext, List<String>> variantAndAssemblyFiles = variantAndAssemblyIds
                .mapValues(v -> v.stream().map(id -> String.format("%s/%s.fastq", fastqDir, id)).collect(Collectors.toList()));
        final JavaPairRDD<StructuralVariantContext, List<SVFastqUtils.FastqRead>> variantReads = variantAndAssemblyFiles
                .mapValues(v -> v.stream().flatMap(file -> SVFastqUtils.readFastqFile(file).stream()).collect(Collectors.toList()));
        final JavaPairRDD<StructuralVariantContext, Map<String, List<SVFastqUtils.FastqRead>>> variantReadsByName = variantReads
                .mapValues(v -> v.stream().collect(Collectors.groupingBy(SVFastqUtils.FastqRead::getName)));
        final JavaPairRDD<StructuralVariantContext, List<Template>> variantTemplates = variantReadsByName.mapValues(map ->
            map.entrySet().stream().map(entry ->
                Template.create(entry.getKey(), removeRepeatedReads(entry.getValue()), read -> new Template.Fragment(read.getName(), read.getFragmentNumber().orElse(0), read.getBases(), ArrayUtils.toInts(read.getQuals(), false)))
            ).collect(Collectors.toList()));
        final BwaVariantTemplateScoreCalculator calculator = new BwaVariantTemplateScoreCalculator(ctx, dist);
        final int padding = this.padding;
        final ReferenceMultiSource reference = getReference();
        variantTemplates.collect().forEach(vt -> {
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
        return variants;
    }

    private static final List<SVFastqUtils.FastqRead> removeRepeatedReads(final List<SVFastqUtils.FastqRead> in) {
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

    private static boolean structuralVariantAlleleIsSupported(final StructuralVariantAllele structuralVariantAllele) {
        switch (structuralVariantAllele) {
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
