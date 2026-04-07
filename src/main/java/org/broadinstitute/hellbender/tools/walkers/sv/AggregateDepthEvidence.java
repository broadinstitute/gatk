package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.arrow.util.VisibleForTesting;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.*;
import org.broadinstitute.hellbender.tools.sv.aggregation.DepthEvidenceTest;
import org.broadinstitute.hellbender.tools.sv.cluster.PloidyTable;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterWalker;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.tsv.*;

import java.io.IOException;
import java.nio.charset.Charset;
import java.util.*;
import java.util.function.Function;

/**
 * <p>This tool assesses depth evidence for structural variants (SVs),annotating records with statistical metrics that
 * can be used to assess a variant's quality. The input VCF should
 * contain multiple samples with GT fields populated. Note that this tool only considers carrier status and does not
 * differentiate heterozygous from homozygous variant genotypes. For paired-end, split-read, and B-allele frequency
 * evidence metrics, see the AggregateSVEvidence tool.</p>
 *
 * <p>Detailed methodology can be found in the supplement of <a href="https://doi.org/10.1038/s41586-020-2287-8">Collins et al. 2020</a>.</p>
 *
 * <p>Briefly, for each variant the median binned read counts are aggreagated across samples. Carrier and control
 * samples are then compared using a permuted t-test. However, if the test is underpowered, </p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 *     <li>
 *         Depth evidence file
 *     </li>
 *     <li>
 *         Median binned read counts table
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk AggregateDepthEvidence \
 *      -V input.vcf.gz \
 *      -O output.vcf.gz \
 *      --median-coverage median_coverage.tsv \
 *      --rd-file all_samples.rd.txt.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        summary = "Read depth evidence assessment tool for copy number variants",
        oneLineSummary = "Read depth evidence assessment tool for copy number variants",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public final class AggregateDepthEvidence extends VariantWalker {

    public static final String DEPTH_EVIDENCE_FILE_PATH_LONG_NAME = "rd-file";
    public static final String MEDIAN_COUNTS_FILE_PATH_LONG_NAME = AggregateSVEvidence.MEDIAN_COVERAGE_LONG_NAME;
    public static final String LARGE_VARIANT_SIZE_LONG_NAME = "large-variant-size";
    public static final String LARGE_VARIANT_POINTS_LONG_NAME = "large-variant-points";
    public static final String LARGE_VARIANT_WINDOW_LONG_NAME = "large-variant-window";
    public static final String NUM_BINS_LONG_NAME = "num-bins";
    public static final String MAX_QUALITY_LONG_NAME = "max-qual";
    public static final String POWER_THRESHOLD_LONG_NAME = "power-threshold";

    @Argument(
            fullName = DEPTH_EVIDENCE_FILE_PATH_LONG_NAME,
            doc = "Indexed read counts file ending with " + DepthEvidenceCodec.FORMAT_SUFFIX + ".gz"
    )
    public GATKPath evidenceFile;

    @Argument(
            fullName = MEDIAN_COUNTS_FILE_PATH_LONG_NAME,
            doc = "Median counts file"
    )
    public GATKPath medianFile;

    /**
     * Expected format is tab-delimited and contains a header with the first column SAMPLE and remaining columns
     * contig names. Each row corresponds to a sample, with the sample ID in the first column and contig ploidy
     * integers in their respective columns.
     */
    @Argument(
            doc = "Sample ploidy table (.tsv), required if using a sites-only VCF",
            fullName = SVClusterWalker.PLOIDY_TABLE_LONG_NAME,
            optional = true
    )
    protected GATKPath ploidyTablePath;

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF"
    )
    public GATKPath outputVcf;

    @Argument(
            fullName = LARGE_VARIANT_SIZE_LONG_NAME,
            doc = "Large variant size",
            minValue = 0
    )
    public long largeVariantSize = 2500000L;

    @Argument(
            fullName = LARGE_VARIANT_POINTS_LONG_NAME,
            doc = "Large variant points",
            minValue = 1
    )
    public int largeVariantPoints = 500;

    @Argument(
            fullName = LARGE_VARIANT_WINDOW_LONG_NAME,
            doc = "Large variant window",
            minValue = 1
    )
    public int largeVariantWindow = 2000;

    @Argument(
            fullName = NUM_BINS_LONG_NAME,
            doc = "Number of bins to resample to",
            minValue = 1
    )
    public int numBins = 10;

    @Argument(
            fullName = MAX_QUALITY_LONG_NAME,
            doc = "Max quality score",
            minValue = 1
    )
    public int maxQual = 20;

    @Argument(
            fullName = POWER_THRESHOLD_LONG_NAME,
            doc = "Power threshold for permuted T tests",
            minValue = 0,
            maxValue = 1
    )
    public double powerThreshold = 0.8;

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private Map<String, Double> sampleMedians;
    private PloidyTable ploidyTable;
    private DepthMatrixLoader loader;
    private DepthEvidenceTest depthEvidenceTest;
    private FeatureDataSource<DepthEvidence> evidenceDataSource;

    @Override
    public void onTraversalStart() {
        evidenceDataSource = new FeatureDataSource<>(evidenceFile.toString());
        sampleMedians = loadMedianSampleCoverageTable();
        dictionary = getBestAvailableSequenceDictionary();
        loader = new DepthMatrixLoader(evidenceDataSource, numBins, largeVariantSize, largeVariantPoints, largeVariantWindow, null, dictionary);
        depthEvidenceTest = new DepthEvidenceTest(powerThreshold);
        writer = createVCFWriter(outputVcf);
        final VCFHeader header = createHeader(getHeaderForVariants());
        writer.writeHeader(header);
    }

    @VisibleForTesting
    protected Map<String, Double> loadMedianSampleCoverageTable() {
        final List<String> lines = IOUtils.readLines(BucketUtils.openFile(medianFile.toString()), Charset.defaultCharset());
        Utils.validate(lines.size() >= 2, "Median coverage file must contain at least two lines");
        final String[] samples = lines.get(0).split("\t");
        final String[] values = lines.get(1).split("\t");
        Utils.validate(samples.length == values.length,
                "Median file's first two lines must have the same number of columns");
        final Map<String, Double> sampleMedians = new HashMap<>();
        try {
            for (int i = 0; i < samples.length; i++) {
                sampleMedians.put(samples[i], Double.valueOf(values[i]));
            }
        } catch (final NumberFormatException nfe) {
            throw new UserException.BadInput(nfe.getMessage());
        }
        final List<String> vcfSamples = getHeaderForVariants().getSampleNamesInOrder();
        Utils.validate(sampleMedians.keySet().containsAll(vcfSamples), "Median counts table does not contain all samples in the VCF");
        return sampleMedians;
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        // Must be a CNV
        final StructuralVariantType svtype = variant.getStructuralVariantType();
        if (svtype != StructuralVariantType.DEL && svtype != StructuralVariantType.DUP && svtype != StructuralVariantType.CNV) {
            writer.add(variant);
            return;
        }
        SVCallRecord record = SVCallRecordUtils.create(variant, dictionary);
        final DepthMatrix depthMatrix = loader.load(new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionB()), sampleMedians);

        DepthEvidenceTest.DepthTestResult result = null;
        if (svtype == StructuralVariantType.DEL || svtype == StructuralVariantType.DUP) {
            result = depthEvidenceTest.test(record, depthMatrix);
            if (result != null) {
                record = depthEvidenceTest.applyToRecord(record, result, maxQual, dictionary);
            }
        }

        if (result == null) {
            writer.add(variant);
        } else {
            writer.add(SVCallRecordUtils.getVariantBuilder(record).make());
        }
    }

    @Override
    public Object onTraversalSuccess() {
        if (writer != null) {
            writer.close();
        }
        return null;
    }

    private VCFHeader createHeader(VCFHeader header) {
        final List<String> samples = header.getSampleNamesInOrder();
        if (samples.isEmpty()) {
            Utils.nonNull(ploidyTablePath, "Ploidy table required for sites-only VCFs");
            ploidyTable = new PloidyTable(ploidyTablePath.toPath());
            final List<String> featureSamples = ((SVFeaturesHeader) evidenceDataSource.getHeader()).getSampleNames();
            for (final String sample : samples) {
                Utils.validate(ploidyTable.contains(sample), "Ploidy table does not contain sample " + sample + " from the depth file");
            }
            header = new VCFHeader(header.getMetaDataInInputOrder(), featureSamples);
        }
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.READ_DEPTH_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Depth evidence quality (-log10 p-value)"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.READ_DEPTH_SECOND_MAX_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Depth evidence second highest quality across individual bins (-log10 p-value)"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.READ_DEPTH_MEDIAN_SEPARATION_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Median normalized depth difference between carriers and background samples"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 1, VCFHeaderLineType.Integer, "Expected copy number for ref genotype"));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
        return header;
    }
}