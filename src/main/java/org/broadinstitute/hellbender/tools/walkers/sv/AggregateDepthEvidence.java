package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.arrow.util.VisibleForTesting;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.nio.charset.Charset;
import java.util.*;

/**
 * AggregateDepthEvidence
 */
@DocumentedFeature
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Batch CNV Interval Genotyping - GATK VariantWalker Port of RdTest.R",
        oneLineSummary = "Batch CNV Interval Genotyping - GATK Java VariantWalker port",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public final class AggregateDepthEvidence extends VariantWalker {

    public static final String DEPTH_EVIDENCE_FILE_PATH_LONG_NAME = "rd-file";
    public static final String MEDIAN_COUNTS_FILE_PATH_LONG_NAME = "median-counts";
    public static final String LARGE_VARIANT_SIZE_LONG_NAME = "large-variant-size";
    public static final String LARGE_VARIANT_POINTS_LONG_NAME = "large-variant-points";
    public static final String LARGE_VARIANT_WINDOW_LONG_NAME = "large-variant-window";
    public static final String NUM_BINS_LONG_NAME = "num-bins";
    public static final String MAX_QUALITY_LONG_NAME = "max-qual";
    public static final String POWER_THRESHOLD_LONG_NAME = "power-thresh";

    @Argument(
            fullName = DEPTH_EVIDENCE_FILE_PATH_LONG_NAME,
            doc = "Read counts file ending with " + DepthEvidenceCodec.FORMAT_SUFFIX + "' (may be gzipped)."
    )
    public GATKPath evidenceFile;

    @Argument(
            fullName = MEDIAN_COUNTS_FILE_PATH_LONG_NAME,
            doc = "Median counts file"
    )
    public GATKPath medianFile;

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
    public int maxQual = 99;

    @Argument(
            fullName = POWER_THRESHOLD_LONG_NAME,
            doc = "Power threshold for permuted T tests",
            minValue = 0,
            maxValue = 1
    )
    public double powerThreshold = 0.8;

    private VariantContextWriter writer;
    private SAMSequenceDictionary samSequenceDictionary;
    private Map<String, Double> sampleMedians;
    private DepthMatrixLoader loader;
    private DepthEvidenceTest depthEvidenceTest;
    private FeatureDataSource<DepthEvidence> evidenceDataSource;

    @Override
    public void onTraversalStart() {
        samSequenceDictionary = getBestAvailableSequenceDictionary();
        evidenceDataSource = new FeatureDataSource<>(evidenceFile.toString());
        sampleMedians = loadMedianSampleCoverageTable();
        final SAMSequenceDictionary dictionary = getBestAvailableSequenceDictionary();
        loader = new DepthMatrixLoader(evidenceDataSource, numBins, largeVariantSize, largeVariantPoints, largeVariantWindow, dictionary);
        depthEvidenceTest = new DepthEvidenceTest(powerThreshold);
        writer = createVCFWriter(outputVcf);
        writer.writeHeader(createHeader(getHeaderForVariants()));
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
        final Set<String> vcfSamples = new HashSet<>(getHeaderForVariants().getSampleNamesInOrder());
        Utils.validate(vcfSamples.containsAll(sampleMedians.keySet()), "Median counts table does not contain all samples in the VCF");
        return sampleMedians;
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        // Must be a bi-allelic CNV
        final StructuralVariantType svtype = variant.getStructuralVariantType();
        if (svtype != StructuralVariantType.DEL && svtype != StructuralVariantType.DUP) {
            writer.add(variant);
            return;
        }
        final SVCallRecord record = SVCallRecordUtils.create(variant, samSequenceDictionary);
        final DepthMatrix depthMatrix = loader.load(new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionB()), sampleMedians);
        final DepthEvidenceTest.AggregateDepthResult result = depthEvidenceTest.test(record, depthMatrix);
        if (result == null ) {
            writer.add(variant);
        } else {
            final VariantContext outputVariant = SVCallRecordUtils.getVariantBuilder(record).make();
            final VariantContextBuilder builder = new VariantContextBuilder(outputVariant);
            builder.attribute(GATKSVVCFConstants.READ_DEPTH_QUALITY_ATTRIBUTE, Math.min(-10. * Math.log10(result.pValue()), maxQual));
            builder.attribute(GATKSVVCFConstants.READ_DEPTH_SECOND_MAX_QUALITY_ATTRIBUTE, Math.min(-10. * Math.log10(result.secondMaxP()), maxQual));
            builder.attribute(GATKSVVCFConstants.READ_DEPTH_MEDIAN_SEPARATION_ATTRIBUTE, result.medianSeparation());
            writer.add(builder.make());
        }
    }

    @Override
    public Object onTraversalSuccess() {
        if (writer != null) {
            writer.close();
        }
        return null;
    }

    private VCFHeader createHeader(final VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.READ_DEPTH_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Depth evidence quality"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.READ_DEPTH_SECOND_MAX_QUALITY_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Depth evidence second highest quality across individual bins"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.READ_DEPTH_MEDIAN_SEPARATION_ATTRIBUTE, 1, VCFHeaderLineType.Float, "Median normalized depth difference between carriers and background samples"));
        return header;
    }
}