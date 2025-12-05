package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.codecs.*;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;


import java.util.Collections;
/**
 * <p>Converts a read counts file (*.counts.tsv) to a depth file (*.rd.txt). </p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         A read counts file (*.counts.tsv).
 *     </li>
 *     <li>Optional:  Sample label of the input file. If not provided, this tool extracts it from the input file header.</li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         An output file containing converted read counts to depth file format.
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk ConvertCountsToDepthFile \
 *       --counts-filename input.counts.tsv \
 *       --output output.rd.txt.gz
 * </pre>
 */
@DocumentedFeature
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Converts a counts file (*.counts.tsv) to a depth file (*" + DepthEvidenceCodec.FORMAT_SUFFIX + ").",
        oneLineSummary = "Converts *.counts.tsv to " + DepthEvidenceCodec.FORMAT_SUFFIX,
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public class ConvertCountsToDepthFile extends FeatureWalker<SimpleCount> {
    public static final String COMPRESSION_LEVEL_NAME = "compression-level";
    public static final String COUNTS_FILE_PATH_ARG_FULL_NAME = "counts-filename";
    public static final String COUNTS_FILE_PATH_ARG_SHORT_NAME = StandardArgumentDefinitions.FEATURE_SHORT_NAME;
    public static final String SAMPLE_NAME_ARG_FULL_NAME = "sample-name";
    public static final String SAMPLE_NAME_ARG_SHORT_NAME = "n";


    @Argument(
            doc = "Input file of counts (*.counts.tsv).",
            fullName = COUNTS_FILE_PATH_ARG_FULL_NAME,
            shortName = COUNTS_FILE_PATH_ARG_SHORT_NAME
    )
    private FeatureInput<SimpleCount> countsPath;

    @Argument(
            doc = "Output file of depth (*" + DepthEvidenceCodec.FORMAT_SUFFIX + ").",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFilePath;

    @Argument(
            doc = "Sample label of the input file. " +
                  "If not provided, this tool extracts it from the input file header, " +
                  "and will throw an error if the header do not contain a sample name.",
            fullName = SAMPLE_NAME_ARG_FULL_NAME,
            shortName = SAMPLE_NAME_ARG_SHORT_NAME,
            optional = true
    )
    private String sampleName;

    @Argument(
            doc = "Output compression level",
            fullName = COMPRESSION_LEVEL_NAME,
            minValue = 0, maxValue = 9, optional = true
    )
    private int compressionLevel = 4;

    private FeatureOutputStream<DepthEvidence> outputSink;

    @Override
    protected boolean isAcceptableFeatureType(Class<? extends Feature> featureType) {
        return SimpleCount.class.isAssignableFrom(featureType);
    }

    @Override
    public void onTraversalStart() {
        final SampleLocatableMetadata header = (SampleLocatableMetadata)getDrivingFeaturesHeader();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();

        final FeatureOutputCodec<? extends Feature, ? extends FeatureSink<? extends Feature>> codec =
                FeatureOutputCodecFinder.find(outputFilePath);

        if (!DepthEvidence.class.isAssignableFrom(codec.getFeatureType()))
            throw new UserException.BadInput(
                    outputFilePath.getURIString() + " is not a valid format for " +
                    DepthEvidence.class.getSimpleName());

        if(sampleName == null)
        {
            sampleName = header.getSampleName();

            if (sampleName == null)
            {
                throw new UserException.BadInput(
                    "Sample name is not provided and the input file header does not contain sample name." +
                    "You may either provide a label for input sample " +
                    "(using --" + SAMPLE_NAME_ARG_FULL_NAME + " argument) or update input file header.");
            }
        }

        @SuppressWarnings("unchecked")
        final FeatureOutputCodec<DepthEvidence, FeatureOutputStream<DepthEvidence>> depthEvidenceCodec =
                (FeatureOutputCodec<DepthEvidence, FeatureOutputStream<DepthEvidence>>) codec;

        outputSink = depthEvidenceCodec.makeSink(
                outputFilePath,
                sequenceDictionary,
                Collections.singletonList(sampleName),
                compressionLevel
        );
    }

    @Override
    public void apply(SimpleCount feature, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final DepthEvidence depthEvidence = new DepthEvidence(
                feature.getContig(),
                feature.getStart(),
                feature.getEnd(),
                new int[] {feature.getCount()});

        outputSink.write(depthEvidence);
    }

    @Override
    public GATKPath getDrivingFeaturePath() {
        return new GATKPath(countsPath.getFeaturePath());
    }

    @Override
    public Object onTraversalSuccess() {
        outputSink.close();
        return super.onTraversalSuccess();
    }
}
