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
import org.broadinstitute.hellbender.utils.codecs.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Prints SV evidence records. Can be used with -L to retrieve records on a set of intervals. Supports streaming input
 * from GCS buckets.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Coordinate-sorted and indexed evidence file URI
 *     </li>
 *     <li>
 *         Reference sequence dictionary
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Coordinate-sorted evidence file, automatically indexed if ending with ".gz"
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk PrintSVEvidence \
 *       --evidence-file gs://my-bucket/batch_name.SR.txt.gz \
 *       -L intervals.bed \
 *       --sequence-dictionary ref.dict \
 *       -O local.SR.txt.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Prints SV evidence records to a file",
        oneLineSummary = "Prints SV evidence records",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public final class PrintSVEvidence <F extends Feature> extends FeatureWalker<F> {

    public static final String EVIDENCE_FILE_NAME = "evidence-file";
    public static final String COMPRESSION_LEVEL_NAME = "compression-level";

    @Argument(
            doc = "Input file URI with extension '"
                    + SplitReadEvidenceCodec.FORMAT_SUFFIX + "', '"
                    + DiscordantPairEvidenceCodec.FORMAT_SUFFIX + "', '"
                    + LocusDepthCodec.FORMAT_SUFFIX + "', '"
                    + BafEvidenceCodec.FORMAT_SUFFIX + "', or '"
                    + DepthEvidenceCodec.FORMAT_SUFFIX + "' (may be gzipped). "
                    + "Can also handle bci rather than txt files.",
            fullName = EVIDENCE_FILE_NAME
    )
    private GATKPath inputFilePath;

    @Argument(
            doc = "Output file with an evidence extension matching the input. Will be indexed if it has a " +
                    "block-compressed extension (e.g. '.gz' or '.bci').",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFilePath;

    @Argument(
            doc = "Output compression level",
            fullName = COMPRESSION_LEVEL_NAME,
            minValue = 0, maxValue = 9, optional = true
    )
    private int compressionLevel = 4;

    @Argument(doc = "List of sample names", fullName = "sample-names", optional = true)
    private List<String> sampleNames = new ArrayList<>();

    private FeatureSink<F> outputSink;
    private Class<F> evidenceClass;

    private static final List<FeatureOutputCodec<? extends Feature, ? extends FeatureSink<?>>> outputCodecs =
            new ArrayList<>(10);
    static {
        outputCodecs.add(new BafEvidenceCodec());
        outputCodecs.add(new DepthEvidenceCodec());
        outputCodecs.add(new DiscordantPairEvidenceCodec());
        outputCodecs.add(new LocusDepthCodec());
        outputCodecs.add(new SplitReadEvidenceCodec());
        outputCodecs.add(new BafEvidenceBCICodec());
        outputCodecs.add(new DepthEvidenceBCICodec());
        outputCodecs.add(new DiscordantPairEvidenceBCICodec());
        outputCodecs.add(new LocusDepthBCICodec());
        outputCodecs.add(new SplitReadEvidenceBCICodec());
    }

    @Override
    @SuppressWarnings("unchecked")
    protected boolean isAcceptableFeatureType( final Class<? extends Feature> featureType ) {
        for ( final FeatureOutputCodec<?, ?> codec : outputCodecs ) {
            if ( featureType.equals(codec.getFeatureType()) ) {
                evidenceClass = (Class<F>)featureType;
                return true;
            }
        }
        return false;
    }

    @Override
    public GATKPath getDrivingFeaturePath() {
        return inputFilePath;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        initializeOutput();
    }

    @SuppressWarnings("unchecked")
    private void initializeOutput() {
        final FeatureOutputCodec<?, ?> outputCodec = findOutputCodec(outputFilePath);
        final Class<?> outputClass = outputCodec.getFeatureType();
        if ( !evidenceClass.equals(outputClass) ) {
            throw new UserException("the input file contains " + evidenceClass.getSimpleName() +
                    " features, but the output file would be expected to contain " +
                    outputClass.getSimpleName() + " features");
        }
        final SAMSequenceDictionary dict;
        final List<String> samples;
        final Object headerObj = getDrivingFeaturesHeader();
        if ( headerObj instanceof FeaturesHeader ) {
            final FeaturesHeader header = (FeaturesHeader)headerObj;
            dict = header.getDictionary() == null ?
                    getBestAvailableSequenceDictionary() :
                    header.getDictionary();
            samples = header.getSampleNames() == null ? sampleNames : header.getSampleNames();
        } else {
            dict = getBestAvailableSequenceDictionary();
            samples = sampleNames;
        }
        outputSink = (FeatureSink<F>)outputCodec.makeSink(outputFilePath, dict, samples, compressionLevel);
    }

    private static FeatureOutputCodec<?, ?> findOutputCodec( final GATKPath outputFilePath ) {
        final String outputFileName = outputFilePath.toString();
        for ( final FeatureOutputCodec<?, ?> codec : outputCodecs ) {
            if ( codec.canDecode(outputFileName) ) {
                return codec;
            }
        }
        throw new UserException("no codec found for path " + outputFileName);
    }

    @Override
    public void apply(final F feature,
                      final ReadsContext readsContext,
                      final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        outputSink.write(feature);
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();
        outputSink.close();
        return null;
    }
}
