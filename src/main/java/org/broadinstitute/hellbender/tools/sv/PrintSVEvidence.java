package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.BafEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.DiscordantPairEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.SplitReadEvidenceCodec;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStreamFactory;

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
 *         Coordinate-sorted evidence file, indexed if block compressed
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
public final class PrintSVEvidence extends FeatureWalker<Feature> {

    public static final String EVIDENCE_FILE_NAME = "evidence-file";
    public static final String COMPRESSION_LEVEL_NAME = "compression-level";

    @Argument(
            doc = "Input file URI with extension '"
                    + SplitReadEvidenceCodec.FORMAT_SUFFIX + "', '"
                    + DiscordantPairEvidenceCodec.FORMAT_SUFFIX + "', '"
                    + BafEvidenceCodec.FORMAT_SUFFIX + "', or '"
                    + DepthEvidenceCodec.FORMAT_SUFFIX + "' (may be gzipped).",
            fullName = EVIDENCE_FILE_NAME
    )
    private GATKPath inputFilePath;

    @Argument(
            doc = "Output file with an evidence extension matching the input. Will be indexed if it has a " +
                    "block-compressed extension (e.g. '.gz').",
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

    private FeatureOutputStream outputStream;

    @Override
    protected boolean isAcceptableFeatureType(final Class<? extends Feature> featureType) {
        return featureType.equals(BafEvidence.class) || featureType.equals(DepthEvidence.class)
                || featureType.equals(DiscordantPairEvidence.class) || featureType.equals(SplitReadEvidence.class);
    }

    @Override
    public GATKPath getDrivingFeaturePath() {
        return inputFilePath;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        validateInputs();
        initializeOutput();
        writeHeader();
    }

    private void validateInputs() {
        final Class inputClass = FeatureManager.getCodecForFile(inputFilePath.toPath()).getClass();
        final Class outputClass = FeatureManager.getCodecForFile(outputFilePath.toPath()).getClass();
        Utils.validate(inputClass == outputClass, "Input and output file types do not match");
    }

    private void initializeOutput() {
        outputStream = new FeatureOutputStreamFactory().create(outputFilePath, SVIOUtils::encodeSVEvidenceFeature,
                getBestAvailableSequenceDictionary(), compressionLevel);
    }

    private void writeHeader() {
        final Object header = getDrivingFeaturesHeader();
        if (header != null) {
            if (header instanceof String) {
                outputStream.writeHeader((String) header);
            } else {
                throw new IllegalArgumentException("Expected header object of type " + String.class.getSimpleName());
            }
        }
    }

    @Override
    public void apply(final Feature feature,
                      final ReadsContext readsContext,
                      final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        // All evidence data types implement an encoding with toString()
        outputStream.add(feature);
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();
        if (outputStream != null) {
            outputStream.close();
        }
        return null;
    }
}
