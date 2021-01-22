package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.codecs.BafEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.DiscordantPairEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.SplitReadEvidenceCodec;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;

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

    private FeatureOutputStream<DiscordantPairEvidence> peStream;
    private FeatureOutputStream<SplitReadEvidence> srStream;
    private FeatureOutputStream<BafEvidence> bafStream;
    private FeatureOutputStream<DepthEvidence> rdStream;
    private FeatureCodec<? extends Feature, ?> featureCodec;
    private Class<? extends Feature> evidenceClass;

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
        featureCodec = FeatureManager.getCodecForFile(inputFilePath.toPath());
        evidenceClass = featureCodec.getFeatureType();
        initializeOutput();
        writeHeader();
    }

    private void initializeOutput() {
        if (evidenceClass.equals(DiscordantPairEvidence.class)) {
            peStream = new FeatureOutputStream<>(outputFilePath, featureCodec, DiscordantPairEvidenceCodec::encode,
                    getBestAvailableSequenceDictionary(), compressionLevel);
        } else if (evidenceClass.equals(SplitReadEvidence.class)) {
            srStream = new FeatureOutputStream<>(outputFilePath, featureCodec, SplitReadEvidenceCodec::encode,
                    getBestAvailableSequenceDictionary(), compressionLevel);
        } else if (evidenceClass.equals(BafEvidence.class)) {
            bafStream = new FeatureOutputStream<>(outputFilePath, featureCodec, BafEvidenceCodec::encode,
                    getBestAvailableSequenceDictionary(), compressionLevel);
        } else if (evidenceClass.equals(DepthEvidence.class)) {
            rdStream = new FeatureOutputStream<>(outputFilePath, featureCodec, DepthEvidenceCodec::encode,
                    getBestAvailableSequenceDictionary(), compressionLevel);
        } else {
            throw new UserException.BadInput("Unsupported evidence type: " + evidenceClass.getSimpleName());
        }
    }

    private void writeHeader() {
        final Object header = getDrivingFeaturesHeader();
        if (header != null) {
            if (header instanceof String) {
                if (peStream != null) {
                    peStream.writeHeader((String) header);
                } else if (srStream != null) {
                    srStream.writeHeader((String) header);
                } else if (bafStream != null) {
                    bafStream.writeHeader((String) header);
                } else {
                    rdStream.writeHeader((String) header);
                }
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
        if (peStream != null) {
            peStream.add((DiscordantPairEvidence) feature);
        } else if (srStream != null) {
            srStream.add((SplitReadEvidence) feature);
        } else if (bafStream != null) {
            bafStream.add((BafEvidence) feature);
        } else {
            rdStream.add((DepthEvidence) feature);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();
        if (peStream != null) {
            peStream.close();
        } else if (srStream != null) {
            srStream.close();
        } else if (bafStream != null) {
            bafStream.close();
        } else {
            rdStream.close();
        }
        return null;
    }
}
