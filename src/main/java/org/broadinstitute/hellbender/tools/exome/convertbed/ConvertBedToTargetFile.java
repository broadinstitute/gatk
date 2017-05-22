package org.broadinstitute.hellbender.tools.exome.convertbed;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.bed.BEDFeature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetWriter;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;
/**
 * Converts 0-based BED format intervals to 1-based GATK-style intervals.
 *
 * <p>
 *     Conversion retains data from the first four columns of a BED file.
 *     The tool adds a line at top to label the four columns.
 *     These columns are the contig, start, stop and name.
 *     For example, the following BED intervals:
 *     <pre>
 *         chr1    11868   12227
 *         chr1    12612   12721
 *     </pre>
 *
 *     become the following GATK-style intervals:
 *     <pre>
 *         contig  start   stop    name
 *         chr1    11869   12227
 *         chr1    12613   12721
 *     </pre>
 * </p>
 *
 * <h3>Example</h3>
 *
 * <pre>
 * java -Xmx4g -jar $gatk_jar ConvertBedToTargetFile
 *   --input targets.bed
 *   --output targets.tsv
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Converts a BED file to the target file format. Drops BED file columns other than the first four.",
        oneLineSummary = "Convert BED to target file",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public class ConvertBedToTargetFile extends CommandLineProgram {

    public static final String BED_INPUT_LONG_NAME = StandardArgumentDefinitions.INPUT_LONG_NAME;
    public static final String BED_INPUT_SHORT_NAME = StandardArgumentDefinitions.INPUT_SHORT_NAME;

    public static final String TARGET_OUTPUT_LONG_NAME = StandardArgumentDefinitions.OUTPUT_LONG_NAME;
    public static final String TARGET_OUTPUT_SHORT_NAME = StandardArgumentDefinitions.OUTPUT_SHORT_NAME;

    @Argument(
            doc = "Input BED file.",
            shortName = BED_INPUT_SHORT_NAME,
            fullName = BED_INPUT_LONG_NAME,
            optional = false
    )
    protected File inputBedFile;

    @Argument(
            doc = "Output target file.",
            shortName = TARGET_OUTPUT_SHORT_NAME,
            fullName = TARGET_OUTPUT_LONG_NAME,
            optional = false
    )
    protected File outFile;

    @Override
    protected Object doWork() {
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(inputBedFile);
        final Class<? extends Feature> featureType = codec.getFeatureType();
        if (BEDFeature.class.isAssignableFrom(featureType)) {
            final FeatureDataSource<? extends BEDFeature> source = new FeatureDataSource<>(inputBedFile);
            try {
                final List<Target> targets = StreamSupport.stream(source.spliterator(), false).map(ConvertBedToTargetFile::createTargetFromBEDFeature)
                        .collect(Collectors.toList());
                TargetWriter.writeTargetsToFile(outFile, targets);
            } catch (final TribbleException e) {
                throw new UserException.BadInput(String.format("'%s' has a .bed extension but does not seem to be a valid BED file.", inputBedFile.getAbsolutePath()));
            }
        } else {
            throw new UserException.BadInput(String.format("'%s' does not seem to be a BED file.", inputBedFile.getAbsolutePath()));
        }
        return "SUCCESS";
    }

    protected static Target createTargetFromBEDFeature(final BEDFeature bedFeature) {
        return new Target(bedFeature.getName(), new SimpleInterval(bedFeature.getContig(), bedFeature.getStart(), bedFeature.getEnd()));
    }
}
