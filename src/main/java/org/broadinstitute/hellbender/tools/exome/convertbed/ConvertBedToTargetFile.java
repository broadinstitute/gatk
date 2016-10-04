package org.broadinstitute.hellbender.tools.exome.convertbed;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
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

@CommandLineProgramProperties(
        summary = "Converts a target bed file to the target file format. Empty files will probably fail. Drops bed file columns other than the first four.",
        oneLineSummary = "Convert bed to target file",
        programGroup = CopyNumberProgramGroup.class
)
public class ConvertBedToTargetFile extends CommandLineProgram {

    public static final String BED_INPUT_LONG_NAME = StandardArgumentDefinitions.INPUT_LONG_NAME;
    public static final String BED_INPUT_SHORT_NAME = StandardArgumentDefinitions.INPUT_SHORT_NAME;

    public static final String TARGET_OUTPUT_LONG_NAME = StandardArgumentDefinitions.OUTPUT_LONG_NAME;
    public static final String TARGET_OUTPUT_SHORT_NAME = StandardArgumentDefinitions.OUTPUT_SHORT_NAME;

    @Argument(
            doc = "target bed file",
            shortName = BED_INPUT_SHORT_NAME,
            fullName = BED_INPUT_LONG_NAME,
            optional = false
    )
    protected File inputBedFile;

    @Argument(
            doc = "output target file",
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
            final List<Target> targets = StreamSupport.stream(source.spliterator(), false).map(ConvertBedToTargetFile::createTargetFromBEDFeature)
                    .collect(Collectors.toList());
            TargetWriter.writeTargetsToFile(outFile, targets);
        } else {
            throw new UserException.BadInput(String.format("currently only BED formatted exome file are supported. '%s' does not seem to be a BED file", inputBedFile.getAbsolutePath()));
        }
        return "SUCCESS";
    }

    protected static Target createTargetFromBEDFeature(final BEDFeature bedFeature) {
        return new Target(bedFeature.getName(), new SimpleInterval(bedFeature.getContig(), bedFeature.getStart(), bedFeature.getEnd()));
    }
}
