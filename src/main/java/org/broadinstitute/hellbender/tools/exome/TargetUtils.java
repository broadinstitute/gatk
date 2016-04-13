package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class TargetUtils {

    protected final static Logger logger = LogManager.getLogger(TargetUtils.class);

    private TargetUtils() {}

    /**
     * Create a TargetCollection from a BED file.
     *
     * @param exomeFile BED file of targets.
     * @return never {@code null}
     */
    public static TargetCollection<? extends BEDFeature> readTargetFile(final File exomeFile) {
        Utils.regularReadableUserFile(exomeFile);
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(exomeFile);
        logger.log(Level.INFO, String.format("Reading target intervals from exome file '%s' ...", exomeFile.getAbsolutePath()));
        final Class<? extends Feature> featureType = codec.getFeatureType();
        final TargetCollection<? extends BEDFeature> resultAsBEDFeatures;
        if (BEDFeature.class.isAssignableFrom(featureType)) {
            @SuppressWarnings("unchecked")
            final FeatureCodec<? extends BEDFeature, ?> bedFeatureCodec = (FeatureCodec<? extends BEDFeature, ?>) codec;
            resultAsBEDFeatures = TargetCollectionUtils.fromBEDFeatureFile(exomeFile, bedFeatureCodec);
        } else {
            throw new UserException.BadInput(String.format("currently only BED formatted exome file are supported. '%s' does not seem to be a BED file", exomeFile.getAbsolutePath()));
        }
        logger.log(Level.INFO, String.format("Found %d targets to analyze.", resultAsBEDFeatures.targetCount()));

        return resultAsBEDFeatures;
    }

    /**
     * write a list of targets to BED file
     *
     */
    public static void writeTargetsAsBed(final File outFile, final List<Target> targets) {

        Utils.nonNull(outFile, "Output file cannot be null.");
        Utils.nonNull(targets, "Targets cannot be null.");

        final boolean areTargetIntervalsAllPopulated = targets.stream().allMatch(t -> t.getInterval() != null);
        if (!areTargetIntervalsAllPopulated) {
            throw new UserException("Cannot write target file with any null intervals.");
        }

        try (final FileWriter fw = new FileWriter(outFile)) {
            fw.write("##CONTIG  START\t\tEND\tNAME\n");
            for (Target t: targets) {
                final List<String> fields = Arrays.asList(t.getContig(), String.valueOf(t.getInterval().getGA4GHStart()),
                        String.valueOf(t.getInterval().getGA4GHEnd()), t.getName());

                final String lineString = StringUtils.join(fields, '\t');
                fw.write(lineString + "\n");
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e);
        }
    }

    /**
     * Creates a string for a locatable that can be used when creating dummy target names
     * @param locatable The genome region to create a unique dummy target name. Never {@code null}
     * @return never {@code null}
     */
    public static String createDummyTargetName(final Locatable locatable){
        Utils.nonNull(locatable, "Output file cannot be null.");
        return "target_" + locatable.getContig() + "_" + String.valueOf(locatable.getStart()) + "_" + String.valueOf(locatable.getEnd());
    }
}
