package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
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
            resultAsBEDFeatures = TargetCollections.fromBEDFeatureFile(exomeFile, bedFeatureCodec);
        } else {
            throw new UserException.BadInput(String.format("currently only BED formatted exome file are supported. '%s' does not seem to be a BED file", exomeFile.getAbsolutePath()));
        }
        logger.log(Level.INFO, String.format("Found %d targets to analyze.", resultAsBEDFeatures.targetCount()));

        return resultAsBEDFeatures;
    }
}
