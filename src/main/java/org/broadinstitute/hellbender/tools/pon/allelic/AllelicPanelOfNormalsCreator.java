package org.broadinstitute.hellbender.tools.pon.allelic;

import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Creates the panel of normals used for allele-bias correction.  See docs/CNVs/CNV-methods.pdf.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicPanelOfNormalsCreator {
    private static final Logger logger = LogManager.getLogger(AllelicPanelOfNormalsCreator.class);

    private final List<File> pulldownFiles;

    /**
     * Constuctor given pulldown files for all normals in the panel.
     * @param pulldownFiles pulldown files for all normals in the panel
     */
    public AllelicPanelOfNormalsCreator(final List<File> pulldownFiles) {
        Utils.nonNull(pulldownFiles, "List of pulldown files cannot be null.");
        ParamUtils.isPositive(pulldownFiles.size(), "List of pulldown files should contain at least one file.");
        pulldownFiles.stream().forEach(Utils::regularReadableUserFile);
        this.pulldownFiles = new ArrayList<>(pulldownFiles);
    }

    /**
     * Creates an {@link AllelicPanelOfNormals} given a site-frequency threshold;
     * sites appearing in strictly less than this fraction of samples will not be included in the panel of normals.
     * @param siteFrequencyThreshold    site-frequency threshold
     * @return                          an {@link AllelicPanelOfNormals} containing sites
     *                                  above the site-frequency threshold
     */
    public AllelicPanelOfNormals create(final double siteFrequencyThreshold) {
        logger.info("Creating allelic panel of normals...");
        final Map<SimpleInterval, MutableInt> numberOfSamplesMap = new HashMap<>(); //used to filter on frequency
        final Map<SimpleInterval, AllelicCount> totalCountsMap = new HashMap<>();   //store only the total counts (smaller memory footprint)
        int pulldownFileCounter = 1;
        final int totalNumberOfSamples = pulldownFiles.size();
        for (final File pulldownFile : pulldownFiles) {
            logger.info("Processing pulldown file " + pulldownFileCounter++ + "/" + totalNumberOfSamples + " (" + pulldownFile + ")...");
            final AllelicCountCollection pulldownCounts = new AllelicCountCollection(pulldownFile);
            for (final AllelicCount count : pulldownCounts.getCounts()) {
                //update the sum of ref and alt counts at each site
                final SimpleInterval site = count.getInterval();
                final AllelicCount currentCountAtSite = totalCountsMap.getOrDefault(site, new AllelicCount(site, 0, 0));
                final AllelicCount updatedCountAtSite = new AllelicCount(
                        site,
                        currentCountAtSite.getRefReadCount() + count.getRefReadCount(),
                        currentCountAtSite.getAltReadCount() + count.getAltReadCount());
                totalCountsMap.put(site, updatedCountAtSite);
                //update the number of samples seen possessing each site
                final MutableInt numberOfSamplesAtSite = numberOfSamplesMap.get(site);
                if (numberOfSamplesAtSite == null) {
                    numberOfSamplesMap.put(site, new MutableInt(1));
                } else {
                    numberOfSamplesAtSite.increment();
                }
            }
        }
        logger.info("Total number of unique sites present in samples: " + totalCountsMap.size());
        //filter out sites that appear at a frequency strictly less than the provided threshold
        final AllelicCountCollection totalCounts = new AllelicCountCollection();
        numberOfSamplesMap.entrySet().stream()
                .filter(e -> e.getValue().doubleValue() / totalNumberOfSamples >= siteFrequencyThreshold)
                .map(e -> totalCountsMap.get(e.getKey()))
                .forEach(totalCounts::add);
        logger.info(String.format("Number of unique sites present in samples above site frequency = %4.3f: %d", siteFrequencyThreshold, totalCounts.getCounts().size()));
        return new AllelicPanelOfNormals(totalCounts);
    }
}
