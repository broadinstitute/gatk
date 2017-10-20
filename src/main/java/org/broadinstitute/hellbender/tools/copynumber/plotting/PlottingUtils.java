package org.broadinstitute.hellbender.tools.copynumber.plotting;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SampleLocatableCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class PlottingUtils {
    static final String CNV_PLOTTING_R_LIBRARY = "CNVPlottingLibrary.R";

    static final String MINIMUM_CONTIG_LENGTH_LONG_NAME = "minimumContigLength";
    static final String MINIMUM_CONTIG_LENGTH_SHORT_NAME = "minContigLength";

    static final String CONTIG_DELIMITER = "CONTIG_DELIMITER";  //used to delimit contig names and lengths passed to the R script
    static final int DEFAULT_MINIMUM_CONTIG_LENGTH = 1000000;   //can be used to filter out mitochondrial contigs, unlocalized contigs, etc.

    private PlottingUtils() {}

    static Map<String, Integer> getContigLengthMap(final File sequenceDictionaryFile,
                                                   final int minContigLength,
                                                   final Logger logger) {
        final SAMSequenceDictionary sequenceDictionary = ReferenceUtils.loadFastaDictionary(sequenceDictionaryFile);
        Utils.validateArg(sequenceDictionary.getSequences().stream().map(SAMSequenceRecord::getSequenceName).noneMatch(n -> n.contains(CONTIG_DELIMITER)),
                String.format("Contig names cannot contain \"%s\".", CONTIG_DELIMITER));
        final Map<String, Integer> contigLengthMap = sequenceDictionary.getSequences().stream()
                .filter(s -> s.getSequenceLength() >= minContigLength)
                .collect(Collectors.toMap(SAMSequenceRecord::getSequenceName, SAMSequenceRecord::getSequenceLength,
                        (c, l) -> {
                            throw new IllegalArgumentException(String.format("Duplicate contig in sequence dictionary: %s", c));
                        },
                        LinkedHashMap::new));
        Utils.validateArg(contigLengthMap.size() > 0,
                "There must be at least one contig above the threshold length in the sequence dictionary.");
        logger.info("Contigs above length threshold: " + contigLengthMap.toString());
        return contigLengthMap;
    }

    //validate contig names and lengths
    static <T extends Locatable> void validateContigs(final Map<String, Integer> contigLengthMap,
                                                      final SampleLocatableCollection<T> locatableCollection,
                                                      final File file,
                                                      final Logger logger) {
        if (locatableCollection == null) {
            return;
        }
        final Set<String> contigNames = contigLengthMap.keySet();
        final Set<String> fileContigNames = locatableCollection.getRecords().stream().map(T::getContig).collect(Collectors.toSet());
        if (!contigNames.containsAll(fileContigNames)) {
            logger.warn(String.format("Contigs present in the file %s are missing from the sequence dictionary and will not be plotted.", file));
        }
        final Map<String, Integer> fileContigMaxPositionMap = locatableCollection.getIntervals().stream().filter(i -> contigNames.contains(i.getContig()))
                .collect(Collectors.toMap(SimpleInterval::getContig, SimpleInterval::getEnd, Integer::max));
        fileContigMaxPositionMap.keySet().forEach(c -> Utils.validateArg(fileContigMaxPositionMap.get(c) <= contigLengthMap.get(c),
                String.format("Position present in the file %s exceeds contig length in the sequence dictionary.", file)));
    }

    static String addTrailingSlashIfNecessary(final String outputDir) {
        return outputDir.endsWith(File.separator) ? outputDir : outputDir + File.separator;
    }
}
