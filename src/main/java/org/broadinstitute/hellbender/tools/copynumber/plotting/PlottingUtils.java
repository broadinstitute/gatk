package org.broadinstitute.hellbender.tools.copynumber.plotting;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AbstractSampleLocatableCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class PlottingUtils {
    static final String CNV_PLOTTING_R_LIBRARY = "CNVPlottingLibrary.R";

    static final String MINIMUM_CONTIG_LENGTH_LONG_NAME = "minimum-contig-length";

    static final String CONTIG_DELIMITER = "CONTIG_DELIMITER";  //used to delimit contig names and lengths passed to the R script
    static final int DEFAULT_MINIMUM_CONTIG_LENGTH = 1000000;   //can be used to filter out mitochondrial contigs, unlocalized contigs, etc.

    static final String SEQUENCE_DICTIONARY_DOC_STRING = "File containing a sequence dictionary, which specifies the contigs to be plotted and their relative lengths. " +
            "The sequence dictionary must be a subset of those contained in other input files. " +
            "Contigs will be plotted in the order given. " +
            "Contig names should not include the string \"" + PlottingUtils.CONTIG_DELIMITER + "\". " +
            "The tool only considers contigs in the given dictionary for plotting, and " +
            "data for contigs absent in the dictionary generate only a warning. In other words, you may " +
            "modify a reference dictionary for use with this tool to include only contigs for which plotting is desired, " +
            "and sort the contigs to the order in which the plots should display the contigs.";

    static final String MINIMUM_CONTIG_LENGTH_DOC_STRING = "Threshold length (in bp) for contigs to be plotted. " +
            "Contigs with lengths less than this threshold will not be plotted. " +
            "This can be used to filter out mitochondrial contigs, unlocalized contigs, etc.";

    private PlottingUtils() {}

    static Map<String, Integer> getContigLengthMap(final SAMSequenceDictionary sequenceDictionary,
                                                   final int minContigLength,
                                                   final Logger logger) {
        Utils.nonNull(sequenceDictionary);
        ParamUtils.isPositiveOrZero(minContigLength, "Minimum contig length must be non-negative.");
        Utils.nonNull(logger);
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
                                                      final AbstractSampleLocatableCollection<T> locatableCollection,
                                                      final File file,
                                                      final Logger logger) {
        Utils.nonNull(contigLengthMap);
        Utils.nonNull(logger);
        if (locatableCollection == null) {
            Utils.validateArg(file == null, "File can only be null if collection is also null.");
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

    //validate sequence dictionary contains another as a subset, using only sequence names and lengths
    static void validateSequenceDictionarySubset(final SAMSequenceDictionary sequenceDictionary,
                                                 final SAMSequenceDictionary sequenceDictionarySubset) {
        Utils.nonNull(sequenceDictionary);
        Utils.nonNull(sequenceDictionarySubset);
        final Set<SAMSequenceRecord> sequences = sequenceDictionary.getSequences().stream()
                .map(s -> new SAMSequenceRecord(s.getMd5(), s.getSequenceLength()))
                .collect(Collectors.toCollection(LinkedHashSet::new));
        final Set<SAMSequenceRecord> sequencesSubset = sequenceDictionary.getSequences().stream()
                .map(s -> new SAMSequenceRecord(s.getMd5(), s.getSequenceLength()))
                .collect(Collectors.toCollection(LinkedHashSet::new));
        Utils.validateArg(sequences.containsAll(sequencesSubset),
                String.format("Sequence dictionary (%s) must be a subset of those contained in other input files (%s).",
                        sequencesSubset, sequences));
    }

    static String addTrailingSlashIfNecessary(final String outputDir) {
        Utils.nonEmpty(outputDir);
        return outputDir.endsWith(File.separator) ? outputDir : outputDir + File.separator;
    }
}
