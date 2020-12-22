package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.ParsingUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

public final class IndexUtils {
    private IndexUtils(){}

    private static final Logger logger = LogManager.getLogger(IndexUtils.class);

    /**
     * Load a Tribble .idx index from disk, checking for out of date indexes and old versions
     * @return an Index, or null if we're unable to load
     */
    public static Index loadTribbleIndex(final Path featureFile) {
        Utils.nonNull(featureFile);
        final String indexFile = Tribble.indexFile(featureFile.toString());
        final Path indexPath = IOUtils.getPath(indexFile);
        if (! Files.isReadable(indexPath)) {
            return null;
        }
        logger.debug("Loading Tribble index from disk for file " + featureFile);
        try {
            final Index index = IndexFactory.loadIndex(indexFile);
            checkIndexVersionAndModificationTime(featureFile, indexPath, index);
            return index;
        } catch (final RuntimeException e){
            return null;
        }
    }

    /**
     * Try to load the tabix index from disk, checking for out of date indexes and old versions
     * @return an Index, or null if we're unable to load
     */
    public static Index loadTabixIndex(final Path featureFile) {
        Utils.nonNull(featureFile);
        try {
            final String path = featureFile.toUri().toString();
            final boolean isTabix = new AbstractFeatureReader.ComponentMethods().isTabix(path, null);
            if (! isTabix){
                return null;
            }
            final String indexPath = ParsingUtils.appendToPath(path, FileExtensions.TABIX_INDEX);
            logger.debug("Loading tabix index from disk for file " + featureFile);
            final Index index = IndexFactory.loadIndex(indexPath);
            final Path indexFile = IOUtils.getPath(indexPath);
            checkIndexVersionAndModificationTime(featureFile, indexFile, index);
            return index;
        } catch (final IOException | RuntimeException e) {
            return null;
        }
    }

    /**
     * @throws UserException the index is not the current version.
     * Prints a warning if the index is not up-to date (ie older than the feature file)
     */
    public static void checkIndexVersionAndModificationTime(final Path featureFile, final Path indexFile, final Index index) {
        Utils.nonNull(featureFile, "feature-file");
        Utils.nonNull(indexFile, "indexFile");
        Utils.nonNull(index, "index");
        if (! index.isCurrentVersion()) {
            // we've loaded an old version of the index, we want to remove it
            throw new UserException("Index file " + indexFile + " is out of date (old version). Use " + IndexFeatureFile.class.getSimpleName() +  " to make an index.");
        } else try {
            if (Files.getLastModifiedTime(indexFile).compareTo(Files.getLastModifiedTime(featureFile)) < 0) {
                logger.warn("Index file " + indexFile + " is out of date (index older than input file). Use " + IndexFeatureFile.class.getSimpleName() +  " to make a new index.");
            }
        } catch (IOException e) {
            logger.warn("Could not validate index modification time.  Assuming index is up-to-date.");
        }
    }

    /**
     * get the sequence dictionary contig list that is in the index or null if there is no index or no contigs
     * Note: the dictionary returned will not have the contig lengths filled in {@link SAMSequenceRecord#UNKNOWN_SEQUENCE_LENGTH} is used.
     * Note: this method is specifically designed for getting sequence dictionaries from indices on Feature files (tribble or tabix indices)
     * @return a SAMSequenceDictionary or null if the index cannot be loaded or there are no contigs in the index
     */
    public static SAMSequenceDictionary createSequenceDictionaryFromFeatureIndex(final Path featureFile) {
        Utils.nonNull(featureFile);
        logger.warn(
                String.format(
                    "Feature file \"%s\" appears to contain no sequence dictionary. " +
                    "Attempting to retrieve a sequence dictionary from the associated index file",
                    featureFile.toAbsolutePath().toUri().toString())
        );
        final Index index = loadIndex(featureFile);
        return index == null ? null : getSamSequenceDictionaryFromIndex(index);
    }

    /**
     * Determine if <code>sequenceDictionary</code> has the form of a dictionary created from an index file by
     * {@link #createSequenceDictionaryFromFeatureIndex}).
     * @param sequenceDictionary
     * @return
     */
    public static boolean isSequenceDictionaryFromIndex(final SAMSequenceDictionary sequenceDictionary) {
        Utils.nonNull(sequenceDictionary);
        return sequenceDictionary.getSequences()
                .stream()
                .allMatch(seqRec -> seqRec.getSequenceLength() == SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH);
    }

    /**
     * Loads and returns the index of the feature file or null if there is no index.
     * First it tries to load the tribble index and then to load the tabix index.
     */
    public static Index loadIndex(final Path featureFile){
        Utils.nonNull(featureFile);
        final Index tribbleIndex = loadTribbleIndex(featureFile);
        if (tribbleIndex != null) {
            return tribbleIndex;
        }

        final Index tabixIndex = loadTabixIndex(featureFile);
        if (tabixIndex != null){
            return tabixIndex;
        }
        return null;
    }

    private static SAMSequenceDictionary getSamSequenceDictionaryFromIndex(final Index index) {
        final List<String> seqNames = index.getSequenceNames();
        if (seqNames == null || seqNames.isEmpty()) {
            return null;
        }
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        //use UNKNOWN_SEQUENCE_LENGTH to indicate contigs that will not be compared by length (see SequenceDictionaryUtils.sequenceRecordsAreEquivalent)
        seqNames.forEach(seqName -> dict.addSequence(new SAMSequenceRecord(seqName, SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH)));
        return dict;
    }

}
