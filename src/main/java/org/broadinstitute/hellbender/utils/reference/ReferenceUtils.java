package org.broadinstitute.hellbender.utils.reference;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.repackaged.com.google.common.io.CountingOutputStream;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.engine.datasources.ReferenceAPISource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceFileSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceHadoopSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.seqdoop.hadoop_bam.FastaInputFormat;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.Writer;
import java.nio.file.Path;

/**
 * A collection of static methods for dealing with references.
 */
public final class ReferenceUtils {

    // Private so that no one will instantiate this class.
    private ReferenceUtils() {
    }

    /**
     * Given a fasta filename, return the name of the corresponding index file.
     * (This also works if the file is in gs://)
     */
    public static String getFastaIndexFileName(String fastaFilename) {
        return fastaFilename + ".fai";
    }

    /**
     * Given a fasta filename, return the name of the corresponding dictionary file.
     * (This also works if the file is in gs://)
     */
    public static String getFastaDictionaryFileName(String fastaFilename) {
        int lastDot = fastaFilename.lastIndexOf('.');
        return fastaFilename.substring(0, lastDot) + ".dict";
    }

    /**
     * Given a fasta dictionary file, returns its sequence dictionary
     *
     * @param fastaDictionaryFile fasta dictionary file
     * @return the SAMSequenceDictionary from fastaDictionaryFile
     */
    public static SAMSequenceDictionary loadFastaDictionary(final File fastaDictionaryFile) {
        try (final FileInputStream fastaDictionaryStream = new FileInputStream(fastaDictionaryFile)) {
            return loadFastaDictionary(fastaDictionaryStream);
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile("Error loading fasta dictionary file " + fastaDictionaryFile, e);
        } catch (UserException.MalformedFile e) {
            throw new UserException.MalformedFile(
                    "Could not read sequence dictionary from given fasta file " +
                            fastaDictionaryFile
            );
        }
    }

    /**
     * Given an InputStream connected to a fasta dictionary, returns its sequence dictionary
     * <p>
     * Note: does not close the InputStream it's passed
     *
     * @param fastaDictionaryStream InputStream connected to a fasta dictionary
     * @return the SAMSequenceDictionary from the fastaDictionaryStream
     */
    public static SAMSequenceDictionary loadFastaDictionary(final InputStream fastaDictionaryStream) {
        // Don't close the reader when we're done, since we don't want to close the client's InputStream for them
        final BufferedLineReader reader = new BufferedLineReader(fastaDictionaryStream);

        final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
        final SAMFileHeader header = codec.decode(reader, fastaDictionaryStream.toString());

        // Make sure we have a valid sequence dictionary before continuing:
        if (header.getSequenceDictionary() == null || header.getSequenceDictionary().isEmpty()) {
            throw new UserException.MalformedFile(
                    "Could not read sequence dictionary from given fasta stream " +
                            fastaDictionaryStream
            );
        }

        return header.getSequenceDictionary();
    }

    /**
     * Dump a reference into a path.
     * <p>
     * Notice that this will take long with large reference input sources.
     * </p>
     *
     * @param input   the source reference.
     * @param whereTo where to dump the reference.
     * @return a new reference source pointing to the dump.
     */
    public static ReferenceSource dumpReferenceSource(final ReferenceSource input, final Path whereTo, final int basesPerLine)
            throws IOException {
        // get approx ~ 2^16 (64Kb) bases from the source each time.
        // adjusted up to the closest multiple of basesPerLine.
        final int requestSize = (int) Math.ceil((1 << 16) / (double) basesPerLine) * basesPerLine;
        Utils.nonNull(input, "the input reference source cannot be negative");
        Utils.nonNull(whereTo, "the output path cannot be negative");
        ParamUtils.isPositive(basesPerLine, "bases per line must be greater han 0");
        final Path fastaPath = whereTo;
        final Path indexPath = fastaPath.getParent().resolve(whereTo.getFileName() + ".fai");
        final CountingOutputStream fastaOutStream = new CountingOutputStream(BucketUtils.createFile(fastaPath.toString()));
        final OutputStream indexOutStream = BucketUtils.createFile(indexPath.toString());
        final Writer fastaWriter = new PrintWriter(fastaOutStream);
        final Writer indexWriter = new PrintWriter(indexOutStream);
        final SAMSequenceDictionary dictionary = input.getReferenceSequenceDictionary();
        for (final SAMSequenceRecord sequence : dictionary.getSequences()) {
            final String sequenceName = sequence.getSequenceName();
            final int sequenceLength = sequence.getSequenceLength();
            fastaWriter.append('>').append(sequenceName).append(System.lineSeparator()).flush();
            indexWriter.append(sequenceName).append('\t')
                    .append(String.valueOf(sequenceLength)).append('\t')
                    .append(String.valueOf(fastaOutStream.getCount())).append('\t')
                    .append(String.valueOf(basesPerLine)).append(System.lineSeparator());

            int nextIdx = 0;
            while (nextIdx < sequenceLength) {
                final SimpleInterval request = new SimpleInterval(sequenceName, nextIdx + 1, Math.min(nextIdx + requestSize, sequenceLength));
                final byte[] bases = input.getReferenceBases(request).getBases();
                final int stop = bases.length - basesPerLine;
                int offset = 0;
                for (offset = 0; offset <= stop; offset += basesPerLine) {
                    fastaOutStream.write(bases, offset, basesPerLine);
                    fastaWriter.append(System.lineSeparator()).flush();
                }
                // a partial last line:
                if (offset < bases.length) {
                    fastaOutStream.write(bases, offset, bases.length - offset);
                    fastaWriter.append(System.lineSeparator()).flush();
                }
                nextIdx += bases.length;
            }
        }
        fastaWriter.close();
        indexWriter.close();
        if (BucketUtils.isCloudStorageUrl(whereTo)) {
            return new ReferenceAPISource(null, whereTo.toString());
        } else if (BucketUtils.isHadoopUrl(whereTo.toString())) {
            return new ReferenceHadoopSource(whereTo.toString());
        } else if (BucketUtils.isFileUrl(whereTo.toString())) {
            return new ReferenceFileSource(whereTo.toString());
        } else {
            throw new GATKException.ShouldNeverReachHereException("this should never have been reached.");
        }
    }
}
