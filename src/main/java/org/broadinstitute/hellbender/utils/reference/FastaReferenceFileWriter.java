package org.broadinstitute.hellbender.utils.reference;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.io.CountingOutputStream;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.codehaus.plexus.util.StringUtils;

import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Utility class to compose FASTA formatted reference files.
 */
public class FastaReferenceFileWriter implements AutoCloseable {

    public static final int DEFAULT_BASES_PER_LINE = 60;

    public static final char HEADER_START_CHAR = '>';
    public static final char HEADER_NAME_AND_DESCRIPTION_SEPARATOR = ' ';

    private static final byte[] LINE_SEPARATOR = System.lineSeparator().getBytes();

    // we use this counter for basic byte array outputs.
    private final CountingOutputStream fastaOutput;
    // we use the writer for String / number and other formatted outputs.
    private final Writer fastaWriter;
    private final Writer faiWriter;
    private final OutputStream dictOutput;
    private final int defaultBasePerLine;
    private final Map<String, Long> sequenceNamesAndSizes;
    private int currentBasePerLine;
    private int lineBasesCount;
    private long currentBasesCount;
    private long currentSequenceOffset;
    private String currentSequenceName;
    private boolean closed;

    /**
     * Creates a reference FASTA file writer.
     * <p>
     *     The default bases-per-line is set to {@link #DEFAULT_BASES_PER_LINE}.
     * </p>
     * <p>
     *     Names for the fasta index and dictionary are constructed from the FASTA output file using common practices
     *     as resolved by {@link ReferenceSequenceFileFactory#getFastaIndexFileName(Path)}
     *     and {@link ReferenceSequenceFileFactory#getDefaultDictionaryForReferenceSequence(Path)}
     *     respectively.
     * </p>
     *
     * @param fastaFile the output fasta file path.
     * @param makeFaiOutput whether an index must be generated.
     * @param makeDictOutput whether a dictionary must be generated.
     * @throws IllegalArgumentException if {@code fastaFile} is {@code null}.
     * @throws IOException if such exception is thrown when accessing the output path resources.
     */
    public FastaReferenceFileWriter(final Path fastaFile, final boolean makeFaiOutput, final boolean makeDictOutput)
        throws IOException
    {
        this(fastaFile, DEFAULT_BASES_PER_LINE, makeFaiOutput, makeDictOutput);
    }

    /**
     * Creates a reference FASTA file writer.
     * <p>
     *     Names for the fasta index and dictionary are constructed from the FASTA output file using common practices
     *     as resolved by {@link ReferenceSequenceFileFactory#getFastaIndexFileName(Path)}
     *     and {@link ReferenceSequenceFileFactory#getDefaultDictionaryForReferenceSequence(Path)}
     *     respectively.
     * </p>
     *
     * @param fastaFile the output fasta file path.
     * @param basesPerLine default bases per line.
     * @param makeFaiOutput whether an index must be generated.
     * @param makeDictOutput whether a dictionary must be generated.
     * @throws IllegalArgumentException if {@code fastaFile} is {@code null} or {@code basesPerLine} is 0 or negative.
     * @throws IOException if such exception is thrown when accessing the output path resources.
     */
    public FastaReferenceFileWriter(final Path fastaFile, final int basesPerLine, final boolean makeFaiOutput,
                                    final boolean makeDictOutput)
        throws IOException
    {
        this(Utils.nonNull(fastaFile, "the input fasta-file cannot be null"),
                basesPerLine,
                defaultFaiFile(makeFaiOutput, fastaFile),
                defaultDictFile(makeDictOutput, fastaFile));
    }

    /**
     * Creates a reference FASTA file writer.
     * <p>
     *     The default bases-per-line is set to {@link #DEFAULT_BASES_PER_LINE}.
     * </p>
     * <p>
     *     You can specify a specific path for the index and dictionary file. If either set to {@code null} such
     *     a file won't be generated.
     * </p>
     *
     * @param fastaFile the output fasta file path.
     * @param faiFile the path of the index file, if requested, {@code null} if none should be generated.
     * @param dictFile the path of the dictFile, if requested, {@code null} if nono should be generated.
     * @throws IllegalArgumentException if {@code fastaFile} is {@code null}.
     * @throws IOException if such exception is thrown when accessing the output path resources.
     */
    public FastaReferenceFileWriter(final Path fastaFile, final Path faiFile, final Path dictFile)
        throws IOException
    {
        this(fastaFile, DEFAULT_BASES_PER_LINE, faiFile, dictFile);
    }

    /**
     * Creates a reference FASTA file writer.
     * <p>
     *     You can specify a specific path for the index and dictionary file. If either set to {@code null} such
     *     a file won't be generated.
     * </p>
     *
     * @param fastaFile the output fasta file path.
     * @param faiFile the path of the index file, if requested, {@code null} if none should be generated.
     * @param dictFile the path of the dictFile, if requested, {@code null} if nono should be generated.
     * @throws IllegalArgumentException if {@code fastaFile} is {@code null} or {@code basesPerLine} is 0 or negative.
     * @throws IOException if such exception is thrown when accessing the output path resources.
     */
    public FastaReferenceFileWriter(final Path fastaFile, final int basesPerLine, final Path faiFile, final Path dictFile)
        throws IOException
    {
        this(BucketUtils.createFile(fastaFile.toString()),
             basesPerLine,
             faiFile == null ? null : BucketUtils.createFile(faiFile.toString()),
             dictFile == null ? null : BucketUtils.createFile(dictFile.toString()));
    }

    // The actual common constructor, all other constructor delegate on this one:
    private FastaReferenceFileWriter(final OutputStream fastaOutput,
                                     final int basesPerLine,
                                     final OutputStream faiOutput,
                                     final OutputStream dictOutput) {
        this.defaultBasePerLine = checkBasesPerLine(basesPerLine);
        this.fastaOutput = new CountingOutputStream(fastaOutput);
        this.fastaWriter = new OutputStreamWriter(fastaOutput);
        this.faiWriter = faiOutput == null ? null : new OutputStreamWriter(faiOutput);
        this.dictOutput = dictOutput;
        this.sequenceNamesAndSizes = new LinkedHashMap<>();
    }

    private static Path defaultFaiFile(final boolean makeFaiFile, final Path fastaFile) {
        if (makeFaiFile) {
            return ReferenceSequenceFileFactory.getFastaIndexFileName(fastaFile);
        } else {
            return null;
        }
    }

    private static Path defaultDictFile(final boolean makeDictFile, final Path fastaFile) {
        if (makeDictFile) {
            return ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(fastaFile);
        } else {
            return null;
        }
    }

    // checks that a sequence name is valid.
    private static String checkSequenceName(final String name) {
        Utils.nonNull(name, "the sequence name cannot be null");
        Utils.validateArg(!name.isEmpty(), "the input sequence name cannot be null");
        for (int i = 0; i < name.length(); i++) {
            final char ch = name.charAt(i);
            Utils.validateArg(!Character.isWhitespace(ch), "the input name contains a blank characters: '" + StringUtils.escape(name) + "'");
            Utils.validateArg(!Character.isISOControl(ch), "the input name contains control characters: '" + StringUtils.escape(name) + "'");
            Utils.validateArg(ch != HEADER_START_CHAR, "the input name contains '" + HEADER_START_CHAR + "' character: '" + StringUtils.escape(name) + "'");
        }
        return name;
    }

    private void checkSequenceBases(final byte[] bases, final int offset, final int length) {
        Utils.nonNull(bases, "the input bases array cannot be null");
        final int to = offset + length;
        for (int i = offset; i < to; i++) {
            Utils.validateArg(Nucleotide.valueOf(bases[i]) != Nucleotide.INVALID,
                    "the input sequence contains invalid base calls like: " + StringUtils.escape(""+ (char) bases[i]));
        }
    }

    private static String checkDescription(final String description) {
        if (description == null || description.isEmpty()) {
            return "";
        } else {
            for (int i = 0; i < description.length(); i++) {
                Utils.validateArg(Character.isISOControl(description.charAt(i)), "the input name contains control characters: '" + StringUtils.escape(description) + "'");
            }
            return description;
        }
    }

    private static int checkBasesPerLine(final int value) {
        return ParamUtils.isPositive(value, "base per line must be 1 or greater");
    }


    /**
     * Starts the input of the bases of a new sequence.
     * <p>
     *     This operation automatically closes the previous sequence base input if any.
     * </p>
     * <p>
     *     The sequence name cannot contain any blank characters (as determined by {@link Character#isWhitespace(char)}),
     *     control characters (as determined by {@link Character#isISOControl(char)}) or the the FASTA header star character
     *     {@value #HEADER_START_CHAR}. It cannot be the empty string either ("").
     * </p>
     * <p>
     *     No description is included in the output.
     * </p>
     * <p>
     *     The input bases-per-line is set to the default provided at construction or {@link #DEFAULT_BASES_PER_LINE}
     *     if none was provided.
     * </p>
     * <p>
     *     This method cannot be called after the writer has been closed.
     * </p>
     * <p>
     *     It also will fail if no base was added to the previous sequence if any.
     * </p>
     * @param sequenceName the name of the new sequence.
     * @return this instance.
     * @throws IllegalArgumentException if any argument does not comply with requirements listed above or if a sequence
     *  with the same name has already been added to the writer.
     * @throws IllegalStateException if no base was added to the previous base or the writer is already closed.
     * @throws IOException if such exception is thrown when writing into the output resources.
     */
    public FastaReferenceFileWriter startSequence(final String sequenceName)
        throws IOException
    {
        return startSequence(sequenceName, "", defaultBasePerLine);
    }

    /**
     * Starts the input of the bases of a new sequence.
     * <p>
     *     This operation automatically closes the previous sequence base input if any.
     * </p>
     * <p>
     *     The sequence name cannot contain any blank characters (as determined by {@link Character#isWhitespace(char)}),
     *     control characters (as determined by {@link Character#isISOControl(char)}) or the the FASTA header star character
     *     {@value #HEADER_START_CHAR}. It cannot be the empty string either ("").
     * </p>
     * <p>
     *     The input bases-per-line must be 1 or greater.
     * </p>
     * <p>
     *     This method cannot be called after the writer has been closed.
     * </p>
     * <p>
     *     It also will fail if no base was added to the previous sequence if any.
     * </p>
     * @param sequenceName the name of the new sequence.
     * @param basesPerLine number of bases per line for this sequence.
     * @return this instance.
     * @throws IllegalArgumentException if any argument does not comply with requirements listed above or if a sequence
     *  with the same name has already been added to the writer.
     * @throws IllegalStateException if no base was added to the previous base or the writer is already closed.
     * @throws IOException if such exception is thrown when writing into the output resources.
     */
    public FastaReferenceFileWriter startSequence(final String sequenceName, final int basesPerLine)
        throws IOException
    {
        return startSequence(sequenceName, "", checkBasesPerLine(basesPerLine));
    }

    /**
     * Starts the input of the bases of a new sequence.
     * <p>
     *     This operation automatically closes the previous sequence base input if any.
     * </p>
     * <p>
     *     The sequence name cannot contain any blank characters (as determined by {@link Character#isWhitespace(char)}),
     *     control characters (as determined by {@link Character#isISOControl(char)}) or the the FASTA header star character
     *     {@value #HEADER_START_CHAR}. It cannot be the empty string either ("").
     * </p>
     * <p>
     *     The description cannot contain {@link Character#isISOControl(char)}. If set to {@code null} or the empty
     *     string ("") no description will be outputted.
     * </p>
     * <p>
     *     The input bases-per-line is set to the default provided at construction or {@link #DEFAULT_BASES_PER_LINE}
     *     if none was provided.
     * </p>
     * <p>
     *     This method cannot be called after the writer has been closed.
     * </p>
     * <p>
     *     It also will fail if no base was added to the previous sequence if any.
     * </p>
     * @param sequenceName the name of the new sequence.
     * @param description optional description for that sequence.
     * @return this instance.
     * @throws IllegalArgumentException if any argument does not comply with requirements listed above or if a sequence
     *  with the same name has already been added to the writer.
     * @throws IllegalStateException if no base was added to the previous base or the writer is already closed.
     * @throws IOException if such exception is thrown when writing into the output resources.
     */
    public FastaReferenceFileWriter startSequence(final String sequenceName, final String description)
        throws IOException
    {
        return startSequence(sequenceName, description, defaultBasePerLine);
    }

    /**
     * Starts the input of the bases of a new sequence.
     * <p>
     *     This operation automatically closes the previous sequence base input if any.
     * </p>
     * <p>
     *     The sequence name cannot contain any blank characters (as determined by {@link Character#isWhitespace(char)}),
     *     control characters (as determined by {@link Character#isISOControl(char)}) or the the FASTA header star character
     *     {@value #HEADER_START_CHAR}. It cannot be the empty string either ("").
     * </p>
     * <p>
     *     The description cannot contain {@link Character#isISOControl(char)}. If set to {@code null} or the empty
     *     string ("") no description will be outputted.
     * </p>
     * <p>
     *     The input bases-per-line must be 1 or greater.
     * </p>
     * <p>
     *     This method cannot be called after the writer has been closed.
     * </p>
     * <p>
     *     It also will fail if no base was added to the previous sequence if any.
     * </p>
     * @param sequenceName the name of the new sequence.
     * @param description optional description for that sequence.
     * @param basesPerLine number of bases per line for this sequence.
     * @return this instance.
     * @throws IllegalArgumentException if any argument does not comply with requirements listed above or if a sequence
     *  with the same name has already been added to the writer.
     * @throws IllegalStateException if no base was added to the previous base or the writer is already closed.
     * @throws IOException if such exception is thrown when writing into the output resources.
     */
    public FastaReferenceFileWriter startSequence(final String sequenceName, final String description, final int basesPerLine)
        throws IOException
    {
        assertIsNotClosed();
        checkSequenceName(sequenceName);
        final String nonNullDescription = checkDescription(description);
        checkBasesPerLine(basesPerLine);
        closeSequence();
        if (sequenceNamesAndSizes.containsKey(currentSequenceName)) {
            throw new IllegalArgumentException("the input sequence name '" + sequenceName + "' has already been added");
        }
        currentSequenceName = sequenceName;
        currentBasePerLine = basesPerLine;
        fastaWriter.append(HEADER_START_CHAR).append(sequenceName);
        if (!nonNullDescription.isEmpty()) {
            fastaWriter.append(HEADER_NAME_AND_DESCRIPTION_SEPARATOR).append(nonNullDescription);
        }
        fastaWriter.append(System.lineSeparator()).flush();
        currentSequenceOffset = fastaOutput.getCount();
        return this;
    }

    private void closeSequence()
        throws IOException
    {
        if (currentSequenceName != null) {
            if (currentBasesCount == 0) {
                throw new IllegalStateException("no base was added");
            }
            closeFaiSequence();
            sequenceNamesAndSizes.put(currentSequenceName, currentBasesCount);
            fastaOutput.write(LINE_SEPARATOR);
            currentBasesCount = 0;
        }
    }

    private void closeFaiSequence()
        throws IOException
    {
        if (faiWriter != null) {
            faiWriter.append(currentSequenceName).append('\t')
                     .append(String.valueOf(currentBasesCount)).append('\t')
                     .append(String.valueOf(currentSequenceOffset)).append('\t')
                     .append(String.valueOf(currentBasePerLine)).append(System.lineSeparator());
        }
    }

    /**
     * Adds bases to current sequence from a {@code byte} array.
     *
     * @param bases array containing the bases to be added.
     * @return this instance.
     * @throws IllegalArgumentException if {@bases} is {@code null} or
     *              the input array contains invalid bases (as assessed by: {@link Nucleotide#valueOf(byte)}).
     * @throws IllegalStateException if no sequence was started or the writer is already closed.
     * @throws IOException if such exception is throw when writing in any of the outputs.
     */
    public FastaReferenceFileWriter appendBases(final byte[] bases)
        throws IOException
    {
        return appendBases(bases, 0, bases.length);
    }

    /**
     * Adds bases to current sequence from a range in a {@code byte} array.
     *
     * @param bases array containing the bases to be added.
     * @param offset the position of the first base to add.
     * @param length how many bases to be added starting from position {@code offset}.
     * @return this instance.
     * @throws IllegalArgumentException if {@bases} is {@code null} or
     *              {@code offset} and {@code length} do not entail a valid range in {@code bases} or
     *              that range in {@base} contain invalid bases (as assessed by: {@link Nucleotide#valueOf(byte)}).
     * @throws IllegalStateException if no sequence was started or the writer is already closed.
     * @throws IOException if such exception is throw when writing in any of the outputs.
     */
    public FastaReferenceFileWriter appendBases(final byte[] bases, final int offset, final int length)
        throws IOException
    {
        assertIsNotClosed();
        assertSequenceOpen();
        checkSequenceBases(bases, offset, length);
        Utils.nonNull(bases, "the input bases cannot be negative");
        ParamUtils.isPositiveOrZero(offset, "the input offset cannot be negative");
        ParamUtils.isPositiveOrZero(length, "the input length must not be negative");
        final int to = offset + length;
        Utils.validateArg(to >= bases.length, "the length + offset goes beyond the end of " +
                "the input base array: '" + to + "' >= '" + bases.length + "'");
        int next = offset;
        while (next < to) {
            if (lineBasesCount == currentBasePerLine) {
                fastaOutput.write(LINE_SEPARATOR);
                lineBasesCount = 0;
            }
            int nextLength = Math.min(to - next, currentBasePerLine);
            fastaOutput.write(bases, next, nextLength);
            lineBasesCount += nextLength;
            next += nextLength;
        }
        currentBasesCount += length;
        return this;
    }

    private void assertSequenceOpen() {
        if (currentSequenceName == null) {
            throw new IllegalStateException("trying to add bases without starting a sequence");
        }
    }

    private void assertIsNotClosed() {
        if (closed) {
            throw new IllegalStateException("already closed");
        }
    }

    /**
     * Closes this writer flushing all remaining writing operation input the output resources.
     * <p>
     *     Further calls to {@link #appendBases} or {@link #startSequence} will result in an exception.
     * </p>
     *
     * @throws IOException if such exception is thrown when closing output writers and output streams.
     */
    public void close() throws IOException
    {
        if (!closed) {
            closeSequence();
            fastaWriter.close();
            //Implied by the one above:
            //fastaOutput.close();
            if (faiWriter != null) faiWriter.close();
            if (dictOutput != null) {
                outputDictionary();
                dictOutput.close();
            }
            closed = true;
        }
    }

    // print out the dictionary
    private void outputDictionary() {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(
                sequenceNamesAndSizes.entrySet().stream()
                    .map(e -> new SAMSequenceRecord(e.getKey(),  e.getValue().intValue()))
                    .collect(Collectors.toList())
        );
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(dictionary);
        try (@SuppressWarnings("unused") final SAMFileWriter samWriter =
                     new SAMFileWriterFactory().makeSAMWriter(header, false, dictOutput)) {}
    }
}
