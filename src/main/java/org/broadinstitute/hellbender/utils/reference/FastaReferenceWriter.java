package org.broadinstitute.hellbender.utils.reference;

import com.google.common.io.CountingOutputStream;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.apache.commons.io.output.NullWriter;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.codehaus.plexus.util.StringUtils;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Writes a FASTA formatted reference file.
 * <p>
 *     In addition it can also compose the index and dictionary files for the newly written reference file.
 * </p>
 * <p>
 *    Example:
 *    <code>
 *        String[] seqNames = ...;
 *        byte[][] seqBases = ...;
 *        ...
 *        try (final FastaReferenceWriter writer = new FastaReferenceFileWriter(outputFile)) {
 *            for (int i = 0; i < seqNames.length; i++) {
 *              writer.startSequence(seqNames[i]).appendBases(seqBases[i]);
 *            }
 *        }
 *    </code>
 * </p>
 * <p>
 *     The two main operations that one can invoke on a opened writer is {@link #startSequence} and {@link #appendBases}.
 *     The former indicates that we are going to append a new sequence to the output and is invoked once per sequence.
 *     The latter adds bases to the current sequence and can be called as many times as is needed.
 * </p>
 * <p>
 *     The writer will make sure that the output adheres to the FASTA reference sequence file format restrictions:
 *     <ul>
 *         <li>Sequence names are valid (non-empty, without space/blank, control characters),</li>
 *         <li>sequence description are valid (without control characters),</li>
 *         <li>bases are valid nucleotides ore IUPAC redundancy codes and X [ACGTNX...] (lower or uppercase are accepted),</li>
 *         <li>sequence cannot have 0 length,</li>
 *         <li>and that each sequence can only appear once in the output</li>
 *     </ul>
 * </p>
 */
public final class FastaReferenceWriter implements AutoCloseable {

    /**
     * Default number of bases per line.
     */
    public static final int DEFAULT_BASES_PER_LINE = 60;

    /**
     * Sequence header start character.
     */
    public static final char HEADER_START_CHAR = '>';

    /**
     * Character used to separate the sequence name and the description if any.
     */
    public static final char HEADER_NAME_AND_DESCRIPTION_SEPARATOR = ' ';

    /**
     * Charset used for all outputs; fixed to UTF-8.
     */
    private static final Charset CHARSET = Charset.forName("UTF-8");

    /**
     * The line separator string.
     */
    private static final char LINE_SEPARATOR_CHR = '\n';

    /**
     * Character used to separate the fields in a index file line.
     */
    private static final char INDEX_FIELD_SEPARATOR_CHR = '\t';

    /**
     * Convenient cached {@code byte[]} representation of the line separator.
     */
    private static final byte[] LINE_SEPARATOR = String.valueOf(LINE_SEPARATOR_CHR).getBytes(CHARSET);

    /**
     * Output stream to the main FASTA output.
     * <p>
     *     We use it also to count the number of bytes so far outputted thus the offset included in
     *     the index file entry.
     * </p>
     */
    private final CountingOutputStream fastaStream;

    /**
     * Writer for the index file.
     */
    private final Writer indexWriter;

    /**
     * Output writer to the output dictionary.
     */
    private final Writer dictWriter;

    /**
     * Output codec for the dictionary.
     */
    private final SAMSequenceDictionaryCodec dictCodec;

    /**
     * Default number of bases per line to be applied unless one is
     */
    private final int defaultBasePerLine;

    /**
     * Records the sequences that have been already fully appended to this writer.
     * <p>
     *     The key is the sequence name.
     * </p>
     * <p>
     *     The value is the sequence length in bases.
     * </p>
     */
    private final Map<String, Long> sequenceNamesAndSizes;

    /**
     * Bases per line to be applied to the sequence that is been currently appended to the output.
     */
    private int currentBasesPerLine;

    /**
     * Holds the number of bases in the current output line.
     */
    private int currentLineBasesCount;

    /**
     * Holds the number of bases so far appended for the current sequence.
     */
    private long currentBasesCount;

    /**
     * Holds the FASTA output file offset for the current sequence.
     */
    private long currentSequenceOffset;

    /**
     * Holds the name of the sequence that is been appended currently.
     */
    private String currentSequenceName;

    /**
     * Flag indicating whether this writer has been already closed.
     */
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
    public FastaReferenceWriter(final Path fastaFile, final boolean makeFaiOutput, final boolean makeDictOutput)
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
    public FastaReferenceWriter(final Path fastaFile, final int basesPerLine, final boolean makeFaiOutput,
                                final boolean makeDictOutput)
        throws IOException
    {
        this(Utils.nonNull(fastaFile, "the output fasta-file cannot be null"),
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
     * @param indexFile the path of the index file, if requested, {@code null} if none should be generated.
     * @param dictFile the path of the dictFile, if requested, {@code null} if nono should be generated.
     * @throws IllegalArgumentException if {@code fastaFile} is {@code null}.
     * @throws IOException if such exception is thrown when accessing the output path resources.
     */
    public FastaReferenceWriter(final Path fastaFile, final Path indexFile, final Path dictFile)
        throws IOException
    {
        this(fastaFile, DEFAULT_BASES_PER_LINE, indexFile, dictFile);
    }

    /**
     * Creates a reference FASTA file writer.
     * <p>
     *     You can specify a specific path for the index and dictionary file. If either set to {@code null} such
     *     a file won't be generated.
     * </p>
     *
     * @param fastaFile the output fasta file path.
     * @param indexFile the path of the index file, if requested, {@code null} if none should be generated.
     * @param dictFile the path of the dictFile, if requested, {@code null} if nono should be generated.
     * @throws IllegalArgumentException if {@code fastaFile} is {@code null} or {@code basesPerLine} is 0 or negative.
     * @throws IOException if such exception is thrown when accessing the output path resources.
     */
    public FastaReferenceWriter(final Path fastaFile, final int basesPerLine, final Path indexFile, final Path dictFile)
        throws IOException
    {
        // This code is a slight repeat of {@link #FastaReferenceWriter(OutputStream,int,OutputStream,OutputStream)
        // for the sake of avoiding creating output if basesPerLine is invalid.
        this.defaultBasePerLine = checkBasesPerLine(basesPerLine);

        this.fastaStream = new CountingOutputStream(Files.newOutputStream(Utils.nonNull(fastaFile)));
        this.indexWriter = indexFile == null ? new NullWriter() : new OutputStreamWriter(Files.newOutputStream(indexFile), CHARSET);
        final BufferedWriter dictWriter = new BufferedWriter(dictFile == null ? new NullWriter() : new OutputStreamWriter(Files.newOutputStream(dictFile), CHARSET));
        this.dictWriter = dictWriter;
        this.dictCodec = new SAMSequenceDictionaryCodec(dictWriter);
        this.dictCodec.encodeHeaderLine(false);
        this.sequenceNamesAndSizes = new LinkedHashMap<>();
    }

    /**
     * Creates a reference FASTA file writer.
     * <p>
     *     You can specify a specific output stream to each file: the main fasta output, its index and its dictionary.
     * </p>
     *
     * @param fastaOutput the output fasta file path.
     * @param indexOutput the output stream to the index file, if requested, {@code null} if none should be generated.
     * @param dictOutput the output stream to the dictFile, if requested, {@code null} if none should be generated.
     * @throws IllegalArgumentException if {@code fastaFile} is {@code null} or {@code basesPerLine} is 0 or negative.
     */
    public FastaReferenceWriter(final OutputStream fastaOutput,
                                 final int basesPerLine,
                                 final OutputStream indexOutput,
                                 final OutputStream dictOutput) {
        this.defaultBasePerLine = checkBasesPerLine(basesPerLine);
        this.fastaStream = new CountingOutputStream(Utils.nonNull(fastaOutput));
        this.indexWriter = indexOutput == null ? new NullWriter() : new OutputStreamWriter(indexOutput, CHARSET);
        final BufferedWriter dictWriter = new BufferedWriter(dictOutput == null ? new NullWriter() : new OutputStreamWriter(dictOutput, CHARSET));
        this.dictWriter = dictWriter;
        this.dictCodec = new SAMSequenceDictionaryCodec(dictWriter);
        this.dictCodec.encodeHeaderLine(false);
        this.sequenceNamesAndSizes = new LinkedHashMap<>();
    }

    private static Path defaultFaiFile(final boolean makeFaiFile, final Path fastaFile) {
        return makeFaiFile ? ReferenceSequenceFileFactory.getFastaIndexFileName(fastaFile) : null;
    }

    private static Path defaultDictFile(final boolean makeDictFile, final Path fastaFile) {
        return makeDictFile ? ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(fastaFile) : null;
    }

    // checks that a sequence name is valid.
    private static void checkSequenceName(final String name) {
        Utils.nonNull(name, "the sequence name cannot be null");
        Utils.validateArg(!name.isEmpty(), "the input sequence name cannot be null");
        for (int i = 0; i < name.length(); i++) {
            final char ch = name.charAt(i);
            if (Character.isWhitespace(ch)) {
                throw new IllegalArgumentException("the input name contains blank characters: '" + StringUtils.escape(name) + "'");
            } else if (Character.isISOControl(ch)) {
                throw new IllegalArgumentException("the input name contains control characters: '" + StringUtils.escape(name) + "'");
            }
        }
    }

    private static void checkSequenceBases(final byte[] bases, final int offset, final int length) {
        Utils.nonNull(bases, "the input bases array cannot be null");
        final int to = offset + length;
        for (int i = offset; i < to; i++) {
            final byte b = bases[i];
            if (!Nucleotide.decode(b).isValid()) {
                throw new IllegalArgumentException( "the input sequence contains invalid base calls like: "
                        + StringUtils.escape(""+ (char) b));
            }
        }
    }

    private static String checkDescription(final String description) {
        if (description == null || description.isEmpty()) {
            return "";
        } else {
            for (int i = 0; i < description.length(); i++) {
                final char c = description.charAt(i);
                if (Character.isISOControl(c) && c != '\t') { // tab is the only valid control char in the description.
                    throw new IllegalArgumentException("the input name contains non-tab control characters: '"
                            + StringUtils.escape(description) + "'");
                }
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
     * @throws IllegalStateException if no base was added to the previous sequence or the writer is already closed.
     * @throws IOException if such exception is thrown when writing into the output resources.
     */
    public FastaReferenceWriter startSequence(final String sequenceName)
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
     * @throws IllegalStateException if no base was added to the previous sequence or the writer is already closed.
     * @throws IOException if such exception is thrown when writing into the output resources.
     */
    public FastaReferenceWriter startSequence(final String sequenceName, final int basesPerLine)
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
     * @throws IllegalStateException if no base was added to the previous sequence or the writer is already closed.
     * @throws IOException if such exception is thrown when writing into the output resources.
     */
    public FastaReferenceWriter startSequence(final String sequenceName, final String description)
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
     * @throws IllegalArgumentException if any argument does not comply with requirements listed above.
     * @throws IllegalStateException if no base was added to the previous sequence or the writer is already closed of
     *         the sequence has been already added.
     * @throws IOException if such exception is thrown when writing into the output resources.
     */
    public FastaReferenceWriter startSequence(final String sequenceName, final String description, final int basesPerLine)
        throws IOException
    {
        assertIsNotClosed();
        checkSequenceName(sequenceName);
        final String nonNullDescription = checkDescription(description);
        checkBasesPerLine(basesPerLine);
        closeSequence();
        if (sequenceNamesAndSizes.containsKey(sequenceName)) {
            throw new IllegalStateException("the input sequence name '" + sequenceName + "' has already been added");
        }
        currentSequenceName = sequenceName;
        currentBasesPerLine = basesPerLine;
        final StringBuilder builder = new StringBuilder(sequenceName.length() + nonNullDescription.length() + 10);
        builder.append(HEADER_START_CHAR).append(sequenceName);
        if (!nonNullDescription.isEmpty()) {
            builder.append(HEADER_NAME_AND_DESCRIPTION_SEPARATOR).append(nonNullDescription);
        }
        fastaStream.write(builder.toString().getBytes(CHARSET));
        fastaStream.write(LINE_SEPARATOR);
        currentSequenceOffset = fastaStream.getCount();
        return this;
    }

    private void closeSequence()
        throws IOException
    {
        if (currentSequenceName != null) {
            if (currentBasesCount == 0) {
                throw new IllegalStateException("no base was added");
            }
            sequenceNamesAndSizes.put(currentSequenceName, currentBasesCount);
            writeIndexEntry();
            writeDictEntry();
            fastaStream.write(LINE_SEPARATOR);
            currentBasesCount = 0;
            currentLineBasesCount = 0;
            currentSequenceName = null;
        }
    }

    private void writeIndexEntry()
        throws IOException
    {
        indexWriter.append(currentSequenceName).append(INDEX_FIELD_SEPARATOR_CHR)
                     .append(String.valueOf(currentBasesCount)).append(INDEX_FIELD_SEPARATOR_CHR)
                     .append(String.valueOf(currentSequenceOffset)).append(INDEX_FIELD_SEPARATOR_CHR)
                     .append(String.valueOf(currentBasesPerLine)).append(INDEX_FIELD_SEPARATOR_CHR)
                     .append(String.valueOf(currentBasesPerLine + LINE_SEPARATOR.length)).append(LINE_SEPARATOR_CHR);
    }

    private void writeDictEntry() {
        dictCodec.encodeSequenceRecord(new SAMSequenceRecord(currentSequenceName, (int) currentBasesCount));
    }

    /**
     * Adds bases to current sequence from a {@code byte} array.
     *
     * @param bases array containing the bases to be added.
     * @return this instance.
     * @throws IllegalArgumentException if {@bases} is {@code null} or
     *              the input array contains invalid bases (as assessed by: {@link Nucleotide#decode(byte)}).
     * @throws IllegalStateException if no sequence was started or the writer is already closed.
     * @throws IOException if such exception is throw when writing in any of the outputs.
     */
    public FastaReferenceWriter appendBases(final byte[] bases)
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
     *              that range in {@base} contain invalid bases (as assessed by: {@link Nucleotide#decode(byte)}).
     * @throws IllegalStateException if no sequence was started or the writer is already closed.
     * @throws IOException if such exception is throw when writing in any of the outputs.
     */
    public FastaReferenceWriter appendBases(final byte[] bases, final int offset, final int length)
        throws IOException
    {
        assertIsNotClosed();
        assertSequenceOpen();
        checkSequenceBases(bases, offset, length);
        ParamUtils.isPositiveOrZero(offset, "the input offset cannot be negative");
        ParamUtils.isPositiveOrZero(length, "the input length must not be negative");
        final int to = offset + length;
        Utils.validateArg(to <= bases.length, "the length + offset goes beyond the end of " +
                "the input base array: '" + to + "' > '" + bases.length + "'");

        int next = offset;
        while (next < to) {
            if (currentLineBasesCount == currentBasesPerLine) {
                fastaStream.write(LINE_SEPARATOR);
                currentLineBasesCount = 0;
            }
            final int nextLength = Math.min(to - next, currentBasesPerLine - currentLineBasesCount);
            fastaStream.write(bases, next, nextLength);
            currentLineBasesCount += nextLength;
            next += nextLength;
        }
        currentBasesCount += length;
        return this;
    }

    /**
     * Appends a new sequence to the output.
     * <p>
     *     This is a convenient short handle for {@code startSequence(name).appendBases(bases)}.
     * </p>
     * <p>
     *     The new sequence remains open meaning that additional bases for that sequence can be added with additional calls to {@link #appendBases}.
     * </p>
     * @param name the name of the new sequence.
     * @param bases the (first) bases of the sequence.
     * @return a reference to this very same writer.
     * @throws IOException if such an exception is thrown when actually writing into the output streams/channels.
     * @throws IllegalArgumentException if either {@code name} or {@code bases} is {@code null} or contains an invalid value (e.g. unsupported bases or sequence names).
     * @throws IllegalStateException if the writer is already closed, a previous sequence (if any was opened) has no base appended to it or a sequence
     *   with such name was already appended to this writer.
     */
    public FastaReferenceWriter appendSequence(final String name, final byte[] bases) throws IOException {
        return startSequence(name).appendBases(bases);
    }

    /**
     * Appends a new sequence to the output with or without a description.
     * <p>
     *     This is a convenient short handle for {@code startSequence(name, description).appendBases(bases)}.
     * </p>
     * <p>
     *     A {@code null} or empty ("") description will be ignored (no description will be output).
     * </p>
     * <p>
     *     The new sequence remains open meaning that additional bases for that sequence can be added with additional calls to {@link #appendBases}.
     * </p>
     * @param name the name of the new sequence.
     * @param bases the (first) bases of the sequence.
     * @param description the description for the new sequence.
     * @return a reference to this very same writer.
     * @throws IOException if such an exception is thrown when actually writing into the output streams/channels.
     * @throws IllegalArgumentException if either {@code name} or {@code bases} is {@code null} or contains an invalid value (e.g. unsupported bases or sequence names). Also when
     *  the {@code description} contains unsupported characters.
     * @throws IllegalStateException if the writer is already closed, a previous sequence (if any was opened) has no base appended to it or a sequence
     *   with such name was already appended to this writer.
     */
    public FastaReferenceWriter appendSequence(final String name, final String description, final byte[] bases) throws IOException {
        return startSequence(name, description).appendBases(bases);
    }

    /**
     * Appends a new sequence to the output with or without a description and an alternative number of bases-per-line.
     * <p>
     *     This is a convenient short handle for {@code startSequence(name, description, bpl).appendBases(bases)}.
     * </p>
     * <p>
     *     A {@code null} or empty ("") description will be ignored (no description will be output).
     * </p>
     * <p>
     *     The new sequence remains open meaning that additional bases for that sequence can be added with additional calls to {@link #appendBases}.
     * </p>
     * @param name the name of the new sequence.
     * @param bases the (first) bases of the sequence.
     * @param description the description for the sequence.
     * @param basesPerLine alternative number of bases per line to be used for the sequence.
     * @return a reference to this very same writer.
     * @throws IOException if such an exception is thrown when actually writing into the output streams/channels.
     * @throws IllegalArgumentException if either {@code name} or {@code bases} is {@code null} or contains an invalid value (e.g. unsupported bases or sequence names). Also when the
     *  {@code description} contains unsupported characters or {@code basesPerLine} is 0 or negative.
     * @throws IllegalStateException if the writer is already closed, a previous sequence (if any was opened) has no base appended to it or a sequence
     *   with such name was already appended to this writer.
     */
    public FastaReferenceWriter appendSequence(final String name, final String description, final int basesPerLine, final byte[] bases) throws IOException {
        return startSequence(name, description, basesPerLine).appendBases(bases);
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
            if (sequenceNamesAndSizes.isEmpty()) {
                throw new IllegalStateException("no sequences where added to the reference");
            }
            fastaStream.close();
            indexWriter.close();
            dictWriter.close();
            closed = true;
        }
    }

    /**
     * Convenient method to write a FASTA file with a single sequence.
     *
     * @param whereTo the path to. must not be null.
     * @param makeIndex whether the index file should be written at its standard location.
     * @param makeDict whether the dictionary file should be written at it standard location.
     * @param name the sequence name, cannot contain white space, or control chracter or the header start character.
     * @param description the sequence description, "" if no description.
     * @param bases the sequence bases, cannot be {@code null}.
     * @throws IOException if such exception is thrown when writing in the output resources.
     */
    public static void writeSingleSequenceReference(final Path whereTo, final boolean makeIndex,
                                                    final boolean makeDict, final String name,
                                                    final String description, final byte[] bases)
            throws IOException
    {
        try (final FastaReferenceWriter writer = new FastaReferenceWriter(whereTo, makeIndex, makeDict)) {
            writer.startSequence(name, description);
            writer.appendBases(bases);
        }
    }

    /**
     * Convenient method to write a FASTA file with a single sequence.
     *
     * @param whereTo the path to. must not be null.
     * @param basesPerLine number of bases per line. must be 1 or greater.
     * @param makeIndex whether the index file should be written at its standard location.
     * @param makeDict whether the dictionary file should be written at it standard location.
     * @param name the sequence name, cannot contain white space, or control chracter or the header start character.
     * @param description the sequence description, "" if no description.
     * @param bases the sequence bases, cannot be {@code null}.
     * @throws IOException if such exception is thrown when writing in the output resources.
     */
    public static void writeSingleSequenceReference(final Path whereTo, final int basesPerLine, final boolean makeIndex,
                                                    final boolean makeDict, final String name,
                                                    final String description, final byte[] bases)
        throws IOException
    {
        try (final FastaReferenceWriter writer = new FastaReferenceWriter(whereTo, basesPerLine, makeIndex, makeDict)) {
            writer.startSequence(name, description);
            writer.appendBases(bases);
        }
    }
}
