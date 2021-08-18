package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.utils.ValidationUtils;
import org.apache.commons.compress.utils.CountingOutputStream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HashSet;
import java.util.Set;


/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/**
 * Writes a FASTA formatted reference file.
 * In addition it can also compose the index and dictionary files for the newly written reference file.
 * </p>
 * <p>
 * Example:
 * <pre>
 * String[] seqNames = ...;
 * byte[][] seqBases = ...;
 * ...
 * try (final FastaReferenceWriter writer = new FastaReferenceFileWriter(outputFile)) {
 *      for (int i = 0; i < seqNames.length; i++) {
 *          writer.startSequence(seqNames[i]).appendBases(seqBases[i]);
 *      }
 * }
 * </pre>
 * </p>
 * <p>
 * The two main operations that one can invoke on a opened writer is {@link #startSequence} and {@link #appendBases}.
 * The former indicates that we are going to append a new sequence to the output and is invoked once per sequence.
 * The latter adds bases to the current sequence and can be called as many times as is needed.
 * </p>
 * <p>
 * The writer will make sure that the output adheres to the FASTA reference sequence file format restrictions:
 * <ul>
 * <li>Sequence names are valid (non-empty, without space/blank, control characters),</li>
 * <li>Sequence description are valid (without control characters),</li>
 * <li>Bases are valid nucleotides or IUPAC redundancy codes and X [ACGTNX...] (lower or uppercase are accepted),</li>
 * <li>Sequence cannot have 0 length,</li>
 * <li>And that each sequence can only appear once in the output</li>
 * </ul>
 * </p>
 */
public class FastaIndexDictWriter implements AutoCloseable{

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
     * We use it also to count the number of bytes so far outputted thus the offset included in
     * the index file entry.
     * </p>
     */
    private Writer fastaStream;

    /**
     * Writer for the FAI index file.
     *
     * NOTE: GZI index writing (if necessary) is handled to BlockCompressedOutputStream class and thus is not controlled here.
     */
    private Writer faiIndexWriter;

    /**
     * Output writer to the output dictionary.
     */
    private Writer dictWriter;


    /**
     * the md5 digester (or null if not adding md5)
     */
    private MessageDigest md5Digester;

    /**
     * Output codec for the dictionary.
     */
    private SAMSequenceDictionaryCodec dictCodec;

    /**
     * Default number of bases per line to be applied unless one is
     */
    private int defaultBasePerLine;

    /**
     * Records the sequence names that have been already fully appended to this writer.
     */
    private final Set<String> sequenceNames = new HashSet<>();

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
     * Creates a reference FASTA file writer (private...use the builder: {@link FastaIndexDictWriter}.
     * <p>
     * You can specify a specific output stream to each file: the main fasta output, its index and its dictionary.
     * You can only provide a compressed stream to the fastaOutput, and only in the case that an index isn't written.
     * <p>
     * <p>
     * </p>
     *
     * @param fastaOutput the (uncompressed) output fasta file path.
     * @param indexOutput the (uncompressed) output stream to the index file, if requested, {@code null} if none should be generated.
     * @param dictOutput  the (uncompressed) output stream to the dictFile, if requested, {@code null} if none should be generated.
     * @throws IllegalArgumentException if {@code fastaFile} is {@code null} or {@code basesPerLine} is 0 or negative.
     */
    FastaIndexDictWriter(final int basesPerLine, final boolean addMd5,
                         final OutputStream fastaOutput,
                         final OutputStream indexOutput,
                         final OutputStream dictOutput) {

        try {
            this.md5Digester = addMd5 ? MessageDigest.getInstance("MD5") : null;
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("Couldn't get md5 algorithm!", e);
        }

        this.defaultBasePerLine = basesPerLine;
        this.fastaStream = NullWriter.NULL_WRITER;
        this.faiIndexWriter = indexOutput == null ? NullWriter.NULL_WRITER : new OutputStreamWriter(indexOutput, CHARSET);
        this.dictWriter = dictOutput == null ? NullWriter.NULL_WRITER : new OutputStreamWriter(dictOutput, CHARSET);
        this.dictCodec = new SAMSequenceDictionaryCodec(dictWriter);
        this.dictCodec.encodeHeaderLine(false);
    }

    // checks that a sequence name is valid.
    private static void checkSequenceName(final String name) {
        ValidationUtils.nonEmpty(name, "Sequence name");

        for (int i = 0; i < name.length(); i++) {
            final char ch = name.charAt(i);
            if (Character.isWhitespace(ch)) {
                throw new IllegalArgumentException("the input name contains blank characters: '" + name + "'");
            } else if (Character.isISOControl(ch)) {
                throw new IllegalArgumentException("the input name contains control characters: '" + name + "'");
            }
        }
    }

    private static void checkSequenceBases(final byte[] bases, final int offset, final int length) {
        ValidationUtils.nonNull(bases, "input bases");
        ValidationUtils.validateArg(bases.length >= offset + length, "Cannot validate bases beyond end of array.");
        final int to = offset + length;
        for (int i = offset; i < to; i++) {
            final byte b = bases[i];
            if (!SequenceUtil.isIUPAC(b)) {
                throw new IllegalArgumentException("the input sequence contains invalid base calls like: " + (char) b);
            }
        }
    }

    private static String checkDescription(final String description) {
        if (description == null || description.isEmpty()) {
            return "";
        }
        for (int i = 0; i < description.length(); i++) {
            final char c = description.charAt(i);
            if (Character.isISOControl(c) && c != '\t') { // tab is the only valid control char in the description.
                throw new IllegalArgumentException("the input name contains non-tab control characters: '" +
                        description + "'");
            }
        }
        return description;
    }

    /**
     * Starts the input of the bases of a new sequence.
     * <p>
     * This operation automatically closes the previous sequence base input if any.
     * </p>
     * <p>
     * The sequence name cannot contain any blank characters (as determined by {@link Character#isWhitespace(char)}),
     * control characters (as determined by {@link Character#isISOControl(char)}) or the the FASTA header start character
     * {@value #HEADER_START_CHAR}. It cannot be the empty string either ("").
     * </p>
     * <p>
     * No description is included in the output.
     * </p>
     * <p>
     * The input bases-per-line is set to the default provided at construction or {@link #DEFAULT_BASES_PER_LINE}
     * if none was provided.
     * </p>
     * <p>
     * This method cannot be called after the writer has been closed.
     * </p>
     * <p>
     * It also will fail if no base was added to the previous sequence if any.
     * </p>
     *
     * @param sequenceName the name of the new sequence.
     * @return this instance.
     * @throws IllegalArgumentException if any argument does not comply with requirements listed above or if a sequence
     *                                  with the same name has already been added to the writer.
     * @throws IllegalStateException    if no base was added to the previous sequence or the writer is already closed.
     * @throws IOException              if such exception is thrown when writing into the output resources.
     */
    public FastaIndexDictWriter startSequence(final String sequenceName)
            throws IOException {
        return startSequence(sequenceName, "", defaultBasePerLine);
    }

    /**
     * Starts the input of the bases of a new sequence.
     * <p>
     * This operation automatically closes the previous sequence base input if any.
     * </p>
     * <p>
     * The sequence name cannot contain any blank characters (as determined by {@link Character#isWhitespace(char)}),
     * control characters (as determined by {@link Character#isISOControl(char)}) or the the FASTA header start character
     * {@value #HEADER_START_CHAR}. It cannot be the empty string either ("").
     * </p>
     * <p>
     * The input bases-per-line must be 1 or greater.
     * </p>
     * <p>
     * This method cannot be called after the writer has been closed.
     * </p>
     * <p>
     * It also will fail if no base was added to the previous sequence if any.
     * </p>
     *
     * @param sequenceName the name of the new sequence.
     * @param basesPerLine number of bases per line for this sequence.
     * @return this instance.
     * @throws IllegalArgumentException if any argument does not comply with requirements listed above or if a sequence
     *                                  with the same name has already been added to the writer.
     * @throws IllegalStateException    if no base was added to the previous sequence or the writer is already closed.
     * @throws IOException              if such exception is thrown when writing into the output resources.
     */
    public FastaIndexDictWriter startSequence(final String sequenceName, final int basesPerLine)
            throws IOException {
        return startSequence(sequenceName, "", checkBasesPerLine(basesPerLine));
    }

    /**
     * Starts the input of the bases of a new sequence.
     * <p>
     * This operation automatically closes the previous sequence base input if any.
     * </p>
     * <p>
     * The sequence name cannot contain any blank characters (as determined by {@link Character#isWhitespace(char)}),
     * control characters (as determined by {@link Character#isISOControl(char)}) or the the FASTA header start character
     * {@value #HEADER_START_CHAR}. It cannot be the empty string either ("").
     * </p>
     * <p>
     * The description cannot contain {@link Character#isISOControl(char)}. If set to {@code null} or the empty
     * string ("") no description will be outputted.
     * </p>
     * <p>
     * The input bases-per-line is set to the default provided at construction or {@link #DEFAULT_BASES_PER_LINE}
     * if none was provided.
     * </p>
     * <p>
     * This method cannot be called after the writer has been closed.
     * </p>
     * <p>
     * It also will fail if no base was added to the previous sequence if any.
     * </p>
     *
     * @param sequenceName the name of the new sequence.
     * @param description  optional description for that sequence.
     * @return this instance.
     * @throws IllegalArgumentException if any argument does not comply with requirements listed above or if a sequence
     *                                  with the same name has already been added to the writer.
     * @throws IllegalStateException    if no base was added to the previous sequence or the writer is already closed.
     * @throws IOException              if such exception is thrown when writing into the output resources.
     */
    public FastaIndexDictWriter startSequence(final String sequenceName, final String description)
            throws IOException {
        return startSequence(sequenceName, description, defaultBasePerLine);
    }

    /**
     * Starts the input of the bases of a new sequence.
     * <p>
     * This operation automatically closes the previous sequence base input if any.
     * </p>
     * <p>
     * The sequence name cannot contain any blank characters (as determined by {@link Character#isWhitespace(char)}),
     * control characters (as determined by {@link Character#isISOControl(char)}) or the the FASTA header start character
     * {@value #HEADER_START_CHAR}. It cannot be the empty string either ("").
     * </p>
     * <p>
     * The description cannot contain {@link Character#isISOControl(char)}. If set to {@code null} or the empty
     * string ("") no description will be outputted.
     * </p>
     * <p>
     * The input bases-per-line must be 1 or greater.
     * </p>
     * <p>
     * This method cannot be called after the writer has been closed.
     * </p>
     * <p>
     * It also will fail if no base was added to the previous sequence if any.
     * </p>
     *
     * @param sequenceName the name of the new sequence.
     * @param description  optional description for that sequence.
     * @param basesPerLine number of bases per line for this sequence.
     * @return this instance.
     * @throws IllegalArgumentException if any argument does not comply with requirements listed above.
     * @throws IllegalStateException    if no base was added to the previous sequence or the writer is already closed of
     *                                  the sequence has been already added.
     * @throws IOException              if such exception is thrown when writing into the output resources.
     */
    public FastaIndexDictWriter startSequence(final String sequenceName, final String description, final int basesPerLine)
            throws IOException {
        assertIsNotClosed();
        //checkSequenceName(sequenceName);
        final String nonNullDescription = checkDescription(description);
        checkBasesPerLine(basesPerLine);
        currentSequenceName = sequenceName;
        currentBasesCount = basesPerLine;
        closeSequence();
        if (sequenceNames.contains(sequenceName)) {
            throw new IllegalStateException("the input sequence name '" + sequenceName + "' has already been added");
        }
        currentBasesPerLine = basesPerLine;
        final StringBuilder builder = new StringBuilder(sequenceName.length() + nonNullDescription.length() + 2);
        builder.append(HEADER_START_CHAR).append(sequenceName);
        if (!nonNullDescription.isEmpty()) {
            builder.append(HEADER_NAME_AND_DESCRIPTION_SEPARATOR).append(nonNullDescription);
        }

        if (md5Digester != null) {
            md5Digester.reset();
        }
        return this;
    }

    private void closeSequence()
            throws IOException {
        if (currentSequenceName != null) {
            if (currentBasesCount == 0) {
                throw new IllegalStateException("no base was added");
            }
            sequenceNames.add(currentSequenceName);
            writeIndexEntry();
            writeDictEntry();
            currentBasesCount = 0;
            currentLineBasesCount = 0;
            currentSequenceName = null;
        }
    }

    private void writeIndexEntry()
            throws IOException {
        faiIndexWriter.append(currentSequenceName).append(INDEX_FIELD_SEPARATOR_CHR)
                .append(String.valueOf(currentBasesCount)).append(INDEX_FIELD_SEPARATOR_CHR)
                .append(String.valueOf(currentSequenceOffset)).append(INDEX_FIELD_SEPARATOR_CHR)
                .append(String.valueOf(currentBasesPerLine)).append(INDEX_FIELD_SEPARATOR_CHR)
                .append(String.valueOf(currentBasesPerLine + LINE_SEPARATOR.length)).append(LINE_SEPARATOR_CHR);
    }

    private void writeDictEntry() {
        final SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord(currentSequenceName, (int) currentBasesCount);
        if (md5Digester != null) {
            samSequenceRecord.setMd5(SequenceUtil.md5DigestToString(md5Digester.digest()));
        }
        dictCodec.encodeSequenceRecord(samSequenceRecord);
    }

    /**
     * Adds bases to current sequence from a {@code byte} array.
     *
     * @param basesBases String containing the bases to be added.
     *                   string will be interpreted using ascii and will throw if any character is >= 127.
     * @return this instance.
     * @throws IllegalArgumentException if {@code bases} is {@code null} or
     *                                  the input array contains invalid bases (as assessed by: {@link SequenceUtil#isIUPAC(byte)}).
     * @throws IllegalStateException    if no sequence was started or the writer is already closed.
     * @throws IOException              if such exception is throw when writing in any of the outputs.
     */
    public FastaIndexDictWriter appendBases(final String basesBases)
            throws IOException {
        return appendBases(basesBases.getBytes(StandardCharsets.US_ASCII));
    }

    /**
     * Adds bases to current sequence from a {@code byte} array.
     * Will throw if any character is >= 127.
     *
     * @param bases array containing the bases to be added.
     * @return this instance.
     * @throws IllegalArgumentException if {@code bases} is {@code null} or
     *                                  the input array contains invalid bases (as assessed by: {@link SequenceUtil#isIUPAC(byte)}).
     * @throws IllegalStateException    if no sequence was started or the writer is already closed.
     * @throws IOException              if such exception is throw when writing in any of the outputs.
     */
    public FastaIndexDictWriter appendBases(final byte[] bases)
            throws IOException {
        return appendBases(bases, 0, bases.length);
    }

    /**
     * Adds bases to current sequence from a range in a {@code byte} array.
     * Will throw if any character is >= 127.
     *
     * @param bases  array containing the bases to be added.
     * @param offset the position of the first base to add.
     * @param length how many bases to be added starting from position {@code offset}.
     * @return this instance.
     * @throws IllegalArgumentException if {@code bases} is {@code null} or
     *                                  {@code offset} and {@code length} do not entail a valid range in {@code bases} or
     *                                  that range in {@code base} contain invalid bases (as assessed by: {@link SequenceUtil#isIUPAC(byte)}).
     * @throws IllegalStateException    if no sequence was started or the writer is already closed.
     * @throws IOException              if such exception is throw when writing in any of the outputs.
     */
    public FastaIndexDictWriter appendBases(final byte[] bases, final int offset, final int length)
            throws IOException {
        assertIsNotClosed();
        assertSequenceOpen();
        checkSequenceBases(bases, offset, length);
        ValidationUtils.validateArg(offset >= 0, "the input offset cannot be negative");
        ValidationUtils.validateArg(length >= 0, "the input length must not be negative");
        final int to = offset + length;
        ValidationUtils.validateArg(to <= bases.length, "the length + offset goes beyond the end of " +
                "the input base array: '" + to + "' > '" + bases.length + "'");

        int next = offset;
        while (next < to) {
            if (currentLineBasesCount == currentBasesPerLine) {
                currentLineBasesCount = 0;
            }
            final int nextLength = Math.min(to - next, currentBasesPerLine - currentLineBasesCount);
            if (md5Digester != null) {
                md5Digester.update(new String(bases, next, nextLength).toUpperCase().getBytes());
            }
            currentLineBasesCount += nextLength;
            next += nextLength;
        }
        currentBasesCount += length;
        return this;
    }

    /**
     * Appends a new sequence to the output.
     * <p>
     * This is a convenient short handle for {@code startSequence(name).appendBases(bases)}.
     * </p>
     * <p>
     * The new sequence remains open meaning that additional bases for that sequence can be added with additional calls to {@link #appendBases}.
     * </p>
     *
     * @param sequence a {@link ReferenceSequence} to add.
     * @return a reference to this very same writer.
     * @throws IOException              if such an exception is thrown when actually writing into the output streams/channels.
     * @throws IllegalArgumentException if either {@code name} or {@code bases} is {@code null} or contains an invalid value (e.g. unsupported bases or sequence names).
     * @throws IllegalStateException    if the writer is already closed, a previous sequence (if any was opened) has no base appended to it or a sequence
     *                                  with such name was already appended to this writer.
     */
    public FastaIndexDictWriter addSequence(ReferenceSequence sequence) throws IOException {
        return startSequence(sequence.getName()).appendBases(sequence.getBases());
    }

    /**
     * Appends a new sequence to the output with or without a description.
     * <p>
     * This is a convenient short handle for {@code startSequence(name, description).appendBases(bases)}.
     * </p>
     * <p>
     * A {@code null} or empty ("") description will be ignored (no description will be output).
     * </p>
     * <p>
     * The new sequence remains open meaning that additional bases for that sequence can be added with additional calls to {@link #appendBases}.
     * </p>
     *
     * @param name        the name of the new sequence.
     * @param bases       the (first) bases of the sequence.
     * @param description the description for the new sequence.
     * @return a reference to this very same writer.
     * @throws IOException              if such an exception is thrown when actually writing into the output streams/channels.
     * @throws IllegalArgumentException if either {@code name} or {@code bases} is {@code null} or contains an invalid value (e.g. unsupported bases or sequence names). Also when
     *                                  the {@code description} contains unsupported characters.
     * @throws IllegalStateException    if the writer is already closed, a previous sequence (if any was opened) has no base appended to it or a sequence
     *                                  with such name was already appended to this writer.
     */
    public FastaIndexDictWriter appendSequence(final String name, final String description, final byte[] bases) throws IOException {
        return startSequence(name, description).appendBases(bases);
    }

    /**
     * Appends a new sequence to the output with or without a description and an alternative number of bases-per-line.
     * <p>
     * This is a convenient short handle for {@code startSequence(name, description, bpl).appendBases(bases)}.
     * </p>
     * <p>
     * A {@code null} or empty ("") description will be ignored (no description will be output).
     * </p>
     * <p>
     * The new sequence remains open meaning that additional bases for that sequence can be added with additional calls to {@link #appendBases}.
     * </p>
     *
     * @param name         the name of the new sequence.
     * @param bases        the (first) bases of the sequence.
     * @param description  the description for the sequence.
     * @param basesPerLine alternative number of bases per line to be used for the sequence.
     * @return a reference to this very same writer.
     * @throws IOException              if such an exception is thrown when actually writing into the output streams/channels.
     * @throws IllegalArgumentException if either {@code name} or {@code bases} is {@code null} or contains an invalid value (e.g. unsupported bases or sequence names). Also when the
     *                                  {@code description} contains unsupported characters or {@code basesPerLine} is 0 or negative.
     * @throws IllegalStateException    if the writer is already closed, a previous sequence (if any was opened) has no base appended to it or a sequence
     *                                  with such name was already appended to this writer.
     */
    public FastaIndexDictWriter appendSequence(final String name, final String description, final int basesPerLine, final byte[] bases) throws IOException {
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
     * Further calls to {@link #appendBases} or {@link #startSequence} will result in an exception.
     * </p>
     *
     * @throws IOException if such exception is thrown when closing output writers and output streams.
     * @throws IllegalStateException if closing without writing any sequences or closing when writing a sequence is in progress
     */
    @Override
    public void close() throws IOException {
        if (!closed) {
            try {
                closeSequence();
                if (sequenceNames.isEmpty()) {
                    throw new IllegalStateException("no sequences were added to the reference");
                }
            } finally {
                closed = true;
//                fastaStream.close();
                faiIndexWriter.close();
                dictWriter.close();
            }
        }
    }

    /**
     * Convenient method to write a FASTA file with a single sequence.
     *
     * @param whereTo     the path to. must not be null.
     * @param makeIndex   whether the index file should be written at its standard location.
     * @param makeDict    whether the dictionary file should be written at it standard location.
     * @param name        the sequence name, cannot contain white space, or control chracter or the header start character.
     * @param description the sequence description, can be null or "" if no description.
     * @param bases       the sequence bases, cannot be {@code null}.
     * @throws IOException if such exception is thrown when writing in the output resources.
     */
    public static void writeSingleSequenceReference(final Path whereTo, final boolean makeIndex,
                                                    final boolean makeDict, final String name,
                                                    final String description, final byte[] bases)
            throws IOException {
        try (final FastaIndexDictWriter writer = new FastaIndexDictWriterBuilder().setFastaFile(whereTo).setMakeFaiOutput(makeIndex).setMakeDictOutput(makeDict).build()) {
            writer.startSequence(name, description);
            writer.appendBases(bases);
        }
    }

    /**
     * Convenient method to write a FASTA file with a single sequence.
     *
     * @param whereTo      the path to. must not be null.
     * @param basesPerLine number of bases per line. must be 1 or greater.
     * @param makeIndex    whether the index file should be written at its standard location.
     * @param makeDict     whether the dictionary file should be written at it standard location.
     * @param name         the sequence name, cannot contain white space, or control chracter or the header start character.
     * @param description  the sequence description, can be null or "" if no description.
     * @param bases        the sequence bases, cannot be {@code null}.
     * @throws IOException if such exception is thrown when writing in the output resources.
     */
    public static void writeSingleSequenceReference(final Path whereTo, final int basesPerLine, final boolean makeIndex,
                                                    final boolean makeDict, final String name,
                                                    final String description, final byte[] bases)
            throws IOException {
        try (final FastaIndexDictWriter writer = new FastaIndexDictWriterBuilder().setBasesPerLine(basesPerLine).setFastaFile(whereTo).setMakeFaiOutput(makeIndex).setMakeDictOutput(makeDict).build()) {
            writer.startSequence(name, description);
            writer.appendBases(bases);
        }
    }

    private static class NullWriter extends Writer {

        @Override
        public void write(char[] cbuf, int off, int len) throws IOException {
            // no op
        }

        @Override
        public void flush() throws IOException {
            // no op
        }

        @Override
        public void close() throws IOException {
            // no op
        }

        private NullWriter() {
        }

        /**
         * The only singleton instance of this class (no need for more!)
         */
        public final static FastaIndexDictWriter.NullWriter NULL_WRITER = new FastaIndexDictWriter.NullWriter();
    }

    protected static int checkBasesPerLine(final int value) {
        ValidationUtils.validateArg(value > 0, "bases per line must be 1 or greater");
        return value;
    }

}


