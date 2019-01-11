package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMFlag;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.reference.FastaReferenceWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * Encompasses an aligner to a single-sequence reference.
 * <p>
 *     This is an {@link AutoCloseable} mean to be use in try-with-resource constructs.
 * </p>
 * <p>
 *     If you don't do so, please remember to call {@link #close} at the end to free
 *     resources.
 * </p>
 */
public class SingleSequenceReferenceAligner<T, U> implements AutoCloseable {

    private final File image;
    private final BwaMemAligner aligner;
    private final BwaMemIndex index;
    private final List<String> refNames;
    private boolean closed = false;

    @FunctionalInterface
    public interface TriFunction<T, U, V, W>{
        W apply(T t, U u, V v);
    }

    public static final Predicate<? super BwaMemAlignment> NO_FILTER = x -> true;

    private final Function<? super T, List<byte[]>> basesOf;
    private final TriFunction<? super T, List<List<BwaMemAlignment>>, List<String>, ? extends U> alignmentOf;
    private final Predicate<? super BwaMemAlignment> alignmentFilter;

    public SingleSequenceReferenceAligner(final String name, final byte[] bases,
                                          final Function<? super T, List<byte[]>> basesOf,
                                          final TriFunction<? super T, List<List<BwaMemAlignment>>, List<String>, ? extends U> alignmentOf) {
        this(name, bases, basesOf, alignmentOf, NO_FILTER);
    }

    public SingleSequenceReferenceAligner(final String name, final byte[] bases,
                                          final Function<? super T, List<byte[]>> basesOf,
                                          final TriFunction<? super T, List<List<BwaMemAlignment>>, List<String>, ? extends U> alignmentOf,
                                          final Predicate<? super BwaMemAlignment> alignmentFilter) {
        Utils.nonNull(name, "the input reference name cannot be null");
        Utils.nonNull(bases, "the input bases cannot be null");
        this.basesOf = Utils.nonNull(basesOf);
        this.alignmentOf = Utils.nonNull(alignmentOf);
        this.alignmentFilter = Utils.nonNull(alignmentFilter);
        Utils.validate(bases.length > 0, "the reference contig bases sequence must have at least one base");
        try {
            final File fasta = File.createTempFile("ssvh-temp", ".fasta");
            fasta.deleteOnExit();
            image = new File(fasta.getParentFile(), fasta.getName().replace(".fasta", ".img"));
            image.deleteOnExit();
            FastaReferenceWriter.writeSingleSequenceReference(fasta.toPath(), false, false, name, null, bases);
            BwaMemIndex.createIndexImageFromFastaFile(fasta.toString(), image.toString());
            fasta.delete(); // we don't need the fasta around.
            index = new BwaMemIndex(image.toString());
            aligner = new BwaMemAligner(index);
        } catch (final IOException ex) {
            throw new GATKException("could not create index files", ex);
        }
        refNames = Collections.singletonList(name);
    }

    /**
     * Gives access to the underlying aligner so that you can modify its options.
     * <p>
     *     You could align sequences directly thru the return object but in that case you will lose the
     *     {@link BwaMemAlignment} to {@link AlignmentInterval} translation.
     * </p>
     *
     * @return never {@code null}.
     */
    public BwaMemAligner getAligner() {
        return aligner;
    }

    public List<U> align(final Iterable<? extends T> inputs) {
        Utils.nonNull(inputs);
        return align(Utils.stream(inputs).collect(Collectors.toList()));
    }

    /**
     * Aligns the input object returning a list of the outputs in the corresponding order.
     * @param inputs
     * @return
     */
    public List<U> align(final List<? extends T> inputs) {
        checkNotClosed();
        Utils.nonNull(inputs, "the input sequence array cannot be null");
        final List<List<byte[]>> seqs = inputs.stream().map(basesOf).collect(Collectors.toList());
        final List<byte[]> flattenSeqs = seqs.stream().flatMap(Collection::stream).collect(Collectors.toList());
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(flattenSeqs);
        if (alignments.size() != flattenSeqs.size()) { // paranoiah??
            throw new IllegalStateException("something went terribly wrong and the number of returned alignment list does " +
                    "not correspond to the number of input sequences: " + alignments.size() + " != " + flattenSeqs.size());
        }
        final List<U> result = new ArrayList<>(inputs.size());
        int nextAlignmentIndex = 0;
        for (int i = 0; i < inputs.size(); i++) {
            final T inputObject = inputs.get(i);
            final List<byte[]> sequences = seqs.get(i);
            final List<List<BwaMemAlignment>> relevantAlignments =
                    alignments.subList(nextAlignmentIndex, nextAlignmentIndex += sequences.size());
            final List<List<BwaMemAlignment>> filteredAlignments;
            if (alignmentFilter == NO_FILTER) {
                filteredAlignments = relevantAlignments;
            } else {
                filteredAlignments = new ArrayList<>(relevantAlignments.size());
                for (int j = 0; j < relevantAlignments.size(); j++) {
                    filteredAlignments.add(relevantAlignments.get(j).stream().filter(alignmentFilter).collect(Collectors.toList()));
                }
            }
            final U outputObject = alignmentOf.apply(inputObject, filteredAlignments, refNames);
            result.add(outputObject);
        }
        return result;
    }

    /**
     * Composes a map of the aligned sequences.
     * <p>
     *     The key of such a map would be determined by the input object and the output alignment.
     * </p>
     * <p>
     *     Iterations over the entries, keys and values of the resulting tree will follow the order
     *     of the input objects. In case of key collisions (more than one input object, alignment result in the same key)
     *     then we only keep the first occurrence of such key.
     * </p>
     * @param inputs the input objects to align.
     * @param keyOf function that composed the key given the input and output objects.
     * @param <V>
     * @return never {@code null}.
     */
    public <V> Map<V, U> align(final List<? extends T> inputs, final BiFunction<? super T, ? super U, ? extends V> keyOf) {
        final List<U> outs = align(inputs);
        final LinkedHashMap<V, U> result = new LinkedHashMap<>(outs.size());
        for (int i = 0; i < outs.size(); i++) {
            result.putIfAbsent(keyOf.apply(inputs.get(i), outs.get(i)), outs.get(i));
        }
        return result;
    }

    /**
     * Composes a contig aligner from an arbitrary input type given contig name and base sequence generation functions.
     * <p>
     *     Supplementary alignments will be ignored.
     * </p>
     *
     * @param refName the name of the only reference sequence.
     * @param refBases the bases for the only reference sequence.
     * @param nameOf function that produces the aligned contig name based on the input object.
     * @param basesOf function that produces the contig bases sequence based on the input object.
     * @param <T> the type-parameter of the input object.
     * @return never {@code null}.
     */
    public static <T> SingleSequenceReferenceAligner<T, AlignedContig> contigsAligner(final String refName, final byte[] refBases,
                                                                                      final Function<? super T, String> nameOf,
                                                                                      final Function<? super T, byte[]> basesOf) {
        return new SingleSequenceReferenceAligner<>(refName, refBases,
                t -> Collections.singletonList(basesOf.apply(t)),
                (t, bma, refNames) -> {
                        final String name = nameOf.apply(t);
                        final byte[] bases = basesOf.apply(t);
                        final List<AlignmentInterval> intervals = bma.get(0).stream().map(b -> new AlignmentInterval(b, refNames, bases.length)).collect(Collectors.toList());
                        return new AlignedContig(name, bases, intervals);
                },
                bma -> bma.getRefId() >= 0 && SAMFlag.SECONDARY_ALIGNMENT.isUnset(bma.getSamFlag()));
    }

    private void checkNotClosed() {
        if (closed) {
            throw new IllegalStateException("operation not allowed once the aligner is closed");
        }
    }

    @Override
    public void close() throws IOException {
        if (!closed) {
            aligner.close();
            closed = true;
            try {
                index.close();
            } finally {
                image.delete();
            }
        }
    }
}
