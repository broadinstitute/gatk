package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFlag;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.FilterLongReadAlignmentsSAMSpark;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.reference.FastaReferenceWriter;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Created by valentin on 10/11/17.
 */
@DefaultSerializer(ShortSVHaplotype.Serializer.class)
public class ShortSVHaplotype extends SVHaplotype {

    private final String name;
    private final List<AlignmentInterval> intervals;
    private final byte[] bases;

    protected ShortSVHaplotype(final String name, final List<AlignmentInterval> intervals, final byte[] bases) {
        this.name = name;
        this.intervals = intervals;
        this.bases = bases;
    }

    private ShortSVHaplotype(final Kryo kryo, final Input input) {
        this.name = input.readString();
        final int length = input.readInt();
        this.bases = input.readBytes(length);
        final int numberOfIntervals = input.readInt();
        final List<AlignmentInterval> intervals = new ArrayList<>(numberOfIntervals);
        for (int i = 0; i < numberOfIntervals; i++) {
            intervals.add(kryo.readObject(input, AlignmentInterval.class));
        }
        this.intervals = Collections.unmodifiableList(intervals);
    }

    @Override
    public List<AlignmentInterval> getReferenceAlignmentIntervals() {
        return intervals;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public int getLength() {
        return bases.length;
    }

    protected void unsafeCopyBases(final int offset, final byte[] dest, final int destOffset, final int length) {
        System.arraycopy(bases, offset, dest, destOffset, length);
    }

    public <T> List<List<AlignmentInterval>> align(final Iterable<T> input, Function<T, byte[]> basesOf) {
        try {
            final Path tempFasta = Files.createTempFile("ssvh-ref", ".fasta");
            try {
                final Path tempImg = Files.createTempFile("ssvh-ref", ".img");
                FastaReferenceWriter.writeSingleSequenceReference(tempFasta, false, false, name, null, bases);
                BwaMemIndex.createIndexImageFromFastaFile(tempFasta.toString(), tempImg.toString());
                try (final BwaMemIndex index = new BwaMemIndex(tempImg.toString());
                     final BwaMemAligner aligner = new BwaMemAligner(index)) {
                    final List<String> haplotypeNames = Collections.singletonList(name);
                    final List<byte[]> seqs = Utils.stream(input).map(basesOf).collect(Collectors.toList());
                    final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(seqs);
                    final List<List<AlignmentInterval>> result = new ArrayList<>(alignments.size());

                    for (int i = 0; i < alignments.size(); i++) {
                         final int queryLength = seqs.get(i).length;
                         final List<AlignmentInterval> intervals = alignments.get(i).stream()
                                 .filter(bwa -> bwa.getRefId() >= 0)
                                 .filter(bwa -> SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(bwa.getSamFlag()))
                                 .map(bma -> new AlignmentInterval(bma, haplotypeNames, queryLength))
                                 .collect(Collectors.toList()); // ignore secondary alignments.
                         result.add(intervals);
                    }
                    return result;
                } finally {
                    tempImg.getFileSystem().provider().deleteIfExists(tempImg);
                }
            } catch (final IOException ex) {
                throw new GATKException("could not create temporal image file");
            } finally {
                tempFasta.getFileSystem().provider().deleteIfExists(tempFasta);
            }
        } catch (final IOException ex) {
            throw new GATKException("could not create temporal fasta file.");
        }
    }

    public static class Serializer extends com.esotericsoftware.kryo.Serializer<ShortSVHaplotype> {

        @Override
        public void write(final Kryo kryo, final Output output, final ShortSVHaplotype object) {
            output.writeString(object.name);
            output.writeInt(object.bases.length);
            output.write(object.bases);
            output.writeInt(object.intervals.size());
            for (final AlignmentInterval interval : object.intervals) {
                kryo.writeObject(output, interval);
            }
        }

        @Override
        public ShortSVHaplotype read(Kryo kryo, Input input, Class<ShortSVHaplotype> type) {
            return new ShortSVHaplotype(kryo, input);
        }
    }
}
