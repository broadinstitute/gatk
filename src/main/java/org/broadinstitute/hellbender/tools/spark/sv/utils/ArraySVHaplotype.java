package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFlag;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.reference.FastaReferenceWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Created by valentin on 10/11/17.
 */
@DefaultSerializer(ArraySVHaplotype.Serializer.class)
public class ArraySVHaplotype extends AbstractSVHaplotype {

    private final byte[] bases;

    public ArraySVHaplotype(final String name, final List<AlignmentInterval> intervals, final byte[] bases) {
        super(name, intervals);
        this.bases = bases;
    }

    private ArraySVHaplotype(final Kryo kryo, final Input input) {
        super(kryo, input);
        final int length = input.readInt();
        this.bases = input.readBytes(length);
    }

    @Override
    public int getLength() {
        return bases.length;
    }

    protected void unsafeCopyBases(final int offset, final byte[] dest, final int destOffset, final int length) {
        System.arraycopy(bases, offset, dest, destOffset, length);
    }

    public <T> List<List<AlignmentInterval>> align(final Iterable<T> input, Function<T, byte[]> basesOf) {
        try (final SingleReferenceSequenceAligner aligner = new SingleReferenceSequenceAligner(name, bases)) {
            final List<byte[]> seqs = Utils.stream(input).map(basesOf).collect(Collectors.toList());
            return aligner.align(seqs);
        } catch (final IOException ex) {
            throw new GATKException("could not create aligner", ex);
        }
    }

    public static class Serializer extends com.esotericsoftware.kryo.Serializer<ArraySVHaplotype> {

        @Override
        public void write(final Kryo kryo, final Output output, final ArraySVHaplotype object) {
            output.writeString(object.name);
            output.writeInt(object.intervals.size());
            for (final AlignmentInterval interval : object.intervals) {
                kryo.writeObject(output, interval);
            }
            output.writeInt(object.bases.length);
            output.write(object.bases);
        }

        @Override
        public ArraySVHaplotype read(Kryo kryo, Input input, Class<ArraySVHaplotype> type) {
            return new ArraySVHaplotype(kryo, input);
        }
    }

    private static class SingleReferenceSequenceAligner implements AutoCloseable {

        private final File fasta;
        private final File image;
        private final BwaMemAligner aligner;
        private final BwaMemIndex index;
        private final String name;

        public SingleReferenceSequenceAligner(final String name, final byte[] bases) {
            try {
                this.name = name;
                fasta = File.createTempFile("ssvh-temp", ".fasta");
                fasta.deleteOnExit();
                image = new File(fasta.getParentFile(), fasta.getName().replace(".fasta", ".img"));
                image.deleteOnExit();
                FastaReferenceWriter.writeSingleSequenceReference(fasta.toPath(), false, false, name, null, bases);
                BwaMemIndex.createIndexImageFromFastaFile(fasta.toString(), image.toString());
                index = new BwaMemIndex(image.toString());
                aligner = new BwaMemAligner(index);
            } catch (final IOException ex) {
                throw new GATKException("could not create index files", ex);
            }
        }

        public final List<List<AlignmentInterval>> align(final List<byte[]> seqs) {
            final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(seqs);
            final List<List<AlignmentInterval>> result = new ArrayList<>(alignments.size());
            final List<String> refNames = Collections.singletonList(name);
            for (int i = 0; i < alignments.size(); i++) {
                final int queryLength = seqs.get(i).length;
                final List<AlignmentInterval> intervals = alignments.get(i).stream()
                        .filter(bwa -> bwa.getRefId() >= 0)
                        .filter(bwa -> SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(bwa.getSamFlag()))
                        .map(bma -> new AlignmentInterval(bma, refNames, queryLength))
                        .collect(Collectors.toList()); // ignore secondary alignments.
                result.add(intervals);
            }
            return result;
        }

        @Override
        public void close() throws IOException {
            aligner.close();
            index.close();
            image.delete();
            fasta.delete();
        }
    }
}
