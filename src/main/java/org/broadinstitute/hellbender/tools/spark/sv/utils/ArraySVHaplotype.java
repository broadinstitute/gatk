package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFlag;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.FastaReferenceWriter;
import org.broadinstitute.hellbender.utils.report.GATKReportColumnFormat;

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

    protected final byte[] bases;
    protected final boolean isContig;

    public ArraySVHaplotype(final String name, final List<AlignmentInterval> intervals, final byte[] bases, final String variantId, final SimpleInterval variantLocation, final boolean isContig) {
        super(name, intervals, variantId, variantLocation);
        this.bases = bases;
        this.isContig = isContig;
    }

    public static ArraySVHaplotype of(final GATKRead read) {
        final byte[] bases = read.getBases();
        final boolean isContig = "CTG".equals(read.getReadGroup());
        final String variantId = read.getAttributeAsString("VC");
        final List<AlignmentInterval> saIntervals = read.hasAttribute("SA") ? AlignmentInterval.decodeList(read.getAttributeAsString("SA")) : Collections.emptyList();
        final List<AlignmentInterval> allIntervals;
        if (read.isUnmapped()) {
            allIntervals = saIntervals;
        } else {
            allIntervals = new ArrayList<>();
            allIntervals.add(new AlignmentInterval(read));
            allIntervals.addAll(saIntervals);
        }
        return new ArraySVHaplotype(read.getName(), allIntervals, bases, variantId,
                new SimpleInterval(read.getAssignedContig(), read.getAssignedStart()), isContig);
    }

    protected ArraySVHaplotype(final Kryo kryo, final Input input) {
        super(kryo, input);
        final int length = input.readInt();
        this.bases = input.readBytes(length);
        this.isContig = input.readBoolean();
    }

    @Override
    public int getLength() {
        return bases.length;
    }

    @Override
    public boolean isContig() {
        return isContig;
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

    public static class Serializer<S extends ArraySVHaplotype> extends AbstractSVHaplotype.Serializer<S> {

        @Override
        public void write(final Kryo kryo, final Output output, final S object) {
            super.write(kryo, output, object);
            output.writeInt(object.bases.length);
            output.write(object.bases);
            output.writeBoolean(object.isContig);
        }

        @Override
        @SuppressWarnings("unchecked")
        public S read(Kryo kryo, Input input, Class<S> type) {
            if (type != ArraySVHaplotype.class) {
                throw new IllegalArgumentException("the input class must be: " + ArraySVHaplotype.class);
            } else {
                return (S) new ArraySVHaplotype(kryo, input);
            }
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
                        .filter(bwa -> SAMFlag.SECONDARY_ALIGNMENT.isUnset(bwa.getSamFlag()))
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
