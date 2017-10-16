package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.FilterLongReadAlignmentsSAMSpark;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Created by valentin on 10/11/17.
 */
public abstract class SVHaplotype  {

    public static final String REF_HAPLOTYPE_NAME = "ref";
    public static final String ALT_HAPLOTYPE_NAME = "alt";

    /**
     * Returns the cigar of this haplotype versus the reference
     * @return
     */
    public abstract List<AlignmentInterval> getReferenceAlignmentIntervals();

    public abstract String getName();

    public abstract int getLength();

    /**
     * Copies haplotypes bases into an array assuming that the argument values passed are correct.
     *
     * @param offset
     * @param dest
     * @param destOffset
     * @param length
     */
    protected abstract void unsafeCopyBases(final int offset, final byte[] dest, final int destOffset, final int length);

    /**
     * Copies bases from the haplotype into an array.
     * @param offset
     * @param whereTo
     * @param destOffset
     * @param length
     * @throws IllegalArgumentException if {@code whereTo} is {@code null} or the indeces and length passed are not correct.
     */
    public void copyBases(final int offset, final byte[] whereTo, final int destOffset, final int length) {
        ParamUtils.isPositiveOrZero(length, "length must be positive or zero");
        ParamUtils.isPositiveOrZero(offset, "dest offset");
        ParamUtils.isPositiveOrZero(destOffset, "the offset must be positive or zero");
        Utils.validateArg(destOffset + length <= whereTo.length, "the to index must be less than the length of the haplotype");
        Utils.validateArg(offset + length <= getLength(), "the from index cannot be larger than the to index" );
        Utils.nonNull(whereTo);
        unsafeCopyBases(offset, whereTo, destOffset, length);
    }

    public byte[] getBases(final int from, final int length) {
        Utils.validateArg(from + length <= getLength(), "the to index must be less than the length of the haplotype");
        ParamUtils.isPositiveOrZero(from, "the from index cannot be negative");
        final byte[] result = new byte[length];
        unsafeCopyBases(from, result, 0, length);
        return result;
    }

    public byte[] getBases() {
        return getBases(0, getLength());
    }

    public abstract <T> List<List<AlignmentInterval>> align(final Iterable<T> input, Function<T, byte[]> basesOf);

    public <T> List<AlignedContig> alignContigs(final Iterable<AlignedContig> contigs) {
        final List<String> names = Utils.stream(contigs).map(c -> c.contigName).collect(Collectors.toList());
        final List<byte[]> bases = Utils.stream(contigs).map(c -> c.contigSequence).collect(Collectors.toList());
        final List<List<AlignmentInterval>> intervals = align(bases, Function.identity());
        final List<AlignedContig> result = new ArrayList<>(intervals.size());
        final Set<String> haplotypeName = Collections.singleton(getName());
        for (int i = 0; i < intervals.size(); i++) {
            final List<List<AlignmentInterval>> bestCombos = FilterLongReadAlignmentsSAMSpark.pickBestConfigurations(names.get(i), intervals.get(i), haplotypeName);
            final AlignedContig alignedContig = new AlignedContig(names.get(i), bases.get(i), bestCombos.get(0), bestCombos.size() > 1);
            result.add(alignedContig);
        }
        return result;
    }

}
