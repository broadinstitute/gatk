package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.ArraySVHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.OptionalDouble;

/**
 * Created by valentin on 9/16/17.
 */
@DefaultSerializer(SVContig.Serializer.class)
public class SVContig extends ArraySVHaplotype {

    /**
     * @return {@link Double#NaN} if there is no alignment against the reference haplotype.
     */
    public double getReferenceHaplotypeScore() {
        return referenceAlignmentScore == null ? Double.NaN : referenceAlignmentScore.getLog10Prob();
    }

    /**
     * @return {@link Double#NaN} if there is no alignment againts the alternative haplotype.
     */
    public double getAlternativeHaplotypeScore() {
        return alternativeAlignmentScore == null ? Double.NaN  : alternativeAlignmentScore.getLog10Prob();
    }

    public String getAlternativeScoreString() {
        return alternativeAlignmentScore == null ? "." : alternativeAlignmentScore.toString();
    }

    public String getReferenceScoreString() {
        return referenceAlignmentScore == null ? "." : referenceAlignmentScore.toString();
    }

    private List<AlignmentInterval> referenceAlignment;
    private List<AlignmentInterval> alternativeAlignment;
    private RealignmentScore referenceAlignmentScore;
    private RealignmentScore alternativeAlignmentScore;
    private OptionalDouble callQuality;

    public static SVContig of(final AlignedContig contig, final SVContext variant, final SAMSequenceDictionary dictionary, final int padding) {
        final AlignmentInterval primary = contig.getAlignments().stream()
                .filter(ai -> !ai.cigarAlong5to3DirectionOfContig.containsOperator(CigarOperator.H))
                .sorted(Comparator.comparingInt(ai -> -ai.alnScore))
                .findFirst().orElse(null);
        if (primary == null) {
            throw new IllegalArgumentException("no primary alignment!");
        }
        final List<SimpleInterval> breakPoints = variant.getBreakPointIntervals(padding, dictionary, true);
        final int mappingQuality = contig.getAlignments().stream()
                .filter(ai -> breakPoints.stream().anyMatch(bp -> bp.overlaps(ai.referenceSpan)))
                .mapToInt(ai -> ai.mapQual)
                .filter(mq -> mq != SAMRecord.UNKNOWN_MAPPING_QUALITY)
                .max().orElse(0);

        final String variantId = variant.getUniqueID();
        final SimpleInterval location = primary.referenceSpan.getStartInterval();

        return new SVContig(contig.getContigName(), location, variantId, contig.getContigSequence(), contig.getAlignments(), null, null, null, null, mappingQuality);
    }

    public static SVContig of(final GATKRead read, final RealignmentScoreParameters scoreParameters) {
        final String variantId = getMandatoryAttribute(read, GenotypeStructuralVariantsSpark.VARIANT_CONTEXT_TAG);
        final List<AlignmentInterval> refAln = getAlignmentIntervalsAttribute(read, GenotypeStructuralVariantsSpark.REFERENCE_ALIGNMENT_TAG);
        final List<AlignmentInterval> altAln = getAlignmentIntervalsAttribute(read, GenotypeStructuralVariantsSpark.ALTERNATIVE_ALIGNMENT_TAG);
        final RealignmentScore refScore = getOptionalAlignmentScore(read, GenotypeStructuralVariantsSpark.REFERENCE_SCORE_TAG, scoreParameters);
        final RealignmentScore altScore = getOptionalAlignmentScore(read, GenotypeStructuralVariantsSpark.ALTERNATIVE_SCORE_TAG, scoreParameters);
        final SimpleInterval location = new SimpleInterval(read.getAssignedContig(), read.getAssignedStart(), read.getAssignedStart());
        final int mappingQuality = read.getMappingQuality();
        return new SVContig(read.getName(), location, variantId, read.getBases(), AlignmentInterval.decodeList(read.getAttributeAsString(SAMTag.SA.name())), refAln, refScore, altAln, altScore, mappingQuality);
    }

    public void setReferenceHaplotypeAlignment(final List<AlignmentInterval> intervals, final RealignmentScore score) {
        referenceAlignment = intervals;
        referenceAlignmentScore = score;
        callQuality = OptionalDouble.empty();
    }

    public void setAlternativeHaplotypeAlignment(final List<AlignmentInterval> intervals, final RealignmentScore score) {
        alternativeAlignment = intervals;
        alternativeAlignmentScore = score;
        callQuality = OptionalDouble.empty();
    }

    public List<AlignmentInterval> geReferenceHaplotypeAlignment() {
        return referenceAlignment == null ? Collections.emptyList() : referenceAlignment;
    }

    public List<AlignmentInterval> getAlternativeHaplotypeAlignment() {

        return alternativeAlignment == null ? Collections.emptyList() : alternativeAlignment;
    }

    private static RealignmentScore getOptionalAlignmentScore(final GATKRead read, final String tag, final RealignmentScoreParameters parameters) {
        final String str = read.getAttributeAsString(tag);
        return str == null ? null : RealignmentScore.decode(str, parameters);
    }

    private static String getMandatoryAttribute(final GATKRead read, final String tag) {
        return ReadUtils.getOptionalStringAttribute(read, tag)
                .orElseThrow(() -> new UserException.BadInput("input read missing '" + tag + "' attribute"));
    }

    private static List<AlignmentInterval> getAlignmentIntervalsAttribute(final GATKRead read, final String tag) {
        final Optional<String> str = ReadUtils.getOptionalStringAttribute(read, tag);
        return str.map(AlignmentInterval::decodeList).orElse(null);
    }


    public SVContig(final String name, final Locatable loc, final String variantId,
                    final byte[] bases, final List<AlignmentInterval> originalReferenceAlignment, final List<AlignmentInterval> refAln, final RealignmentScore refScore,
                    final List<AlignmentInterval> altAln, final RealignmentScore altScore, final int mappingQuality) {
        super(name, originalReferenceAlignment, bases, variantId, SimpleInterval.valueOf(loc), mappingQuality, true);
        if (isReference() || isAlternative()) {
            throw new IllegalArgumentException("invalid assembled contig name, must not be reference or alternative like: " + name);
        }
        this.referenceAlignment = refAln;
        this.alternativeAlignment = altAln;
        this.referenceAlignmentScore = refScore;
        this.alternativeAlignmentScore = altScore;
        this.callQuality = OptionalDouble.empty();
    }

    private SVContig(final Kryo kryo, final Input input) {
        super(kryo, input);
        this.referenceAlignment = Serializer.readAlignment(kryo, input);
        this.alternativeAlignment = Serializer.readAlignment(kryo, input);
        this.referenceAlignmentScore = kryo.readObjectOrNull(input, RealignmentScore.class);
        this.alternativeAlignmentScore = kryo.readObjectOrNull(input, RealignmentScore.class);
        this.callQuality = OptionalDouble.empty();
    }

    public static class Serializer extends ArraySVHaplotype.Serializer<SVContig> {

        @Override
        public void write(Kryo kryo, Output output, SVContig object) {
            super.write(kryo, output, object);
            writeAlignment(kryo, output, object.referenceAlignment);
            writeAlignment(kryo, output, object.alternativeAlignment);
            kryo.writeObjectOrNull(output, object.referenceAlignmentScore, RealignmentScore.class);
            kryo.writeObjectOrNull(output, object.alternativeAlignmentScore, RealignmentScore.class);
        }

        @Override
        public SVContig read(Kryo kryo, Input input, Class<SVContig> type) {
            return new SVContig(kryo, input);
        }

        private void writeAlignment(final Kryo kryo, final Output output, final List<AlignmentInterval> alignment) {
            output.writeInt(alignment == null ? 0 : alignment.size());
            if (alignment != null) {
                for (final AlignmentInterval interval : alignment) {
                    kryo.writeObject(output, interval);
                }
            }
        }

        private static List<AlignmentInterval> readAlignment(final Kryo kryo, final Input input) {
            final int length = input.readInt();
            if (length <= 0) {
                return null;
            } else {
                final AlignmentInterval[] intervals = new AlignmentInterval[length];
                for (int i = 0; i < length; i++) {
                    intervals[i] = kryo.readObject(input, AlignmentInterval.class);
                }
                return Collections.unmodifiableList(Arrays.asList(intervals));
            }
        }
    }

    private int getMappingQuality() {
        return getReferenceAlignment().stream().mapToInt(ai -> ai.mapQual).filter(mq -> mq != SAMRecord.UNKNOWN_MAPPING_QUALITY).max().orElse(0);
    }

    public int getCallQuality() {
        if (!callQuality.isPresent()) {
            if (referenceAlignmentScore == null || alternativeAlignmentScore == null) {
                if (referenceAlignmentScore == alternativeAlignmentScore) {
                    callQuality = OptionalDouble.of(0);
                } else {
                    callQuality = OptionalDouble.of(getMappingQuality());
                }
            }
            final double refScore = getReferenceHaplotypeScore();
            final double altScore = getAlternativeHaplotypeScore();
            final double haplotypeQualDiff = Math.abs(refScore - altScore);
            final double originalMappingQuality = getMappingQuality();
            callQuality = OptionalDouble.of(Math.min(haplotypeQualDiff, originalMappingQuality));
        }
        return (int) Math.round(callQuality.getAsDouble());
    }
}
