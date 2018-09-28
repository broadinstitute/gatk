package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;

@DefaultSerializer(SVKmerHuge.Serializer.class)
public final class SVKmerHuge extends SVKmer implements Comparable<SVKmerHuge> {
    private final long[] vals; // little-endian order
    private int hashVal;

    public SVKmerHuge( final int kSize ) {
        Utils.validateArg( kSize > 0, "Kmer length must be positive.");
        vals = new long[(kSize + 31)/32];
    }

    public SVKmerHuge( final SVKmerHuge that ) {
        this.vals = Arrays.copyOf( that.vals, that.vals.length );
        this.hashVal = that.hashVal;
    }

    private SVKmerHuge( final long[] vals ) {
        this.vals = vals;
    }

    private SVKmerHuge( final Input input ) {
        final int len = input.readInt();
        vals = input.readLongs(len);
    }

    private void serialize( final Output output ) {
        output.writeInt(vals.length);
        output.writeLongs(vals);
    }

    @Override
    public SVKmerHuge successor( final Base base, final int kSize ) {
        final long[] newVals = Arrays.copyOf(vals, vals.length);
        int idx = newVals.length;
        newVals[idx - 1] &= (1L << (((kSize - 1) & 0x1F) << 1)) - 1L;
        while ( --idx > 0 ) {
            newVals[idx] = (newVals[idx] << 2) | (newVals[idx-1] >>> 62);
        }
        newVals[0] = (newVals[0] << 2) | base.value;
        return new SVKmerHuge(newVals);
    }

    @Override
    public SVKmerHuge predecessor( final Base base, final int kSize ) {
        final long[] newVals = Arrays.copyOf(vals, vals.length);
        final int nnn = newVals.length - 1;
        for ( int idx = 0; idx != nnn; ++idx ) {
            newVals[idx] = (newVals[idx] >>> 2) | (newVals[idx+1] << 62);
        }
        newVals[nnn] = (newVals[nnn] >>> 2) | (base.value << (((kSize - 1) & 0x1F) << 1));
        return new SVKmerHuge(newVals);
    }

    @Override
    public boolean isCanonical( final int kSize ) {
        Utils.validateArg( (kSize & 1) != 0, "Kmer length must be odd to canonicalize.");
        return (vals[kSize >>> 6] & (1L << (kSize & 0x3f))) == 0;
    }

    @Override
    public SVKmerHuge reverseComplement( final int kSize ) {
        final long[] newVals = Arrays.copyOf(vals, vals.length);
        for ( int idx1 = 0, idx2 = newVals.length - 1; idx1 < idx2; ++idx1, --idx2 ) {
            long tmp = reverseComplement(newVals[idx1]);
            newVals[idx1] = reverseComplement(newVals[idx2]);
            newVals[idx2] = tmp;
        }
        if ( (newVals.length & 1) != 0 ) {
            final int idx = newVals.length >> 1;
            newVals[idx] = reverseComplement(newVals[idx]);
        }
        final int nnn = newVals.length - 1;
        final int lShift = (kSize & 0x1F) << 1;
        if ( lShift != 0 ) {
            final int rShift = 64 - lShift;
            for ( int idx = 0; idx != nnn; ++idx ) {
                newVals[idx] = (newVals[idx] >>> rShift) | (newVals[idx + 1] << lShift);
            }
            newVals[nnn] >>>= rShift;
        }
        return new SVKmerHuge(newVals);
    }

    @Override
    public SVKmerHuge canonical( final int kSize ) { return isCanonical(kSize) ? this : reverseComplement(kSize); }

    @Override
    public Base firstBase( final int kSize ) {
        return baseValues.get((int)((vals[vals.length-1] >>> (((kSize - 1) & 0x1f) << 1)) & 0x3L));
    }

    @Override
    public Base lastBase() {
        return baseValues.get((int)(vals[0] & 0x3L));
    }

    @Override
    public int firstTrimer( final int kSize ) {
        final int idx = vals.length - 1;
        final int basesInLastWord = kSize & 0x1F;
        switch ( basesInLastWord ) {
            case 0: return (int)(vals[idx] >>> 58);
            case 1: return (int)((vals[idx] << 4) | (vals[idx-1] >>> 60));
            case 2: return (int)((vals[idx] << 2) | (vals[idx-1] >>> 62));
            case 3: return (int)(vals[idx]);
            default: return (int)(vals[idx] >>> ((basesInLastWord - 3) << 1));
        }
    }

    @Override
    public int lastTrimer() {
        return (int)(vals[0] & 0x3FL);
    }

    @Override
    public SVKmerHuge removeFirstAndLastBase( final int kSize ) {
        final SVKmerHuge result = predecessor(Base.T, kSize);
        result.vals[vals.length-1] &= (1L << (((kSize - 2) & 0x1f) << 1)) - 1L;
        return result;
    }

    @Override
    public boolean equals( final Object obj ) {
        return obj instanceof SVKmerHuge && equals((SVKmerHuge)obj);
    }

    public boolean equals( final SVKmerHuge that ) {
        return Arrays.equals(this.vals, that.vals);
    }

    @Override
    public int hashCode() {
        if ( hashVal == 0 ) {
            long longHashVal = SVUtils.FNV64_DEFAULT_SEED;
            for ( final long val : vals ) {
                longHashVal = SVUtils.fnvLong64(hashVal, val);
            }
            hashVal = (int)(longHashVal ^ (longHashVal >>> 32));
        }
        return hashVal;
    }

    // results could be bogus if the kmers are of unequal sizes
    @Override
    public int compareTo( final SVKmerHuge that ) {
        int idx = Math.min(this.vals.length, that.vals.length);
        while ( idx-- > 0 ) {
            int cmp = Long.compare(this.vals[idx], that.vals[idx]);
            if ( cmp != 0 ) return cmp;
        }
        return Integer.compare(this.vals.length, that.vals.length);
    }

    @Override
    public String toString( final int kSize ) {
        final StringBuilder sb = new StringBuilder(kSize);
        int count = 0;
        for ( long val : vals ) {
            int nnn = 32;
            while ( nnn-- > 0 ) {
                sb.append(baseChars.charAt((int)val & 3));
                val >>>= 2;
                if ( ++count == kSize ) return sb.reverse().toString();
            }
        }
        throw new GATKException("kSize doesn't match array length");
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SVKmerHuge> {
        @Override
        public void write( final Kryo kryo, final Output output, final SVKmerHuge svKmer ) {
            svKmer.serialize(output);
        }

        @Override
        public SVKmerHuge read( final Kryo kryo, final Input input, final Class<SVKmerHuge> klass ) {
            return new SVKmerHuge(input);
        }
    }
}
