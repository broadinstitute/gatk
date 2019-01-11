package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

public class ReadTransformingIteratorUnitTest extends GATKBaseTest {

    private static final List<GATKRead> READS = Arrays.asList(
            ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10S101M10S")),
            ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("50M")),
            ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M5I10M")),
            ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("30M20S")),
            ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("101M"))
    );

    private static final List<Integer> READ_CIGAR_SIZES = Arrays.asList(3,1,3,2,1);



    private static List<GATKRead> transformedReads(final ReadTransformer readTransformer) {
        // note that this is a deep copy -- transformation is done in place and we don't want to modify {@code READS}
        final Iterator<GATKRead> iterator = READS.stream().map(GATKRead::deepCopy).collect(Collectors.toList()).iterator();
        final ReadTransformingIterator readTransformingIterator = new ReadTransformingIterator(iterator, readTransformer);
        return Utils.stream(readTransformingIterator.iterator()).collect(Collectors.toList());
    }

    @Test
    public void testIdentityTransformation() {
        final List<GATKRead> transformedReads = transformedReads(ReadTransformer.identity());
        Assert.assertEquals(transformedReads.size(), READS.size());
        for (int n = 0; n < READS.size(); n++) {
            Assert.assertEquals(transformedReads.get(n).getSAMString(), READS.get(n).getSAMString());
        }
    }

    @Test
    public void testModifyingTransformation() {
        final String cigarString = "10S11I54M92S";

        final ReadTransformer transformer = read -> {
            read.setCigar(cigarString);
            return read;
        };

        final List<GATKRead> transformedReads = transformedReads(transformer);
        Assert.assertEquals(transformedReads.size(), READS.size());
        for (int n = 0; n < READS.size(); n++) {
            Assert.assertEquals(transformedReads.get(n).getCigar().toString(), cigarString);
        }
    }

    @Test
    public void testAnnotatingTransformation() {
        final String cigarLengthAttributeName = "CL";

        final ReadTransformer transformer = read -> {
            read.setAttribute(cigarLengthAttributeName, read.getCigarElements().size());
            return read;
        };

        final List<GATKRead> transformedReads = transformedReads(transformer);
        Assert.assertEquals(transformedReads.size(), READS.size());
        for (int n = 0; n < READS.size(); n++) {
            Assert.assertEquals(transformedReads.get(n).getAttributeAsInteger(cigarLengthAttributeName), READ_CIGAR_SIZES.get(n));
        }
    }
}