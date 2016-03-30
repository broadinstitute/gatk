package org.broadinstitute.hellbender.tools.walkers.variantutils;

import com.google.cloud.dataflow.sdk.transforms.SerializableComparator;
import com.google.common.hash.BloomFilter;
import com.google.common.hash.Funnels;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.Hidden;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Encodes sites form the given variant file as a Bloom filter or a int[] or List<GATKVariant> and serializes",
        oneLineSummary = "Encodes and serializes sites form the given variant file",
        programGroup = VariantProgramGroup.class
)
@Hidden //hidden and unsupported
public final class SerializeVCFSites extends VariantWalker {
    static final Logger logger = LogManager.getLogger(SerializeVCFSites.class);

    /*
     * BloomFilters need to know the number of elements a priori so we
     * store data in the bitset first and then convert to BloomFilters after the traversal
     *
     * This map stores them sorted by contig name.
     */
    private SortedMap<String, BitSet> bitSets;

    @Argument(shortName = "fpp", fullName = "desired_false_positive_probability", optional = true)
    private double fpp = 0.001;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "File name prefix to which serialized encodings should be written (defaults to null which does not serialize)",
            optional = true)
    private String outPrefix = null;

    private List<GATKVariant> variantList;

    private SerializableComparator<String> byContigIndex(final SAMSequenceDictionary dict){
        return (contig1, contig2) -> Integer.compare(dict.getSequenceIndex(contig1), dict.getSequenceIndex(contig2));
    }

    @Override
    public void onTraversalStart() {
        final SAMSequenceDictionary dict = getSequenceDictionary();

        bitSets = makeBitSets(dict);

        variantList = new ArrayList<>(1000000);//pre allocate a bunch of space
    }

    private final SortedMap<String, BitSet> makeBitSets(final SAMSequenceDictionary dict) {
        final SortedMap<String, BitSet> bitSets = new TreeMap<>(byContigIndex(dict));
        for(final SAMSequenceRecord seq : dict.getSequences()){
            bitSets.put(seq.getSequenceName(), new BitSet(seq.getSequenceLength()+1));
        }
        return bitSets;
    }

    private SAMSequenceDictionary getSequenceDictionary() {
        return getHeaderForVariants().getSequenceDictionary();
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        bitSets.get(variant.getContig()).set(variant.getStart(), variant.getEnd()+1);
        variantList.add(VariantContextVariantAdapter.sparkVariantAdapter(variant));
    }

    @Override
    public Object onTraversalSuccess() {
        final long genomeLength = getSequenceDictionary().getSequences().stream().mapToLong(s->s.getSequenceLength()).sum();
        logger.info("Total genome length: " + genomeLength);
        final long variantLength = bitSets.values().stream().mapToLong(bs -> bs.cardinality()).sum();
        logger.info("Total variant length: " + variantLength + " "+ ((double)variantLength/genomeLength) + " of genome");

        final SortedMap<String, BloomFilter<Integer>> filters = makeFilters();
        final SortedMap<String, int[]> positionsArrays = makePositionArrays();

        if (outPrefix == null){
            logger.info("Not serializing results");
        } else {
            serialize((Serializable) filters, outPrefix + ".bloom_fpp" + fpp);
            serialize((Serializable) positionsArrays, outPrefix + ".positions");
            serialize((Serializable) this.variantList, outPrefix + ".gatvariants");
        }
        return bitSets.values().stream().mapToInt(BitSet::cardinality).sum();
    }

    //Contig -> sorted array of ints
    private SortedMap<String, int[]> makePositionArrays() {
        final SortedMap<String, int[]> arrays = new TreeMap<>(byContigIndex(getSequenceDictionary()));
        for (final Map.Entry<String, BitSet> kv : bitSets.entrySet()){
            final String contig = kv.getKey();
            final BitSet bs = kv.getValue();
            final int cardinality = bs.cardinality();
            final int[] positions = new int[cardinality];
            if (cardinality > 0) {
                logger.info("Position arrays adding " + cardinality + " from contig:" + contig);
            }

            int pos = 0;//current position in the arrat
            for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i+1)) {
                positions[pos++] = i;
            }
            arrays.put(contig, positions);
        }
        return arrays;
    }

    private SortedMap<String, BloomFilter<Integer>> makeFilters() {
        final SortedMap<String, BloomFilter<Integer>> filters = new TreeMap<>(byContigIndex(getSequenceDictionary()));
        for (final Map.Entry<String, BitSet> kv : bitSets.entrySet()){
            final String contig = kv.getKey();
            final BitSet bs = kv.getValue();
            final int cardinality = bs.cardinality();
            final BloomFilter<Integer> filter = BloomFilter.create(Funnels.integerFunnel(), cardinality, fpp);
            if (cardinality > 0) {
                logger.info("Bloom filter adding " + cardinality + " from contig:" + contig);
            }

            for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i+1)) {
              filter.put(i);
            }
            filters.put(contig, filter);
        }
        return filters;
    }

    private void serialize(final Serializable result, final String outFile) {
        logger.info("Serializing to " + outFile);

        try(final FileOutputStream fout = new FileOutputStream(outFile);
            final ObjectOutputStream oos = new ObjectOutputStream(fout)){
            oos.writeObject(result);
        } catch (final IOException e) {
            throw new GATKException("Problem while serializing the results", e);
        }
    }
}
