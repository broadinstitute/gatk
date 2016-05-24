package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class AssemblyRegionTestDataSetUnitTest extends BaseTest {

    @Test(dataProvider="activeRegionTestDataSets")
    @SuppressWarnings("fallthrough")
    public void testAssemblyRegionsDataSet(final AssemblyRegionTestDataSet as, final int kmerSize, final int readLength, final String variation, final int readCount, final int regionSize, final byte bq, final byte iq, final byte dq) {
        Assert.assertNotNull(as);
        Assert.assertEquals(as.assemblyResultSet().getMaximumKmerSize(), kmerSize);
        final List<GATKRead> reads = as.readList();
        Assert.assertEquals(reads.size(), readCount);
        for (final GATKRead r : reads) {
            Assert.assertEquals(r.getLength(), readLength);
        }

        final List<Haplotype> haplotypes = as.haplotypeList();
        final List<Civar> haplotypeCivars = Civar.fromCharSequence(variation).optionalizeAll().unroll();

        Assert.assertEquals(haplotypes.size(), haplotypeCivars.size());
        Assert.assertTrue(haplotypeCivars.size() > 1);
        int variants = 0;
        for (int i = 0; i < variation.length(); i++) {
           final char c = variation.charAt(i);
           switch (c) {
               case 'W':
               case 'T':
               case 'C':
               case 'D':
               case 'I':
                    variants++;
               default:

           }
        }

        Assert.assertEquals(haplotypes.size(), (int) Math.pow(2, variants));

        final Map<String,Integer> haplotypeNumberByString = new HashMap<>();
        for (int i = 0; i < haplotypes.size(); i++) {
            final Haplotype hap = haplotypes.get(i);
            final Civar civar = haplotypeCivars.get(i);
            Assert.assertEquals(hap.getBaseString(), civar.applyTo(as.getReference()));
            if (i == 0) {
                Assert.assertEquals(hap.getBaseString(), as.getReference());
            } else {
                Assert.assertNotEquals(hap.getBaseString(), as.getReference());
            }
            Assert.assertFalse(haplotypeNumberByString.containsKey(hap.getBaseString()));
            haplotypeNumberByString.put(hap.getBaseString(), i);
        }

        final int[] hapReadsNotInReference = new int[haplotypes.size()];

        for (int i = 0; i < readCount; i++) {
            final GATKRead r = as.readList().get(i);

            final int hapNumber = i % haplotypes.size();
            final int offset = i % (haplotypes.get(hapNumber).length() - readLength);
            Assert.assertEquals(getReadString(r), haplotypes.get(hapNumber).getBaseString().substring(offset, offset + readLength));
            if (!as.getReference().contains(getReadString(r))) {
                hapReadsNotInReference[hapNumber]++;
            }
        }

        Assert.assertEquals(hapReadsNotInReference[0], 0);

        for (int i = 1; i < hapReadsNotInReference.length; i++) {
            Assert.assertNotEquals(hapReadsNotInReference[i], 0);
        }

    }

    private static String getReadString(final GATKRead read){
        final byte[] readBases = read.getBases();
        if (readBases.length == 0) {
            return SAMRecord.NULL_SEQUENCE_STRING;
        }
        return new String(readBases);
    }

    /**
     * Constructs a test data-set based on the given parameters.
     * @param kmerSize length of the kmer.
     * @param readLength length of the read.
     * @param variation variation in that active region.
     * @param readCount number of reads in the active region
     * @param regionSize Active region size (~ size of the haplotype(s))
     * @param bq Base quality value common for all base-calls.
     * @param iq Insertion quality based for all read positions.
     * @param dq Deletion quality based for all read positions.
     * @return never null.
     */
    public static AssemblyRegionTestDataSet createAssemblyRegionTestDataSet(final int kmerSize, final int readLength, final String variation, final int readCount, final int regionSize, final byte bq, final byte iq, final byte dq) {

        final String reference = REF.substring(0, regionSize);

        final AssemblyRegionTestDataSet result = new AssemblyRegionTestDataSet(kmerSize, reference,
                new String[]{"Civar:" + variation},
                new String[]{"*:" + readCount + ":" + readLength}, byteRepeat(bq, readLength), byteRepeat(dq, readLength), byteRepeat(iq, readLength));


        return result;
    }

    @DataProvider(name="activeRegionTestDataSets")
    public Iterator<Object[]> activeRegionTestDataSets() {
        return new java.util.Iterator<Object[]>() {

            private int i = 0;

            @Override
            public boolean hasNext() {
                return i < ACTIVE_REGION_TEST_DATA_SET_PARAMETERS.length;
            }

            @Override
            public Object[] next() {

                final Object[] params = ACTIVE_REGION_TEST_DATA_SET_PARAMETERS[i++];
                final int kmerSize = (Integer) params[0];
                final int readLength = (Integer) params[1];
                final String variation = (String) params[2];
                final int readCount = (Integer) params[3];
                final int regionSize = (Integer) params[4];
                final AssemblyRegionTestDataSet dataSet = createAssemblyRegionTestDataSet(kmerSize, readLength, variation, readCount, regionSize, (byte) 20, (byte) 35, (byte) 35);
                return new Object[] { dataSet , kmerSize, readLength, variation, readCount, regionSize, (byte)20, (byte) 35, (byte) 35};
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    private static int[] KMER_SIZES = new int[] { 10 };
    private static int[] READ_COUNTS = new int[] { 1000 };
    private static int[] READ_LENGTHS = new int[] { 100 };
    private static String[] VARIATION_CIVARS = new String[] {
            "*1T*",
            "*3Iacg*",
            "*30Igctcggatgccttgcggggctccagagtcc*",
            "*3D*",
            "*30D*",
            "*1T3=3Iacg*",
            "*1T*3Iacg*",
            "*1T8=1T8=1T8=1T8=1T*",
            "*1T*1T*1T*1T*1T*"
     };

    private static int[] REGION_SIZE = new int[] { 300 };


    private static int[] intValues(final String[] kmerSizes) {
        final int[] result = new int[kmerSizes.length];
        for (int i = 0; i < result.length; i++)
            result[i] = Integer.parseInt(kmerSizes[i]);
        return result;
    }

    private static final Object[][] ACTIVE_REGION_TEST_DATA_SET_PARAMETERS;

    static {
        final int totalLength = KMER_SIZES.length * READ_COUNTS.length * READ_LENGTHS.length * VARIATION_CIVARS.length * REGION_SIZE.length;
        ACTIVE_REGION_TEST_DATA_SET_PARAMETERS = new Object[totalLength][];
        int next = 0;
        for (final int ks : KMER_SIZES)
            for (final int rc : READ_COUNTS)
                for (final int rl : READ_LENGTHS)
                    for (final String v : VARIATION_CIVARS)
                        for (final int rs : REGION_SIZE)
                            ACTIVE_REGION_TEST_DATA_SET_PARAMETERS[next++] = new Object[] { ks, rl, v, rc, rs };

    }

    private static byte[] byteRepeat(final byte bq, final int readLength) {
        final byte[] result = new byte[readLength];
        Arrays.fill(result, bq);
        return result;
    }



    private static final String REF =
            "TCGAGAAATTTGTATCCCGCCCCCGCAGCTTGCCAGCTCTTTCAGTATCATGGAGCCCAT" +
                    "GGTTGAATGAGTCCAATAACGAACTTCGACATGATAAAATCCCCCCCTCGCGACTTCCAG" +
                    "AGAAGAAGACTACTGACTTGAGCGTTCCCAGCACTTCAGCCAAGGAAGTTACCAATTTTT" +
                    "TGTTTCCGAATGACACGCGTCTCCTTGCGGGTAGATCGCCGACCGCAGAACTTACGAGCC" +
                    "AGGGGAAACAGTAAGGCCTAATTAGGTAAAGGGAGTAAGTGCTCGAACGCTTCAGATGTA" +
                    "ACCATATACTTACGCTGGATCTTCTCCCGCGAATTTTAACCCTCACCAACTACGAGATTT" +
                    "GAGGTAAACCAAATAAGCACGTAGTGGCGCTATCCGACTGTTCCCAAATTGTAACTTATC" +
                    "GTTCCGTGAAGGCCAGAGTTACTTCCCGGCCCTTTCCATGCGCGCACCATACCCTCCTAG" +
                    "TTCCCCGGTTATCTCTCCGAGGAGGGAGTGAGCGATCCTCCGTTTACGTTTTGTTACCAA" +
                    "TGACGTAGCTATGTATTTTGTACAGGTTGCCAACGGGTTTCACAATTCACAGATAGTGGG" +
                    "GATCCCGGCAAAGGGCCTATATTTGCGGTCCAACTTAGGCGTAAACTACGATGGTACCTA" +
                    "CTCAGACCCAGCTCGCGCGGCGTAAATAACGCACTCATCCCAGCTGATTCTCGGCGATCT" +
                    "ACGCAGCGACATGATTATCAACAGCTGTCTGGCAGCTCTAATCTTTTACCATGGTCGTAA" +
                    "AAGCCTCCAAGAGTTAGATCATACCTAACGCCACAAAAGTGACACGACGCCGATGGGTAC" +
                    "CGGACTTTAGGTCGACCACAGTTCGGTAAGGGAGAGGCCCTGCGGCGTACTTCATTTTGT" +
                    "ATATGCAACGTGCCCAAGTGGCGCCAGGCAAGTCTCAGCTGGTTCCTGTGTTAGCTCGAG" +
                    "GCTAGGCATGGGAGCTGATTGAACATGGGTTGGGGGCCTCGAACCGTCGAGGACCCCATA" +
                    "GTACCTCGGACACCAAGTAGGGCAGCCTATAGTTTGAAGCAGTACTATTTCAGGGGGGGA" +
                    "GCCCTCATGGTCTCTTCTACTGATGACTCAACACGCTAGGGACGTGAAGTCGATTCCTTC" +
                    "GATGGTTATAAATCAAAGGCTCAGAGTGCAGTCTGGAGCGCCCATCTAACGGTACGCATC" +
                    "TCGATTGCTCGGTCGCCTTTCACACTCCGCGAAAATTCATACCGCTCATTCACTAGGTTG" +
                    "CGAAGCCTACACTGATATATGAATCCAAGCTAGAGCAGGGCTCTTAAAATTCGGAGTTGT" +
                    "AGATGCTCAATACTCCAATCGGTTTTTTCGTGCACCACCGCGGGTGGCTGACAAGGGTTT" +
                    "GACATCGAGAAACAAGGCAGTTCCGGGCTGAAAGTAGCGCCGGGTAAGGTACGCGCCTGG" +
                    "TATGGCAGGACTATGAAGCCAATACAAAGGCTACATCCTCACTCGGGTGGACGGAAACGC" +
                    "AGAATTATGGTTACTTTTTGGATACGTGAAACATGTCCCATGGTAGCCCAAAGACTTGGG" +
                    "AGTCTATCACCCCTAGGGCCCATTTCTGGATATAGACGCCAGGTTGAATCCGTATTTGGA" +
                    "GGTACGATGGATCAGTCTGGGTGGGACGTGCTCCATTTATACCCTGCGCAGGCTGGACCG" +
                    "AGGACCGCAAGATGCGACGGTGCACAAGTAATTGACAACAAACCATCGTGTTTTCATTAT" +
                    "GGTACCAGGATCTTCAAGCCGAGTCAATCAAGCTCGGATTACAGTGTTTACCGCGTCTTG" +
                    "CGGTTACTCACAAACTGTAATCCACCACAAGTCAAGCCATTGCCTCTCTGAGACGCCGTA" +
                    "TGAATTAATATGTAAACTTTGCGCGGGTTCACTGCGATCCGTTCAGTCTCGTCCAAGGGC" +
                    "ACAATCGAATTCCCATTTGTATGTTCGGCTAACTTCTACCCATCCCCCGAAGTTTAGCAG" +
                    "GTCGTGAGGTGTCATGGAGGCTCTCGTTCATCCCGTGGGACATCAAGCTTCGCCTTGATA" +
                    "AAGCACCCCGCTCGGGTGTAGCAGAGAAGACGCCTACTGAATTGTGCGATCCCTCCACCT" +
                    "CAGCTAAGGTAGCTACCAATATTTAGTTTTTTAGCCTTGCGACAGACCTCCTACTTAGAT" +
                    "TGCCACGCATTGAGCTAGCGAGTCAGCGATAAGCATGACGCGCTTTCAAGCGTCGCGAGT" +
                    "ATGTGAACCAAGGCTCCGGACAGGACTATATACTTGGGTTTGATCTCGCCCCGACAACTG" +
                    "CAAACCTCAACATTTATAGATTATAAGGTTAGCCGAAATTGCACGTGGTGGCGCCCGCCG" +
                    "ACTGCTCCCCGAGTGTGGCTCTTTGATCTGACAACGCGCGACCTCCATCGCGGCCGATTG" +
                    "TTTCTGCGGACCATGTCGTCCTCATAGTTTGGGCATGTTTCCGTTGTAGGAGTGAAGCCA" +
                    "CTTAGCTTTGCGCCGTAGTCCCAATGAAAAACCTATGGACTTTGTTTTGGGTAGCATCAG" +
                    "GAATCTGAACCCTGTGAATGTGGGGGTCGCGCGCATAGACCTTTATCTCCGGTTCAAGTT" +
                    "AGGCATGAGGCTGCATGCTACGTTGTCACACCTACACTGCTCGAAGTAAATATGGGAAGC" +
                    "GCGCGGCCTGGCCCGAGGCGTTCCGCGCCGCCACGTGTTCGTTAACTGTTGATTGGTGGC" +
                    "ACATAAGCAATACCGTAGTCCCTCAAATTCAGCTCTGTTATCTCGAGCGTTATGTGTCAA" +
                    "ATGGCGTAGAACGGGATTGACTGTTTGACACTAGCTGGTGTTCGGTTCGGTAACGGAGAA" +
                    "TCTGTGGGGCTATGTCACTAATACTTTCGAAACGCCCCGTACCGATGCTGAACAAGTCGA" +
                    "TGCAGGCTCCCGTCTTTGAATAGGGGTAAACATACAAGTCGATAGAAGATGGGTAGGGGC" +
                    "CTCCAATTCATCCAACACTCTACGCCTTCTCCAAGAGCTAGTAGGGCACCCTGCAGTTGG" +
                    "AAAGGGAACTATTTCGTAGGGCGAGCCCATACCGTCTCTCTTGCGGAAGACTTAACACGA" +
                    "TAGGAAGCTGGAATAGTTTCGAACGATGGTTATTAATCCTAATAACGGAACGCTGTCTGG" +
                    "AGGATGAGTGTGACGGAGTGTAACTCGATGAGTTACCCGCTAATCGAACTGGGCGAGAGA" +
                    "TCCCAGCGCTGATGCACTCGATCCCGAGGCCTGACCCGACATATCAGCTCAGACTAGAGC" +
                    "GGGGCTGTTGACGTTTGGGGTTGAAAAAATCTATTGTACCAATCGGCTTCAACGTGCTCC" +
                    "ACGGCTGGCGCCTGAGGAGGGGCCCACACCGAGGAAGTAGACTGTTGCACGTTGGCGATG" +
                    "GCGGTAGCTAACTAAGTCGCCTGCCACAACAACAGTATCAAAGCCGTATAAAGGGAACAT" +
                    "CCACACTTTAGTGAATCGAAGCGCGGCATCAGAATTTCCTTTTGGATACCTGATACAAAG" +
                    "CCCATCGTGGTCCTTAGACTTCGTGCACATACAGCTGCACCGCACGCATGTGGAATTAGA" +
                    "GGCGAAGTACGATTCCTAGACCGACGTACGATACAACTATGTGGATGTGACGAGCTTCTT" +
                    "TTATATGCTTCGCCCGCCGGACCGGCCTCGCGATGGCGTAGCTGCGCATAAGCAAATGAC" +
                    "AATTAACCACTGTGTACTCGTTATAACATCTGGCAGTTAAAGTCGGGAGAATAGGAGCCG" +
                    "CAATACACAGTTTACCGCATCTAGACCTAACTGAGATACTGCCATAGACGACTAGCCATC" +
                    "CCTCTGGCTCTTAGATAGCCCGATACAGTGATTTTGAAAGGTTTGCGGGGCACAGCTATG" +
                    "ACTTGCTTAGCTACGTGTGAGGGAAGGAACTTTTGCGTATTTGTATGTTCACCCGTCTAC" +
                    "TACGCATGCGGGCAGATTATGTAGGTTGAGAGATGCGGGAGAAGTTCTCGACCTTCCCGT" +
                    "GGGACGTGAACCTATCCCCTAATAGAGCATTCCGTTCGAGCATGGCAGTAAGTACGCCTT" +
                    "CTCAATTGTGCTAACCTTCATCCCTATCAAAGCTTGGAGCCAATGATCAGGGTTATTCCC" +
                    "TTGGGACAGACTTCCTACTCACAGTCGGTCACATTGGGCTACTCCATGGGTCTTCAGCTT" +
                    "GACCCGGTCTGTTGGGCCGCGATTACGTGAGTTAGGGCCCCGGACTGCGCTGTATAGTCG" +
                    "ATTCTCATCCGGCCCCCACATCTGGAAACCCCAACTTATTTAGATAACATGATTAGCCGA" +
                    "AGTTGCACGGCGTGTCCACCGTGGAGTCCTCCCCGGGTGTCCCTCCTTCATTTGACGATA" +
                    "AGCAGCGGCTACCACCATTGATTAACACAAGGAACGGTGATGTTAACATAGATTCGGCAC" +
                    "ATTACTCTTGTAGGTGTGGAATCACTTAGCTACGCGGCGAAGCCTTATGGCAAAACCGAT" +
                    "GGGCAATGATTCGGGTAGCGCTAAAAGTCCATAGCACGTGCATCCCAACGTGGCGTGCGT" +
                    "ACAGCTTGACCACCGCTTCACGCTAAGGTGCTGGCCACATGCTAAATTGATGCGCCTGCA" +
                    "CTGCTCAAAGGATAATTACGAAGCGGGCGGCCTGGCGGGAGCACTACCCCATCGACGCGT" +
                    "ACTCGAATACTGTTTATTGCTCACACATGAACAAATTAGTAGAGTGCCACTTTCAGCCCT" +
                    "CTTGTCGTCGGCGATGTGTGTAAAATGGCGTTGATGTGGATCGACTCTATAAAGGTATCT" +
                    "ACTGATGCGTAGGGAGATCCGGAATCTATTGGCCTATGTCACTGAAACTATCCAAACACC" +
                    "CCATGTCGATACTGAACGTATCGACGCATACCTCCTTCCTTGAAAACGCACAATCATACA" +
                    "ACTGGGCACATAATGCGTACGCCCATCTAGTACACCCATCTCTGTGGGTCCAGTTCAAGA" +
                    "GCTGGAAGAGCACCCTCCACAAGGTCAAGTGGTATCCTGGTAAGGTAAGCTCGTACCGTG" +
                    "ATTCATGCGACAGGGGTAAGACCATCAGTAGTAGGGATAGTGCCAAACCTCACTCACCAC" +
                    "TGCCAATAAGGGGTCCTTACCTGAAGAATAAGTGTCAGCCAGTGTAACCCGATGAGGAAC" +
                    "CCAAAAGGCGAACCGGGCCAGACAACCCGGCGGTATCGCACTCAAAGCCGGGACACGACG" +
                    "CGTCACAGCCGGTAAGAGTAACCCCGGAGTGAAGACCTATGGGGCTGGATAAAACTGCCG" +
                    "TGGTAACCGCCTTCAACAACCCGAATACGTGGCACTTCAGGAGGCGCCCGGAGGGGGGAT" +
                    "GTTTTCTACTATTCGAGGCCGTTCGTTATAACTAGTTGCGTTCCTAGCCGCTATAATTGT" +
                    "CTCTTTGCCGACTAATGAGAACAACCACACCATAGCGATTTGACGCGGCGCCTCGGAATA" +
                    "CCGTTTCAGCAGGCGCTTGGTAAGGCCATCGCGAATACCAGGTATCGTGTAAGTAGCGTA" +
                    "GGCCCGCACGCAAGATAAACTGCTAGGGAACCGCGTTTCCACGACCGGTGCACGATTTAA" +
                    "TTTCGCCGACGTGATGACATTCCAGGCAGTGCCTCTGCCGCCGGACCCCTCTCGTGATTG" +
                    "GGTAGCTGGACATGCCCTTGTAAGATATAACAAGAGCCTGCCTGTCTAATGATCTCACGG" +
                    "CGAAAGTCGGGGAGACAGCAGCGGCTGCAGACATTATACCGCAACAACACTAAGGTGAGA" +
                    "TAACTCCGTAATTGACTACGCGTTCCTCTAGACCTTACTTGACCGGATACAGTGTCTTTG" +
                    "ACACGTTTATGGGTTACAGCAATCACATCCAAGACTGGCTATGCACGAAGCAACTCTTGA" +
                    "GTGTTAAAATGTTGACCCCTGTATTTGGGATGCGGGTAGTAGATGAGTGCAGGGACTCCG" +
                    "AGGTCAAGTACATTACCCTCTCATAGGGGGCGTTCTAGATCACGTTACCACCATATCATT" +
                    "CGAGCATGACACCATCTCCGCTGTGCCCATCCTAGTAGTCATTATTCCTATCACGCTTTC" +
                    "GAGTGTCTGGTGGCGGATATCCCCCACGAATGAAAATGTTTTTCGCTGACAGTCATATTG" +
                    "GGGTGCTCCTAAGCTTTTCCACTTGGCTGGGTCAGCTAGGCCTCCGTGCCCGGAGTTTCG" +
                    "GCGCAGTGCTGCCGACAGCCGGCCATTGTCTTTGGGGCCTCATTCGAGGGTACCCGGACC" +
                    "TATCTTGTCGGGACCACCCGGGGTAGTCGTTGGGCTTATGCACCGAAAAGCCCTGCGCCG" +
                    "GCCTCCCCGCTACGGAAGGTGATAAGCTCCGGCAAGCAATTATGAACAACGCAAGGATCG" +
                    "CGGATATAAACAGAGAAACGGCTGATTACACCTGTTCGTGTGGTATCGGTAAATAGCCTC" +
                    "GCGGAGCCTTATGCCATACTCGTCCGCGGAGCACTCTGGTAATGCATATGGTCCACAGGA" +
                    "CATTCGTCGCTTCCGGGTATGCGCTCTATTTGACGGTCCTTTGGCGCACAGATGCTGGCC" +
                    "ACCATTTAAATTAGAGCGACTCCACATCTGTAAGGTCCGCCACGCAGACGACAGCCCAGG" +
                    "GAGACCACTGACCGATCTACCTGAACGGCAACCTTCTGTATCGTACTGGGGCGGAGAGAT" +
                    "AACTACAGTGCCGCTTACAGCCCCTCTGTCGTCGCCGACGTCTGTAGTCTAGCCTCATTA" +
                    "TGATTGCACGCTATTGAGGCATTGACTGATGCCGGAAGACATCTGAAATGAACTGGTCTA" +
                    "TGCGACAGAAACCGTGCACCTACCAAATCTCCTTAGTGTAGGTTCTGACCGATTCGTGCT" +
                    "TCGTTGAGAACTCACATTTTAACAACAGAGGACATATGCCCTACCTCCATGATCTACTGA" +
                    "CGTCCCTGAGGCTGCAATTCATGTAATGGGGCAGTATCCGCGGCAAGTCCTAGTGCAATG" +
                    "GCGGTTTTTTACCCTCGTTCTGAAGAAGAGGCGACGCGGGTGCGGTCATCACTAATGTGG" +
                    "AAATTGGGAAGACTCTCGGGCCTCCGCCTTTAGGCGGTGCTTACTCTTTCATAAAGGGGC" +
                    "TGTTAGTTATGCCCCGCGAGGATTCGAAAAGGTGAGCCAACTCGGCCGATCCGGAGAGAC" +
                    "GGGCTTCAAAGCTGCCTGACGACGGTTGTGGGCCCGTAACAAAATCCTCCCAATAAGCCC" +
                    "CCGTGAGCGTCGGTTGAACAGCCCTGGTCGGCCCGACCAGAAGCCCGAATATATCGCTTT" +
                    "ACGGCTCTTGGGCCGGGGTGCGTTACCTTGCAGAAATCGAGGCCGTCCGTTAATTCCTGT" +
                    "TGCATTCATACCGCGTATATTTGTCTCTTTACCCGCTTACTTGGATAAGCATGGCATAGC" +
                    "TTTTTATCGGAGCGCCTCCGTACACGGTACGATCGCACGCCTCGTGAGATCAATACGTAT" +
                    "ACCAGGTGTCCTGTGAGCAGCGAAAGCCTAAACGGGAGATACACCGCCAAAAGTCCGTGT" +
                    "GAATACGAGTCGTGGCAAATTTGGTCTGGCTGTGATCTAGATATTCCAGGCGGTACGTCT" +
                    "GCTCTCGCGTGCCTCTAGTGGCTCGCTAGATAGTCTAGCCGCTGGTAAACACTCCATGAC" +
                    "CCCGGCTCTCCATTGATGCCACGGCGATTGTTGGAGAGCCAGCAGCGACTGCAAACGTCA" +
                    "GATCAGAGTAATACTAGCAAGCGATAAGTCCCTAACTGGTTGTGGCCTTCTGTAGAGTGA" +
                    "ACTTCACCACATATGCTGTCTCTGGCACGTGGATGGTTTGGAGAAATCAGATTCAAGTCT" +
                    "GATCAACCTTCAAACAGATCTAGAGTCTAAAACAGTGATCTCCTGCGTGCGAGATAGAAA" +
                    "TACTAGGTAACTACAGGGACTGCGACGTTTTAAACGTTGGTCCGTCAGAAGCGCCATTCA" +
                    "GGATCACGTTACCCCGAAAAAAAGGTACCAGGAGCTCTTCTCCTCTGCAGTCAGGTCTAT" +
                    "AGAAACTACACCATTAACCTTCCTGAGAACCGGGAGGTGGGAATCCGTCACATATGAGAA" +
                    "GGTATTTGCCCGATAATCAATACTCCAGGCTTCTAACTTTTTCCACTCGCTTGAGCCGGC" +
                    "TTGGCCTTTCTGCCTGAAGATTCGTTGGACTGGTGCCAACGCGCAGGCATAGTTCCAGGA" +
                    "GAATTATCCGGGGGCAGTGACAACCAACATCTCGGGTCTTGCCCAACCGGTCTACACGCT" +
                    "GATATAGCGAATCACCGAGAACCCGGCGCCACGCAATGGAACGTCCTTAACTCTGGCAGG" +
                    "CAATTAAAGGGAACGTATATATAACGCAAAAAAACTGGAAAATTGGCGAGAGAATCTTCT" +
                    "CTCTGTCTATCGAAGAATGGCCACGCGGAGGCATGCGTCATGCTAGCGTGCGGGGTACTC" +
                    "TTGCTATCCATTTGGGTCACAGGACACTCGCTGTTTTCGAATTTACCCTTTATGCGCCGG" +
                    "TATTGAACCACGCTTATGCCCAGCATCGTTACAACCAGACTGATACTAGATGTATAATGT" +
                    "CCGCCATGCAGACGAAACCAGTCGGAGATTACCGAGCATTCTATCACGTCGGCGACCACT" +
                    "AGTGAGCTACTGGAGCCGAGGGGTAACGATGATGCCCCTAAGAACCTCTCGGTCGACGCA" +
                    "AGCGATTACACTCCTGTCACATCATAATCGTTTGCTATTCAGGGGTTGACCAACACCGGA" +
                    "TAGCTTTTCACTTGAAGTATTATGCACGACAGGGTGCGTGTACCAACTAAACCTGTTTTA" +
                    "ACTTACCTCAGACTAGTTGGAAGTGTGGCTAGATCTTAGCTTTCGTCACTAGAGGGCCCA" +
                    "CGCTTAGTTTTTATGATCCATTGATCTCCTAGACGCTGCAAGATTTGCAACCAGGCAGAC" +
                    "TTAGCGGTAGGTCCTAGTGCAGCGGGACTTTTTTTCTATAGTCGTTGAGAGGAGGAGTCG" +
                    "TCAGACCAGATACCTTTGATGTCCTGATTGGAAGGACCGTTGGCCCCCGACCCTTAGACA" +
                    "GTGTACTCAGTTCTATAAACGAGCTATTAGATATGAGATCCGTAGATTGAAAAGGGTGAC" +
                    "GGAATTCGCCCGGACGCAAAAGACGGACAGCTAGGTATCCTGAGCACGGTTGCGCGTCCG" +
                    "AATCAAGCTCCTCTTTACAGGCCCCGGTTTCTGTTGGTCGTAGAGCGCAGAACGGATTGG" +
                    "GGGGATGTACGACAATATCTCTTAGTCACCTTTGGGTCACGGTCTGCTACCTTACAGGAA" +
                    "TTCAGACCGTCCTTTAATTTCCCTTGCATATATGTTGCGTTTCTTCGACCTTCTAACCGC" +
                    "ACCCTTAGGACGAAGACAGATACGTTCTTACCCATACTCCACCGTTGGCAGCGGGATCGC" +
                    "ATGTCCCACGTGAAACATTGCTAAACCCTCAGGCCTCTGAGCGACAAAAGCTTTAAAGGG" +
                    "AAATTCGCGCCCATAACTTGGTCCGAATACGGGTTCTAGCATCGTTCGTCTGAGTTTGTT" +
                    "CTATATAAAACGGGCGCAATGTCTGCTTTGATCAACCTCCAATACCTCGTATGATTGTGC" +
                    "ACCCGCCGGTGACCACTCAATGATGTGGGGTCCCCGTTGCAACTACGAGGATTTATTGAG" +
                    "ACCGACCTACGTTCGGCATTGTGGGCAGAGTGAAGTATTGGCAAACGTTAAGTGCCGAAC" +
                    "TAGATCTGACCTAACGGTAAGAGAGTTTCATAATACGTCCAGCCGCATGCGCAGGGTACA" +
                    "TTTGGACAGTATTGAATGGACTCTGATCAACCTTCACACCGATCTAGAAACGAGTGCGTA" +
                    "GATCAGCCAGGTGCAAACCAAAAATTCTAGGTTACTAGAAGTTTTGCGACGTTCTAAGAG" +
                    "TTGGACGAAATGTTTCGCGACCTAGGATGAGGTCGCCCTAGAAAATAGATTTCTGCTACT" +
                    "CTCCTCATAAGCAGTCCGGTGTATCGAAAGTACAAGACTAGCCTTGCTAGCAACCGCGGG" +
                    "CTGGGAGCCTAAGGCATCACTCAAGATACAGGCTCGGTAACGTACGCTCTAGCCATCTAA" +
                    "CTATCCCCTATGTCTTATAGGGACCTACGTTATCTGCCTG";

    protected final static String REF_MD5 = Utils.calcMD5(REF);

}