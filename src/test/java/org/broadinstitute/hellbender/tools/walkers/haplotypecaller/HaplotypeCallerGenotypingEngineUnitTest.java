/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2015 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.tools.walkers.haplotypecaller;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 3/15/12
 */

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.*;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.haplotype.EventMap;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.smithwaterman.Parameters;
import org.broadinstitute.gatk.utils.smithwaterman.SWPairwiseAlignment;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Unit tests for {@link HaplotypeCallerGenotypingEngine}.
 */
public class HaplotypeCallerGenotypingEngineUnitTest extends BaseTest {

    private static ReferenceSequenceFile seq;
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        genomeLocParser = new GenomeLocParser(seq);
    }

    private class BasicGenotypingTestProvider extends TestDataProvider {
        byte[] ref;
        byte[] hap;
        Map<Integer,Byte> expected;

        public BasicGenotypingTestProvider(String refString, String hapString, Map<Integer, Byte> expected) {
            super(BasicGenotypingTestProvider.class, String.format("Haplotype to VCF test: ref = %s, alignment = %s", refString,hapString));
            ref = refString.getBytes();
            hap = hapString.getBytes();
            this.expected = expected;
        }
        
        public Map<Integer,VariantContext> calcAlignment() {
            final SWPairwiseAlignment alignment = new SWPairwiseAlignment(ref, hap, new Parameters(3,-1,-4, -1));
            final Haplotype h = new Haplotype(hap, false, alignment.getAlignmentStart2wrt1(), alignment.getCigar());
            return HaplotypeCallerGenotypingEngine.generateVCsFromAlignment(h, ref, genomeLocParser.createGenomeLoc("4", 1, 1 + ref.length), "name");
        }

        public String toString() {
            return "REF:" + new String(ref) + ",ALT:" + new String(hap);
        }
    }

    @DataProvider(name = "BasicGenotypingTestProvider")
    public Object[][] makeBasicGenotypingTests() {

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGACTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(2 + contextSize, (byte)'M');
            map.put(21 + contextSize, (byte)'M');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG", "ATCTCGCATCGCGAGCATCGCCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'I');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGACACTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'D');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1, (byte)'M');
            map.put(20, (byte)'D');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider("AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(2 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'I');
            map.put(30 + contextSize, (byte)'D');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "ACCTCGCATCGCGAGCATCGTTACTAGCCGATG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'D');
            map.put(28 + contextSize, (byte)'M');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCCATAG", map);
        }

        return BasicGenotypingTestProvider.getTests(BasicGenotypingTestProvider.class);
    }
    
    @Test(dataProvider = "BasicGenotypingTestProvider", enabled = true)
    public void testHaplotypeToVCF(BasicGenotypingTestProvider cfg) {
        Map<Integer,VariantContext> calculatedMap = cfg.calcAlignment();
        Map<Integer,Byte> expectedMap = cfg.expected;
        logger.warn(String.format("Test: %s", cfg.toString()));
        if(!compareVCMaps(calculatedMap, expectedMap)) {
            logger.warn("calc map = " + calculatedMap);
            logger.warn("expected map = " + expectedMap);
        }
        Assert.assertTrue(compareVCMaps(calculatedMap, expectedMap),"" + cfg);
    }

    @Test(dataProvider="AddMiscellaneousDataProvider", enabled=false)
    public void testAddMiscellaneousAllele(final String readBases, final int readOffset,
                                           final String ref, final int refOffset,
                                           final String referenceAllele, final String[] alternatives, final double[] likelihoods, final double[] expected) {
        final byte baseQual = (byte)30;

        final byte[] baseQuals = Utils.dupBytes(baseQual, readBases.length());
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(readBases.getBytes(), baseQuals, readBases.length() + "M");
        final GenomeLoc loc = new UnvalidatingGenomeLoc("20",0,refOffset,refOffset);
        final ReadBackedPileup pileup = new ReadBackedPileupImpl(loc,Collections.singletonList(read),readOffset);
        final VariantContextBuilder vcb = new VariantContextBuilder();
        final GenotypeBuilder gb = new GenotypeBuilder();
        final List<String> alleleStrings = new ArrayList<>( 1 + alternatives.length);
        alleleStrings.add(referenceAllele);
        alleleStrings.addAll(Arrays.asList(alternatives));

        gb.AD(new int[] { 1 });
        gb.DP(1);
        gb.PL(likelihoods);

        vcb.alleles(alleleStrings);
        vcb.loc("20",refOffset,refOffset + referenceAllele.length() -1);

        vcb.genotypes(gb.make());

        final VariantContext vc = vcb.make();

        final VariantContext updatedVc = null; // GenotypingEngine.addMiscellaneousAllele(vc,pileup,ref.getBytes(),0);
        final GenotypeLikelihoods updatedLikelihoods = updatedVc.getGenotype(0).getLikelihoods();
        Assert.assertEquals(updatedLikelihoods.getAsVector().length, expected.length);
        final double[] updatedLikelihoodsArray = updatedVc.getGenotype(0).getLikelihoods().getAsVector();
        for (int i = 0; i < updatedLikelihoodsArray.length; i++) {
            Assert.assertEquals(updatedLikelihoodsArray[i],expected[i],0.0001);
        }
        Allele altAllele = null;
        for (final Allele allele : updatedVc.getAlleles())
            if (allele.isSymbolic() && allele.getBaseString().equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE_NAME))
                altAllele = allele;
        Assert.assertNotNull(altAllele);
    }

    @DataProvider(name="AddMiscellaneousDataProvider")
    public Iterator<Object[]> addMiscellaneousAlleleDataProvider() {
        return Arrays.asList(ADD_MISCELLANEOUS_ALLELE_DATA).iterator();
    }

    private static final double MATCH_LnLK = QualityUtils.qualToProbLog10((byte)30);
    private static final double MISS_LnLK = QualityUtils.qualToErrorProbLog10((byte)30);

    private static final Object[][] ADD_MISCELLANEOUS_ALLELE_DATA = new Object[][] {
            new Object[] {"ACTG", 0,"ACTGTGAGTATTCC",0,"A",new String[]{}, new double[] {MATCH_LnLK * MATCH_LnLK}, 6,
                    new double[] {MATCH_LnLK * MATCH_LnLK,MATCH_LnLK * MISS_LnLK, MISS_LnLK * MISS_LnLK}}
    };

    /**
     * Private function to compare Map of VCs, it only checks the types and start locations of the VariantContext
     */
    private boolean compareVCMaps(Map<Integer, VariantContext> calc, Map<Integer, Byte> expected) {
        if( !calc.keySet().equals(expected.keySet()) ) { return false; } // sanity check
        for( Integer loc : expected.keySet() ) {
            Byte type = expected.get(loc);
            switch( type ) {
                case 'I':
                    if( !calc.get(loc).isSimpleInsertion() ) { return false; }
                    break;
                case 'D':
                    if( !calc.get(loc).isSimpleDeletion() ) { return false; }
                    break;
                case 'M':
                    if( !(calc.get(loc).isMNP() || calc.get(loc).isSNP()) ) { return false; }
                    break;
                default:
                    return false;
            }
        }
        return true;
    }


    @DataProvider(name = "CreateHaplotypeMappingProvider")
    public Object[][] makeCreateHaplotypeMappingData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final Set<Haplotype> haplotypes = new HashSet<>();
        final Allele ref = Allele.create("A", true);
        final Allele altC = Allele.create("C", false);
        final Allele altT = Allele.create("T", false);

        final Haplotype AtoC1 = new Haplotype("AACAA".getBytes());
        final VariantContext vc1 = new VariantContextBuilder().chr("20").start(3).stop(3).alleles(Arrays.asList(ref, altC)).make();
        AtoC1.setEventMap(new EventMap(Arrays.asList(vc1)));
        AtoC1.getEventMap().put(3, vc1);
        haplotypes.add(AtoC1);

        final Haplotype AtoC2 = new Haplotype("AAACA".getBytes());
        final VariantContext vc2 = new VariantContextBuilder().chr("20").start(4).stop(4).alleles(Arrays.asList(ref, altT)).make();
        AtoC2.setEventMap(new EventMap(Arrays.asList(vc2)));
        AtoC2.getEventMap().put(4, vc2);
        haplotypes.add(AtoC2);

        tests.add(new Object[]{vc1, haplotypes, AtoC1});
        tests.add(new Object[]{vc2, haplotypes, AtoC2});
        tests.add(new Object[]{new VariantContextBuilder().chr("20").start(1).stop(1).alleles(Arrays.asList(ref, altT)).make(), haplotypes, null});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider="CreateHaplotypeMappingProvider")
    public void testCreateHaplotypeMapping(final VariantContext vc, final Set<Haplotype> haplotypes, final Haplotype expected) {
        final Map<VariantContext, Set<Haplotype>> mapping = HaplotypeCallerGenotypingEngine.constructHaplotypeMapping(Arrays.asList(vc), haplotypes);
        final Set<Haplotype> actual = mapping.get(vc);
        if ( expected == null )
            Assert.assertTrue(actual.isEmpty(), actual.toString());
        else {
            Assert.assertEquals(actual.size(), 1);
            Assert.assertEquals(actual.iterator().next(), expected);
        }
    }

    @DataProvider(name = "ConstructPhaseSetMappingProvider")
    public Object[][] makeConstructPhaseSetMappingData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final Allele ref = Allele.create("A", true);
        final Allele altC = Allele.create("C", false);
        final Allele altT = Allele.create("T", false);

        final VariantContext vc1 = new VariantContextBuilder().chr("20").start(1).stop(1).alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc2 = new VariantContextBuilder().chr("20").start(2).stop(2).alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc3 = new VariantContextBuilder().chr("20").start(3).stop(3).alleles(Arrays.asList(ref, altT)).make();
        final VariantContext vc4 = new VariantContextBuilder().chr("20").start(4).stop(4).alleles(Arrays.asList(ref, altC)).make();
        final List<VariantContext> calls = Arrays.asList(vc2, vc3, vc4);

        final Haplotype pos1 = new Haplotype("CAAAA".getBytes());
        pos1.setEventMap(new EventMap(Arrays.asList(vc1)));
        pos1.getEventMap().put(1, vc1);
        final Haplotype pos2 = new Haplotype("ACAAA".getBytes());
        pos2.setEventMap(new EventMap(Arrays.asList(vc2)));
        pos2.getEventMap().put(2, vc2);
        final Haplotype pos3 = new Haplotype("AACAA".getBytes());
        pos3.setEventMap(new EventMap(Arrays.asList(vc3)));
        pos3.getEventMap().put(3, vc3);
        final Haplotype pos4 = new Haplotype("AAACA".getBytes());
        pos4.setEventMap(new EventMap(Arrays.asList(vc4)));
        pos4.getEventMap().put(4, vc4);
        final Haplotype pos24 = new Haplotype("ACACA".getBytes());
        pos24.setEventMap(new EventMap(Arrays.asList(vc2, vc4)));
        pos24.getEventMap().put(2, vc2);
        pos24.getEventMap().put(4, vc4);
        final Haplotype pos34 = new Haplotype("AACCA".getBytes());
        pos34.setEventMap(new EventMap(Arrays.asList(vc3, vc4)));
        pos34.getEventMap().put(3, vc3);
        pos34.getEventMap().put(4, vc4);
        final Haplotype pos234 = new Haplotype("ACCCA".getBytes());
        pos234.setEventMap(new EventMap(Arrays.asList(vc2, vc3, vc4)));
        pos234.getEventMap().put(2, vc2);
        pos234.getEventMap().put(3, vc3);
        pos234.getEventMap().put(4, vc4);

        final Map<VariantContext, Set<Haplotype>> haplotypeMap = new HashMap<>();

        // test no phased variants #1
        final Set<Haplotype> haplotypes2 = new HashSet<>();
        haplotypes2.add(pos2);
        haplotypeMap.put(vc2, haplotypes2);
        tests.add(new Object[]{Arrays.asList(vc2), new HashMap<>(haplotypeMap), 2, 0, 0, 0, 0});

        // test no phased variants #2
        final Set<Haplotype> haplotypes3 = new HashSet<>();
        haplotypes3.add(pos3);
        haplotypeMap.put(vc3, haplotypes3);
        tests.add(new Object[]{Arrays.asList(vc2, vc3), new HashMap<>(haplotypeMap), 3, 0, 0, 0, 0});

        // test opposite phase
        tests.add(new Object[]{Arrays.asList(vc2, vc3), new HashMap<>(haplotypeMap), 2, 2, 1, 1, 1});

        // test no phased variants #3
        final Set<Haplotype> haplotypes4 = new HashSet<>();
        haplotypes4.add(pos4);
        haplotypeMap.put(vc4, haplotypes4);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 3, 0, 0, 0, 0});

        // test mixture
        final Set<Haplotype> haplotypes24 = new HashSet<>();
        haplotypes24.add(pos24);
        haplotypeMap.put(vc2, haplotypes24);
        haplotypeMap.put(vc4, haplotypes24);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 2, 3, 1, 2, 1});

        // test 2 hets
        haplotypeMap.remove(vc3);
        tests.add(new Object[]{Arrays.asList(vc2, vc4), new HashMap<>(haplotypeMap), 1, 2, 1, 2, 0});

        // test 2 with opposite phase
        final Set<Haplotype> haplotypes1 = new HashSet<>();
        haplotypes1.add(pos1);
        haplotypeMap.put(vc1, haplotypes1);
        tests.add(new Object[]{Arrays.asList(vc1, vc2, vc4), new HashMap<>(haplotypeMap), 2, 3, 1, 1, 2});

        // test homs around a het
        final Set<Haplotype> haplotypes2hom = new HashSet<>();
        haplotypes2hom.add(pos24);
        haplotypes2hom.add(pos234);
        final Set<Haplotype> haplotypes4hom = new HashSet<>();
        haplotypes4hom.add(pos24);
        haplotypes4hom.add(pos234);
        final Set<Haplotype> haplotypes3het = new HashSet<>();
        haplotypes3het.add(pos234);
        haplotypeMap.put(vc2, haplotypes2hom);
        haplotypeMap.put(vc3, haplotypes3het);
        haplotypeMap.put(vc4, haplotypes4hom);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 2, 3, 1, 3, 0});

        // test hets around a hom
        final Set<Haplotype> haplotypes2het = new HashSet<>();
        haplotypes2het.add(pos234);
        final Set<Haplotype> haplotypes4het = new HashSet<>();
        haplotypes4het.add(pos234);
        final Set<Haplotype> haplotypes3hom = new HashSet<>();
        haplotypes3hom.add(pos3);
        haplotypes3hom.add(pos234);
        haplotypeMap.put(vc2, haplotypes2het);
        haplotypeMap.put(vc3, haplotypes3hom);
        haplotypeMap.put(vc4, haplotypes4het);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 2, 3, 1, 3, 0});

        // test no phased variants around a hom
        final Set<Haplotype> haplotypes2incomplete = new HashSet<>();
        haplotypes2incomplete.add(pos24);
        final Set<Haplotype> haplotypes3incomplete = new HashSet<>();
        haplotypes3incomplete.add(pos34);
        final Set<Haplotype> haplotypes4complete = new HashSet<>();
        haplotypes4complete.add(pos24);
        haplotypes4complete.add(pos34);
        haplotypes4complete.add(pos234);
        haplotypeMap.put(vc2, haplotypes2incomplete);
        haplotypeMap.put(vc3, haplotypes3incomplete);
        haplotypeMap.put(vc4, haplotypes4complete);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 0, 0, 0, 0, 0});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider="ConstructPhaseSetMappingProvider")
    public void testConstructPhaseSetMapping(final List<VariantContext> calls,
                                         final Map<VariantContext, Set<Haplotype>> haplotypeMap,
                                         final int totalHaplotypes,
                                         final int expectedMapSize,
                                         final int expectedNumGroups,
                                         final int expectedNum01,
                                         final int expectedNum10) {
        final Map<VariantContext, Pair<Integer, String>> actualPhaseSetMapping = new HashMap<>();
        final int actualNumGroups = HaplotypeCallerGenotypingEngine.constructPhaseSetMapping(calls, haplotypeMap, totalHaplotypes, actualPhaseSetMapping);
        Assert.assertEquals(actualNumGroups, expectedNumGroups);
        Assert.assertEquals(actualPhaseSetMapping.size(), expectedMapSize);

        int num01 = 0, num10 = 0;
        for ( final Pair<Integer, String> phase : actualPhaseSetMapping.values() ) {
            if ( phase.second.equals("0|1") )
                num01++;
            else if ( phase.second.equals("1|0") )
                num10++;
        }
        Assert.assertEquals(num01, expectedNum01);
        Assert.assertEquals(num10, expectedNum10);
    }

    @DataProvider(name = "ConstructPhaseGroupsProvider")
    public Object[][] makeConstructPhaseGroupsData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final Allele ref = Allele.create("A", true);
        final Allele altC = Allele.create("C", false);

        final Genotype g1 = new GenotypeBuilder().alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc1 = new VariantContextBuilder().chr("20").start(1).stop(1).alleles(Arrays.asList(ref, altC)).genotypes(g1).make();
        final Genotype g2 = new GenotypeBuilder().alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc2 = new VariantContextBuilder().chr("20").start(2).stop(2).alleles(Arrays.asList(ref, altC)).genotypes(g2).make();
        final Genotype g3 = new GenotypeBuilder().alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc3 = new VariantContextBuilder().chr("20").start(3).stop(3).alleles(Arrays.asList(ref, altC)).genotypes(g3).make();
        final List<VariantContext> calls = Arrays.asList(vc1, vc2, vc3);

        // test no phased variants, empty map
        final Map<VariantContext, Pair<Integer, String>> nonePhased1 = new HashMap<>();
        tests.add(new Object[]{calls, nonePhased1, 0, 0, 0});

        // test no phased variants, full map, exception expected
        final Map<VariantContext, Pair<Integer, String>> nonePhased2 = new HashMap<>();
        nonePhased2.put(vc1, new Pair<>(0, "0/1"));
        nonePhased2.put(vc2, new Pair<>(1, "0/1"));
        nonePhased2.put(vc3, new Pair<>(2, "0/1"));
        tests.add(new Object[]{calls, nonePhased2, 3, -1, -1});

        // test 2 phased variants
        final Map<VariantContext, Pair<Integer, String>> twoPhased = new HashMap<>();
        twoPhased.put(vc1, new Pair<>(0, "0/1"));
        twoPhased.put(vc2, new Pair<>(0, "0/1"));
        tests.add(new Object[]{calls, twoPhased, 1, 1, 2});

        // test all phased variants
        final Map<VariantContext, Pair<Integer, String>> allPhased = new HashMap<>();
        allPhased.put(vc1, new Pair<>(0, "0/1"));
        allPhased.put(vc2, new Pair<>(0, "0/1"));
        allPhased.put(vc3, new Pair<>(0, "0/1"));
        tests.add(new Object[]{calls, allPhased, 1, 1, 3});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider="ConstructPhaseGroupsProvider")
    public void testConstructPhaseGroups(final List<VariantContext> calls,
                                         final Map<VariantContext, Pair<Integer, String>> phaseMap,
                                         final int endIndex,
                                         final int expectedNumGroups,
                                         final int expectedGroupSize) {
        final List<VariantContext> actualPhasedCalls;
        try {
            actualPhasedCalls = HaplotypeCallerGenotypingEngine.constructPhaseGroups(calls, phaseMap, endIndex);
        } catch (IllegalStateException e) {
            Assert.assertEquals(-1, expectedNumGroups);
            return;
        }

        final Set<String> uniqueGroups = new HashSet<>();
        int counter = 0;
        for ( final VariantContext call : actualPhasedCalls ) {
            for ( final Genotype g : call.getGenotypes() ) {
                if ( g.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY) ) {
                    uniqueGroups.add(g.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY).toString());
                    counter++;
                }
            }
        }

        Assert.assertEquals(uniqueGroups.size(), expectedNumGroups);
        Assert.assertEquals(counter, expectedGroupSize);
    }
}
