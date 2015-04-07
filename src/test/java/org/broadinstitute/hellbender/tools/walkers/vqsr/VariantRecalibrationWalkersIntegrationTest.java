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
* Copyright 2012-2014 Broad Institute, Inc.
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

package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.variant.VCIterable;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class VariantRecalibrationWalkersIntegrationTest extends WalkerTest {
    private static class VRTest {
        String inVCF;
        String aggregateVCF;
        String tranchesMD5;
        String recalMD5;
        String cutVCFMD5;

        public VRTest(String inVCF, String tranchesMD5, String recalMD5, String cutVCFMD5) {
            this.inVCF = inVCF;
            this.tranchesMD5 = tranchesMD5;
            this.recalMD5 = recalMD5;
            this.cutVCFMD5 = cutVCFMD5;
        }

        public VRTest(String inVCF, String aggregateVCF, String tranchesMD5, String recalMD5, String cutVCFMD5) {
            this.inVCF = inVCF;
            this.aggregateVCF = aggregateVCF;
            this.tranchesMD5 = tranchesMD5;
            this.recalMD5 = recalMD5;
            this.cutVCFMD5 = cutVCFMD5;
        }

        @Override
        public String toString() {
            return "VRTest{inVCF='" + inVCF +"'}";
        }
    }

    VRTest lowPass = new VRTest(validationDataLocation + "phase1.projectConsensus.chr20.raw.snps.vcf",
            "41e2d951a17de433fe378bb3d9ec75d4",  // tranches
            "04336b2453202f286da05b69e57f66ed",  // recal file
            "d29fd0bdc1c8c3a171e10d29f7ffeaec"); // cut VCF

    VRTest lowPassPlusExomes = new VRTest(validationDataLocation + "phase1.projectConsensus.chr20.raw.snps.vcf",
            validationDataLocation + "1kg_exomes_unfiltered.AFR.unfiltered.vcf",
            "ce4bfc6619147fe7ce1f8331bbeb86ce",  // tranches
            "1b33c10be7d8bf8e9accd11113835262",  // recal file
            "4700d52a06f2ef3a5882719b86911e51"); // cut VCF

    @DataProvider(name = "VRTest")
    public Object[][] createData1() {
        return new Object[][]{ {lowPass} };
    }

    @DataProvider(name = "VRAggregateTest")
    public Object[][] createData2() {
        return new Object[][]{ {lowPassPlusExomes} };
    }

    @Test(dataProvider = "VRTest")
    public void testVariantRecalibrator(VRTest params) {
        //System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -resource:known=true,prior=10.0 " + GATKDataLocation + "dbsnp_132_b37.leftAligned.vcf" +
                        " -resource:truth=true,training=true,prior=15.0 " + comparisonDataLocation + "Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf" +
                        " -resource:training=true,truth=true,prior=12.0 " + comparisonDataLocation + "Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf" +
                        " -T VariantRecalibrator" +
                        " -input " + params.inVCF +
                        " -L 20:1,000,000-40,000,000" +
                        " --no_cmdline_in_header" +
                        " -an QD -an HaplotypeScore -an HRun" +
                        " --trustAllPolymorphic" + // for speed
                        " -recalFile %s" +
                        " -mode SNP" +
                        " -tranchesFile %s",
                Arrays.asList(params.recalMD5, params.tranchesMD5));
        final List<File> outputFiles = executeTest("testVariantRecalibrator-"+params.inVCF, spec).getFirst();
        setPDFsForDeletion(outputFiles);
    }

    @Test(dataProvider = "VRTest",dependsOnMethods="testVariantRecalibrator")
    public void testApplyRecalibration(VRTest params) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -T ApplyRecalibration" +
                        " -L 20:12,000,000-30,000,000" +
                        " --no_cmdline_in_header" +
                        " -input " + params.inVCF +
                        " -U LENIENT_VCF_PROCESSING -o %s" +
                        " -mode SNP" +
                        " -tranchesFile " + getMd5DB().getMD5FilePath(params.tranchesMD5, null) +
                        " -recalFile " + getMd5DB().getMD5FilePath(params.recalMD5, null),
                Arrays.asList(params.cutVCFMD5));
        spec.disableShadowBCF(); // TODO -- enable when we support symbolic alleles
        final List<File> outputFiles = executeTest("testApplyRecalibration-"+params.inVCF, spec).getFirst();
        setPDFsForDeletion(outputFiles);
    }

    @Test(dataProvider = "VRAggregateTest")
    public void testVariantRecalibratorAggregate(VRTest params) {
        //System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -resource:known=true,prior=10.0 " + GATKDataLocation + "dbsnp_132_b37.leftAligned.vcf" +
                        " -resource:truth=true,training=true,prior=15.0 " + comparisonDataLocation + "Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf" +
                        " -resource:training=true,truth=true,prior=12.0 " + comparisonDataLocation + "Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf" +
                        " -T VariantRecalibrator" +
                        " -input " + params.inVCF +
                        " -aggregate " + params.aggregateVCF +
                        " -L 20:1,000,000-40,000,000" +
                        " --no_cmdline_in_header" +
                        " -an QD -an HaplotypeScore -an MQ" +
                        " --trustAllPolymorphic" + // for speed
                        " -recalFile %s" +
                        " -mode SNP" +
                        " -tranchesFile %s",
                Arrays.asList(params.recalMD5, params.tranchesMD5));
        final List<File> outputFiles = executeTest("testVariantRecalibratorAggregate-"+params.inVCF, spec).getFirst();
        setPDFsForDeletion(outputFiles);
    }

    @Test(dataProvider = "VRAggregateTest",dependsOnMethods="testVariantRecalibratorAggregate")
    public void testApplyRecalibrationAggregate(VRTest params) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -T ApplyRecalibration" +
                        " -L 20:12,000,000-30,000,000" +
                        " --no_cmdline_in_header" +
                        " -input " + params.inVCF +
                        " -U LENIENT_VCF_PROCESSING -o %s" +
                        " -mode SNP" +
                        " -tranchesFile " + getMd5DB().getMD5FilePath(params.tranchesMD5, null) +
                        " -recalFile " + getMd5DB().getMD5FilePath(params.recalMD5, null),
                Arrays.asList(params.cutVCFMD5));
        spec.disableShadowBCF(); // TODO -- enable when we support symbolic alleles
        final List<File> outputFiles = executeTest("testApplyRecalibrationAggregate-"+params.inVCF, spec).getFirst();
        setPDFsForDeletion(outputFiles);
    }

    VRTest bcfTest = new VRTest(privateTestDir + "vqsr.bcf_test.snps.unfiltered.bcf",
            "3ad7f55fb3b072f373cbce0b32b66df4",  // tranches
            "e747c08131d58d9a4800720f6ca80e0c",  // recal file
            "e5808af3af0f2611ba5a3d172ab2557b"); // cut VCF

    @DataProvider(name = "VRBCFTest")
    public Object[][] createVRBCFTest() {
        return new Object[][]{ {bcfTest} };
    }

    @Test(dataProvider = "VRBCFTest")
    public void testVariantRecalibratorWithBCF(VRTest params) {
        //System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -resource:known=true,prior=10.0 " + GATKDataLocation + "dbsnp_132_b37.leftAligned.vcf" +
                        " -resource:truth=true,training=true,prior=15.0 " + comparisonDataLocation + "Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf" +
                        " -resource:training=true,truth=true,prior=12.0 " + comparisonDataLocation + "Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf" +
                        " -T VariantRecalibrator" +
                        " -input " + params.inVCF +
                        " -L 20:10,000,000-20,000,000" +
                        " --no_cmdline_in_header" +
                        " -an AC " + // integer value
                        " -an QD -an ReadPosRankSum -an FS -an InbreedingCoeff " + // floats value
                        " -mG 2 "+
                        " -recalFile %s" +
                        " -mode SNP" +
                        " -tranchesFile %s",
                2,
                Arrays.asList("bcf", "txt"),
                Arrays.asList(params.recalMD5, params.tranchesMD5));
        final List<File> outputFiles = executeTest("testVariantRecalibrator-"+params.inVCF, spec).getFirst();
        setPDFsForDeletion(outputFiles);
    }

    @Test(dataProvider = "VRBCFTest", dependsOnMethods="testVariantRecalibratorWithBCF")
    public void testApplyRecalibrationWithBCF(VRTest params) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -T ApplyRecalibration" +
                        " -L 20:10,000,000-20,000,000" +
                        " --no_cmdline_in_header" +
                        " -input " + params.inVCF +
                        " -U LENIENT_VCF_PROCESSING -o %s" +
                        " -mode SNP" +
                        " -tranchesFile " + getMd5DB().getMD5FilePath(params.tranchesMD5, null) +
                        " -recalFile " + getMd5DB().getMD5FilePath(params.recalMD5, null),
                Arrays.asList(params.cutVCFMD5));
        spec.disableShadowBCF();
        final List<File> outputFiles = executeTest("testApplyRecalibration-"+params.inVCF, spec).getFirst();
        setPDFsForDeletion(outputFiles);
    }


    VRTest indelUnfiltered = new VRTest(
            validationDataLocation + "combined.phase1.chr20.raw.indels.unfiltered.sites.vcf", // all FILTERs as .
            "9a331328370889168a7aa3a625f73620",  // tranches
            "2cbbd146d68c40200b782e0226f71976",  // recal file
            "64dd98a5ab80cf5fd9a36eb66b38268e"); // cut VCF

    VRTest indelFiltered = new VRTest(
            validationDataLocation + "combined.phase1.chr20.raw.indels.filtered.sites.vcf", // all FILTERs as PASS
            "9a331328370889168a7aa3a625f73620",  // tranches
            "2cbbd146d68c40200b782e0226f71976",  // recal file
            "c0ec662001e829f5779a9d13b1d77d80"); // cut VCF

    @DataProvider(name = "VRIndelTest")
    public Object[][] createTestVariantRecalibratorIndel() {
        return new Object[][]{ {indelUnfiltered}, {indelFiltered} };
    }

    @Test(dataProvider = "VRIndelTest")
    public void testVariantRecalibratorIndel(VRTest params) {
        //System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -resource:known=true,prior=10.0 " + GATKDataLocation + "dbsnp_132_b37.leftAligned.vcf" +
                        " -resource:training=true,truth=true,prior=15.0 " + comparisonDataLocation + "Validated/Mills_Devine_Indels_2011/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.vcf" +
                        " -T VariantRecalibrator" +
                        " -input " + params.inVCF +
                        " -L 20:1,000,000-40,000,000" +
                        " --no_cmdline_in_header" +
                        " -an QD -an ReadPosRankSum -an HaplotypeScore" +
                        " -mode INDEL -mG 3" +
                        " --trustAllPolymorphic" + // for speed
                        " -recalFile %s" +
                        " -tranchesFile %s",
                Arrays.asList(params.recalMD5, params.tranchesMD5));
        final List<File> outputFiles = executeTest("testVariantRecalibratorIndel-"+params.inVCF, spec).getFirst();
        setPDFsForDeletion(outputFiles);
    }

    @Test(dataProvider = "VRIndelTest",dependsOnMethods="testVariantRecalibratorIndel")
    public void testApplyRecalibrationIndel(VRTest params) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -T ApplyRecalibration" +
                        " -L 20:12,000,000-30,000,000" +
                        " -mode INDEL" +
                        " -U LENIENT_VCF_PROCESSING --no_cmdline_in_header" +
                        " -input " + params.inVCF +
                        " -o %s" +
                        " -tranchesFile " + getMd5DB().getMD5FilePath(params.tranchesMD5, null) +
                        " -recalFile " + getMd5DB().getMD5FilePath(params.recalMD5, null),
                Arrays.asList(params.cutVCFMD5));
        spec.disableShadowBCF(); // has to be disabled because the input VCF is missing LowQual annotation
        final List<File> outputFiles = executeTest("testApplyRecalibrationIndel-" + params.inVCF, spec).getFirst();
        setPDFsForDeletion(outputFiles);
    }

    @Test
    public void testApplyRecalibrationSnpAndIndelTogether() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -T ApplyRecalibration" +
                        " -L 20:1000100-1000500" +
                        " -mode BOTH" +
                        " --no_cmdline_in_header" +
                        " -input " + privateTestDir + "VQSR.mixedTest.input" +
                        " -o %s" +
                        " -tranchesFile " + privateTestDir + "VQSR.mixedTest.tranches" +
                        " -recalFile " + privateTestDir + "VQSR.mixedTest.recal",
                Arrays.asList("03a0ed00af6aac76d39e569f90594a02"));
        final List<File> outputFiles = executeTest("testApplyRecalibrationSnpAndIndelTogether", spec).getFirst();
        setPDFsForDeletion(outputFiles);
    }

    @Test(enabled = true)
    public void testApplyRecalibrationSnpAndIndelTogetherExcludeFiltered() throws Exception {
        final String base = "-R " + b37KGReference +
                " -T ApplyRecalibration" +
                " -L 20:1000100-1000500" +
                " -mode BOTH" +
                " --excludeFiltered -ts_filter_level 90.0" +
                " --no_cmdline_in_header" +
                " -input " + privateTestDir + "VQSR.mixedTest.input" +
                " -o %s" +
                " -tranchesFile " + privateTestDir + "VQSR.mixedTest.tranches" +
                " -recalFile " + privateTestDir + "VQSR.mixedTest.recal";

        final WalkerTestSpec spec = new WalkerTestSpec(base, 1, Arrays.asList(""));
        spec.disableShadowBCF();
        final List<File> outputFiles = executeTest("testApplyRecalibrationSnpAndIndelTogether", spec).getFirst();
        setPDFsForDeletion(outputFiles);
        final File VCF = outputFiles.get(0);
        for( final VariantContext VC : VCIterable.readAllVCs(VCF, new VCFCodec()).getSecond() ) {
            if( VC != null ) {
                Assert.assertTrue(VC.isNotFiltered()); // there should only be unfiltered records in the output VCF file
            }
        }
    }

    private void setPDFsForDeletion( final List<File> walkerOutputFiles ) {
        for ( final File outputFile : walkerOutputFiles ) {
            new File(outputFile.getAbsolutePath() + ".pdf").deleteOnExit();
        }
    }
}

