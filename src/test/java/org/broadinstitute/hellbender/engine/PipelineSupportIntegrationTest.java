//package org.broadinstitute.hellbender.engine;
//
//import org.broadinstitute.hellbender.GATKBaseTest;
//import org.broadinstitute.hellbender.utils.runtime.ProcessController;
//import org.broadinstitute.hellbender.testutils.BaseTest;
//import org.testng.Assert;
//import org.testng.annotations.Test;
//
//import java.io.File;
//
//public class PipelineSupportIntegrationTest extends GATKBaseTest {
//
//    private static final String INPUT_BAM = packageRootTestDir + "engine/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10000000-10000020.with.unmapped.bam";
//
//    @Test
//    // Asserting that Picard programs are able to pipe their output without errors
//    public void testPipeForPicardTools() {
//        File output = createTempFile("testOutput",".bam");
//        String gatkLaunchScript = System.getenv("GATK_LAUNCH_SCRIPT");
//        String gatkLaunchScriptSubsitution = gatkLaunchScript == null ? "./gatk" : gatkLaunchScript;
//        BaseTest.runProcess(ProcessController.getThreadLocal(), new String[]{"/bin/sh","-c"," "+gatkLaunchScriptSubsitution+" SortSam -I "+INPUT_BAM+" -O /dev/stdout -SO coordinate " +
//                "| "+gatkLaunchScriptSubsitution+" SetNmMdAndUqTags -I /dev/stdin -O "+ output.getAbsolutePath()+ " --CREATE_INDEX true -R "+b37_reference_20_21});//.split(" "));
//        try (ReadsDataSource inputReads = new ReadsPathDataSource(new File(INPUT_BAM).toPath());
//             ReadsDataSource outputReads = new ReadsPathDataSource(output.toPath())) {
//            Assert.assertTrue(inputReads.iterator().hasNext());
//            Assert.assertTrue(outputReads.iterator().hasNext());
//            final int[] count = {0};
//            inputReads.forEach(r -> {
//                count[0]++;});
//            outputReads.forEach(r -> {
//                count[0]--;});
//            Assert.assertEquals(count[0],0);
//        }
//    }
//}