package org.broadinstitute.hellbender.testutils;

import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public final class ArgumentsBuilderTest{
    @Test
    public void testArgumentsBuilder() {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--value1");
        args.add(1);
        args.add("Input=")
            .add(new File("path/to/somewhere"))
            .add("Value2=2")
            .add(" Value3= 3 Value4=4");

        Assert.assertEquals(args.getArgsArray(), new String[]{"--value1", "1",
                "--Input", "path/to/somewhere","--Value2","2","--Value3","3","--Value4","4"});
    }
    @Test
    public void testOneBigString(){
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add(" Value=1 Value=2 Value=3 Value= 4 ");
        Assert.assertEquals(args.getArgsList().toArray(), new String[]{"--Value", "1", "--Value", "2",
                "--Value", "3","--Value", "4"});

    }

    @Test
    public void testFromArray(){
        ArgumentsBuilder args = new ArgumentsBuilder(new Object[]{"Option=" + new File("path/to"), "OtherOption=" + -1, new File("somewhere")});
        Assert.assertEquals(args.getArgsArray(), new String[]{"--Option","path/to", "--OtherOption", "-1", "somewhere"});
    }

    @Test
    public void testAddArgument() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument("foo", "bar");
        Assert.assertEquals(args.getArgsArray(), new String[]{"--foo","bar"});
        Assert.assertEquals(args.getString(), "--foo bar");
    }

    @Test
    public void testAddPositionalArgument() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addPositionalArgument("foo");
        args.addPositionalArgument("bar");
        Assert.assertEquals(args.getArgsArray(), new String[]{"foo", "bar"});
    }

    @Test
    public void testAddBooleanArgument() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addBooleanArgument("foo", true);
        args.addBooleanArgument("bar", false);
        Assert.assertEquals(args.getArgsArray(), new String[]{"--foo", "true", "--bar", "false"});
    }

    @Test
    public void testAddFileArgument() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();
        File file = new File("path/to/somewhere");
        args.addFileArgument("foo", file);
        Assert.assertEquals(args.getArgsArray(), new String[]{"--foo", file.getAbsolutePath()});
    }

    @Test
    public void testInput() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();
        File file = new File("path/to/somewhere");
        args.addInput(file);
        Assert.assertEquals(args.getArgsArray(), new String[]{"--" + StandardArgumentDefinitions.INPUT_LONG_NAME, file.getAbsolutePath()});
    }

    @Test
    public void testOutput() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();
        File file = new File("path/to/somewhere");
        args.addOutput(file);
        Assert.assertEquals(args.getArgsArray(), new String[]{"--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, file.getAbsolutePath()});
    }

    @Test
    public void testReference() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();
        File file = new File("path/to/somewhere");
        args.addReference(file);
        Assert.assertEquals(args.getArgsArray(), new String[]{"--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME, file.getAbsolutePath()});
    }

    @Test
    public void testVcf() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();
        File file = new File("path/to/somewhere");
        args.addVCF(file);
        Assert.assertEquals(args.getArgsArray(), new String[]{"--" + StandardArgumentDefinitions.VARIANT_LONG_NAME, file.getAbsolutePath()});
    }

}
