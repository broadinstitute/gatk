package org.broadinstitute.hellbender.utils.test;

import org.apache.commons.lang3.StringUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ArgumentsBuilderTest{
    @Test
    public void testArgumentsBuilder() {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--value1");
        args.add(1);
        args.add("Input=")
            .add(new File("/path/to/somewhere"))
            .add("Value2=2")
            .add(" Value3= 3 Value4=4");

        Assert.assertEquals(args.getArgsArray(), new String[]{"--value1", "1",
                "--Input", "/path/to/somewhere","--Value2","2","--Value3","3","--Value4","4"});
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
        ArgumentsBuilder args = new ArgumentsBuilder(new Object[]{"Option=" + new File("/path/to"), "OtherOption=" + -1});
        Assert.assertEquals(args.getArgsArray(), new String[]{"--Option","/path/to", "--OtherOption", "-1"});
    }


}
