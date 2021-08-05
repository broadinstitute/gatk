package org.broadinstitute.hellbender.tools.funcotator;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorDataSourceBundlerUtils;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorDataSourceBundler;

/**
 * A unit test suite for the {@link FuncotatorDataSourceBundlerUtils} class.
 * Created by Hailey on 8/5/21.
 */
public class FuncotatorDataSourceBundlerUtilsUnitTest extends CommandLineProgramTest{

    //==================================================================================================================
    // Static Variables:
    private static final String TEST_WRONG_ORGNAME = "virus";
    private static final String TEST_WRONG_SPECIESNAME = "ashbya_gossypii";

}
