package org.broadinstitute.hellbender.cmdline;

import java.io.File;
import java.util.ArrayList;
import java.util.List;


/**
 * This collection is for arguments that require special treatment by the arguments parser itself.
 * It should not grow beyond a very short list.
 */
public class SpecialArgumentsCollection implements ArgumentCollectionDefinition {
    public static final String HELP_FULLNAME = "help";
    public static final String VERSION_FULLNAME = "version";
    public static final String ARGUMENTS_FILE_FULLNAME = "arguments_file";

    @Option(shortName = "h", fullName = HELP_FULLNAME, doc= "display the help message", special = true)
    public boolean HELP = false;

    @Option(fullName = VERSION_FULLNAME, doc="display the version number for this tool", special = true)
    public boolean VERSION = false;

    @Option(fullName = ARGUMENTS_FILE_FULLNAME, doc="read one or more arguments files and add them to the command line", special = true)
    public List<File> ARGUMENTS_FILE = new ArrayList<>();

}
