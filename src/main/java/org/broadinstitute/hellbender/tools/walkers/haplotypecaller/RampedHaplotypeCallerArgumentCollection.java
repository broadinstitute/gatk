package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;

public class RampedHaplotypeCallerArgumentCollection {

    public static final String RAMPS_DEBUG_READS_LONG_NAME = "ramps-debug-reads";
    public static final String RAMPS_DEBUG_POST_ASSEMBLER_ON_LONG_NAME = "ramps-debug-post-assembler-on";

    /**
     * The following definition related to an haplotype caller feature under development (experimental),
     * colloquially called 'ramps'. See additional information below, before the related parameters.
     */
    public final String OFF_RAMP_TYPE = "off-ramp-type";
    public final String OFF_RAMP_FILE = "off-ramp-file";
    public final String ON_RAMP_TYPE = "on-ramp-type";
    public final String ON_RAMP_FILE = "on-ramp-file";

    public enum OffRampTypeEnum {
        NONE,                           // no off ramp
        PRE_FILTER_OFF,                 // off ramp before the filtering step
        PRE_ASSEMBLER_OFF,              // off ramp before the assemnbler
        POST_ASSEMBLER_OFF,             // off ramp after the assembler
    };

    public enum OnRampTypeEnum {
        NONE,                           // no on ramp
        POST_FILTER_ON,                 // on ramp after the filter
        POST_ASSEMBLER_ON               // on tamp after the aseembler
    };

    /**
     * The following parameters related to an haplotype caller feature under development (experimental),
     * colloquially called 'ramps'. In its essence, it attempts to break the monilithic
     * haplotype calling process into a step-wise process, which can execute up until
     * a specific step (an off ramp) and later be restored to run from that point (an on ramp)).
     *
     * When running the ramped haplotype caller one normally specifies either an off ramp or an on ramp.
     *
     * Specifying an off ramp of a specific type results in the process halting at the step
     * associated with that ramp type and a state file (zip) being saved.
     *
     * Specifying an on ramp results in the process being restarted at the step associated
     * with the ramp type using a state file previously created using an off ramp.
     *
     * So summarize, a ramp is defined by its type and the name of a state file. Off ramps
     * create state files where on ramps read them.
     *
     * Ramp points (steps) are defined by the enums OnRampTypeEnum and OffRampTypeEnum.
     *
     */

    @Advanced
    @Hidden
    @Argument(fullName = OFF_RAMP_TYPE, doc = "ramps: Type of off-ramp, i.e. step in haplotype caller where the process should halt and a ramp state file be created", optional=true)
    public OffRampTypeEnum offRampType=null;

    @Advanced
    @Hidden
    @Argument(fullName = OFF_RAMP_FILE, doc = "ramps: File to use for writing the off-ramp", optional=true)
    public String offRampFile=null;

    @Advanced
    @Hidden
    @Argument(fullName = ON_RAMP_TYPE, doc = "ramp: Type of on-ramp, i.e. step in haplotype-caller where the process should be restored from using a given state file", optional=true)
    public OnRampTypeEnum onRampType=null;

    @Advanced
    @Hidden
    @Argument(fullName = ON_RAMP_FILE, doc = "ramps: File to use for reading on-ramp", optional=true)
    public String onRampFile=null;

    @Advanced
    @Hidden
    @Argument(fullName = RAMPS_DEBUG_READS_LONG_NAME, doc = "ramp: reads to print debug messages for", optional=true)
    public String rampsDebugReads =null;

    @Advanced
    @Hidden
    @Argument(fullName = RAMPS_DEBUG_POST_ASSEMBLER_ON_LONG_NAME, doc = "ramp: debug post assembler on ramp", optional=true)
    public boolean rampsDebugPostAssemblerOn =false;
}
