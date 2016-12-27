package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.GenomicsDBProgramGroup;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;

@CommandLineProgramProperties(
  summary = "Reads variants from GenomicsDB into Spark RDDs of multi-sample VariantContexs.",
  oneLineSummary = "Reads variants from GenomicsDB into Spark",
  programGroup = GenomicsDBProgramGroup.class
)

public final class ReadGenomicsDB {
}
