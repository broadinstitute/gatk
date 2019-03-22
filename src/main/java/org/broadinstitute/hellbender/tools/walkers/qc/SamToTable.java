package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.tribble.Feature;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 blah
 */

@CommandLineProgramProperties(
        summary = org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator.USAGE_SUMMARY,
        oneLineSummary = org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator.USAGE_ONE_LINE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
public final class SamToTable extends ReadWalker {
    public static final String USAGE_ONE_LINE_SUMMARY = "converts a samfile into a TSV";
    public static final String USAGE_SUMMARY = "SAM -> TSV";

    public static final String KNOWN_SITES_ARG_FULL_NAME = "known-sites";

    protected static final Logger logger = LogManager.getLogger(org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator.class);


    @Argument(fullName = KNOWN_SITES_ARG_FULL_NAME, doc = "One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.", optional = false)
    private List<FeatureInput<Feature>> knownSites;

    /**
     * After the header, data records occur one per line until the end of the file. The first several items on a line are the
     * values of the individual covariates and will change depending on which covariates were specified at runtime. The last
     * three items are the data- that is, number of observations for this combination of covariates, number of reference mismatches,
     * and the raw empirical quality score calculated by phred-scaling the mismatch rate.   Use '/dev/stdout' to print to standard out.
     */
    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, doc = "The output table file to create", optional = false)
    private File tableOutput = null;


    @Argument()
    private String tileRegEx = ".*";

    @Argument()
    private String readgroupRegEx = ".*";

    private PrintStream tableOutputStream;

    private ReferenceDataSource referenceDataSource; // datasource for the reference. We're using a different one from the engine itself to avoid messing with its caches.

    final private List<ReadElementExtractor> extractors = getExtractors();

    public boolean requiresReference() {
        return true;
    }

    /**
     * Parse the -cov arguments and create a list of covariates to be used here
     * Based on the covariates' estimates for initial capacity allocate the data hashmap
     */
    @Override
    public void onTraversalStart() {
        Utils.warnOnNonIlluminaReadGroups(getHeaderForReads(), logger);
        referenceDataSource = ReferenceDataSource.of(referenceArguments.getReferencePath());

        try {
            tableOutputStream = new PrintStream(tableOutput);
        } catch (FileNotFoundException e) {
            throw new GATKException("problem trying to open stream against " + tableOutput, e);
        }

        tableOutputStream.println(extractors.stream().map(ReadElementExtractor::header).collect(Collectors.joining()));

    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        List<ReadFilter> filters = new ArrayList<>(6);
        filters.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_SUPPLEMENTARY_ALIGNMENT);
        return filters;
    }

    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     * whether or not the base matches the reference at this particular location
     */
    @Override
    public void apply(GATKRead read, ReferenceContext ref, FeatureContext featureContext ) {
        final String collect = extractors.stream().map(e -> e.extractElement(read)).collect(Collectors.joining());
        tableOutputStream.println(collect);

    }
    private List<ReadElementExtractor> getExtractors() {

        return Arrays.asList(
                new ReadElementExtractorImpl.Tile(),
                new ReadElementExtractorImpl.XCoord(),
                new ReadElementExtractorImpl.YCoord(),
                new ReadElementExtractorImpl.BaseQual(),
                new ReadElementExtractorImpl.Duplicate(),
                new ReadElementExtractorImpl.Errors(),
                new ReadElementExtractorImpl.InsertSize(),
                new ReadElementExtractorImpl.Length(),
                new ReadElementExtractorImpl.Length2(),
                new ReadElementExtractorImpl.Mapped(),
                new ReadElementExtractorImpl.MappingQ()
        );
    }


    @Override
    public Object onTraversalSuccess() {
        tableOutputStream.close();
        return 0;
    }



}