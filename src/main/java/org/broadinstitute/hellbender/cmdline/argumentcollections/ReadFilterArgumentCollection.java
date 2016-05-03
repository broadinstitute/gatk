package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.engine.filters.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.List;
import java.util.ArrayList;
import java.util.function.BiFunction;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * {@link org.broadinstitute.hellbender.cmdline.ArgumentCollection} for allowing read filters to
 * be enabled/disabled on the command line.
 */
public class ReadFilterArgumentCollection implements ArgumentCollectionDefinition {
    private static final Logger logger = LogManager.getLogger(ReadFilterArgumentCollection.class);
    private static final long serialVersionUID = 1L;

    /**
     * Enum of read filters that can be enabled/disabled from the command line
     */
    public enum CommandLineReadFilter implements Serializable {
        /**
         * No-arg filters: these filters have no arguments and one static instance,
         * so the supplier function just returns the static instance and ignores
         * the hdr and args arguments
         */
        ALLOW_ALL                   ((hdr, args) -> ReadFilterLibrary.ALLOW_ALL_READS),
        CIGAR_IS_SUPPORTED          ((hdr, args) -> ReadFilterLibrary.CIGAR_IS_SUPPORTED),
        GOOD_CIGAR                  ((hdr, args) -> ReadFilterLibrary.GOOD_CIGAR),
        HAS_READ_GROUP              ((hdr, args) -> ReadFilterLibrary.HAS_READ_GROUP),
        HAS_MATCHING_BASES_AND_QUALS((hdr, args) -> ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS),
        MAPPED                      ((hdr, args) -> ReadFilterLibrary.MAPPED),
        MAPPING_QUALITY_AVAILABLE   ((hdr, args) -> ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE),
        MAPPING_QUALITY_NOT_ZERO    ((hdr, args) -> ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO),
        MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE((hdr, args) -> ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE),
        MATE_DIFFERENT_STRAND       ((hdr, args) -> ReadFilterLibrary.MATE_DIFFERENT_STRAND),
        NOT_DUPLICATE               ((hdr, args) -> ReadFilterLibrary.NOT_DUPLICATE),
        PRIMARY_ALIGNMENT           ((hdr, args) -> ReadFilterLibrary.PRIMARY_ALIGNMENT),
        PASSES_VENDOR_QUALITY_CHECK ((hdr, args) -> ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK),
        READLENGTH_EQUALS_CIGARLENGTH ((hdr, args) -> ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH),
        SEQ_IS_STORED               ((hdr, args) -> ReadFilterLibrary.SEQ_IS_STORED),
        VALID_ALIGNMENT_START       ((hdr, args) -> ReadFilterLibrary.VALID_ALIGNMENT_START),
        VALID_ALIGNMENT_END         ((hdr, args) -> ReadFilterLibrary.VALID_ALIGNMENT_END),

        /**
         * Read filters with user-defined arguments. The supplier function validates the
         * user provided values and returns the already-populated filter embedded in the
         * argument collection.
         *
         * For read filters with user-defined arguments that also require a SAMFileHeader,
         * the supplier function creates a new instance of the read filter using the
         * given header, populating the other args from the instance embedded in
         * the read filter argument collection; updates the embedded filter with the
         * new one and returns it.
         */
        ALIGNMENT_AGREES_WITH_HEADER((hdr, args) -> new AlignmentAgreesWithHeaderReadFilter(hdr)),
        FRAGMENT_LENGTH((hdr, args) -> {
            validateArgs(CommandLineReadFilter.valueOf("FRAGMENT_LENGTH"), args.fragmentLengthFilter::validate);
            return args.fragmentLengthFilter;
        }),
        LIBRARY((hdr, args) -> {
            validateArgs(CommandLineReadFilter.valueOf("LIBRARY"), args.libraryFilter::validate);
            LibraryReadFilter lrf = new LibraryReadFilter(hdr);
            lrf.libraryToKeep = args.libraryFilter.libraryToKeep;
            // update the actual argument collection filter just in case...
            args.libraryFilter = lrf;
            return lrf;
        }),
        MAPPING_QUALITY((hdr, args) -> {
            validateArgs(CommandLineReadFilter.valueOf("MAPPING_QUALITY"), args.mappingQualityFilter::validate);
            MappingQualityReadFilter mqf = new MappingQualityReadFilter(args.mappingQualityFilter.minMappingQualityScore);
            // update the actual argument collection filter just in case...
            args.mappingQualityFilter = mqf;
            return mqf;
        }),
        PLATFORM((hdr, args) -> {
            validateArgs(CommandLineReadFilter.valueOf("PLATFORM"), args.platformFilter::validate);
            PlatformReadFilter pf = new PlatformReadFilter(hdr);
            pf.PLFilterNames = args.platformFilter.PLFilterNames;
            // update the actual argument collection filter just in case...
            args.platformFilter = pf;
            return pf;
        }),
        PLATFORM_UNIT((hdr, args) -> {
            validateArgs(CommandLineReadFilter.valueOf("PLATFORM_UNIT"), args.platformUnitFilter::validate);
            PlatformUnitReadFilter pfuf = new PlatformUnitReadFilter(hdr);
            pfuf.blackListedLanes = args.platformUnitFilter.blackListedLanes;
            // update the actual argument collection filter just in case...
            args.platformUnitFilter = pfuf;
            return pfuf;
        }),
        READ_GROUP((hdr, args) -> {
            validateArgs(CommandLineReadFilter.valueOf("READ_GROUP"), args.readGroupFilter::validate);
            return args.readGroupFilter;
        }),
        READ_GROUP_BLACK_LISTED((hdr, args) -> {
            validateArgs(CommandLineReadFilter.valueOf("READ_GROUP_BLACK_LISTED"), args.readGroupBlackListFilter::validate);
            ReadGroupBlackListReadFilter rgblf = new ReadGroupBlackListReadFilter(
                    args.readGroupBlackListFilter.blackList, hdr);
            // update the actual argument collection filter just in case...
            args.readGroupBlackListFilter = rgblf;
            return rgblf;
        }),
        READ_NAME((hdr, args) -> {
            validateArgs(CommandLineReadFilter.valueOf("READ_NAME"), args.readNameFilter::validate);
            return args.readNameFilter;
        }),
        READ_LENGTH((hdr, args) -> {
            validateArgs(CommandLineReadFilter.valueOf("READ_LENGTH"), args.readLengthFilter::validate);
            return args.readLengthFilter;
        }),
        READ_STRAND((hdr, args) -> { return args.readStrandFilter; }),
        SAMPLE((hdr, args) -> {
            validateArgs(CommandLineReadFilter.valueOf("SAMPLE"), args.sampleFilter::validate);
            SampleReadFilter sf = new SampleReadFilter(hdr);
            sf.samplesToKeep = args.sampleFilter.samplesToKeep;
            // update the actual argument collection filter just in case...
            args.sampleFilter = sf;
            return sf;
        }),
        WELLFORMED((hdr, args) -> new WellformedReadFilter(hdr));

        /**
         * Each filter in the enum has an accompanying filter supplier function that accepts
         * a ReadFilterArgumentCollection object, extracts the filter's required arguments
         * from the argument collection, validates the arguments, and returns the initialized
         * read filter.
         */
        private BiFunction<SAMFileHeader, ReadFilterArgumentCollection, ReadFilter> filterSupplier = null;
        CommandLineReadFilter(final BiFunction<SAMFileHeader, ReadFilterArgumentCollection,ReadFilter> filterSupplier) {
            this.filterSupplier = filterSupplier;
        }

        /**
         * Execute the filter supplier and return the initialized filter
         *
         */
        public ReadFilter getFilterInstance(
                final SAMFileHeader samHeader,
                final ReadFilterArgumentCollection readFilterCollection)
        {
            return filterSupplier.apply(samHeader, readFilterCollection);
        };

        public static boolean validateArgs(
                final CommandLineReadFilter filter,
                final Supplier<String> validator)
        {
            String errorMessage = validator.get();
            if (null != errorMessage) {
                throw new UserException.CommandLineException("The " + filter.name() + " read filter " + errorMessage);
            }
            return true;
        }
    }

    @Argument(fullName = StandardArgumentDefinitions.READ_FILTER_LONG_NAME,
            shortName = StandardArgumentDefinitions.READ_FILTER_SHORT_NAME,
            doc = "Filters to apply to reads before analysis",
            common = true,
            optional = true)
    public List<CommandLineReadFilter> enabledCommandLineFilters = new ArrayList<>();

    @Argument(fullName = "disableReadFilter",
            shortName = "df",
            doc = "Filters to be disabled before analysis",
            common = true,
            optional = true)
    public List<CommandLineReadFilter> disabledCommandLineFilters = new ArrayList<>();

    @Argument(fullName = "disableAllReadFilters", shortName = "f", doc = "Disable all read filters", common = false, optional = true)
    public boolean disableAllReadFilters = false;

    //************************************************************
    // Argument collections for read filters that accept arguments

    @ArgumentCollection(doc="Read filter for fragment length", dependsOnArgument="RF", dependsOnValue="FRAGMENT_LENGTH")
    FragmentLengthReadFilter fragmentLengthFilter = new FragmentLengthReadFilter();

    @ArgumentCollection(doc="Read filter for library name", dependsOnArgument="RF", dependsOnValue="LIBRARY")
    LibraryReadFilter libraryFilter = new LibraryReadFilter(null);

    @ArgumentCollection(doc="Read filter for mapping quality", dependsOnArgument="RF", dependsOnValue="MAPPING_QUALITY")
    MappingQualityReadFilter mappingQualityFilter = new MappingQualityReadFilter(0);

    @ArgumentCollection(doc="Read filter for platform", dependsOnArgument="RF", dependsOnValue="PLATFORM")
    PlatformReadFilter platformFilter = new PlatformReadFilter(null);

    @ArgumentCollection(doc="Read filter for filtering lanes", dependsOnArgument="RF", dependsOnValue="PLATFORM_UNIT")
    PlatformUnitReadFilter platformUnitFilter = new PlatformUnitReadFilter(null);

    @ArgumentCollection(doc="Read filter for read name", dependsOnArgument="RF", dependsOnValue="READ_NAME")
    ReadNameReadFilter readNameFilter = new ReadNameReadFilter();

    @ArgumentCollection(doc="Read filter for read length", dependsOnArgument="RF", dependsOnValue="READ_LENGTH")
    ReadLengthReadFilter readLengthFilter = new ReadLengthReadFilter();

    @ArgumentCollection(doc="Read filter for read group", dependsOnArgument="RF", dependsOnValue="READ_GROUP")
    ReadGroupReadFilter readGroupFilter = new ReadGroupReadFilter();

    @ArgumentCollection(doc="Read filter for read group tag", dependsOnArgument="RF", dependsOnValue="READ_GROUP_BLACK_LISTED")
    ReadGroupBlackListReadFilter readGroupBlackListFilter = new ReadGroupBlackListReadFilter(new ArrayList<String>(), null);

    @ArgumentCollection(doc="Read filter for strand", dependsOnArgument="RF", dependsOnValue="READ_STRAND")
    ReadStrandFilter readStrandFilter = new ReadStrandFilter();

    @ArgumentCollection(doc="Read filter for sample name", dependsOnArgument="RF", dependsOnValue="SAMPLE")
    SampleReadFilter sampleFilter = new SampleReadFilter(null);

    /**
     * Returns the list of read filters that result from combining the provided default
     * read filters with command line arguments for enabling/disabling individual
     * filters.
     *
     * @param defaultFilters List of default command line read filters to be merged
     *                       with command line arguments. May not be null.
     * @return List of CommandLineReadFilters, resulting from combining the default
     * filters with all filters enabled/disabled on the command line.
     */
    public List<CommandLineReadFilter> getMergedCommandLineFilters(final List<CommandLineReadFilter> defaultFilters) {
        Utils.nonNull(defaultFilters, "A default read filter must be provided");

        // warn if any filters that were requested to be disabled are not enabled
        disabledCommandLineFilters.stream()
                .filter(f -> !defaultFilters.contains(f))
                .forEach(f -> logger.warn("Disabled filter: " + f.name() + " is already disabled"));

        final List<CommandLineReadFilter> curatedFilters = new ArrayList<>();
        if (!disableAllReadFilters) {
            // add in the defaults first (order matters), removing any that were explicitly
            // disabled on the command line, followed by the ones that were explicitly enabled
            // on the command line
            curatedFilters.addAll(defaultFilters.stream().filter(
                    f -> !disabledCommandLineFilters.contains(f)).collect(Collectors.toList())
            );
            // if you disabled AND enabled, enabled wins
            curatedFilters.addAll(enabledCommandLineFilters);
        }

        return curatedFilters;
    }

}
