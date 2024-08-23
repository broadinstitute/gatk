package org.broadinstitute.hellbender.tools.sv.stratify;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

// Groups variants by SVTYPE, SVLEN, and overlap with one or more interval sets
public class SVStatificationEngine {

    // Configuration table column names
    public static final String NAME_COLUMN = "NAME";
    public static final String SVTYPE_COLUMN = "SVTYPE";
    public static final String MIN_SIZE_COLUMN = "MIN_SIZE";
    public static final String MAX_SIZE_COLUMN = "MAX_SIZE";
    public static final String CONTEXT_COLUMN = "CONTEXT";
    protected static final Set<String> COLUMN_NAMES = ImmutableSet.of(NAME_COLUMN, SVTYPE_COLUMN, MIN_SIZE_COLUMN, MAX_SIZE_COLUMN, CONTEXT_COLUMN);
    public static final String CONTEXT_COLUMN_DELIMITER = ",";

    public static final Set<String> NULL_TABLE_VALUES = Set.of("-1", "", "NULL", "NA");

    protected final Map<String, OverlapDetector<Locatable>> contextMap;
    protected final Map<String, Stratum> strata;
    protected final SAMSequenceDictionary dictionary;

    public SVStatificationEngine(final SAMSequenceDictionary dictionary) {
        contextMap = new HashMap<>();
        strata = new HashMap<>();
        this.dictionary = Utils.nonNull(dictionary);
    }

    public void addContext(final String name, final List<Locatable> intervals) {
        Utils.nonNull(name);
        Utils.nonNull(intervals);
        Utils.validateArg(!contextMap.containsKey(name), "Context with name " + name + " already exists");
        contextMap.put(name, OverlapDetector.create(intervals));
    }

    /**
     * Adds a new stratification group
     * @param name a unique ID
     * @param svType SV type, may be null
     * @param minSize minimum size in bp (inclusive), may be null
     * @param maxSize maximum size in bp (exclusive), may be null
     * @param contextNames reference context names
     */
    public void addStratification(final String name, final GATKSVVCFConstants.StructuralVariantAnnotationType svType,
                                  final Integer minSize, final Integer maxSize, final Set<String> contextNames) {
        addStratification(new Stratum(name, svType, minSize, maxSize, contextNames));
    }

    protected void addStratification(final Stratum stratification) {
        Utils.validateArg(!strata.containsKey(stratification.getName()), "Encountered duplicate name " + stratification.getName());
        strata.put(stratification.getName(), stratification);
    }

    /**
     * Retrieves intervals for the given context
     * @param name context ID
     * @return searchable interval set
     */
    public OverlapDetector<Locatable> getContextIntervals(final String name) {
        return contextMap.get(name);
    }

    /**
     * Factory method for creating a new engine from a config file and set of reference contexts. The config file
     * is a table parsable by {@link TableReader}, with mandatory columns defined in {@link #COLUMN_NAMES}.
     * @param contextMap map from reference context name to interval set
     * @param configFilePath path to stratification config table
     * @param dictionary reference dict
     * @return new engine
     */
    public static SVStatificationEngine create(final Map<String, List<Locatable>> contextMap,
                                               final GATKPath configFilePath,
                                               final SAMSequenceDictionary dictionary) {
        Utils.nonNull(contextMap);
        Utils.nonNull(configFilePath);
        final SVStatificationEngine engine = new SVStatificationEngine(dictionary);
        for (final Map.Entry<String, List<Locatable>> entry : contextMap.entrySet()) {
            engine.addContext(entry.getKey(), entry.getValue());
        }
        try (final TableReader<Stratum> tableReader = TableUtils.reader(configFilePath.toPath(), engine::tableParser)) {
            for (final Stratum stratification : tableReader) {
                engine.addStratification(stratification);
            }
        } catch (final IOException e) {
            throw new GATKException("IO error while reading config table", e);
        }
        return engine;
    }

    /**
     * Get all stratification groups matching a given query record.
     * @param record query record
     * @param overlapFraction minimum overlap fraction (0 to 1)
     * @param numBreakpointOverlaps minimum number of breakpoint ends that must lie in the reference context(s) (0, 1, 2)
     * @param numBreakpointOverlapsInterchrom minimum breakpoint ends for interchromosomal variants (1, 2)
     * @return all matching strata
     */
    public Collection<Stratum> getMatches(final SVCallRecord record, final double overlapFraction, final int numBreakpointOverlaps, final int numBreakpointOverlapsInterchrom) {
        Utils.nonNull(record);
        final List<Stratum> result = new ArrayList<>();
        for (final Stratum stratification : strata.values()) {
            if (stratification.matches(record, overlapFraction, numBreakpointOverlaps, numBreakpointOverlapsInterchrom)) {
                result.add(stratification);
            }
        }
        return result;
    }

    protected Function<DataLine, Stratum> tableParser(TableColumnCollection columns, Function<String, RuntimeException> exceptionFactory) {
        // Check for expected columns
        for (final String column : COLUMN_NAMES) {
            if (!columns.contains(column)) {
                throw exceptionFactory.apply("Missing column " + column);
            }
        }
        // Check there are no extra columns
        if (columns.columnCount() != COLUMN_NAMES.size()) {
            throw exceptionFactory.apply("Expected " + columns.columnCount() + " columns but found " + columns.columnCount());
        }
        return this::parseTableLine;
    }

    protected Stratum parseTableLine(final DataLine dataLine) {
        final GATKSVVCFConstants.StructuralVariantAnnotationType svType = GATKSVVCFConstants.StructuralVariantAnnotationType.valueOf(dataLine.get(SVTYPE_COLUMN));
        final String name = dataLine.get(NAME_COLUMN);
        final Integer minSize = parseIntegerMaybeNull(dataLine.get(MIN_SIZE_COLUMN));
        final Integer maxSize = parseIntegerMaybeNull(dataLine.get(MAX_SIZE_COLUMN));
        final Set<String> contextNames = parseContextString(dataLine.get(CONTEXT_COLUMN));
        return new Stratum(name, svType, minSize, maxSize, contextNames);
    }

    protected Set<String> parseContextString(final String val) {
        if (NULL_TABLE_VALUES.contains(val)) {
            return Collections.emptySet();
        } else {
            final String[] contextArray = val.split(CONTEXT_COLUMN_DELIMITER);
            for (final String context : contextArray) {
                if (!contextMap.containsKey(context)) {
                    throw new GATKException("Could not find context with name " + context);
                }
            }
            return Lists.newArrayList(contextArray).stream().collect(Collectors.toUnmodifiableSet());
        }
    }

    protected Integer parseIntegerMaybeNull(final String val) {
        if (NULL_TABLE_VALUES.contains(val)) {
            return null;
        } else {
            return Integer.valueOf(val);
        }
    }

    public Collection<Stratum> getStrata() {
        return strata.values();
    }

    public class Stratum {

        final GATKSVVCFConstants.StructuralVariantAnnotationType svType;
        final int minSize;  // inclusive
        final int maxSize;  // exclusive
        final List<String> contextNames;
        final String name;

        Stratum(final String name, final GATKSVVCFConstants.StructuralVariantAnnotationType svType,
                final Integer minSize, final Integer maxSize, final Set<String> contextNames) {
            this.name = Utils.nonNull(name);
            for (final String contextName : contextNames) {
                if (contextName != null && !contextMap.containsKey(contextName)) {
                    throw new IllegalArgumentException("Unregistered context name " + contextName);
                }
            }
            if (maxSize != null && minSize != null && maxSize <= minSize) {
                throw new IllegalArgumentException("Min size must be strictly less than max size");
            }
            if (maxSize != null && maxSize < 0) {
                throw new IllegalArgumentException("Max size cannot be less than 0");
            }
            if (maxSize != null && maxSize == Integer.MAX_VALUE) {
                throw new IllegalArgumentException("Max size " + Integer.MAX_VALUE + " is reserved");
            }
            if (minSize != null && minSize < 0) {
                throw new IllegalArgumentException("Min size cannot be less than 0");
            }
            if ((svType == GATKSVVCFConstants.StructuralVariantAnnotationType.BND || svType == GATKSVVCFConstants.StructuralVariantAnnotationType.CTX) && (minSize != null || maxSize != null)) {
                throw new IllegalArgumentException("BND/CTX categories cannot have min or max size (" + name + ")");
            }
            this.svType = svType;
            // Map min from any negative number to negative infinity
            if (minSize == null) {
                this.minSize = Integer.MIN_VALUE;
            } else {
                this.minSize = minSize;
            }
            // Map max from any negative number to infinity
            if (maxSize == null) {
                this.maxSize = Integer.MAX_VALUE;
            } else {
                this.maxSize = maxSize;
            }
            this.contextNames = contextNames.stream().sorted().collect(Collectors.toList());
        }

        protected boolean matches(final SVCallRecord record, final double overlapFraction,
                                  final int numBreakpointOverlaps, final int numBreakpointOverlapsInterchrom) {
            return matchesType(record) && matchesSize(record) && matchesContext(record, overlapFraction, numBreakpointOverlaps, numBreakpointOverlapsInterchrom);
        }

        protected boolean matchesType(final SVCallRecord record) {
            return record.getType() == svType;
        }

        protected boolean matchesSize(final SVCallRecord record) {
            final Integer length = record.getLength();
            if (length == null) {
                // Undefined length requires null min/max boundaries
                return minSize == Integer.MIN_VALUE && maxSize == Integer.MAX_VALUE;
            } else {
                return length >= minSize && length < maxSize;
            }
        }

        /**
         * Determines whether a given query record belongs to this context.
         * @param record query record
         * @param overlapFraction minimum variant interval overlap fraction
         * @param numBreakpointOverlaps minimum number of breakpoint ends that must lie in the context
         * @param numBreakpointOverlapsInterchrom minimum breakpoint ends if the variant is intermchromosomal
         * @return true if the SV matches this context
         */
        public boolean matchesContext(final SVCallRecord record,
                                      final double overlapFraction,
                                      final int numBreakpointOverlaps,
                                      final int numBreakpointOverlapsInterchrom) {
            Utils.nonNull(record);
            Utils.validate(overlapFraction >= 0 && overlapFraction <= 1,
                    "Overlap fraction threshold " + overlapFraction + " must be on [0, 1]");
            Utils.validate(numBreakpointOverlaps >= 0 && numBreakpointOverlaps <= 2,
                    "Breakpoint overlaps threshold " + numBreakpointOverlaps + " must be 0, 1, or 2");
            Utils.validate(numBreakpointOverlapsInterchrom == 1 || numBreakpointOverlapsInterchrom == 2,
                    "Interchromosomal breakpoint overlaps threshold " + numBreakpointOverlapsInterchrom + " must be 1 or 2");
            Utils.validate(!(overlapFraction == 0 && numBreakpointOverlaps == 0),
                    "Overlap fraction and overlapping breakpoints thresholds cannot both be 0");
            if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.INS) {
                // Just require the insertion locus to fall in an interval
                return matchesContextBreakpointOverlap(record, 1);
            } else if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.BND || record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CTX) {
                // Interchromosomal variants
                return matchesContextBreakpointOverlap(record, numBreakpointOverlapsInterchrom);
            } else {
                return matchesContextIntrachromosomal(record, overlapFraction, numBreakpointOverlaps);
            }
        }

        protected boolean matchesContextIntrachromosomal(final SVCallRecord record,
                                                         final double overlapFraction,
                                                         final int numBreakpointOverlaps) {
            return matchesContextOverlapFraction(record, overlapFraction) && matchesContextBreakpointOverlap(record, numBreakpointOverlaps);
        }

        protected boolean matchesContextOverlapFraction(final SVCallRecord record, final double overlapFraction) {
            if (overlapFraction > 0 && !contextNames.isEmpty()) {
                if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CPX) {
                    throw new GATKException("Context overlap for CPX types not currently supported (" + name + ")");
                }
                final SimpleInterval interval = new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionB());
                final List<SimpleInterval> overlaps = new ArrayList<>();
                for (final String context : contextNames) {
                    overlaps.addAll(contextMap.get(context).getOverlaps(interval).stream().map(SimpleInterval::new).collect(Collectors.toList()));
                }
                final List<SimpleInterval> mergedOverlaps = IntervalUtils.sortAndMergeIntervals(overlaps, dictionary, IntervalMergingRule.ALL)
                        .values().stream().flatMap(List::stream).collect(Collectors.toList());
                long overlapLength = 0;
                for (final Locatable overlap : mergedOverlaps) {
                    overlapLength += interval.intersect(overlap).size();
                }
                return overlapLength / (double) interval.getLengthOnReference() >= overlapFraction;
            } else {
                return true;
            }
        }

        protected boolean matchesContextBreakpointOverlap(final SVCallRecord record, final int numBreakpointOverlaps) {
            if (numBreakpointOverlaps > 0 && !contextNames.isEmpty()) {
                if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CPX) {
                    throw new GATKException("Context overlap for CPX types not currently supported (" + name + ")");
                }
                final SimpleInterval intervalA = new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionA());
                final SimpleInterval intervalB = new SimpleInterval(record.getContigB(), record.getPositionB(), record.getPositionB());
                return countAnyContextOverlap(intervalA) + countAnyContextOverlap(intervalB) >= numBreakpointOverlaps;
            } else {
                return true;
            }
        }

        protected int countAnyContextOverlap(final SimpleInterval interval) {
            for (final String context : contextNames) {
                if (contextMap.get(context).overlapsAny(interval)) {
                    return 1;
                }
            }
            return 0;
        }

        public GATKSVVCFConstants.StructuralVariantAnnotationType getSvType() {
            return svType;
        }

        public Integer getMinSize() {
            return minSize;
        }

        public Integer getMaxSize() {
            return maxSize;
        }

        public List<String> getContextNames() {
            return contextNames;
        }

        public String getName() {
            return name;
        }
    }
}
