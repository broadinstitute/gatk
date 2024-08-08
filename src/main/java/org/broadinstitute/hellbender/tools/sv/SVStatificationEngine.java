package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.IOException;
import java.util.*;
import java.util.function.Function;

public class SVStatificationEngine {

    // Configuration table column names
    public static final String SVTYPE_COLUMN = "SVTYPE";
    public static final String MIN_SIZE_COLUMN = "MIN_SIZE";
    public static final String MAX_SIZE_COLUMN = "MAX_SIZE";
    public static final String CONTEXT_COLUMN = "CONTEXT";

    public static final Set<String> NULL_TABLE_VALUES = Set.of("-1", "", "NULL", "NA");

    final Map<String, OverlapDetector<Locatable>> contextMap;
    final Collection<SVStratification> stratifications;

    public SVStatificationEngine() {
        contextMap = new HashMap<>();
        stratifications = new ArrayList<>();
    }

    public void addContext(final String name, final List<Locatable> intervals) {
        Utils.validateArg(!contextMap.containsKey(name), "Context with name " + name + " already exists");
        contextMap.put(name, OverlapDetector.create(intervals));
    }

    public void addStratification(final GATKSVVCFConstants.StructuralVariantAnnotationType svType,
                                  final int minSize, final int maxSize, final String contextName) {
        Utils.validateArg(contextMap.containsKey(contextName), "Unknown context name: " + contextName);
        addStratification(new SVStratification(svType, minSize, maxSize, contextName, contextMap.get(contextName)));
    }

    public OverlapDetector<Locatable> getContextIntervals(final String name) {
        return contextMap.get(name);
    }

    public static SVStatificationEngine create(final Map<String, List<Locatable>> contextMap,
                                               final GATKPath configFilePath) {
        final SVStatificationEngine engine = new SVStatificationEngine();
        for (final Map.Entry<String, List<Locatable>> entry : contextMap.entrySet()) {
            engine.addContext(entry.getKey(), entry.getValue());
        }
        try (final TableReader<SVStratification> tableReader = TableUtils.reader(configFilePath.toPath(), engine::tableParser)) {
            for (final SVStatificationEngine.SVStratification stratification : tableReader) {
                engine.addStratification(stratification);
            }
        } catch (final IOException e) {
            throw new GATKException("IO error while reading config table", e);
        }
        return engine;
    }


    public SVStratification getMatch(final SVCallRecord record, final double overlapFraction,
                                     final int numBreakpointOverlaps, final int numBreakpointOverlapsInterchrom) {
        for (final SVStratification stratification : stratifications) {
            if (stratification.matches(record, overlapFraction, numBreakpointOverlaps, numBreakpointOverlapsInterchrom)) {
                return stratification;
            }
        }
        return null;
    }

    protected void addStratification(final SVStratification stratification) {
        for (final SVStratification other : stratifications) {
            if (!stratification.isMutuallyExclusive(other)) {
                throw new IllegalArgumentException("Attempted to add stratification " + stratification.getName() +
                        " but it is not mutually exclusive with other " + other.getName());
            }
        }
        stratifications.add(stratification);
    }

    protected Function<DataLine, SVStratification> tableParser(TableColumnCollection columns, Function<String, RuntimeException> exceptionFactory) {
        if (columns.columnCount() != 4) {
            throw exceptionFactory.apply("Expected 4 columns but found " + columns.columnCount());
        }
        if (!columns.contains(SVTYPE_COLUMN)) {
            throw exceptionFactory.apply("Missing column " + SVTYPE_COLUMN);
        }
        if (!columns.contains(MIN_SIZE_COLUMN)) {
            throw exceptionFactory.apply("Missing column " + MIN_SIZE_COLUMN);
        }
        if (!columns.contains(MAX_SIZE_COLUMN)) {
            throw exceptionFactory.apply("Missing column " + MAX_SIZE_COLUMN);
        }
        if (!columns.contains(CONTEXT_COLUMN)) {
            throw exceptionFactory.apply("Missing column " + CONTEXT_COLUMN);
        }
        return this::parseTableLine;
    }

    protected SVStatificationEngine.SVStratification parseTableLine(final DataLine dataLine) {
        final GATKSVVCFConstants.StructuralVariantAnnotationType svType = GATKSVVCFConstants.StructuralVariantAnnotationType.valueOf(dataLine.get(SVTYPE_COLUMN));
        final Integer minSize = parseIntegerMaybeNull(dataLine.get(MIN_SIZE_COLUMN));
        final Integer maxSize = parseIntegerMaybeNull(dataLine.get(MAX_SIZE_COLUMN));
        final String context = parseStringMaybeNull(dataLine.get(CONTEXT_COLUMN));
        if (!contextMap.containsKey(context) && context != null) {
            throw new GATKException("Could not find context with name " + context);
        }
        return new SVStatificationEngine.SVStratification(svType, minSize, maxSize, context, contextMap.get(context));
    }

    protected String parseStringMaybeNull(final String val) {
        if (NULL_TABLE_VALUES.contains(val)) {
            return null;
        } else {
            return val;
        }
    }

    protected Integer parseIntegerMaybeNull(final String val) {
        if (NULL_TABLE_VALUES.contains(val)) {
            return null;
        } else {
            return Integer.valueOf(val);
        }
    }

    public Collection<SVStratification> getStratifications() {
        return stratifications;
    }

    public class SVStratification {

        final GATKSVVCFConstants.StructuralVariantAnnotationType svType;
        final int minSize;  // inclusive
        final int maxSize;  // exclusive
        final String contextName;
        final OverlapDetector<Locatable> contextIntervals;
        final String name;

        SVStratification(final GATKSVVCFConstants.StructuralVariantAnnotationType svType,
                         final Integer minSize, final Integer maxSize, final String contextName,
                         final OverlapDetector<Locatable> contextIntervals) {
            if ((contextName == null && contextIntervals != null) || (contextIntervals == null && contextName != null)) {
                throw new IllegalArgumentException("Context name and intervals must either both or neither be null");
            }
            if (maxSize != null && minSize != null && maxSize >= 0 && minSize >= 0 && maxSize < minSize) {
                throw new IllegalArgumentException("Max size cannot be less than min size");
            }
            this.svType = svType;
            // Map min from any negative number to -1
            final String minSizeString;
            if (minSize == null || minSize < 0) {
                this.minSize = -1;
                minSizeString = "";
            } else {
                this.minSize = minSize;
                minSizeString = "_" + minSize;
            }
            // Map max from any negative number to infinity
            final String maxSizeString;
            if (maxSize == null || maxSize < 0) {
                this.maxSize = Integer.MAX_VALUE;
                maxSizeString = "";
            } else {
                this.maxSize = maxSize;
                maxSizeString = "_" + maxSize;
            }
            final String contextNameString = contextName == null ? "" : "_" + contextName;
            this.contextName = contextName;
            this.contextIntervals = contextIntervals;
            this.name = svType.name() + minSizeString + maxSizeString + contextNameString;
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
            return length == null || (length >= minSize && length < maxSize);
        }

        public boolean matchesContext(final SVCallRecord record,
                                      final double overlapFraction,
                                      final int numBreakpointOverlaps,
                                      final int numBreakpointOverlapsInterchrom) {
            Utils.nonNull(record);
            Utils.validate(overlapFraction >= 0 && overlapFraction <= 1,
                    "Invalid overlap fraction threshold " + overlapFraction);
            Utils.validate(numBreakpointOverlaps >= 0 && numBreakpointOverlaps <= 2,
                    "Invalid breakpoint overlaps threshold " + numBreakpointOverlaps);
            Utils.validate(numBreakpointOverlapsInterchrom >= 0 && numBreakpointOverlapsInterchrom <= 2,
                    "Invalid interchromosomal breakpoint overlaps threshold " + numBreakpointOverlapsInterchrom);
            Utils.validate(!(overlapFraction == 0 && numBreakpointOverlaps == 0),
                    "Overlap fraction and overlapping breakpoints thresholds cannot both be 0");
            if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.INS) {
                // Just require the insertion locus to fall in an interval
                return matchesContextBreakpointOverlap(record, 1);
            } else if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.BND) {
                // TODO handle cpx
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
            // TODO handle cpx
            if (overlapFraction > 0) {
                final SimpleInterval interval = new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionB());
                final Set<Locatable> overlaps = contextIntervals.getOverlaps(interval);
                long overlapLength = 0;
                for (final Locatable overlap : overlaps) {
                    overlapLength += interval.intersect(overlap).size();
                }
                return overlapLength / (double) interval.getLengthOnReference() >= overlapFraction;
            } else {
                return true;
            }
        }

        protected boolean matchesContextBreakpointOverlap(final SVCallRecord record, final int numBreakpointOverlaps) {
            if (numBreakpointOverlaps > 0) {
                final SimpleInterval intervalA = new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionA());
                final SimpleInterval intervalB = new SimpleInterval(record.getContigB(), record.getPositionB(), record.getPositionB());
                return countAnyContextOverlap(intervalA) + countAnyContextOverlap(intervalB) >= numBreakpointOverlaps;
            } else {
                return true;
            }
        }

        protected int countAnyContextOverlap(final SimpleInterval interval) {
            if (contextIntervals.overlapsAny(interval)) {
                return 1;
            } else {
                return 0;
            }
        }

        protected boolean isMutuallyExclusive(final SVStratification other) {
            Utils.nonNull(other);
            if (svType != other.getSvType()) {
                return true;
            } else if (!contextName.equals(other.getContextName())) {
                return true;
            } else if (minSize >= other.getMaxSize() || other.getMinSize() >= maxSize) {
                return true;
            } else {
                return false;
            }
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

        public String getContextName() {
            return contextName;
        }

        public String getName() {
            return name;
        }

        @Override
        public String toString() {
            return "SVStratification{" +
                    "svType=" + svType +
                    ", minSize=" + minSize +
                    ", maxSize=" + maxSize +
                    ", contextName='" + contextName + '\'' +
                    ", name='" + name + '\'' +
                    '}';
        }
    }
}
