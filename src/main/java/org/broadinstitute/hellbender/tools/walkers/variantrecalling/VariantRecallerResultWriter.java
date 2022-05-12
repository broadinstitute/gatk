package org.broadinstitute.hellbender.tools.walkers.variantrecalling;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.Tuple;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class VariantRecallerResultWriter {
    final PrintWriter     pw;
    boolean               first = true;
    final boolean         debugFormat = false;

    protected VariantRecallerResultWriter(final File file) {
        try {
            pw = new PrintWriter(file);
        } catch (IOException e) {
            throw new GATKException("Unable to open output file: " + file.getAbsolutePath(), e);
        }
    }

    protected void close() {
        pw.close();
    }

    protected void add(final Locatable loc, final List<Map<Integer, AlleleLikelihoods<GATKRead, Allele>>> genotypeLikelihoodsList,
                    final List<VariantContext> variants, final List<AssemblyResultSet> assemblyResultList,
                    final List<SAMFileHeader> fileHeaderList) {

        // build a map of vcs by startPos
        final Map<Integer, VariantContext>    vcStartPos = new LinkedHashMap<>();
        variants.forEach(vc -> {
            vcStartPos.put(vc.getStart(), vc);
        });

        // print location (as a separator)
        if ( debugFormat ) {
            if (first) {
                first = false;
            } else {
                pw.println("");
            }
            pw.println("loc: " + loc);
            pw.println("ref: " + new String(assemblyResultList.get(0).getFullReferenceWithPadding()));
        }

        // loop on result
        final List<Tuple<Double, String>> vcLines = new LinkedList<>();
        for ( int index = 0 ; index < genotypeLikelihoodsList.size() ; index++ ) {
            final Map<Integer, AlleleLikelihoods<GATKRead, Allele>>   genotypeLikelihoods = genotypeLikelihoodsList.get(index);
            final SAMFileHeader                                       fileHeader = fileHeaderList.get(index);

            genotypeLikelihoods.forEach((startPos, likelihoods) -> {

                // DK: map to vc? ignore unmapped?
                final VariantContext vc = vcStartPos.get(startPos);
                if (vc != null) {

                    if (debugFormat) {
                        pw.println("");
                        pw.println("variant: " + vc.getContig() + ":" + vc.getStart());
                        pw.println("variant-info: " + vc);

                        // reads
                        pw.println("");
                        pw.println("reads: " + likelihoods.evidenceCount());
                        likelihoods.sampleEvidence(0).forEach(read -> {
                            pw.println("read: " + read);
                        });

                        // alleles
                        pw.println("");
                        pw.println("alleles: " + likelihoods.alleles().size());
                        likelihoods.alleles().forEach(allele -> {
                            pw.println("allele: " + allele);
                        });
                    } else {
                        pw.print("#" + vc.getContig() + ":" + vc.getStart());
                        if (vc.getType() != VariantContext.Type.MIXED && (vc.getEnd() != vc.getStart()))
                            pw.print("-" + vc.getEnd());
                        pw.print(" " + loc);
                        likelihoods.alleles().forEach(allele -> {
                            pw.print(" " + allele);
                        });
                        pw.println("");
                    }

                    final SimpleInterval vcSpan = new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd());

                    // matrix
                    if (debugFormat) {
                        pw.println("");
                        pw.println("matrix:");
                    }
                    final LikelihoodMatrix<GATKRead, Allele> matrix = likelihoods.sampleMatrix(0);
                    final double[][] values = new double[matrix.numberOfAlleles()][matrix.evidenceCount()];
                    for (int alleleIndex = 0; alleleIndex < matrix.numberOfAlleles(); alleleIndex++) {
                        matrix.copyAlleleLikelihoods(alleleIndex, values[alleleIndex], 0);
                    }
                    final double[] lineValues = new double[matrix.numberOfAlleles()];
                    for (int evidenceIndex = 0; evidenceIndex < matrix.evidenceCount(); evidenceIndex++) {

                        // determine matrix values
                        boolean allValuesNegativeInfinity = true;
                        double sortKey = Double.NEGATIVE_INFINITY;
                        for (int alleleIndex = 0; alleleIndex < matrix.numberOfAlleles(); alleleIndex++) {
                            lineValues[alleleIndex] = values[alleleIndex][evidenceIndex];
                            if (lineValues[alleleIndex] != Double.NEGATIVE_INFINITY) {
                                allValuesNegativeInfinity = false;
                            }
                            sortKey = lineValues[alleleIndex];
                        }

                        // lines which have all values of -Inf are complete alignment failures. Ignore them
                        if (allValuesNegativeInfinity) {
                            continue;
                        }

                        // determine length in key space
                        final GATKRead read = matrix.evidence().get(evidenceIndex);
                        int keyspaceLength = 0;
                        if (read instanceof FlowBasedRead)
                            keyspaceLength = ((FlowBasedRead) read).getKeyLength();

                        // build basic matrix line
                        final StringBuilder line = new StringBuilder();
                        line.append(read.getName());
                        line.append(' ');
                        line.append(keyspaceLength);
                        line.append(' ');
                        line.append(read.isDuplicate() ? 1 : 0);
                        line.append(' ');
                        line.append(read.isReverseStrand() ? 1 : 0);
                        line.append(' ');
                        line.append(read.getMappingQuality());
                        line.append(' ');
                        line.append(StringUtils.join(ArrayUtils.toObject(lineValues), " "));

                        // add bytes at variant location?
                        final StringBuilder bases = new StringBuilder();
                        int firstBaseUnclippedOfs = 0;
                        if (read.getContig() != null) {
                            SimpleInterval readSpan = new SimpleInterval(read.getContig(), read.getStart(), read.getEnd());
                            if (readSpan.contains(vcSpan)) {
                                int ofs = vcSpan.getStart() - readSpan.getStart();
                                int vcLength = vcSpan.getEnd() - vc.getStart() + 1;
                                for (int i = 0; i < vcLength; i++) {
                                    int readOfs = getOffsetOnRead(read, ofs + i);
                                    if (readOfs >= 0) {
                                        bases.append((char) read.getBase(readOfs));
                                        if (!read.isReverseStrand()) {
                                            firstBaseUnclippedOfs = readOfs + (read.getStart() - read.getUnclippedStart());
                                        } else {
                                            firstBaseUnclippedOfs = (read.getLength() - readOfs - 1) + (read.getUnclippedEnd() - read.getEnd());
                                        }
                                    } else {
                                        // we don't like '?' bases anymore!
                                        bases.setLength(0);
                                        break;
                                    }
                                }
                            }
                        }

                        if (bases.length() > 0) {
                            line.append(' ');
                            line.append(bases);
                            line.append(' ');
                            line.append(firstBaseUnclippedOfs);
                            line.append(' ');
                            line.append(fileHeader.getReadGroup(read.getReadGroup()).getSample());
                            vcLines.add(new Tuple<>(sortKey, line.toString()));
                        }
                    }

                }
                // optional: sort vcLines on second column if present
                vcLines.sort((o1, o2) -> {
                    return -Double.compare(o1.a, o2.a);
                });

                // pour into output file
                vcLines.forEach(doubleStringTuple -> pw.println(doubleStringTuple.b));
                vcLines.clear();
            });
        }
    }

    public static int getOffsetOnRead(final GATKRead read, final int intialOfs) {
        int     readOfs = 0;
        int     ofs = intialOfs;

        Iterator<CigarElement> iter = read.getCigar().iterator();
        while ( iter.hasNext() ) {
            CigarElement    elem = iter.next();
            CigarOperator   op = elem.getOperator();
            if ( op.consumesReadBases() ) {
                if (ofs < elem.getLength() ) {
                    return readOfs + ofs;
                } else {
                    readOfs += elem.getLength();
                }
            }
            ofs -= (op.consumesReferenceBases() ? elem.getLength() : 0);
        }

        // if here, not found
        return -1;
    }
}
