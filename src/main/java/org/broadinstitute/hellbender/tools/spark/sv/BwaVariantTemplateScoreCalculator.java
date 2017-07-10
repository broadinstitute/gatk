package org.broadinstitute.hellbender.tools.spark.sv;

import avro.shaded.com.google.common.collect.Iterables;
import breeze.generic.UFunc;
import breeze.linalg.*;
import breeze.linalg.Vector;
import breeze.linalg.operators.*;
import breeze.linalg.support.*;
import breeze.math.Semiring;
import breeze.storage.Zero;
import htsjdk.samtools.*;
import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.util.DoubleArray;
import org.apache.hadoop.mapred.JobConf;
import org.apache.logging.log4j.LogManager;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Function1;
import scala.Function2;
import scala.Tuple2;
import java.util.Map;

import java.io.*;
import java.lang.Iterable;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;

import org.apache.logging.log4j.Logger;
import scala.collection.*;
import scala.collection.Iterator;
import scala.collection.immutable.IndexedSeq;
import scala.collection.immutable.Set;
import scala.math.Numeric;
import scala.math.Ordering;
import scala.reflect.ClassTag;

/**
 * Created by valentin on 5/19/17.
 */
public class BwaVariantTemplateScoreCalculator implements StructuralVariantTemplateHaplotypeScoreCalculator, Serializable {

    private static final long serialVersionUID = 1L;

    private transient final JavaSparkContext ctx;

    private InsertSizeDistribution insertSizeDistribution;

    private Logger logger = LogManager.getLogger(BwaVariantTemplateScoreCalculator.class);

    private final TemplateAlignmentPenaltiesArgumentCollection penalties;

    public BwaVariantTemplateScoreCalculator(final JavaSparkContext ctx, final InsertSizeDistribution insertSizeDistribution) {
        this.ctx = Utils.nonNull(ctx);
        this.insertSizeDistribution = Utils.nonNull(insertSizeDistribution);
        this.penalties = new TemplateAlignmentPenaltiesArgumentCollection();
    }

    @Override
    public void calculate(final TemplateHaplotypeScoreTable table) {
        final SAMFileHeader samHeader = new SAMFileHeader();
        samHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        final List<GATKRead> unmappedReadList = table.templates().stream()
                .flatMap(t -> {
                    final List<Template.Fragment> fragments = t.fragments();
                    return fragments.stream().map(f -> f.toUnmappedRead(samHeader, fragments.size() > 1));
                })
                .collect(Collectors.toList());
        final JavaPairRDD<String, GATKRead> unmappedReads = ctx.parallelize(unmappedReadList).mapToPair(r -> new Tuple2<>(r.getName(), r));
        final Broadcast<Map<String, Template>> templatesByName = ctx.broadcast(table.templates()
            .stream().collect(Collectors.toMap(Template::name, t -> t)));

        for (int i = 0; i < table.numberOfHaplotypes(); i++) {
            final int alleleIndex = i;
            final Haplotype haplotype = table.haplotypes().get(i);
            final SAMSequenceDictionary referenceDictionary = new SAMSequenceDictionary();
            final SAMSequenceRecord contig = new SAMSequenceRecord("1", table.haplotypes().get(i).length());
            referenceDictionary.addSequence(contig);
            final File referenceImage = composeReference(table.haplotypes().get(i));
            final BwaSparkEngine bwa = prepareBwaEngine(referenceImage, samHeader, referenceDictionary);
            final JavaRDD<GATKRead> bwaInputReads = unmappedReads.partitionBy(new HashPartitioner(1)).values().sortBy(r -> r.getName(), true, 1);
            final JavaRDD<GATKRead> mappedReads = bwa.align(bwaInputReads).filter(r -> !r.isSecondaryAlignment() && !r.isSupplementaryAlignment());
            //final Set<String> withSupplementary = mappedReads.filter(GATKRead::isSupplementaryAlignment).map(GATKRead::getName).collect().stream()
            //        .collect(Collectors.toSet());
            //final Broadcast<Set<String>> withSupplementaryBC = ctx.broadcast(withSupplementary);
            // mappedReads.filter(r -> withSupplementaryBC.getValue().contains(r.getName()))
            //         .mapToPair(r -> new Tuple2<>(r.getName(), r))
            //         .groupByKey()
            //         .foreach(t -> {
            //             final Collection c = CollectionUtils.collect(t._2().iterator(), o -> o);
            //             c.size();
            //         });

            final JavaPairRDD<Template, Iterable<GATKRead>> templatesAndMappedReads = mappedReads.mapToPair(r -> new Tuple2<String, GATKRead>(r.getName(), r)).groupByKey()
                    .mapToPair(p -> new Tuple2<>(templatesByName.getValue().get(p._1()),  p._2()));

            templatesAndMappedReads.toLocalIterator().forEachRemaining(p ->
                table.setMappingInfo(alleleIndex, table.indexOf(p._1()), calculateMappingInfo(haplotype, p._1(), p._2())));
            printOutCounters();
            resetCounters();
            disposeReference(referenceImage);
        }
        final int maximumInsertSize = table.maximumInsertSize();
        final double minimumAlignmentScore = table.minimumAlignmentScore();


        for (int h = 0; h < table.numberOfHaplotypes(); h++) {
            for (int j = 0; j < table.numberOfTemplates(); j++) {
                final TemplateMappingInformation mappingInfo = table.getMappingInfo(h, j);
                final TemplateMappingInformation otherMappingInfo = table.getMappingInfo(h == 0 ? 1 : 0, j);
                table.set(h, j, calculateScore(mappingInfo, maximumInsertSize, minimumAlignmentScore, otherMappingInfo));
            }
        }

    }

    private double calculateScore(TemplateMappingInformation mappingInfo, final int maximumInsertSize, final double minimumAlignmentScore, TemplateMappingInformation other) {
        if (mappingInfo == null)
            return Double.NaN;
        double result = 0;
        if (mappingInfo.insertSize.isPresent()) {
            result += -10.0 * (insertSizeDistribution.logProbability(mappingInfo.insertSize.getAsInt()) / Math.log(10));
        } else {
            result += -10.0 * (insertSizeDistribution.logProbability(maximumInsertSize) / Math.log(10)) + penalties.improperPairPenaltyFactor;
        }
        if (mappingInfo.firstAlignmentScore.isPresent()) {
            result += mappingInfo.firstAlignmentScore.getAsDouble();
        } else {
            result += other.firstAlignmentScore.orElse(minimumAlignmentScore) + penalties.unmappedReadPenaltyFactor;
        }
        if (mappingInfo.secondAlignmentScore.isPresent()) {
            result += mappingInfo.secondAlignmentScore.getAsDouble();
        } else {
            result += other.secondAlignmentScore.orElse(minimumAlignmentScore) + penalties.unmappedReadPenaltyFactor;
        }
        return result * -.1;
    }

    private TemplateMappingInformation calculateMappingInfo(final Haplotype haplotype, final Template template, final Iterable<GATKRead> gatkReads) {
        final List<GATKRead> readList = new ArrayList<>();
        Iterables.addAll(readList, gatkReads);
        if (readList.size() == 0) {
            return new TemplateMappingInformation();
        } else if (readList.size() == 1) {
            final GATKRead singleton = readList.get(0);
            if (singleton.isUnmapped()) {
                return new TemplateMappingInformation();
            } else {
                return new TemplateMappingInformation(calculateAlignmentLikelihood(singleton, haplotype, penalties));
            }
        } else if (readList.size() == 2) {
            final GATKRead first = readList.get(0);
            final GATKRead second = readList.get(1);
            if (first.isUnmapped() && second.isUnmapped()) {
                return new TemplateMappingInformation();
            } else if (first.isUnmapped() != second.isUnmapped()) {
                return new TemplateMappingInformation(calculateAlignmentLikelihood(first.isUnmapped() ? first : second, haplotype, penalties));
            } else {
                final GATKRead forward = first.getUnclippedStart() <= second.getUnclippedStart() ? first : second;
                final GATKRead reverse = forward == first ? second : first;
                return new TemplateMappingInformation(calculateAlignmentLikelihood(forward, haplotype, penalties), calculateAlignmentLikelihood(reverse, haplotype, penalties),
                        reverse.getUnclippedEnd() - forward.getUnclippedStart());
            }
        } else {
            return new TemplateMappingInformation();
        }
    }

    private double calculateAlignmentLikelihood(final GATKRead read, final Haplotype haplotype, final TemplateAlignmentPenaltiesArgumentCollection penalties) {
        final byte[] readBases = read.getBases();
        final byte[] quals = read.getBaseQualities();
        final byte[] haplotypeBases = haplotype.getBases();
        double result = 0;
        int nextReadBase = 0;
        int nextHaplotypeBase = read.getStart() - 1;
        for (final CigarElement ce : read.getCigar().getCigarElements()) {
            switch (ce.getOperator()) {
                case EQ:
                case X:
                case M:
                    for (int i = 0; i < ce.getLength(); i++) {
                        if (readBases[nextReadBase + i] == haplotypeBases[nextHaplotypeBase + i]) {
                            result += QualityUtils.qualToProbLog10(quals[nextReadBase + i]) * -10.0;
                        } else {
                            result += Math.min(penalties.maximumMismatchPenalty, QualityUtils.qualToErrorProbLog10(quals[nextReadBase + i]) * -10.0);
                        }
                    }
                    break;
                case I:
                case D:
                    result += Math.min(penalties.maximumIndelPenalty, penalties.gapOpenPenalty + penalties.gapExtensionPenalty * (ce.getLength() - 1));
                    break;
                case S:
                case N:
                case H:
                    for (int i = 0; i < ce.getLength(); i++) {
                        result += Math.min(penalties.maximumMismatchPenalty, QualityUtils.qualToErrorProbLog10(quals[nextReadBase + i]) * -10.0);
                    }
                    break;
                case P:
            }
            if (ce.getOperator().consumesReferenceBases()) {
                nextHaplotypeBase += ce.getLength();
            }
            if (ce.getOperator().consumesReadBases()) {
                nextReadBase += ce.getLength();
            }
        }
        return result;
    }

    private int unmapped = 0;
    private int onlyOneMapped = 0;
    private int nonProperlyMapped1 = 0;
    private int nonProperlyMapped2 = 0;
    private int properlyMapped = 0;
    private List<Integer> sizes = new ArrayList<>();

    private double calculateLikelihood(final Haplotype haplotype, final Template template, final Iterable<GATKRead> gatkReads) {
        final List<GATKRead> readList = new ArrayList<>(2);
        int pos = 0;
        final List<Integer> variantPositions = new ArrayList<>();
        for (final CigarElement el : haplotype.getCigar().getCigarElements()) {
            if (el.getOperator() != CigarOperator.M || el.getOperator() != CigarOperator.X) {
                variantPositions.add(pos);
                if (el.getOperator().consumesReadBases()) {
                    variantPositions.add(pos + el.getLength() - 1);
                }
            }
            if (el.getOperator().consumesReadBases()) {
                pos += el.getLength();
            }
        }

        gatkReads.forEach(readList::add);
        final int mappedReads = (int) readList.stream().filter( r -> !r.isUnmapped()).count();
        if (mappedReads == 2) {
            final GATKRead forward = readList.stream().sorted(Comparator.comparingInt(GATKRead::getUnclippedStart)).findFirst().orElseThrow(IllegalStateException::new);
            final GATKRead reverse = readList.stream().sorted(Comparator.comparingInt(GATKRead::getUnclippedStart)).skip(1).findFirst().orElseThrow(IllegalStateException::new);
            if (forward.getMappingQuality() < 5 || reverse.getMappingQuality() < 5) {
                sizes.add(-1);
                return Double.NaN;
            } else if (closeToVariantPosition(forward, variantPositions) || closeToVariantPosition(reverse, variantPositions)) {
                sizes.add(-1);
                return Double.NaN;
            }
            final int start = forward.getUnclippedStart();
            final int end = reverse.getUnclippedEnd();
            if (forward.getUnclippedEnd() > end) {
                nonProperlyMapped1++;
                sizes.add(-1);
                return Double.NaN;
            }
            properlyMapped++;
            final int size = end - start;
            final double result = insertSizeDistribution.logProbability(size);
            sizes.add(size);
            return result;
        } else {
            if (mappedReads == 1) onlyOneMapped++; else unmapped++;
            sizes.add(-1);
            return Double.NaN;
        }
    }

    private boolean closeToVariantPosition(final GATKRead read, final List<Integer> variantPositions) {
        final int start = read.getStart();
        final int end = read.getEnd();
        final int threshold = (read.getLength() * 3) / 4;
        for (int i : variantPositions) {
            if (Math.max(Math.abs(i - start), Math.abs(i - end)) < threshold) {
                return true;
            }
        }
        return false;
    }

    void printOutCounters() {
        if (logger.isDebugEnabled()) {
            logger.debug("DEBUG counts are " + String.format("[u=%d, 1=%d, np1=%d, np2=%d, p=%d, t=%d]", unmapped, onlyOneMapped, nonProperlyMapped1, nonProperlyMapped2, properlyMapped,
                    unmapped + onlyOneMapped + properlyMapped + nonProperlyMapped1 + nonProperlyMapped2));
        }
        System.err.println("Sizes " + Arrays.toString(sizes.toArray()));
    }

    void resetCounters() {
        unmapped = onlyOneMapped = properlyMapped = nonProperlyMapped1 = nonProperlyMapped2 = 0;
        sizes.clear();
    }

    private void disposeReference(final File referenceImage) {
        referenceImage.delete();
    }

    private BwaSparkEngine prepareBwaEngine(final File referenceImage, final SAMFileHeader header, final SAMSequenceDictionary dict) {
        return new BwaSparkEngine(ctx, referenceImage.getAbsolutePath(), header, dict);
    }

    private File composeReference(final Haplotype haplotype) {
        File dir = null;
        File img = null;
        try {
            dir = File.createTempFile("gatk-sv-bwa-tmp-ref-", ".dir");
            dir.delete();
            dir.mkdir();
            img = File.createTempFile("gatk-sv-bwa-ref-img-", ".fa.img");
            final File fasta = new File(dir, "ref.fa");
            final File fai = new File(dir, "ref.fa.fai");

            final PrintWriter fastaWriter = new PrintWriter(new FileWriter(fasta));
            fastaWriter.println(">1");
            fastaWriter.flush();
            final long offset = fasta.length();
            int nextIdx = 0;
            while (nextIdx < haplotype.length()) {
                fastaWriter.println(haplotype.getBaseString(nextIdx, Math.min(haplotype.length(), nextIdx += 60)));
            }
            fastaWriter.close();
            final PrintWriter faiWriter = new PrintWriter(new FileWriter(fai));
            faiWriter.println(String.join("\t", "1", "" + haplotype.length(), "" + offset, "60", "61"));
            faiWriter.close();
            while (true) {
                try {
                    final String bwaIndexCommand = String.format("bwa index %s ", fasta.getAbsoluteFile());
                    // need to add build index to the bwa-jni.
                    final Process process = Runtime.getRuntime().exec(bwaIndexCommand);
                    final BufferedReader processInputStream = new BufferedReader(new InputStreamReader(process.getErrorStream()));
                    final StringBuffer stdOutput = new StringBuffer(20000);
                    stdOutput.append("OUTPUT:  \n");
                    final Runnable processOutputProcessor = () -> processInputStream.lines().forEach(l -> stdOutput.append(l).append('\n'));
                    final Thread processOutputThread = new Thread(processOutputProcessor);
                    processOutputThread.setDaemon(true);
                    processOutputThread.start();
                    while (true) {
                        try {
                            processOutputThread.join();
                            break;
                        } catch (final InterruptedException ex) {

                        }
                    }
                    final int exitCode = process.waitFor();

                    if (exitCode != 0) {
                        throw new GATKException(
                                String.format("could not create the haplotype index at '%s' with exit-code '%d' and command '%s' and output/error: %s", dir.getAbsolutePath(), exitCode, "cmd", stdOutput.toString()));
                    } else {
                        //FileUtils.deleteDirectory(dir);
                        break;
                    }
                } catch (final InterruptedException e) {
                    // nothing to do here, just loop again.
                }
            }
            BwaMemIndex.createIndexImage(fasta.getAbsolutePath(), img.getAbsolutePath());
            FileUtils.deleteDirectory(dir);
            return img;
        } catch (final IOException ex) {
            throw new GATKException("problems when creating the haplotype index image at '" + img + "' from '" + dir + "'");
        }
    }


}
