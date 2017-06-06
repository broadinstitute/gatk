package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMHeaderRecordComparator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.hadoop.mapred.JobConf;
import org.apache.logging.log4j.LogManager;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.*;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;

import org.apache.logging.log4j.Logger;

/**
 * Created by valentin on 5/19/17.
 */
public class BwaVariantTemplateScoreCalculator implements StructuralVariantTemplateHaplotypeScoreCalculator, Serializable {

    private static final long serialVersionUID = 1L;

    private transient final JavaSparkContext ctx;

    private InsertSizeDistribution insertSizeDistribution;

    private AlignmentPenalties penalties;

    private Logger logger = LogManager.getLogger(BwaVariantTemplateScoreCalculator.class);

    public BwaVariantTemplateScoreCalculator(final JavaSparkContext ctx, final InsertSizeDistribution insertSizeDistribution) {
        this.ctx = Utils.nonNull(ctx);
        this.insertSizeDistribution = Utils.nonNull(insertSizeDistribution);
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
            final JavaRDD<GATKRead> mappedReads = bwa.align(bwaInputReads);

            final JavaPairRDD<Template, Iterable<GATKRead>> templatesAndMappedReads = mappedReads.mapToPair(r -> new Tuple2<String, GATKRead>(r.getName(), r)).groupByKey()
                    .mapToPair(p -> new Tuple2<>(templatesByName.getValue().get(p._1()),  p._2()));

            templatesAndMappedReads.toLocalIterator().forEachRemaining(p ->
                table.set(alleleIndex, table.indexOf(p._1()), calculateLikelihood(haplotype, p._1(), p._2())));
            printOutCounters();
            resetCounters();
            disposeReference(referenceImage);
        }
    }

    private int unmapped = 0;
    private int onlyOneMapped = 0;
    private int nonProperlyMapped1 = 0;
    private int nonProperlyMapped2 = 0;
    private int properlyMapped = 0;

    private double calculateLikelihood(final Haplotype haplotype, final Template template, final Iterable<GATKRead> gatkReads) {
        final List<GATKRead> readList = new ArrayList<>(2);

        gatkReads.forEach(readList::add);
        final int mappedReads = (int) readList.stream().filter( r -> !r.isUnmapped()).count();
        if (mappedReads == 2) {
            final GATKRead forward = readList.stream().sorted(Comparator.comparingInt(GATKRead::getUnclippedStart)).findFirst().orElseThrow(IllegalStateException::new);
            final GATKRead reverse = readList.stream().sorted(Comparator.comparingInt(GATKRead::getUnclippedStart)).skip(1).findFirst().orElseThrow(IllegalStateException::new);
            final int start = forward.getUnclippedStart();
            final int end = reverse.getUnclippedEnd();
            if (forward.getUnclippedEnd() > end) {
                nonProperlyMapped1++;
                return Double.NaN;
            }
            properlyMapped++;
            final int size = end - start;
            final double result = insertSizeDistribution.logProbability(size);
            return result;
        } else {
            if (mappedReads == 1) onlyOneMapped++; else unmapped++;
            return Double.NaN;
        }
    }

    void printOutCounters() {
        if (logger.isDebugEnabled()) {
            logger.debug("DEBUG counts are " + String.format("[u=%d, 1=%d, np1=%d, np2=%d, p=%d, t=%d]", unmapped, onlyOneMapped, nonProperlyMapped1, nonProperlyMapped2, properlyMapped,
                    unmapped + onlyOneMapped + properlyMapped + nonProperlyMapped1 + nonProperlyMapped2));
        }
    }

    void resetCounters() {
        unmapped = onlyOneMapped = properlyMapped = nonProperlyMapped1 = nonProperlyMapped2 = 0;
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
                    final int exitCode = Runtime.getRuntime().exec(bwaIndexCommand).waitFor();
                    if (exitCode != 0) {
                        throw new GATKException(
                                String.format("could not create the haplotype index at '%s' with exit-code '%d' and command '%s'", dir.getAbsolutePath(), exitCode, ""));
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
