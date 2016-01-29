package org.broadinstitute.hellbender.tools;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.sql.DataFrame;
import org.apache.spark.sql.SQLContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.knowm.xchart.*;
import org.knowm.xchart.internal.style.Styler;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "PrintGVCFDataFromParquet",
        oneLineSummary = "PrintGVCFDataFromParquet",
        programGroup = VariantProgramGroup.class
)
public final class PrintGVCFDataFromParquet extends GATKSparkTool{

    private static final long serialVersionUID = 1L;

    @Argument(shortName = "F", fullName = "file", doc = "File")
    public String file;

    @Argument(shortName = "p", fullName = "position", doc = "Position", optional = true)
    public List<Integer> positionsArg = null;//10026357;

    @Argument(shortName = "c", fullName = "contig", doc = "Contig", optional = true)
    public String contigArg = null;//10026357;

    @Argument(shortName = "g", fullName = "gene", doc = "Gene", optional = false)
    public String geneArg = null;//10026357;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final SQLContext sqlContext = new SQLContext(ctx);
        final DataFrame df = sqlContext.read().parquet(file);

        //This would be nice but it produces a sharded file which is useless for me
//        byPos.coalesce(1).write().format("com.databricks.spark.csv")
//                .option("delimiter", "\t")
//                .option("header", "true")
//                .save("pos10026357.tsv");

        final int size= 800;
        final int maxXY= 100;
        Chart_XY chart = new ChartBuilder_XY().width(size).height(size).theme(Styler.ChartTheme.GGPlot2)
                .build();
        // Customize Chart
        chart.getStyler().setDefaultSeriesRenderStyle(Series_XY.ChartXYSeriesRenderStyle.Scatter);
        chart.getStyler().setChartTitleVisible(false);
        chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideSW);
        chart.getStyler().setMarkerSize(size/maxXY);
        chart.getStyler().setXAxisMax(maxXY);
        chart.getStyler().setYAxisMax(maxXY);
        chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideN);
        chart.getStyler().setSeriesColors(new Color[]{new Color(0,0,100,100)});

        // Create Chart
        Chart_Category histoChart = new ChartBuilder_Category().width(size).height(size).title("Score Histogram")
                .xAxisTitle("AlternativeAlleleFraction").yAxisTitle("Count").theme(Styler.ChartTheme.GGPlot2).build();

        // Customize Chart
        histoChart.getStyler().setLegendPosition(Styler.LegendPosition.InsideN);
        histoChart.getStyler().setBarWidthPercentage(.96);
        histoChart.getStyler().setPlotGridVerticalLinesVisible(false);
        histoChart.getStyler().setXAxisMin(-0.05);
        histoChart.getStyler().setXAxisMax(1.05);
        histoChart.getStyler().setYAxisMin(0);
        histoChart.getStyler().setYAxisMax(1000);

        //make the dir
        final File dir=new File(geneArg);
        dir.mkdirs();

        final List<Pair<String, Integer>> chrPoss;
        if (positionsArg == null || positionsArg.isEmpty()) {
            System.out.println("nothing given, do we'll do all sites");
            chrPoss = df.select("contig", "pos").dropDuplicates().collectAsList().stream().map(row -> Pair.of(row.getString(0), row.getInt(1))).collect(Collectors.toList());
        } else {
            System.out.println("we'll do" + positionsArg.size() + "sites");
            chrPoss= positionsArg.stream().map(p -> Pair.of(contigArg, p)).collect(Collectors.toList());
        }

        for (int i = 0; i < chrPoss.size(); i++){
            final Pair<String, Integer> pair = chrPoss.get(i);
            final String contig = pair.getKey();
            final int position = pair.getValue();
            final String chrPos = contig + "_" + position;

            System.out.println("Processing :"  + chrPos + "  " + i + "/" + chrPoss.size());
            long before = System.currentTimeMillis();

            final DataFrame byPos = df.filter(df.col("contig").equalTo(contig).and(df.col("pos").equalTo(position)));
            final List<String> chrData = new ArrayList<>();
            final List<Integer> posData = new ArrayList<>();
            final List<Integer> refData = new ArrayList<>();
            final List<Integer> altData = new ArrayList<>();
            byPos.select("contig", "pos", "NREF", "NALT").collectAsList().stream().forEach(row-> {
                chrData.add(row.getString(0));
                posData.add(row.getInt(1));
                refData.add(row.getInt(2));
                altData.add(row.getInt(3));
            });

            //2D histogram
            System.out.println("time for query:"+ (System.currentTimeMillis()-before));
            before = System.currentTimeMillis();
            chart.addSeries(chrPos, refData, altData);
            try {
                VectorGraphicsEncoder.saveVectorGraphic(chart, dir.getAbsolutePath() + "/RefAlt_" + chrPos, VectorGraphicsEncoder.VectorGraphicsFormat.PDF);
            } catch (IOException e) {
                throw new GATKException("error", e);
            }
            System.out.println("time for 2D plot:"+ (System.currentTimeMillis()-before));
            before = System.currentTimeMillis();
            chart.removeSeries(chrPos);

            //1D histogram
            final Collection<Double> afData= new ArrayList<>();
            for (int j = 0; j < refData.size(); j++) {
                final int altCounts= altData.get(j);
                final int refCounts= refData.get(j);
                final double af = (altCounts * 1.0) / (refCounts+altCounts);
                afData.add(af);
            }
            final Histogram h1= new Histogram(afData, 201, 0.0, 1.0);
            histoChart.addSeries(chrPos, h1.getxAxisData(), h1.getyAxisData());
            try {
                VectorGraphicsEncoder.saveVectorGraphic(histoChart, dir.getAbsolutePath() + "/Histogram_chr" + chrPos, VectorGraphicsEncoder.VectorGraphicsFormat.PDF);
            } catch (IOException e) {
                throw new GATKException("error", e);
            }
            System.out.println("time for 1D plot:"+ (System.currentTimeMillis()-before));
            before = System.currentTimeMillis();
            histoChart.removeSeries(chrPos);

            //TSV
                try (final PrintStream os = new PrintStream(FileUtils.openOutputStream(new File(dir, chrPos + ".tsv")), true)){
                    os.println("CHR\t" + "POS\t" +"NREF\t" + "NALT");
                    for (int j = 0; j < refData.size(); j++) {
                        os.println(chrData.get(j) + "\t" + posData.get(j) + "\t" + refData.get(j) + "\t" + altData.get(j));
                    }
                } catch (IOException e) {
                    throw new GATKException("tsv writing", e);
                }
            System.out.println("time for TSV plot:"+ (System.currentTimeMillis()-before));
            System.out.println();
        }

    }
}
