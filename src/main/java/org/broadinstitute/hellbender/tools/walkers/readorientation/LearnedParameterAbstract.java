//package org.broadinstitute.hellbender.tools.walkers.readorientation;
//
//import org.broadinstitute.hellbender.utils.tsv.*;
//
//import java.io.IOException;
//import java.nio.file.Path;
//import java.util.Arrays;
//import java.util.List;
//import java.util.stream.Collectors;
//
//public abstract class LearnedParameterAbstract {
//    public abstract String getReferenceContext();
//    public abstract String getRevCompContext();
//    public abstract double getParameter(final ArtifactState artifactState);
//    public abstract int getNumAltExamples();
//    public abstract int getNumExamples();
//
//    public static class LearnedParameterTableWriter extends TableWriter<LearnedParameterAbstract> {
//        public LearnedParameterTableWriter(final Path output, final String sample) throws IOException {
//            super(output, LearnedParameterTableColumn.COLUMNS);
//            writeMetadata(TableUtils.SAMPLE_METADATA_TAG, sample);
//        }
//
//        @Override
//        protected void composeLine(final LearnedParameterAbstract parameters, final DataLine dataLine) {
//            dataLine.set(LearnedParameterTableColumn.CONTEXT.toString(), parameters.getReferenceContext())
//                    .set(LearnedParameterTableColumn.REV_COMP.toString(), parameters.getRevCompContext())
//                    .set(LearnedParameterTableColumn.F1R2_A.toString(),
//                            parameters.getParameter(ArtifactState.F1R2_A))
//                    .set(LearnedParameterTableColumn.F1R2_C.toString(),
//                            parameters.getParameter(ArtifactState.F1R2_C))
//                    .set(LearnedParameterTableColumn.F1R2_G.toString(),
//                            parameters.getParameter(ArtifactState.F1R2_G))
//                    .set(LearnedParameterTableColumn.F1R2_T.toString(),
//                            parameters.getParameter(ArtifactState.F1R2_T))
//                    .set(LearnedParameterTableColumn.F2R1_A.toString(),
//                            parameters.getParameter(ArtifactState.F2R1_A))
//                    .set(LearnedParameterTableColumn.F2R1_C.toString(),
//                            parameters.getParameter(ArtifactState.F2R1_C))
//                    .set(LearnedParameterTableColumn.F2R1_G.toString(),
//                            parameters.getParameter(ArtifactState.F2R1_G))
//                    .set(LearnedParameterTableColumn.F2R1_T.toString(),
//                            parameters.getParameter(ArtifactState.F2R1_T))
//                    .set(LearnedParameterTableColumn.HOM_REF.toString(),
//                            parameters.getParameter(ArtifactState.HOM_REF))
//                    .set(LearnedParameterTableColumn.GERMLINE_HET.toString(),
//                            parameters.getParameter(ArtifactState.GERMLINE_HET))
//                    .set(LearnedParameterTableColumn.SOMATIC_HET.toString(),
//                            parameters.getParameter(ArtifactState.SOMATIC_HET))
//                    .set(LearnedParameterTableColumn.HOM_VAR.toString(),
//                            parameters.getParameter(ArtifactState.HOM_VAR))
//                    .set(LearnedParameterTableColumn.N.toString(), parameters.getNumExamples())
//                    .set(LearnedParameterTableColumn.N_ALT.toString(), parameters.getNumAltExamples());
//        }
//    }
//
//    public static class LearnedParameterTableReader<A extends LearnedParameterAbstract> extends TableReader<A> {
//        public LearnedParameterTableReader(final Path table) throws IOException {
//            USER THE READER WITH LAMBADA (or just use one class....)
//            super(table);
//        }
//
//        @Override
//        protected A createRecord(final DataLine dataLine) {
//            final String referenceContext = dataLine.get(LearnedParameterTableColumn.CONTEXT);
//            final double[] parameters = new double[ArtifactState.values().length];
//            for (LearnedParameterTableColumn column : LearnedParameterTableColumn.getArtifactStateColumns()){
//                parameters[column.getState().ordinal()] = Double.valueOf(dataLine.get(column));
//            }
//
//            final int numExamples = Integer.parseInt(dataLine.get(LearnedParameterTableColumn.N));
//            final int numAltExamples = Integer.parseInt(dataLine.get(LearnedParameterTableColumn.N_ALT));
//
//            return new A(referenceContext, parameters, numExamples, numAltExamples);
//        }
//    }
//
//    private enum LearnedParameterTableColumn {
//        CONTEXT("context"),
//        REV_COMP("rev_comp"),
//        F1R2_A("f1r2_a", ArtifactState.F1R2_A),
//        F1R2_C("f1r2_c", ArtifactState.F1R2_C),
//        F1R2_G("f1r2_g", ArtifactState.F1R2_G),
//        F1R2_T("f1r2_t", ArtifactState.F1R2_T),
//        F2R1_A("f2r1_a", ArtifactState.F2R1_A),
//        F2R1_C("f2r1_c", ArtifactState.F2R1_C),
//        F2R1_G("f2r1_g", ArtifactState.F2R1_G),
//        F2R1_T("f2r1_t", ArtifactState.F2R1_T),
//        HOM_REF("hom_ref", ArtifactState.HOM_REF),
//        GERMLINE_HET("germline_het", ArtifactState.GERMLINE_HET),
//        SOMATIC_HET("somatic_het", ArtifactState.SOMATIC_HET),
//        HOM_VAR("hom_var", ArtifactState.HOM_VAR),
//        N("num_examples"),
//        N_ALT("num_alt_examples");
//
//        private String columnName;
//        private ArtifactState state;
//
//        LearnedParameterTableColumn(final String columnName) {
//            this.columnName = columnName;
//        }
//
//        LearnedParameterTableColumn(final String columnName, final ArtifactState state) {
//            this.columnName = columnName;
//            this.state = state;
//        }
//
//        @Override
//        public String toString() {
//            return columnName;
//        }
//
//        // If used on N or N_ALT, this method will return null
//        public ArtifactState getState() { return state; }
//
//        public static List<LearnedParameterTableColumn> getArtifactStateColumns(){
//            final List<LearnedParameterTableColumn> nonArtifactColumns = Arrays.asList(CONTEXT, REV_COMP, N, N_ALT);
//            return Arrays.stream(LearnedParameterTableColumn.values())
//                    .filter(column -> ! nonArtifactColumns.contains(column))
//                    .collect(Collectors.toList());
//        }
//
//        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
//    }
//}
