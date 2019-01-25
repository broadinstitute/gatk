-- MySQL dump 10.13  Distrib 5.7.19, for osx10.12 (x86_64)
--
-- Host: mysql-prd2.broadinstitute.org    Database: DSPRegressionTesting
-- ------------------------------------------------------
-- Server version	5.7.17-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `Analysis`
--

DROP TABLE IF EXISTS `Analysis`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Analysis` (
  `idAnalysis` int(11) unsigned NOT NULL,
  `version` varchar(1024) COLLATE utf8mb4_unicode_ci NOT NULL,
  `configuration` json DEFAULT NULL,
  `baseData` int(11) unsigned NOT NULL,
  `scenarioOutputForComparison` int(11) unsigned NOT NULL,
  PRIMARY KEY (`idAnalysis`),
  UNIQUE KEY `idAnalysis_UNIQUE` (`idAnalysis`),
  KEY `Analysis_scenarioOutputForComparison_idx` (`scenarioOutputForComparison`),
  KEY `Analysis_baseData_idx` (`baseData`),
  CONSTRAINT `Analysis_baseData` FOREIGN KEY (`baseData`) REFERENCES `OutputFiles` (`idOutputFiles`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `Analysis_scenarioOutputForComparison` FOREIGN KEY (`scenarioOutputForComparison`) REFERENCES `OutputFiles` (`idOutputFiles`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Analysis`
--

LOCK TABLES `Analysis` WRITE;
/*!40000 ALTER TABLE `Analysis` DISABLE KEYS */;
/*!40000 ALTER TABLE `Analysis` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `FileTypes`
--

DROP TABLE IF EXISTS `FileTypes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `FileTypes` (
  `idFileTypes` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `type` varchar(45) COLLATE utf8mb4_unicode_ci NOT NULL COMMENT 'The type of the file.',
  PRIMARY KEY (`idFileTypes`,`type`),
  UNIQUE KEY `idFileTypes_UNIQUE` (`idFileTypes`),
  UNIQUE KEY `FileTypes_type_UNIQUE` (`type`)
) ENGINE=InnoDB AUTO_INCREMENT=9 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `FileTypes`
--

LOCK TABLES `FileTypes` WRITE;
/*!40000 ALTER TABLE `FileTypes` DISABLE KEYS */;
INSERT INTO `FileTypes` VALUES (1,'EXOME_READS'),(2,'GENOME_READS'),(3,'SOMATIC_VCF'),(4,'GERMLINE_VCF'),(5,'EXOME_READS_INDEX'),(6,'GENOME_READS_INDEX'),(7,'SOMATIC_VCF_INDEX'),(8,'GERMLINE_VCF_INDEX');
/*!40000 ALTER TABLE `FileTypes` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `InputFiles`
--

DROP TABLE IF EXISTS `InputFiles`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `InputFiles` (
  `idInputFiles` int(11) unsigned NOT NULL AUTO_INCREMENT COMMENT 'Unique ID of each file.',
  `path` varchar(2048) COLLATE utf8mb4_unicode_ci NOT NULL COMMENT 'Path to the file.',
  `type` varchar(45) COLLATE utf8mb4_unicode_ci NOT NULL COMMENT 'Type of each file.',
  `md5sum` char(32) COLLATE utf8mb4_unicode_ci NOT NULL,
  PRIMARY KEY (`idInputFiles`),
  UNIQUE KEY `idInputFiles_UNIQUE` (`idInputFiles`),
  UNIQUE KEY `InputFiles_md5sum_UNIQUE` (`md5sum`),
  KEY `InputFiles_fileType_idx` (`type`),
  CONSTRAINT `InputFiles_fileType` FOREIGN KEY (`type`) REFERENCES `FileTypes` (`type`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `InputFiles`
--

LOCK TABLES `InputFiles` WRITE;
/*!40000 ALTER TABLE `InputFiles` DISABLE KEYS */;
/*!40000 ALTER TABLE `InputFiles` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `Metrics`
--

DROP TABLE IF EXISTS `Metrics`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Metrics` (
  `idMetrics` int(11) unsigned NOT NULL,
  `metricTableName` varchar(1024) COLLATE utf8mb4_unicode_ci NOT NULL,
  `concreteMetricID` int(11) NOT NULL,
  PRIMARY KEY (`idMetrics`),
  UNIQUE KEY `idMetrics_UNIQUE` (`idMetrics`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Metrics`
--

LOCK TABLES `Metrics` WRITE;
/*!40000 ALTER TABLE `Metrics` DISABLE KEYS */;
/*!40000 ALTER TABLE `Metrics` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `MetricsEvaluator`
--

DROP TABLE IF EXISTS `MetricsEvaluator`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `MetricsEvaluator` (
  `idMetricsEvaluator` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `scenarioInfo` int(11) unsigned NOT NULL,
  `metric` int(11) unsigned NOT NULL,
  `absoluteAllowableMaxValue` double DEFAULT NULL,
  `absoluteAllowableMinValue` double DEFAULT NULL,
  `pairwiseAllowableIncrease` double DEFAULT NULL,
  `pairwiseAllowableDecrease` double DEFAULT NULL,
  `isPairwiseMeasurePercentage` tinyint(4) DEFAULT NULL,
  PRIMARY KEY (`idMetricsEvaluator`),
  UNIQUE KEY `idMetricsEvaluator_UNIQUE` (`idMetricsEvaluator`),
  KEY `MetricsEvaluator_metricID_idx` (`metric`),
  KEY `MetricsEvaluator_scenarioInfoID_idx` (`scenarioInfo`),
  CONSTRAINT `MetricsEvaluator_metricID` FOREIGN KEY (`metric`) REFERENCES `Metrics` (`idMetrics`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `MetricsEvaluator_scenarioInfoID` FOREIGN KEY (`scenarioInfo`) REFERENCES `TestScenarioInfo` (`idTestScenarioInfo`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `MetricsEvaluator`
--

LOCK TABLES `MetricsEvaluator` WRITE;
/*!40000 ALTER TABLE `MetricsEvaluator` DISABLE KEYS */;
/*!40000 ALTER TABLE `MetricsEvaluator` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `MetricsEvaluatorData`
--

DROP TABLE IF EXISTS `MetricsEvaluatorData`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `MetricsEvaluatorData` (
  `metricsEvaluatorID` int(11) unsigned NOT NULL,
  `dataID` int(11) unsigned NOT NULL,
  PRIMARY KEY (`metricsEvaluatorID`,`dataID`),
  KEY `MetricsEvaluatorData_dataID_idx` (`dataID`),
  CONSTRAINT `MetricsEvaluatorData_dataID` FOREIGN KEY (`dataID`) REFERENCES `OutputFiles` (`idOutputFiles`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `MetricsEvaluatorData_metricsEvaluatorID` FOREIGN KEY (`metricsEvaluatorID`) REFERENCES `MetricsEvaluator` (`idMetricsEvaluator`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `MetricsEvaluatorData`
--

LOCK TABLES `MetricsEvaluatorData` WRITE;
/*!40000 ALTER TABLE `MetricsEvaluatorData` DISABLE KEYS */;
/*!40000 ALTER TABLE `MetricsEvaluatorData` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `OutputFiles`
--

DROP TABLE IF EXISTS `OutputFiles`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `OutputFiles` (
  `idOutputFiles` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `type` varchar(45) COLLATE utf8mb4_unicode_ci NOT NULL,
  `path` varchar(2048) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `timeCreated` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `scenario` int(11) unsigned NOT NULL,
  `md5sum` char(32) COLLATE utf8mb4_unicode_ci NOT NULL,
  PRIMARY KEY (`idOutputFiles`),
  UNIQUE KEY `idOutputFiles_UNIQUE` (`idOutputFiles`),
  UNIQUE KEY `OutputFiles_md5sum_UNIQUE` (`md5sum`),
  KEY `OutputFiles_fileType_idx` (`type`),
  KEY `OutputFiles_scenarioID` (`scenario`),
  CONSTRAINT `OutputFiles_fileType` FOREIGN KEY (`type`) REFERENCES `FileTypes` (`type`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `OutputFiles_scenarioID` FOREIGN KEY (`scenario`) REFERENCES `TestScenario` (`idTestScenarioRun`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `OutputFiles`
--

LOCK TABLES `OutputFiles` WRITE;
/*!40000 ALTER TABLE `OutputFiles` DISABLE KEYS */;
/*!40000 ALTER TABLE `OutputFiles` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `TestReport`
--

DROP TABLE IF EXISTS `TestReport`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `TestReport` (
  `idTestReport` int(11) unsigned NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`idTestReport`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `TestReport`
--

LOCK TABLES `TestReport` WRITE;
/*!40000 ALTER TABLE `TestReport` DISABLE KEYS */;
/*!40000 ALTER TABLE `TestReport` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `TestReportData`
--

DROP TABLE IF EXISTS `TestReportData`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `TestReportData` (
  `reportID` int(11) unsigned NOT NULL,
  `dataID` int(11) unsigned NOT NULL,
  PRIMARY KEY (`reportID`,`dataID`),
  KEY `TestReportData_dataID_idx` (`dataID`),
  CONSTRAINT `TestReportData_dataID` FOREIGN KEY (`dataID`) REFERENCES `OutputFiles` (`idOutputFiles`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `TestReportData_reportID` FOREIGN KEY (`reportID`) REFERENCES `TestReport` (`idTestReport`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `TestReportData`
--

LOCK TABLES `TestReportData` WRITE;
/*!40000 ALTER TABLE `TestReportData` DISABLE KEYS */;
/*!40000 ALTER TABLE `TestReportData` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `TestReportMetrics`
--

DROP TABLE IF EXISTS `TestReportMetrics`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `TestReportMetrics` (
  `reportId` int(11) unsigned NOT NULL,
  `metricId` int(11) unsigned NOT NULL,
  KEY `TestReportMetrics_reportID_idx` (`reportId`),
  KEY `TestReportMetrics_metricID_idx` (`metricId`),
  CONSTRAINT `TestReportMetrics_metricID` FOREIGN KEY (`metricId`) REFERENCES `Metrics` (`idMetrics`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `TestReportMetrics_reportID` FOREIGN KEY (`reportId`) REFERENCES `TestReport` (`idTestReport`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `TestReportMetrics`
--

LOCK TABLES `TestReportMetrics` WRITE;
/*!40000 ALTER TABLE `TestReportMetrics` DISABLE KEYS */;
/*!40000 ALTER TABLE `TestReportMetrics` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `TestScenario`
--

DROP TABLE IF EXISTS `TestScenario`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `TestScenario` (
  `idTestScenarioRun` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `scenarioInfo` int(11) unsigned NOT NULL,
  `configuration` json NOT NULL,
  `inputFile` int(11) unsigned NOT NULL,
  `tool` int(11) unsigned NOT NULL,
  PRIMARY KEY (`idTestScenarioRun`),
  UNIQUE KEY `idTestScenario_UNIQUE` (`idTestScenarioRun`),
  KEY `TestScenario_fileID_idx` (`inputFile`),
  KEY `TestScenario_toolID_idx` (`tool`),
  KEY `TestScenario_scenarioInfoID_idx` (`scenarioInfo`),
  CONSTRAINT `TestScenario_fileID` FOREIGN KEY (`inputFile`) REFERENCES `InputFiles` (`idInputFiles`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `TestScenario_scenarioInfoID` FOREIGN KEY (`scenarioInfo`) REFERENCES `TestScenarioInfo` (`idTestScenarioInfo`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `TestScenario_toolID` FOREIGN KEY (`tool`) REFERENCES `Tools` (`idTools`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `TestScenario`
--

LOCK TABLES `TestScenario` WRITE;
/*!40000 ALTER TABLE `TestScenario` DISABLE KEYS */;
/*!40000 ALTER TABLE `TestScenario` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `TestScenarioInfo`
--

DROP TABLE IF EXISTS `TestScenarioInfo`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `TestScenarioInfo` (
  `idTestScenarioInfo` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(512) COLLATE utf8mb4_unicode_ci NOT NULL,
  `description` varchar(2048) COLLATE utf8mb4_unicode_ci NOT NULL,
  PRIMARY KEY (`idTestScenarioInfo`),
  UNIQUE KEY `idTestScenarioInfo_UNIQUE` (`idTestScenarioInfo`),
  UNIQUE KEY `name_UNIQUE` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `TestScenarioInfo`
--

LOCK TABLES `TestScenarioInfo` WRITE;
/*!40000 ALTER TABLE `TestScenarioInfo` DISABLE KEYS */;
/*!40000 ALTER TABLE `TestScenarioInfo` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `Tools`
--

DROP TABLE IF EXISTS `Tools`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Tools` (
  `idTools` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `version` varchar(1024) COLLATE utf8mb4_unicode_ci NOT NULL,
  `name` varchar(1024) COLLATE utf8mb4_unicode_ci NOT NULL,
  `dateReleased` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `wdl` varchar(1024) COLLATE utf8mb4_unicode_ci NOT NULL,
  `wdlChecksum` char(32) COLLATE utf8mb4_unicode_ci NOT NULL,
  PRIMARY KEY (`idTools`),
  UNIQUE KEY `idTools_UNIQUE` (`idTools`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Tools`
--

LOCK TABLES `Tools` WRITE;
/*!40000 ALTER TABLE `Tools` DISABLE KEYS */;
/*!40000 ALTER TABLE `Tools` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2019-03-11 17:13:19
