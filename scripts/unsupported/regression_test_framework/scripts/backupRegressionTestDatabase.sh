#!/usr/bin/env bash

echo "About to dump DSPRegressionTesting database..."
mysqldump -h mysql-prd2.broadinstitute.org -u DSPRegressionTesting -p DSPRegressionTesting > DSPRegressionTesting.$(date +%Y%m%dT%H%M%S).sql

