#!/bin/bash

# Script to validate WDL files using womtool
# Usage: ./validate_wdl.sh [wdl_file]

set -e

WOMTOOL_JAR="/Users/hatcher/Documents/tools/wdl_verification/womtool-83.jar"
WDL_FILE="${1:-GVSRealignVetIndels.wdl}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ ! -f "$WOMTOOL_JAR" ]]; then
    echo "ERROR: womtool jar not found at $WOMTOOL_JAR"
    exit 1
fi

if [[ ! -f "$SCRIPT_DIR/$WDL_FILE" ]]; then
    echo "ERROR: WDL file not found at $SCRIPT_DIR/$WDL_FILE"
    exit 1
fi

echo "Validating WDL file: $WDL_FILE"
echo "Using womtool: $WOMTOOL_JAR"
echo "----------------------------------------"

# Validate the WDL syntax
echo "1. Validating WDL syntax..."
java -jar "$WOMTOOL_JAR" validate "$SCRIPT_DIR/$WDL_FILE"

# Generate inputs template
echo ""
echo "2. Generating inputs template..."
java -jar "$WOMTOOL_JAR" inputs "$SCRIPT_DIR/$WDL_FILE" > "$SCRIPT_DIR/${WDL_FILE%.wdl}.inputs.template.json"
echo "Inputs template saved to: ${WDL_FILE%.wdl}.inputs.template.json"

# If we have example inputs, validate them
for inputs_file in "$SCRIPT_DIR"/${WDL_FILE%.wdl}*.inputs.json; do
    if [[ -f "$inputs_file" ]]; then
        echo ""
        echo "3. Validating inputs file: $(basename "$inputs_file")"
        java -jar "$WOMTOOL_JAR" validate "$SCRIPT_DIR/$WDL_FILE" --inputs "$inputs_file"
    fi
done

echo ""
echo "✅ WDL validation completed successfully!"
