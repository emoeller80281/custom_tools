#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Run everything from the script's own directory so that all downloads and the
# generator's intermediate metadata stay contained here instead of leaking into
# whatever directory the script happens to be invoked from.
cd "${SCRIPT_DIR}"

OUTPUT_DIR="${SCRIPT_DIR}/python_lib"
# get a copy of the Datasets OpenAPI v2 Specification
wget https://www.ncbi.nlm.nih.gov/datasets/docs/v2/openapi3/openapi3.docs.yaml
# Get the OpenAPI library generator (a Java jar file)
wget https://repo1.maven.org/maven2/org/openapitools/openapi-generator-cli/7.17.0/openapi-generator-cli-7.17.0.jar -O openapi-generator-cli.jar
# Create the datasets OpenAPI library for python in the directory ${OUTPUT_DIR}.
# For more info see: https://openapi-generator.tech/docs/usage/#generate
java -jar openapi-generator-cli.jar generate \
    -g python \
    -i openapi3.docs.yaml \
    -o ${OUTPUT_DIR} \
    --package-name "ncbi.datasets.openapi" \
    --additional-properties=pythonAttrNoneIfUnset=true,projectName="ncbi-datasets-pylib"

1