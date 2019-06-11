# CVDI/Disease
CVDI/Disease computes derived phenotypes, with incidence, prevalence, death, and censoring, from tabfiles (Seung Hoan's phenotype definition format). As a brief review, tabfiles contain 3 columns: Field, Coding, and exclude. E.g.:

```
Field	Coding	exclude
20002	1076,1079	0
```

*Field* is the UK Biobank FieldID (e.g., `phenotype.FieldID`)

*Coding* is the UK Biobank value (e.g., `phenotype.value` or `coding.coding`)

*Exclude* is whether the row is an exclusion criterion (1) or not (0)

# Install updated dependencies
`go get -u`

# Build (examples)
`go build -o cvdidisease.osx *.go`

`GOOS=linux go build -o cvdidisease.linux *.go`

# Database dependencies
This requires the materialized tables (defined in the SQL files in this directory) to exist in tables with the same name as their filename (except the suffix).