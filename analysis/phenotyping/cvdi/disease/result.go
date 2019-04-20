package main

import (
	"fmt"
	"os"
	"strconv"
	"strings"
	"time"

	"github.com/carbocation/pfx"

	"cloud.google.com/go/bigquery"
	"google.golang.org/api/iterator"
)

type Result struct {
	SampleID                int64                `bigquery:"sample_id"`
	HasDisease              bigquery.NullInt64   `bigquery:"has_disease"`
	IncidentDisease         bigquery.NullInt64   `bigquery:"incident_disease"`
	PrevalentDisease        bigquery.NullInt64   `bigquery:"prevalent_disease"`
	PhenotypeDateCensor     bigquery.NullDate    `bigquery:"date_censor"`
	RoughPhenotypeAgeCensor bigquery.NullFloat64 `bigquery:"age_censor"` // Note: just uses days/365. Don't use.
	BirthDate               bigquery.NullDate    `bigquery:"birthdate"`
	EnrollDate              bigquery.NullDate    `bigquery:"enroll_date"`
	EnrollAge               bigquery.NullFloat64 `bigquery:"enroll_age"`
	HasDied                 bigquery.NullInt64   `bigquery:"has_died"`
	DeathDate               bigquery.NullDate    `bigquery:"death_date"`
	DeathAge                bigquery.NullFloat64 `bigquery:"death_age"`
	ComputedDate            bigquery.NullDate    `bigquery:"computed_date"`
	MissingFields           bigquery.NullString  `bigquery:"missing_fields"`
}

// PhenotypeAgeCensor computes a proper (leap-year aware) age at phenotype onset
// (or censoring), rather than the heuristic age extracted from the database
// which assumes each year to be the same length.
func (r Result) PhenotypeAgeCensor() (bigquery.NullFloat64, error) {
	if !r.BirthDate.Valid || !r.PhenotypeDateCensor.Valid {
		return bigquery.NullFloat64{}, nil
	}

	birthTime, err := time.Parse("2006-01-02", r.BirthDate.Date.String())
	if err != nil {
		return bigquery.NullFloat64{}, fmt.Errorf("Error parsing BirthDate '%s': %s", r.BirthDate, err.Error())
	}

	phenoTime, err := time.Parse("2006-01-02", r.PhenotypeDateCensor.Date.String())
	if err != nil {
		return bigquery.NullFloat64{}, fmt.Errorf("Error parsing PhenotypeDateCensor '%s': %s", r.PhenotypeDateCensor, err.Error())
	}

	stringYears := TimesToFractionalYears(birthTime, phenoTime)

	floatYears, err := strconv.ParseFloat(stringYears, 64)
	if err != nil {
		return bigquery.NullFloat64{}, fmt.Errorf("Error parsing duration from BirthDate '%s' and PhenotypeDateCensor '%s' : %s. Setting phenotype onset date to birthdate (age 0)", r.BirthDate, r.PhenotypeDateCensor, err.Error())
	}

	return bigquery.NullFloat64{Float64: floatYears, Valid: true}, nil
}

// DeathAgeCensor computes a proper (leap-year aware) age at death (or death
// censoring), rather than the heuristic age extracted from the database which
// assumes each year to be the same length.
func (r Result) DeathAgeCensor() (bigquery.NullFloat64, error) {
	if !r.BirthDate.Valid || !r.DeathDate.Valid {
		return bigquery.NullFloat64{}, nil
	}

	birthTime, err := time.Parse("2006-01-02", r.BirthDate.Date.String())
	if err != nil {
		return bigquery.NullFloat64{}, fmt.Errorf("Error parsing BirthDate '%s': %s", r.BirthDate, err.Error())
	}

	deathTime, err := time.Parse("2006-01-02", r.DeathDate.Date.String())
	if err != nil {
		return bigquery.NullFloat64{}, fmt.Errorf("Error parsing DeathDate '%s': %s", r.DeathDate, err.Error())
	}

	stringYears := TimesToFractionalYears(birthTime, deathTime)

	floatYears, err := strconv.ParseFloat(stringYears, 64)
	if err != nil {
		return bigquery.NullFloat64{}, fmt.Errorf("Error parsing duration from BirthDate '%s' and DeathDate '%s' : %s. Setting death onset date to birthdate (age 0)", r.BirthDate, r.DeathDate, err.Error())
	}

	return bigquery.NullFloat64{Float64: floatYears, Valid: true}, nil
}

func ExecuteQuery(BQ *WrappedBigQuery, query *bigquery.Query, diseaseName string, missingFields []string) error {
	itr, err := query.Read(BQ.Context)
	if err != nil {
		return pfx.Err(fmt.Sprint(err.Error(), query.Parameters))
	}
	todayDate := time.Now().Format("2006-01-02")
	missing := strings.Join(missingFields, ",")
	fmt.Fprintf(STDOUT, "disease\tsample_id\thas_disease\tincident_disease\tprevalent_disease\tcensor_date\tcensor_age\tbirthdate\tenroll_date\tenroll_age\thas_died\tdeath_censor_date\tdeath_censor_age\tcensor_computed_date\tcensor_missing_fields\tcomputed_date\tmissing_fields\n")
	for {
		var r Result
		err := itr.Next(&r)
		if err == iterator.Done {
			break
		}
		if err != nil {
			return pfx.Err(err)
		}

		censoredPhenoAge, err := r.PhenotypeAgeCensor()
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s: setting age_censor to null for %d (birthdate %s phenotype date %s) because of error: %s\n", diseaseName, r.SampleID, r.BirthDate, r.PhenotypeDateCensor, err.Error())

			// UK Biobank uses impossible values (e.g., 1900-01-01) to indicate that
			// the date is not known. See, e.g., FieldID 42000. This does not mean
			// that the value is illegal, so it shouldn't be null. Instead, it should
			// be some legal value. Here, we set the age of incidence to be 0 years,
			// and we set the date of incidence to be the birthdate.
			censoredPhenoAge = bigquery.NullFloat64{Float64: 0, Valid: true}
			r.PhenotypeDateCensor = r.BirthDate
		}

		censoredDeathAge, err := r.DeathAgeCensor()
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s: setting age_censor to null for %d (birthdate %s deat date %s) because of error: %s\n", diseaseName, r.SampleID, r.BirthDate, r.DeathDate, err.Error())

			// UK Biobank uses impossible values (e.g., 1900-01-01) to indicate that
			// the date is not known. See, e.g., FieldID 42000. This does not mean
			// that the value is illegal, so it shouldn't be null. Instead, it should
			// be some legal value. Here, we set the age of incidence to be 0 years,
			// and we set the date of incidence to be the birthdate.
			censoredDeathAge = bigquery.NullFloat64{Float64: 0, Valid: true}
			r.DeathDate = r.BirthDate
		}

		fmt.Fprintf(STDOUT, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
			diseaseName, r.SampleID, NA(r.HasDisease), NA(r.IncidentDisease), NA(r.PrevalentDisease), NA(r.PhenotypeDateCensor), NA(censoredPhenoAge), NA(r.BirthDate), NA(r.EnrollDate), NA(r.EnrollAge), NA(r.HasDied), NA(r.DeathDate), NA(censoredDeathAge), NA(r.ComputedDate), NA(r.MissingFields), todayDate, missing)
	}

	return nil
}

// NA emits an empty string instead of "NULL" since this plays better with
// BigQuery
func NA(input interface{}) interface{} {
	invalid := ""

	switch v := input.(type) {
	case bigquery.NullInt64:
		if !v.Valid {
			return invalid
		}
	case bigquery.NullFloat64:
		if !v.Valid {
			return invalid
		}
	case bigquery.NullString:
		if !v.Valid {
			return invalid
		}
	case bigquery.NullDate:
		if !v.Valid {
			return invalid
		}
	}

	return input
}

func BuildQuery(BQ *WrappedBigQuery, tabs *TabFile, displayQuery bool) (*bigquery.Query, error) {
	query := rawQuery(BQ)

	params := []bigquery.QueryParameter{}

	standardPart := "AND FALSE"
	if len(tabs.Include.Standard)+len(tabs.Exclude.Standard) > 0 {
		standardPart = "AND p.FieldID IN UNNEST(@StandardFieldIDs)"
		fieldIDs := tabs.AllStandardFields()
		params = append(params, bigquery.QueryParameter{Name: "StandardFieldIDs", Value: fieldIDs})
	}

	includePart := ""
	includedValues := tabs.AllIncluded()
	if len(includedValues) > 0 {
		i := 0
		for _, v := range includedValues {
			includePart = includePart + fmt.Sprintf("\nOR (hd.FieldID = @IncludeParts%d AND hd.value IN UNNEST(@IncludeParts%d) )", i, i+1)
			params = append(params, bigquery.QueryParameter{Name: fmt.Sprintf("IncludeParts%d", i), Value: v.FieldID})
			params = append(params, bigquery.QueryParameter{Name: fmt.Sprintf("IncludeParts%d", i+1), Value: v.FormattedValues()})
			i += 2
		}
	}

	excludePart := ""
	exludedValues := tabs.AllExcluded()
	if len(exludedValues) > 0 {
		i := 0
		for _, v := range exludedValues {
			excludePart = excludePart + fmt.Sprintf("\nOR (excl.FieldID = @ExcludeParts%d AND excl.value IN UNNEST(@ExcludeParts%d) )", i, i+1)
			params = append(params, bigquery.QueryParameter{Name: fmt.Sprintf("ExcludeParts%d", i), Value: v.FieldID})
			params = append(params, bigquery.QueryParameter{Name: fmt.Sprintf("ExcludeParts%d", i+1), Value: v.FormattedValues()})
			i += 2
		}
	}

	// Assemble all the parts
	query = fmt.Sprintf(query, standardPart, includePart, excludePart)

	if displayQuery {
		fmt.Println(query)
		for _, v := range params {
			if x, ok := v.Value.([]string); ok {
				fmt.Printf("%v: (\"%s\")\n", v.Name, strings.Join(x, `","`))
				continue
			}
			fmt.Printf("%v: %v\n", v.Name, v.Value)
		}
		return nil, nil
	}

	// Generate the bigquery query object, but don't call it
	bqQuery := BQ.Client.Query(query)
	bqQuery.QueryConfig.Parameters = append(bqQuery.QueryConfig.Parameters, params...)

	return bqQuery, nil
}

// TODO: create a `has_died` field and apply the censor table's
// death_censor_date properly. Rename death_date to death_censor_date. Rename
// death_age to death_censor_age.
//
// TODO: Resolve age_censor vs enroll_age. Choose one or the other (likely the
// latter, so you end up with enroll_age, censor_age, death_censor_age).
func rawQuery(BQ *WrappedBigQuery) string {
	return fmt.Sprintf(`
	WITH undated_fields AS (
		SELECT p.sample_id, p.FieldID, p.value, MIN(SAFE.PARSE_DATE("%%%%E4Y-%%%%m-%%%%d", denroll.value)) first_date
		FROM `+"`%s`"+`.phenotype p
		JOIN `+"`%s`"+`.phenotype denroll ON denroll.FieldID=53 AND denroll.sample_id=p.sample_id AND denroll.instance = 0 AND denroll.array_idx = 0
		WHERE TRUE
		  %%s
		GROUP BY p.FieldID, p.sample_id, p.value
	  ), included_only AS (
	  SELECT sample_id, has_disease, incident_disease, prevalent_disease, date_censor, DATE_DIFF(date_censor,birthdate, DAY)/365.0 age_censor, 
        birthdate, enroll_date, enroll_age, death_date, death_age, computed_date, missing_fields
	  FROM (
		SELECT c.sample_id, 
		  CASE 
			WHEN MIN(hd.first_date) IS NOT NULL THEN 1
			ELSE 0
		  END has_disease,
		  CASE 
			WHEN MIN(hd.first_date) > MIN(c.enroll_date) THEN 1
			WHEN MIN(hd.first_date) IS NOT NULL THEN NULL
			ELSE 0
		  END incident_disease,
		  CASE 
			WHEN MIN(hd.first_date) > MIN(c.enroll_date) THEN 0
			WHEN MIN(hd.first_date) IS NOT NULL THEN 1
			ELSE 0
		  END prevalent_disease,
		  CASE 
			WHEN MIN(hd.first_date) IS NOT NULL THEN MIN(hd.first_date)
			ELSE MIN(c.phenotype_censor_date)
		  END date_censor,
		  MIN(c.birthdate) birthdate,
		  MIN(c.enroll_date) enroll_date,
		  MIN(c.enroll_age) enroll_age,
		  MIN(c.death_date) death_date,
		  MIN(c.death_age) death_age,
		  MIN(c.computed_date) computed_date,
		  MIN(c.missing_fields) missing_fields
		FROM `+"`%s.censor`"+` c
		LEFT OUTER JOIN (
			SELECT * FROM `+"`%s.materialized_hesin_dates`"+`
			UNION DISTINCT
			SELECT * FROM `+"`%s.materialized_special_dates`"+`
			UNION DISTINCT
			SELECT * FROM undated_fields
		  ) hd ON c.sample_id=hd.sample_id
		  AND (
			FALSE
			%%s
		  )

		GROUP BY sample_id
	  )
	  ),excluded_only AS (
		SELECT sample_id, has_disease, incident_disease, prevalent_disease, date_censor, DATE_DIFF(date_censor,birthdate, DAY)/365.0 age_censor, 
		  birthdate, enroll_date, enroll_age, death_date, death_age, computed_date, missing_fields
		FROM (
		  SELECT c.sample_id, 
			CASE 
			  WHEN MIN(excl.first_date) IS NOT NULL THEN 1
			  ELSE 0
			END has_disease,
			CASE 
			  WHEN MIN(excl.first_date) > MIN(c.enroll_date) THEN 1
			  WHEN MIN(excl.first_date) IS NOT NULL THEN NULL
			  ELSE 0
			END incident_disease,
			CASE 
			  WHEN MIN(excl.first_date) > MIN(c.enroll_date) THEN 0
			  WHEN MIN(excl.first_date) IS NOT NULL THEN 1
			  ELSE 0
			END prevalent_disease,
			CASE 
			  WHEN MIN(excl.first_date) IS NOT NULL THEN MIN(excl.first_date)
			  ELSE MIN(c.phenotype_censor_date)
			END date_censor,
			MIN(c.birthdate) birthdate,
			MIN(c.enroll_date) enroll_date,
			MIN(c.enroll_age) enroll_age,
			MIN(c.death_date) death_date,
			MIN(c.death_age) death_age,
			MIN(c.computed_date) computed_date,
			MIN(c.missing_fields) missing_fields
		  FROM `+"`%s.censor`"+` c
		  LEFT OUTER JOIN (
			  SELECT * FROM `+"`%s.materialized_hesin_dates`"+`
			  UNION DISTINCT
			  SELECT * FROM `+"`%s.materialized_special_dates`"+`
			  UNION DISTINCT
			  SELECT * FROM undated_fields
			) excl ON c.sample_id=excl.sample_id
			AND (
			  FALSE
			  %%s
			)
		  GROUP BY sample_id
		)
		)
	  
	  SELECT c.sample_id, 
		CASE 
		  -- Enrollment occurred after exclusion:
		  WHEN eo.has_disease = 1 AND SAFE.DATE_DIFF(c.enroll_date, eo.date_censor, DAY) > 0 THEN NULL
		  -- Exclusion occurred after enrollment and prior to disease onset; we will exclude:
		  WHEN eo.has_disease = 1 AND io.has_disease = 1 AND SAFE.DATE_DIFF(io.date_censor,eo.date_censor, DAY) > 0 THEN NULL
		  -- Exclusion occurred after disease onset; we'll allow it:
		  WHEN eo.has_disease = 1 AND io.has_disease = 1 AND SAFE.DATE_DIFF(eo.date_censor,io.date_censor, DAY) > 0 THEN io.has_disease 
			-- Met exclusion but no inclusion
			WHEN eo.has_disease = 1 AND (io.has_disease = 0 OR io.has_disease IS NULL) THEN NULL
			-- Didn't meet exclusion or inclusion means we censor at the date given by UKBB:
		  WHEN io.has_disease IS NULL THEN 0
		  ELSE io.has_disease
		END has_disease, 
		CASE 
		  -- Enrollment occurred after exclusion:
		  WHEN eo.has_disease = 1 AND SAFE.DATE_DIFF(c.enroll_date, eo.date_censor, DAY) > 0 THEN NULL
		  -- Exclusion occurred after enrollment and prior to disease onset; we will exclude:
		  WHEN eo.has_disease = 1 AND io.has_disease = 1 AND SAFE.DATE_DIFF(io.date_censor,eo.date_censor, DAY) > 0 THEN NULL
		  -- Exclusion occurred after disease onset; we'll allow it:
			WHEN eo.has_disease = 1 AND io.has_disease = 1 AND SAFE.DATE_DIFF(eo.date_censor,io.date_censor, DAY) > 0 THEN io.incident_disease
			-- Met exclusion but no inclusion
			WHEN eo.has_disease = 1 AND (io.has_disease = 0 OR io.has_disease IS NULL) THEN NULL
		  -- Didn't meet exclusion or inclusion means we censor at the date given by UKBB:
		  WHEN io.has_disease IS NULL THEN 0
		  ELSE io.incident_disease
		END incident_disease, 
		CASE 
		  -- Enrollment occurred after exclusion:
		  WHEN eo.has_disease = 1 AND SAFE.DATE_DIFF(c.enroll_date, eo.date_censor, DAY) > 0 THEN NULL
		  -- Exclusion occurred after enrollment and prior to disease onset; we will exclude:
		  WHEN eo.has_disease = 1 AND io.has_disease = 1 AND SAFE.DATE_DIFF(io.date_censor,eo.date_censor, DAY) > 0 THEN NULL
		  -- Exclusion occurred after disease onset; we'll allow it:
			WHEN eo.has_disease = 1 AND io.has_disease = 1 AND SAFE.DATE_DIFF(eo.date_censor,io.date_censor, DAY) > 0 THEN io.prevalent_disease
			-- Met exclusion but no inclusion
			WHEN eo.has_disease = 1 AND (io.has_disease = 0 OR io.has_disease IS NULL) THEN NULL
		  -- Didn't meet exclusion or inclusion means we censor at the date given by UKBB:
		  WHEN io.has_disease IS NULL THEN 0
		  ELSE io.prevalent_disease
		END prevalent_disease, 
		CASE 
		  -- Enrollment occurred after exclusion:
		  WHEN eo.has_disease = 1 AND SAFE.DATE_DIFF(c.enroll_date, eo.date_censor, DAY) > 0 THEN eo.date_censor
		  -- Exclusion occurred after enrollment and prior to disease onset; we will exclude:
		  WHEN eo.has_disease = 1 AND io.has_disease = 1 AND SAFE.DATE_DIFF(io.date_censor,eo.date_censor, DAY) > 0 THEN eo.date_censor
		  -- Exclusion occurred after disease onset; we'll allow it:
			WHEN eo.has_disease = 1 AND io.has_disease = 1 AND SAFE.DATE_DIFF(eo.date_censor,io.date_censor, DAY) > 0 THEN io.date_censor
			-- Met exclusion but no inclusion
			WHEN eo.has_disease = 1 AND (io.has_disease = 0 OR io.has_disease IS NULL) THEN eo.date_censor
		  -- Didn't meet exclusion or inclusion means we censor at the date given by UKBB:
		  WHEN io.has_disease IS NULL THEN c.phenotype_censor_date
		  ELSE io.date_censor
		END date_censor, 
		CASE
		  -- Enrollment occurred after exclusion:
		  WHEN eo.has_disease = 1 AND SAFE.DATE_DIFF(c.enroll_date, eo.date_censor, DAY) > 0 THEN DATE_DIFF(eo.date_censor,c.birthdate, DAY)/365.0
		  -- Exclusion occurred after enrollment and prior to disease onset; we will exclude:
		  WHEN eo.has_disease = 1 AND io.has_disease = 1 AND SAFE.DATE_DIFF(io.date_censor,eo.date_censor, DAY) > 0 THEN DATE_DIFF(eo.date_censor,c.birthdate, DAY)/365.0
		  -- Exclusion occurred after disease onset; we'll allow it:
			WHEN eo.has_disease = 1 AND io.has_disease = 1 AND SAFE.DATE_DIFF(eo.date_censor,io.date_censor, DAY) > 0 THEN DATE_DIFF(io.date_censor,c.birthdate, DAY)/365.0 
			-- Met exclusion but no inclusion
			WHEN eo.has_disease = 1 AND (io.has_disease = 0 OR io.has_disease IS NULL) THEN DATE_DIFF(eo.date_censor,c.birthdate, DAY)/365.0
		  -- Didn't meet exclusion or inclusion means we censor at the date given by UKBB:
		  WHEN io.has_disease IS NULL THEN c.phenotype_censor_age
		  ELSE DATE_DIFF(io.date_censor,c.birthdate, DAY)/365.0 
		END age_censor, 
		c.birthdate, 
		c.enroll_date, 
		c.enroll_age, 
		CASE WHEN c.death_date IS NULL THEN 0 ELSE 1 END has_died,
		CASE WHEN c.death_date IS NULL THEN c.death_censor_date ELSE c.death_date END death_date, 
		CASE WHEN c.death_date IS NULL THEN c.death_censor_age ELSE c.death_age END death_age, 
		c.computed_date, 
		c.missing_fields
	  FROM `+"`%s.censor`"+` c
	  LEFT JOIN included_only io ON io.sample_id=c.sample_id
	  LEFT JOIN excluded_only eo ON eo.sample_id=c.sample_id
	  ORDER BY has_disease DESC, incident_disease DESC, age_censor ASC
	  `,
		BQ.Database, BQ.Database, BQ.Database, materializedDB, materializedDB, BQ.Database, materializedDB, materializedDB, BQ.Database)
}
