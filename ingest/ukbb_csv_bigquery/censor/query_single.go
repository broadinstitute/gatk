package main

import (
	"fmt"

	"cloud.google.com/go/bigquery"
	"google.golang.org/api/iterator"
)

func BigQuerySingleFieldFirst(wbq *WrappedBigQuery, fieldID int64) (map[int64]string, error) {
	out := make(map[int64]string)

	query := wbq.Client.Query(fmt.Sprintf(`SELECT * 
FROM %s.phenotype
WHERE 1=1
AND FieldID=@FieldID
ORDER BY instance ASC, array_idx ASC

-- Uncomment for testing
-- ORDER BY sample_id DESC
-- LIMIT 10
`, wbq.Database))

	query.QueryConfig.Parameters = append(query.QueryConfig.Parameters, []bigquery.QueryParameter{
		{Name: "FieldID", Value: fieldID},
	}...)

	itr, err := query.Read(wbq.Context)
	if err != nil {
		return nil, err
	}
	for {
		var values SamplePheno
		err := itr.Next(&values)
		if err == iterator.Done {
			break
		}
		if err != nil {
			return nil, err
		}

		// Take only the first, since we use this for things like enrollment
		// date. If someone came to a follow-up visit, we don't want to say that
		// they "enrolled" at the time of their follow-up, for example. Relies
		// on sort order specified above in the query.
		if _, exists := out[values.SampleID]; exists {
			continue
		}
		out[values.SampleID] = values.Value
	}

	return out, nil
}
