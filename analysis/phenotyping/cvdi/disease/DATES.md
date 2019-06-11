# Dates in the UK Biobank

## Attended assessment center

* Date FieldID 53

Useful for:

1. Defining threshold date for incidence
1. Defining dates for things that don't otherwise have an associated date

## Birth

* Date FieldID 34: Year of birth
* Date FieldID 52: Month of birth
* Date Field 33: birth date (*Note: this field is restricted due to its precision*)

## Lost to follow-up

* Date FieldID 191

## Died

* Date FieldID 40000

## ICD10

* Date FieldID ==> derived from HESIN (Hospital Episode Statistics data in Showcase)
* Main ICD10: 41202
* Secondary ICD10: 41204
* ICD10 Primary Cause of Death: 40001
* ICD10 Secondary Cause of Death: 40002

## ICD9

* Date FieldID ==> derived from HESIN
* Main ICD9: 41203
* Secondary ICD9: 41205

## Operation (OPCS4)

* Date FieldID ==> derived from HESIN
* Main OPCS4: 41200
* Secondary OPCS4: 41210
* Self-reported:
  * FieldID: 20004
  * *float32 Year*: 20010 (need to truncate and add month/day)

## Special cases

* Myocardial infarction:
  * FieldID: 42001
  * Date: 42000
* Non-cancer illness:
  * FieldID: 20002
  * Date: 20008
* Cancer:
  * FieldID: 20001
  * Date: 20006
