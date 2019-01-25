# Database Naming Conventions:

1. Use lowercase for all names.
2. Separate names by underscores.
3. Do not use quotes in names.
4. Use descriptive names.  Data types are not names (do not use `timestamp`, instead use `file_created_time`).
5. Do not use reserved words as names.
6. Use appropriate abbreviations.  Do not use `first_nm`, instead use `first_name`.  Use `table_idx`, not `table_index`.
7. Use _singular_ names, not plural names, for all fields and tables.
8. All tables should have a _Primary Key_ named `id` that is `auto increment` `unsigned` `int` or similar.


When a table contains a column that is a foreign key, I just copy the column name of that key from whatever table it came from. For example, say table `foo_bar` has the FK `foo_id` (where `foo_id` is the _Primary Key_ of `foo`).
When defining FKs to enforce referencial integrity, I use the following: `tablename_fk_columnname` (e.g. furthering example 3, it would be `foo_bar_foo_id`). Since this is a table name/column name combination, it is guaranteed to be unique within the database.
