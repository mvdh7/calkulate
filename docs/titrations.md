# Individual titrations

To import a single titration dataset ready to work with in Calkulate, just provide the relevant row of your titration table (or the whole table if it only contains one row):

    :::python
    import calkulate as calk
    titration = calk.Titration(titration_table_row)

This imports the data from the titration data file and stores it in a `Titration` object along with the relevant information from the titration table.  Some of these values are available at the top level:

!!! info "Top-level titration properties"

    Some of these come directly from the titration table columns of the same name:
    
    `titration.file_name`, `titration.file_path`

The remaining data are separated into several categories:

  * `titration.analyte` for properties of the analyte.
  * `titration.titrant` for properties of the titrant.
  * `titration.mixture` for properties of the titrant-analyte mixture.
  * `titration.settings` for... settings.

## The analyte

The attributes of `titration.analyte` are all single scalar values.  At first, these include:

  

## The titrant

Some attributes of `titration.titrant` are single scalar values (e.g. its molinity), while others are arrays (e.g. the amount of it added to the analyte).

## The mixture

The attributes of `titration.mixture` are all arrays.

## The settings
