# Cotoneaster_review

Here is a brief overview of what is in the code.

I conducted a very basic cleaning of the gbif records (there is a section dedicated to this in the code) where I checked for records in country capitals, country centroids, equal longitude and latitudes, all zeros, records in the oceans, etc. It is obviously possible to refine it even further if you want to.

I created a translation table using the WCVP using all names from gbif and glonaf. I then standardized the nomenclatures of these databases using the translation table. 

For each record in gbif, I added the standardized taxon name (including wcvp id and authors), the botanical country in which it was found (following the wgsrpd level 3 (Brummitt 2001), and it’s status according to the WCVP and GloNAF. 

In the summary table, I additionally summarised the number of records of each species in each area and the earliest and latest years in which it was recorded.
