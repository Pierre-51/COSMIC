COPY gene FROM './/gene.csv' (FORMAT 'csv', force_not_null ('ID_GENE', 'GENE_SYMBOL'), quote '"', delimiter ',', header 1);
COPY transcript FROM './/transcript.csv' (FORMAT 'csv', force_not_null ('ID_TRANSCRIPT', 'ID_GENE'), quote '"', delimiter ',', header 1);
