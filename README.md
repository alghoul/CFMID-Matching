# CFMID-Matching
A pandas/SQL driven package that parses MS/MS data from MGF files and peforms a spectrum match with CFMID
predicted spectra stored in a SQL database in search for the corresponding chemical. The matches are performed using a simple peak match
algorithm and the scoring is done based on the coine dot product algorithm. 
