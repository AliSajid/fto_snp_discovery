# Download Sequences for both assemblies and extract the UCEs

library(tidyverse)
library(reutils)

options(reutils.email = "Ali.Imami@rockets.utoledo.edu",
        reutils.api.key = "a20864ad0756cdd3c694893bedc036a99708")

retmode <- "text"
rettype <- "fasta"

fasta <- efetch(c("NC_000016.9"), db = "nuccore", rettype = rettype, retmode = retmode)

write(fasta$content, "data/NC_000016.9.fasta")

fasta <- efetch(c("NC_000016.10"), db = "nuccore", rettype = rettype, retmode = retmode)

write(fasta$content, "data/NC_000016.10.fasta")
