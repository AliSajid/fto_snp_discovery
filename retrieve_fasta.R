# Download Sequences for both assemblies and extract the UCEs

library(tidyverse)
library(reutils)

options(reutils.email = "Ali.Imami@rockets.utoledo.edu",
        reutils.api.key = "a20864ad0756cdd3c694893bedc036a99708")

retmode <- "text"
rettype <- "fasta"

fasta <- efetch("NC_000016.9", db = "nuccore", rettype = fasta, retmode = retmode)

write(fasta$content, "data/NC_000016.9.fasta")
