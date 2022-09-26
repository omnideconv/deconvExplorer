# DeconvExplorer

[![R-CMD-check](https://github.com/omnideconv/DeconvExplorer/actions/workflows/R-CMD-check.yaml/badge.svg?branch=master)](https://github.com/omnideconv/DeconvExplorer/actions/workflows/R-CMD-check.yaml)

Shiny App for the omnideconv package. #toBeContinued


## Installation 
```
devtools::install_github("omnideconv/DeconvExplorer")
# install.packages("DeconvExplorer") # for later
```

## Usage 

```
DeconvExplorer()

# start app and upload data in one step 
DeconvExplorer(bulkExpressionData, SingleCellData, CellTypeAnnotations, BatchIDs)
```

The Following formats are allowed for uploading data: 
- txt
- csv 
- tsv 
- rds

For more information about DeconvExplorer and omnideconv please visit https://github.com/omnideconv
