## Generally supported data types:

- rds
- csv
- tsv
- txt

### Data requirements:

#### Single cell RNA-seq data

- **Genes** x **Cells** matrix
- Counts are **not log-transformed**
- Rownames (gene names) are provided in the same format as in the bulk RNA-seq data, for instance HGNC symbols

<img src="../www/sc.png" alt="sc_image" width="800"/>

#### Cell type annotations

- Vector containing cell type annotations
- Annotations are in the same order as the columns of the single cell matrix

<img src="../www/cell_anno.png" alt="anno_image" width="100"/>

#### Batch ids

- Vector containing batch ids, so sample or patient ids
- Ids are in the same order as the columns of the single cell matrix
- This is only necessary for Bisque, MuSiC and SCDC

<img src="../www/batch.png" alt="batch_image" width="90"/>

#### (Marker genes)

- Vector containing gene names
- This is only necessary for BSeq-sc

<img src="../www/markers.png" alt="markers_image" width="90"/>

#### Bulk RNA-seq data

- **Genes** x **Samples** matrix
- Rownames (gene names) are provided in the same format as in the sc RNA-seq data, for instance HGNC symbols

<img src="../www/bulk.png" alt="bulk_image" width="800"/>

#### Signature

Supported data types:

- csv
- tsv
- rds

<span class="help-block">
  For csv and tsv files the first column <strong>must</strong> contain gene identifiers
</span>

<img src="../www/signature.png" alt="signature_image" width="350"/>
