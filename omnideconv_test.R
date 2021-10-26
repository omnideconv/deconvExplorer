library (omnideconv)

single = omnideconv::single_cell_data_1
bulk = omnideconv::bulk
batch = omnideconv::batch_ids_1
cell_type = omnideconv::cell_type_annotations_1

# bisque -running
signature = build_model(single_cell_object = single, cell_type_annotations = cell_type, method = "bisque", batch_ids = batch, bulk_gene_expression = bulk)
deconv_result = omnideconv::deconvolute(bulk_gene_expression = bulk, signature = signature, method = "bisque", single_cell_object = single, cell_type_annotations = cell_type, batch_ids = batch)


# scdc, running, signature returns null but deconvolute is running!
signature = build_model(single_cell_object = single, cell_type_annotations = cell_type, method = "scdc", batch_ids = batch, bulk_gene_expression = bulk)
deconv_result = omnideconv::deconvolute(bulk_gene_expression = bulk, signature = signature, method = "scdc", single_cell_object = single, cell_type_annotations = cell_type, batch_ids = batch)



# momf, running, signature returns null but deconvolute is running!
signature = build_model(single_cell_object = single, cell_type_annotations = cell_type, method = "momf", batch_ids = batch, bulk_gene_expression = bulk)
deconv_result = omnideconv::deconvolute(bulk_gene_expression = bulk, signature = signature, method = "momf", single_cell_object = single, cell_type_annotations = cell_type, batch_ids = batch)

# rcibersort running, need to set credentials before!
signature = build_model(single_cell_object = single, cell_type_annotations = cell_type, method = "bisque", batch_ids = batch, bulk_gene_expression = bulk)
deconv_result = omnideconv::deconvolute(bulk_gene_expression = bulk, signature = signature, method = "cibersortx", single_cell_object = single, cell_type_annotations = cell_type, batch_ids = batch)

# error
signature = build_model(single_cell_object = single, cell_type_annotations = cell_type, method = "music", batch_ids = batch, bulk_gene_expression = bulk)
deconv_result = omnideconv::deconvolute(bulk_gene_expression = bulk, signature = signature, method = "music", single_cell_object = single, cell_type_annotations = cell_type, batch_ids = batch)


# markers missing! error
signature = build_model(single_cell_object = single, cell_type_annotations = cell_type, method = "bseqsc", batch_ids = batch, bulk_gene_expression = bulk)
deconv_result = omnideconv::deconvolute(bulk_gene_expression = bulk, signature = signature, method = "bseqsc", single_cell_object = single, cell_type_annotations = cell_type, batch_ids = batch)


# !
signature = build_model(single_cell_object = single, cell_type_annotations = cell_type, method = "scaden", batch_ids = batch, bulk_gene_expression = bulk)
deconv_result = omnideconv::deconvolute(bulk_gene_expression = bulk, signature = signature, method = "scaden", single_cell_object = single, cell_type_annotations = cell_type, batch_ids = batch)




