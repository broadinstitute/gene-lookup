"""Global constants shared between BigQuery data loading and website generation.

This module contains the BigQuery schema definition with column descriptions that serve
as the single source of truth for column documentation across the project.
"""

# Export group names for organizing columns in the export dialog
GROUP_CORE = "Core"
GROUP_CONSTRAINT = "Constraint"
GROUP_OMIM_CLINVAR = "OMIM & ClinVar"
GROUP_CLINGEN_GENCC = "ClinGen & GenCC"
GROUP_PANEL_APP_UK = "PanelApp UK"
GROUP_PANEL_APP_AU = "PanelApp AU"
GROUP_DECIPHER_OTHER = "Decipher & Other Sources"
GROUP_DBNSFP = "dbNSFP"

GROUP_ORDER = [
    GROUP_CORE,
    GROUP_CONSTRAINT,
    GROUP_OMIM_CLINVAR,
    GROUP_CLINGEN_GENCC,
    GROUP_PANEL_APP_UK,
    GROUP_PANEL_APP_AU,
    GROUP_DECIPHER_OTHER,
    GROUP_DBNSFP,
]

BIGQUERY_COLUMNS = [
    # Core
    {
        "type": "STRING",
        "name": "gene_symbol",
        "description": "HGNC-approved gene symbol (uppercase).",
        "displayName": "Gene Symbol",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CORE,
    },
    {
        "type": "STRING",
        "name": "ensembl_gene_id",
        "description": "Ensembl gene ID (ENSG).",
        "displayName": "Ensembl Gene ID",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CORE,
    },
    {
        "type": "STRING",
        "name": "gene_aliases",
        "description": "Alternative gene symbols from HGNC (uppercase, comma-separated).",
        "displayName": "Gene Aliases",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CORE,
    },
    {
        "type": "STRING",
        "name": "hgnc_gene_id",
        "description": "HGNC gene ID (e.g., HGNC:1234).",
        "displayName": "HGNC Gene ID",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CORE,
    },
    {
        "type": "STRING",
        "name": "gene_chrom",
        "description": "Chromosome name.",
        "displayName": "Chromosome",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CORE,
    },
    {
        "type": "INTEGER",
        "name": "gene_start",
        "description": "Gene start coordinate (GRCh38).",
        "displayName": "Start",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CORE,
    },
    {
        "type": "INTEGER",
        "name": "gene_end",
        "description": "Gene end coordinate (GRCh38).",
        "displayName": "End",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CORE,
    },
    {
        "type": "STRING",
        "name": "inheritance",
        "description": "Summarized inheritance modes across all sources (semicolon-separated).",
        "displayName": "Inheritance",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CORE,
    },

    # AI Summary
    {
        "type": "STRING",
        "name": "LLM_phenotype_summary",
        "description": "LLM-generated concise summary of phenotypes from all sources.",
        "displayName": "Phenotype Summary",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CORE,
    },

    # Sources
    {
        "type": "STRING",
        "name": "sources",
        "description": "Number and list of sources that include this gene (e.g., '3: OMIM, ClinGen, GenCC').",
        "displayName": "Sources",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CORE,
    },

    # Constraint
    {
        "type": "FLOAT",
        "name": "pLI_v2",
        "description": "Probability of loss-of-function intolerance from gnomAD v2.",
        "displayName": "pLI (v2)",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CONSTRAINT,
    },
    {
        "type": "FLOAT",
        "name": "pLI_v4",
        "description": "Probability of loss-of-function intolerance from gnomAD v4.",
        "displayName": "pLI (v4)",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CONSTRAINT,
    },
    {
        "type": "FLOAT",
        "name": "lof_oe_ci_upper_v4",
        "description": "Upper bound of the 90% confidence interval for the loss-of-function observed/expected ratio from gnomAD v4 (LOEUF).",
        "displayName": "LOEUF (v4)",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CONSTRAINT,
    },
    {
        "type": "FLOAT",
        "name": "mis_oe_ci_upper_v4",
        "description": "Upper bound of the 90% confidence interval for the missense observed/expected ratio from gnomAD v4 (MOEUF).",
        "displayName": "MOEUF (v4)",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CONSTRAINT,
    },
    {
        "type": "FLOAT",
        "name": "s_het",
        "description": "Gene constraint score from GeneBayes (Zeng et al. 2024). Scores >0.1 indicate high likelihood of extreme selection.",
        "displayName": "s_het",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CONSTRAINT,
    },

    # OMIM
    {
        "type": "STRING",
        "name": "OMIM_mim_number",
        "description": "OMIM gene MIM number(s) (semicolon-separated if multiple).",
        "displayName": "OMIM MIM Number",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_OMIM_CLINVAR,
    },
    {
        "type": "STRING",
        "name": "OMIM_phenotype_mim_number",
        "description": "OMIM phenotype MIM number(s) (semicolon-separated if multiple).",
        "displayName": "OMIM Phenotype MIM Number",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_OMIM_CLINVAR,
    },
    {
        "type": "STRING",
        "name": "OMIM_inheritance",
        "description": "Inheritance mode(s) from OMIM (semicolon-separated if multiple).",
        "displayName": "OMIM Inheritance",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_OMIM_CLINVAR,
    },
    {
        "type": "STRING",
        "name": "OMIM_phenotype_description",
        "description": "Phenotype description(s) from OMIM (semicolon-separated if multiple).",
        "displayName": "OMIM Phenotype",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_OMIM_CLINVAR,
    },

    # ClinGen
    {
        "type": "STRING",
        "name": "CLINGEN_disease_label",
        "description": "Disease name(s) from ClinGen gene-disease validity curation (semicolon-separated if multiple).",
        "displayName": "ClinGen Disease",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CLINGEN_GENCC,
    },
    {
        "type": "STRING",
        "name": "CLINGEN_disease_mondo_id",
        "description": "MONDO disease ID(s) from ClinGen (semicolon-separated if multiple).",
        "displayName": "ClinGen MONDO ID",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CLINGEN_GENCC,
    },
    {
        "type": "STRING",
        "name": "CLINGEN_inheritance",
        "description": "Mode of inheritance from ClinGen (semicolon-separated if multiple).",
        "displayName": "ClinGen Inheritance",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CLINGEN_GENCC,
    },
    {
        "type": "STRING",
        "name": "CLINGEN_classification",
        "description": "ClinGen gene-disease validity classification (Definitive, Strong, Moderate, or Limited).",
        "displayName": "ClinGen Classification",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CLINGEN_GENCC,
    },
    {
        "type": "STRING",
        "name": "CLINGEN_haploinsufficient",
        "description": "ClinGen haploinsufficiency classification for this gene.",
        "displayName": "ClinGen Haploinsufficiency",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CLINGEN_GENCC,
    },

    # GenCC
    {
        "type": "STRING",
        "name": "GENCC_disease_name",
        "description": "Disease name(s) from GenCC (semicolon-separated if multiple).",
        "displayName": "GenCC Disease",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CLINGEN_GENCC,
    },
    {
        "type": "STRING",
        "name": "GENCC_classification",
        "description": "Gene-disease classification(s) from GenCC (semicolon-separated if multiple).",
        "displayName": "GenCC Classification",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CLINGEN_GENCC,
    },
    {
        "type": "STRING",
        "name": "GENCC_inheritance",
        "description": "Mode of inheritance from GenCC (semicolon-separated if multiple).",
        "displayName": "GenCC Inheritance",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_CLINGEN_GENCC,
    },

    # PanelApp UK
    {
        "type": "STRING",
        "name": "PANEL_APP_UK_confidence",
        "description": "PanelApp UK confidence level(s) for this gene (semicolon-separated if multiple panels).",
        "displayName": "PanelApp UK Confidence",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_UK,
    },
    {
        "type": "STRING",
        "name": "PANEL_APP_UK_penetrance",
        "description": "Penetrance from PanelApp UK (semicolon-separated if multiple panels).",
        "displayName": "PanelApp UK Penetrance",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_UK,
    },
    {
        "type": "STRING",
        "name": "PANEL_APP_UK_mode_of_pathogenicity",
        "description": "Mode of pathogenicity from PanelApp UK (semicolon-separated if multiple panels).",
        "displayName": "PanelApp UK Pathogenicity",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_UK,
    },
    {
        "type": "STRING",
        "name": "PANEL_APP_UK_inheritance",
        "description": "Mode of inheritance from PanelApp UK (semicolon-separated if multiple panels).",
        "displayName": "PanelApp UK Inheritance",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_UK,
    },
    {
        "type": "STRING",
        "name": "PANEL_APP_UK_evidence",
        "description": "Evidence level from PanelApp UK (semicolon-separated if multiple panels).",
        "displayName": "PanelApp UK Evidence",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_UK,
    },
    {
        "type": "STRING",
        "name": "PANEL_APP_UK_phenotypes",
        "description": "Phenotype(s) from PanelApp UK (semicolon-separated if multiple panels).",
        "displayName": "PanelApp UK Phenotypes",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_UK,
    },
    {
        "type": "STRING",
        "name": "PANEL_APP_UK_panel_name",
        "description": "Panel name(s) from PanelApp UK (semicolon-separated if multiple panels).",
        "displayName": "PanelApp UK Panel",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_UK,
    },

    # PanelApp AU
    {
        "type": "STRING",
        "name": "PANEL_APP_AU_confidence",
        "description": "PanelApp Australia confidence level(s) for this gene (semicolon-separated if multiple panels).",
        "displayName": "PanelApp AU Confidence",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_AU,
    },
    {
        "type": "STRING",
        "name": "PANEL_APP_AU_penetrance",
        "description": "Penetrance from PanelApp Australia (semicolon-separated if multiple panels).",
        "displayName": "PanelApp AU Penetrance",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_AU,
    },
    {
        "type": "STRING",
        "name": "PANEL_APP_AU_mode_of_pathogenicity",
        "description": "Mode of pathogenicity from PanelApp Australia (semicolon-separated if multiple panels).",
        "displayName": "PanelApp AU Pathogenicity",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_AU,
    },
    {
        "type": "STRING",
        "name": "PANEL_APP_AU_inheritance",
        "description": "Mode of inheritance from PanelApp Australia (semicolon-separated if multiple panels).",
        "displayName": "PanelApp AU Inheritance",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_AU,
    },
    {
        "type": "STRING",
        "name": "PANEL_APP_AU_evidence",
        "description": "Evidence level from PanelApp Australia (semicolon-separated if multiple panels).",
        "displayName": "PanelApp AU Evidence",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_AU,
    },
    {
        "type": "STRING",
        "name": "PANEL_APP_AU_phenotypes",
        "description": "Phenotype(s) from PanelApp Australia (semicolon-separated if multiple panels).",
        "displayName": "PanelApp AU Phenotypes",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_AU,
    },
    {
        "type": "STRING",
        "name": "PANEL_APP_AU_panel_name",
        "description": "Panel name(s) from PanelApp Australia (semicolon-separated if multiple panels).",
        "displayName": "PanelApp AU Panel",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_PANEL_APP_AU,
    },

    # Decipher
    {
        "type": "STRING",
        "name": "DECIPHER_inheritance",
        "description": "Mode of inheritance from DECIPHER (semicolon-separated if multiple).",
        "displayName": "Decipher Inheritance",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DECIPHER_OTHER,
    },
    {
        "type": "STRING",
        "name": "DECIPHER_disease_names",
        "description": "Disease name(s) from DECIPHER (semicolon-separated if multiple).",
        "displayName": "Decipher Disease",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DECIPHER_OTHER,
    },

    # ClinVar
    {
        "type": "STRING",
        "name": "CLINVAR_phenotypes",
        "description": "Phenotype(s) from ClinVar for pathogenic/likely pathogenic variants in this gene.",
        "displayName": "ClinVar Phenotypes",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_OMIM_CLINVAR,
    },
    {
        "type": "STRING",
        "name": "CLINVAR_clinical_significance",
        "description": "Clinical significance categories from ClinVar for variants in this gene.",
        "displayName": "ClinVar Significance",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_OMIM_CLINVAR,
    },
    {
        "type": "STRING",
        "name": "CLINVAR_stars",
        "description": "ClinVar review status gold stars for variants in this gene.",
        "displayName": "ClinVar Stars",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_OMIM_CLINVAR,
    },
    {
        "type": "STRING",
        "name": "CLINVAR_variant_consequences",
        "description": "Major variant consequence types from ClinVar for pathogenic/likely pathogenic variants in this gene.",
        "displayName": "ClinVar Consequences",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_OMIM_CLINVAR,
    },

    # Fridman
    {
        "type": "STRING",
        "name": "FRIDMAN_omim_phenotype_id",
        "description": "OMIM phenotype ID(s) from Fridman et al. 2025 list of recessive disease genes.",
        "displayName": "Fridman 2025 OMIM ID",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DECIPHER_OTHER,
    },
    {
        "type": "STRING",
        "name": "FRIDMAN_phenotype_category",
        "description": "Disorder group from Fridman et al. 2025 (semicolon-separated if multiple).",
        "displayName": "Fridman 2025 Category",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DECIPHER_OTHER,
    },
    {
        "type": "STRING",
        "name": "FRIDMAN_inheritance",
        "description": "Inheritance mode from Fridman et al. 2025 (AR or AR-AD).",
        "displayName": "Fridman 2025 Inheritance",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DECIPHER_OTHER,
    },

    # dbNSFP
    {
        "type": "STRING",
        "name": "DBNSFP_orphanet_disorder_id",
        "description": "Orphanet disorder ID(s) from dbNSFP (semicolon-separated if multiple).",
        "displayName": "Orphanet Disorder ID",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DECIPHER_OTHER,
    },
    {
        "type": "STRING",
        "name": "DBNSFP_orphanet_disorder",
        "description": "Orphanet disorder name(s) from dbNSFP (semicolon-separated if multiple).",
        "displayName": "Orphanet Disorder",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DECIPHER_OTHER,
    },
    {
        "type": "STRING",
        "name": "DBNSFP_orphanet_association_type",
        "description": "Orphanet gene-disease association type(s) from dbNSFP (semicolon-separated if multiple).",
        "displayName": "Orphanet Association Type",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DECIPHER_OTHER,
    },
    {
        "type": "STRING",
        "name": "DBNSFP_hpo_id",
        "description": "Human Phenotype Ontology (HPO) ID(s) from dbNSFP (semicolon-separated if multiple).",
        "displayName": "HPO ID",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DBNSFP,
    },
    {
        "type": "STRING",
        "name": "DBNSFP_hpo_name",
        "description": "Human Phenotype Ontology (HPO) term name(s) from dbNSFP (semicolon-separated if multiple).",
        "displayName": "HPO Name",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DBNSFP,
    },
    {
        "type": "STRING",
        "name": "DBNSFP_function_description",
        "description": "Gene function description from dbNSFP/UniProt.",
        "displayName": "Function Description",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DBNSFP,
    },
    {
        "type": "STRING",
        "name": "DBNSFP_disease_description",
        "description": "Disease description from dbNSFP/UniProt.",
        "displayName": "Disease Description",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DBNSFP,
    },
    {
        "type": "STRING",
        "name": "DBNSFP_pathway_kegg",
        "description": "KEGG pathway(s) from dbNSFP (semicolon-separated if multiple).",
        "displayName": "Pathway (KEGG)",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DBNSFP,
    },
    {
        "type": "STRING",
        "name": "DBNSFP_pathway_uniprot",
        "description": "UniProt pathway(s) from dbNSFP (semicolon-separated if multiple).",
        "displayName": "Pathway (UniProt)",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DBNSFP,
    },
    {
        "type": "FLOAT",
        "name": "DBNSFP_p_rec",
        "description": "Probability of being a recessive disease gene from dbNSFP.",
        "displayName": "P(rec)",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DBNSFP,
    },
    {
        "type": "STRING",
        "name": "DBNSFP_known_rec_info",
        "description": "Known recessive disease gene information from dbNSFP.",
        "displayName": "Known Recessive Info",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DBNSFP,
    },
    {
        "type": "STRING",
        "name": "DBNSFP_essential_gene",
        "description": "Essential gene annotation from dbNSFP (E=essential, N=non-essential).",
        "displayName": "Essential Gene",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DBNSFP,
    },
    {
        "type": "STRING",
        "name": "DBNSFP_mgi_mouse_phenotype",
        "description": "Mouse phenotype(s) from MGI via dbNSFP (semicolon-separated if multiple).",
        "displayName": "Mouse Phenotype (MGI)",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DBNSFP,
    },
    {
        "type": "STRING",
        "name": "DBNSFP_zfin_zebrafish_phenotype",
        "description": "Zebrafish phenotype tag(s) from ZFIN via dbNSFP (semicolon-separated if multiple).",
        "displayName": "Zebrafish Phenotype (ZFIN)",
        "allowCustomFilter": True,
        "allowExport": True,
        "group": GROUP_DBNSFP,
    },
]


def get_column_descriptions():
    """Return a dictionary mapping column names to their descriptions.

    Returns:
        dict: Mapping of column name to description string.
    """
    return {col["name"]: col.get("description", "") for col in BIGQUERY_COLUMNS}


def get_custom_filter_columns():
    """Return a list of columns that can be used for custom filtering.

    Only columns that have 'allowCustomFilter' set to True are included.

    Returns:
        list: List of dicts with keys: name, displayName, group, description, type
    """
    filterable = []
    for col in BIGQUERY_COLUMNS:
        if col.get("allowCustomFilter", False):
            filterable.append({
                "name": col.get("name", ""),
                "displayName": col.get("displayName", col.get("name", "")),
                "group": col.get("group", ""),
                "description": col.get("description", ""),
                "type": col.get("type", ""),
            })
    return filterable


def get_exportable_columns():
    """Return a list of columns that can be exported, with their display metadata.

    Only columns that have 'allowExport' set to True are included.

    Returns:
        list: List of dicts with keys: name, displayName, group, column, description
              where column is the SQL expression to use (defaults to name).
    """
    exportable = []
    for col in BIGQUERY_COLUMNS:
        if col.get("allowExport", False):
            exportable.append({
                "name": col.get("name", ""),
                "displayName": col.get("displayName", col.get("name", "")),
                "group": col.get("group", ""),
                "column": col.get("exportColumn", col.get("name", "")),
                "description": col.get("description", ""),
            })
    return exportable
