.PHONY: deploy logs test undeploy help load

# This script deploys a Google Cloud Function that proxies requests between the client browser and the BigQuery API.

help:
	@echo "Gene Lookup database proxy deployment Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  deploy    - Deploy or update the Google Cloud Function"
	@echo "  logs      - View recent logs from the deployed function"
	@echo "  test      - Run test queries against the deployed function"
	@echo "  undeploy  - Remove the deployed function and its logs"
	@echo "  load      - Load data into BigQuery by running load_bigquery.py"
	@echo "  help      - Show this help message (default)"

# Deploy the Cloud Function (allowed memory settings: 256Mb, 512M, 1024Mi, 2G, 4Gi, etc.)
deploy:
	gcloud functions deploy query_gene_lookup_db \
		--runtime python310 \
		--trigger-http \
		--allow-unauthenticated \
		--entry-point query_gene_lookup_db \
		--region=us-central1 \
		--memory 256MB \
		--gen2

# View recent logs
logs:
	gcloud functions logs read query_gene_lookup_db --region=us-central1 --limit=100

# Alias for logs target
log: logs

# Run test queries
test:
	@echo "Testing query: SELECT COUNT(*) AS c FROM cmg-analysis.gene_lookup.combined_gene_disease_association_table LIMIT 10"
	curl --silent --show-error --fail -X POST "https://us-central1-cmg-analysis.cloudfunctions.net/query_gene_lookup_db" \
		-H "Content-Type: application/json" \
		-d '{ "sql": "SELECT COUNT(*) AS c FROM cmg-analysis.gene_lookup.combined_gene_disease_association_table LIMIT 10" }' | jq '.'

# Undeploy the function and clean up logs
undeploy:
	gcloud functions delete query_gene_lookup_db --region=us-central1

# Load data into BigQuery
load:
	python3 ../annotation_utils/load_bigquery.py
