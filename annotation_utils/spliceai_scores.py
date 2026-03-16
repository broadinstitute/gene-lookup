import argparse
import requests

def get_spliceai_scores_from_api(chrom, pos, ref, alt):
	variant = f"{chrom}-{pos}-{ref}-{alt}"
	result = requests.get(f"https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=38&distance=500&mask=0&variant={variant}&raw={variant}")
	spliceai_scores = result.json()
	if spliceai_scores.get("error"):
		print(f"WARNING: SpliceAI-lookup returned an error for {variant}: {spliceai_scores['error']}")
		return 0, 0

	spliceai_scores = spliceai_scores.get("scores", [])

	try:
		max_DS_AG = max(float(s["DS_AG"]) for s in spliceai_scores)
	except Exception as e:
		print(f"WARNING: unable to parse DS_AG score for {variant}: {spliceai_scores}")
		max_DS_AG = 0

	try:
		max_DS_DG = max(float(s["DS_DG"]) for s in spliceai_scores)
	except Exception as e:
		print(f"WARNING: unable to parse DS_DG score for {variant}: {spliceai_scores}")
		max_DS_DG = 0

	try:
		max_DS_AL = max(float(s["DS_AL"]) for s in spliceai_scores)
	except Exception as e:
		print(f"WARNING: unable to parse DS_AL score for {variant}: {spliceai_scores}")
		max_DS_AL = 0

	try:
		max_DS_DL = max(float(s["DS_DL"]) for s in spliceai_scores)
	except Exception as e:
		print(f"WARNING: unable to parse DS_DL score for {variant}: {spliceai_scores}")
		max_DS_DL = 0

	spliceai_gain_score = max(max_DS_AG, max_DS_DG)
	spliceai_loss_score = max(max_DS_AL, max_DS_DL)

	return spliceai_gain_score, spliceai_loss_score

#GET https://pangolin-38-xwkwwwxdwq-uc.a.run.app/pangolin/?hg=38&distance=500&mask=0&variant=19-10315746-A-C&raw=chr19:10315746 A>C 200


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("variant")
	args = parser.parse_args()
	chrom, pos, ref, alt = args.variant.split("-")

	spliceai_gain_score, spliceai_loss_score = get_spliceai_scores_from_api(chrom, pos, ref, alt)

	print(f"SpliceAI scores for {chrom}-{pos}-{ref}-{alt}: gain={spliceai_gain_score}, loss={spliceai_loss_score}")

if __name__ == "__main__":
	main()