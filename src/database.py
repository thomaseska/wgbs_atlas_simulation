from urllib import request
import shutil, os

DATA_URL = {
	"hg19_cpgs":"https://www.dropbox.com/scl/fi/k8lk0mbqhxffer0wxzkg9/hg19_cpgs.csv?rlkey=izqwpn87drnr1h63kxqge0njt&st=8oe0ugpo&dl=1",
	"unmethyl_regions":"https://www.dropbox.com/scl/fi/2en8sfz2r7ed95g3kkjul/top250_unmethyl_cell_type.csv?rlkey=4n7itlxigy5vhdlvtq4yhv5dz&st=eab9wfek&dl=1",
	"hg19_genome":"https://www.dropbox.com/scl/fi/egunv56zci8m86zxdbcqf/hg19.pk?rlkey=3ix0bsuqt6lx1418qzhic4t1c&st=c4mqqlr9&dl=1",
	"cell_type_match": "https://github.com/hanyangii/wgbs_atlas_simulation/blob/master/data/cell_type_match.json"
}

DATA_OUT_FILE = {
	"hg19_cpgs":"hg19_cpgs.csv",
	"unmethyl_regions":"top250_unmethyl_cell_type.csv",
	"hg19_genome":"hg19_genome.pk",
	"cell_type_match": "cell_type_match.json"
}

# Find the data dirctory
# assume the code is being run at the designated directory
DATA_DIR = os.path.join(os.getcwd(), "data") 

def download_from_url(url: str, f_path: str) -> None:
	print(f"Downloading {f_path}...")
	with request.urlopen(url) as response, \
		 open(f_path, "wb") as fp:
		shutil.copyfileobj(response, fp)

def download(data: str):
	'''
	Download the given data and return the file path. 
	'''
	if data not in DATA_URL.keys():
		raise ValueError(f"data must be one of {DATA_URL.keys()}. ({data} is given.)")

	f_out_path = os.path.join(DATA_DIR, DATA_OUT_FILE[data])

	if not os.path.exists(f_out_path):
		download_from_url(url = DATA_URL[data], f_path=f_out_path)

	return f_out_path
